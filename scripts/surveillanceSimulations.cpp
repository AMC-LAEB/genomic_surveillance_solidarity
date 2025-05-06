#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix getOnsets(IntegerMatrix& incidence_mat) {
  IntegerMatrix onsets(2, incidence_mat.rows());
  onsets.fill(-1);
  for (int i = 0; i < incidence_mat.rows(); i++ ){
    int idx_1 = -1;
    int idx_10 = -1;
    for (int j = 0; j < incidence_mat.cols(); j++){
      if ((incidence_mat(i,j) > 0) & (idx_1 == -1)){
        idx_1 = j;
        onsets(0,i) = idx_1;}
      if ((incidence_mat(i,j) >= 10) & (idx_10 == -1)){
        idx_10 = j;
        onsets(1,i) = idx_10;
        break;
      }
    }
  }
  return(onsets);
}

// [[Rcpp::export]]
int getOnsetCountry(IntegerMatrix& onsets) {
  int cur_min_onset = 1000;
  int cur_min_onset_idx = -1;
  for (int i = 0; i < onsets.cols(); i++){
    if ((onsets(0,i) < cur_min_onset) & (onsets(0,i) >= 0)){
      cur_min_onset = onsets(0,i);
      cur_min_onset_idx = i;
    }
  }
  return(cur_min_onset_idx);
}

// [[Rcpp::export]]
IntegerVector getInfectedCountries(IntegerMatrix& onsets, int day) {
  IntegerVector infected_countries;
  for (int i = 0; i < onsets.cols(); i++){
    if (onsets(0,i) <= day & onsets(0,i) >= 0){
      infected_countries.push_back(i);
    }
  }
  return(infected_countries);
}

// [[Rcpp::export]]
NumericMatrix getWildtypeCasesAndVariantProportion(int ctr, 
                                                   int nday,
                                                   IntegerMatrix& incidence_mat,
                                                   IntegerMatrix& onsets,
                                                   IntegerMatrix& wildtype_epidemic_mat){
  NumericMatrix wt_and_p_mat(2,nday);
  wt_and_p_mat.fill(0);
  int plus_10_day = onsets(1,ctr);
  IntegerVector wildtype_epidemic = wildtype_epidemic_mat.row(ctr);
  double base_inc = wildtype_epidemic[0];
  if (plus_10_day != -1){
    for (int i = 0; i < plus_10_day; i++){
      wt_and_p_mat(0,i) = base_inc;
    }
    int cnt = 0;
    for (int i = plus_10_day; i < nday; i++){
      wt_and_p_mat(0,i) = wildtype_epidemic[cnt];
      cnt++;
    }
    for (int i = 0; i < nday; i++){
      wt_and_p_mat(1,i) = incidence_mat(ctr,i) / (incidence_mat(ctr,i) +  wt_and_p_mat(0,i));
    }
  } else {
    for (int i = 0; i < nday; i++){
      wt_and_p_mat(0,i) = base_inc;
    }
    for (int i = 0; i < nday; i++){
      wt_and_p_mat(1,i) = incidence_mat(ctr,i) / (incidence_mat(ctr,i) +  wt_and_p_mat(0,i));
    }
  }
  return(wt_and_p_mat);
}

// [[Rcpp::export]]
IntegerMatrix getNumberOfSeqs(int day,
                              int population_size,
                              int onset_day,
                              NumericVector& seqrate_by_tat,
                              NumericMatrix& wt_and_p_mat){
  int time_since_onset = day - onset_day + 1;
  IntegerMatrix n_seq(2,time_since_onset);
  n_seq.fill(0);
  NumericVector true_seqrates(time_since_onset);
  for (int i = 0; i < time_since_onset; i++){
    true_seqrates[i] = seqrate_by_tat[i] * (population_size / 1e6) / 7;
  }
  NumericVector seqrates_rev = rev(true_seqrates);
  for (int i = 0; i < time_since_onset; i++){
    int number_of_sequences = R::rpois(seqrates_rev[i]);
    int n_infections = wt_and_p_mat(0,i + onset_day);
    double variant_proportion = wt_and_p_mat(1,i + onset_day);
    if (number_of_sequences > n_infections){
      number_of_sequences = n_infections;
    }
    int nseq_mt = R::rbinom(number_of_sequences,variant_proportion);
    int nseq_wt = number_of_sequences - nseq_mt;
    n_seq(0,i) = nseq_mt;
    n_seq(1,i) = nseq_wt;
  }
  return(n_seq);
}

// [[Rcpp::export]]
bool checkVariantProportion(int day, 
                            int ctr,
                            double prop_thresh,
                            IntegerMatrix& mt_count_matrix, 
                            IntegerMatrix& wt_count_matrix){
  IntegerVector all_days;
  for (int i = 1; i <= day; i++){
    if (i % 7 == 0){
      all_days.push_back(i);
    }
  }
  for (int test_day_idx = 0; test_day_idx < all_days.size(); test_day_idx++){
    int test_day = all_days[test_day_idx];
    int mt_sum = 0;
    for (int i = test_day - 6; i <= test_day; i++){
      mt_sum += mt_count_matrix(ctr,i);
    }
    if (mt_sum == 0){continue;}
    int wt_sum = 0;
    for (int i = test_day - 6; i <= test_day; i++){
      wt_sum += wt_count_matrix(ctr,i);
    }
    double pvalue = R::pbinom(mt_sum-1, mt_sum+wt_sum, prop_thresh, false, false);
    if (pvalue < 0.05){return(true);}
  }
  return(false);
}


// [[Rcpp::export]]
IntegerVector getStatsSingle(IntegerMatrix& incidence_mat, 
                             bool computeThreshold,
                             NumericVector& popsize_vec,
                             NumericMatrix& seqrate_mat,
                             IntegerMatrix& wildtype_epidemic_mat,
                             IntegerVector& seqrate_changed,
                             int run_freq){
  bool isDetected = false;
  bool isThreshold = false;
  
  int day_of_minimum_threshold = -1;
  int country_of_minimum_threshold = -1;
  int incidence_on_day_of_minimum_threshold = -1;
  int day_of_detection = -1;
  int country_of_detection = -1;
  int incidence_on_day_of_detection = -1;
  
  IntegerMatrix onsets = getOnsets(incidence_mat);
  int index_country = getOnsetCountry(onsets);
  
  int nday = incidence_mat.cols();
  IntegerMatrix mt_count_matrix = IntegerMatrix(196, nday);
  IntegerMatrix wt_count_matrix = IntegerMatrix(196, nday);
  
  for (int day = 0; day < nday; day++){
    IntegerVector infected_countries = getInfectedCountries(onsets, day);
    for (int ctr_idx = 0; ctr_idx < infected_countries.size(); ctr_idx++){
      int ctr = infected_countries[ctr_idx];
      int onset_day = onsets(0,ctr);
      NumericVector seqrate_by_tat = seqrate_mat.row(ctr);
      double population_size = popsize_vec[ctr];
      double prop_thresh = 0.01;
      if (population_size > 1e8){
        prop_thresh = 0.01 * (1e8 / population_size);
      }
      NumericMatrix wt_and_p_mat = getWildtypeCasesAndVariantProportion(ctr, 
                                                                        nday, 
                                                                        incidence_mat,
                                                                        onsets, 
                                                                        wildtype_epidemic_mat);
      IntegerMatrix number_of_seqs_today = getNumberOfSeqs(day,
                                                           population_size,
                                                           onset_day,
                                                           seqrate_by_tat,
                                                           wt_and_p_mat);
      for (int i = 0; i < day - onset_day + 1; i++){
        mt_count_matrix(ctr,onset_day+i) += number_of_seqs_today(0,i);
        wt_count_matrix(ctr,onset_day+i) += number_of_seqs_today(1,i);
      }
      
      if ((seqrate_changed[ctr] == false) | ((seqrate_changed[ctr] == true) & (day % run_freq == 0))){
        if ((isThreshold == false) & computeThreshold){
          bool exceedsProportion = checkVariantProportion(day,
                                                          ctr,
                                                          prop_thresh,
                                                          mt_count_matrix,
                                                          wt_count_matrix);
          if (exceedsProportion){
            day_of_minimum_threshold = day;
            country_of_minimum_threshold = ctr;
            for (int i = 0; i <= day; i++){
              incidence_on_day_of_minimum_threshold += sum(incidence_mat(_,i));
            }
            isThreshold = true;
          }
        }
        int total_seqs = sum(mt_count_matrix(ctr,_));
        if (total_seqs > 0 & isDetected == false){
          day_of_detection = day;
          country_of_detection = ctr;
          for (int i = 0; i <= day; i++){
            incidence_on_day_of_detection += sum(incidence_mat(_,i));
          }
          isDetected = true;
        }
      }
      if ((isDetected & isThreshold) | (isDetected & (computeThreshold == false))){
        IntegerVector out;
        out.push_back(day_of_detection+1);
        out.push_back(country_of_detection+1);
        out.push_back(incidence_on_day_of_detection);
        out.push_back(day_of_minimum_threshold+1);
        out.push_back(country_of_minimum_threshold+1);
        out.push_back(incidence_on_day_of_minimum_threshold);
        out.push_back(index_country+1);
        return(out);
      }
    }
  }
  IntegerVector out;
  out.push_back(day_of_detection+1);
  out.push_back(country_of_detection+1);
  out.push_back(incidence_on_day_of_detection);
  out.push_back(day_of_minimum_threshold+1);
  out.push_back(country_of_minimum_threshold+1);
  out.push_back(incidence_on_day_of_minimum_threshold);
  out.push_back(index_country+1);
  return(out);
}

// [[Rcpp::export]]
IntegerMatrix getTimeToDetection(List simulations,
                                 bool computeThreshold,
                                 NumericVector& popsize_vec,
                                 NumericMatrix& seqrate_mat,
                                 IntegerMatrix& wildtype_epidemic_mat,
                                 IntegerVector& seqrate_changed,
                                 int run_freq){
  IntegerMatrix out(simulations.size(),7);
  for (int i = 0; i < simulations.size(); i++){
    IntegerMatrix incidence_mat = simulations[i];
    IntegerVector time_to_detection = getStatsSingle(incidence_mat,
                                                     computeThreshold,
                                                     popsize_vec,
                                                     seqrate_mat,
                                                     wildtype_epidemic_mat,
                                                     seqrate_changed,
                                                     run_freq);
    out(i,_) = time_to_detection;
  }
  return(out);
}


