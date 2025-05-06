#!/bin/bash

run_func() {
    x=$1
    Rscript -e "source(\"runSurveillanceSimulations.R\"); runAllSurveillanceSimulations($x); flush.console()"
}

export -f run_func

echo {1..144} | tr ' ' '\n' | parallel --ungroup -j80 run_func

wait

Rscript -e 'saveRDS(do.call(rbind,lapply(1:144,function(x)readRDS(paste0("detection_outputs/GS_sim_",x,".rds")))),"combined_outputs.rds")'

wait

Rscript getLeadTimes.R