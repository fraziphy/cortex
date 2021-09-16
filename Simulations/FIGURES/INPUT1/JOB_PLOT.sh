#!/bin/bash
##############################################################
#                              BASH JOB                      #
#                                                            #
##############################################################
# 
# 
# 
# 
SIGNALS="lfp v firing"
POPS="Perturbed Unperturbed"
SUB_POPS="pyr inh"
# 
# 
# 
# 
# 
# 
# 
CONDS="ONE_CORTICAL_COLUMN TWO_CORTICAL_COLUMN"
POP=Perturbed
SUB_POP=pyr
rootMACHINE=$(pwd)
# 
#
# 
# 
export rootMACHINE
# 
#
# 
# 
python3 ../DYNAMICAL_ANALYSIS.py
# 
#
# 
# 
for SUB_POP in $SUB_POPS; do
    COND="ONE_CORTICAL_COLUMN"
    POP=Perturbed
    export COND SUB_POP POP
    python3 ../PRESTIMULUS_PLOT.py
    python3 ../POSTSTIMULUS_PLOT_ONE_CORTICAL_COLUMN.py
    # 
    # 
    # 
    # 
    COND="TWO_CORTICAL_COLUMN"
    python3 ../POSTSTIMULUS_PLOT_TWO_CORTICAL_COLUMN.py
    python3 ../POSTSTIMULUS_PLOT_UNPERTURBED_CORTICAL_COLUMN.py
    python3 ../POSTSTIMULUS_PLOT_UNPERTURBED_CORTICAL_COLUMN_SCREENING.py
    python3 ../POSTSTIMULUS_PLOT_UNPERTURBED_CORTICAL_COLUMN_SCREENING_Bonferroni.py
    for POP in $POPS; do 
        export COND POP
        python3 ../PRESTIMULUS_PLOT.py
    done
done
