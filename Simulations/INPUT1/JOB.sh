#!/bin/bash
##############################################################
#                            BASH JOB                        #
#                                                            #
#                                                            #
#                   Containing SIMULATION sbatch,            #
#                statistical test, and data curation         # 
##############################################################
# 
# 
# 
# 
SIGNALS="lfp v firing"
POPS="Perturbed Unperturbed"
SUB_POP="pyr inh"
# 
# 
# 
# 
# 
# 
# 
# 
rootCLUSTER=home/frazi
# 
# 
# 
# 
root=PROJECTA/INPUT1/
mkdir -p /$rootCLUSTER/scratch/$root
mkdir -p /$rootCLUSTER/FIGURE/$root
jobs_one_group=""
jobs_two_group=""
jobs_one_group_pre=""
jobs_two_group_pre=""
declare -A JOB_G_AMPA_BETA
# 
# 
# 
# 
G_AMPA_STEP=4
G_AMPA_END=10
# 
# 
# 
# 
BETA_STEP=1
BETA_END=5
# 
# 
# 
# 
jid1=$(sbatch -J 1_0 -N 1 -n 20 SIMULATION.sbatch 1_0 $root $rootCLUSTER)
jid2=$(sbatch --dependency=aftercorr:${jid1##* } ../SCRIPT/NEW_G_GABA.sbatch ${G_AMPA_STEP} ${G_AMPA_END} ${BETA_STEP} ${BETA_END} $root $rootCLUSTER)
jid3=$(sbatch --dependency=aftercorr:${jid2##* } -J 2_0 -N 1 -n 20 SIMULATION.sbatch 2_0 $root $rootCLUSTER)
jid4=$(sbatch --dependency=aftercorr:${jid2##* } -J 1_1 -N 1 -n 20 SIMULATION.sbatch 1_1 $root $rootCLUSTER)
# 
# 
# 
# 
for g_ampa in $(seq 2 ${G_AMPA_STEP} ${G_AMPA_END}); do 
    for beta in $(seq 1 ${BETA_STEP} ${BETA_END}); do
        jid5=$(sbatch --dependency=aftercorr:${jid2##* } -J ${g_ampa}_${beta} -N 1 -n 20 SIMULATION.sbatch ${g_ampa}_${beta} $root $rootCLUSTER)
        JOB_G_AMPA_BETA[${g_ampa}_${beta}]="${jid5##* }"
    done
done
# 
# 
# 
# 
for sub_pop in $SUB_POP;do
    for signal in $SIGNALS;do
        jid6=$(sbatch --dependency=aftercorr:${jid1##* } -J 1_0_${sub_pop}_${signal} ../SCRIPT/NON_PARAMETRIC_TEST_WITHIN_POP.sbatch 1_0 Perturbed $sub_pop $signal $root $rootCLUSTER)
        jid7=$(sbatch --dependency=aftercorr:${jid3##* } -J 2_0_${sub_pop}_${signal} ../SCRIPT/NON_PARAMETRIC_TEST_WITHIN_POP.sbatch 2_0 Perturbed $sub_pop $signal $root $rootCLUSTER)
        jid8=$(sbatch --dependency=aftercorr:${jid1##* },${jid3##* } -J 2_0_${sub_pop}_${signal} ../SCRIPT/NON_PARAMETRIC_TEST_ACROSS_STATE.sbatch 2_0 Perturbed $sub_pop $signal $root $rootCLUSTER)
        jobs_one_group="${jobs_one_group}${jid6##* },${jid7##* },${jid8##* },"
        sleep 2
    done
done
# 
# 
# 
# 
for pop in $POPS;do
    for sub_pop in $SUB_POP;do
        for signal in $SIGNALS;do
            jid9=$(sbatch --dependency=aftercorr:${jid4##* } -J 1_1_${pop}_${sub_pop}_${signal} ../SCRIPT/NON_PARAMETRIC_TEST_WITHIN_POP.sbatch 1_1 $pop $sub_pop $signal $root $rootCLUSTER)
            jobs_two_group="${jobs_two_group}${jid9##* },"
            for g_ampa in $(seq 2 ${G_AMPA_STEP} ${G_AMPA_END}); do 
                for beta in $(seq 1 ${BETA_STEP} ${BETA_END}); do
                    jid10=$(sbatch --dependency=aftercorr:${JOB_G_AMPA_BETA[${g_ampa}_${beta}]} -J ${g_ampa}_${beta}_${pop}_${sub_pop}_${signal} ../SCRIPT/NON_PARAMETRIC_TEST_WITHIN_POP.sbatch ${g_ampa}_${beta} $pop $sub_pop $signal $root $rootCLUSTER)
                    jobs_two_group="${jobs_two_group}${jid10##* },"
                    sleep 2
                done
            done
        done
    done
done
# 
# 
# 
# 
for sub_pop in ${SUB_POP}; do
    jid16=$(sbatch --dependency=aftercorr:${jid1##* } -J 1_0_${sub_pop}_Perturbed ../SCRIPT/PRESTIMULUS_PREPARATION_DATA.sbatch 1_0 ${sub_pop} Perturbed $root $rootCLUSTER)
    jid17=$(sbatch --dependency=aftercorr:${jid3##* } -J 2_0_${sub_pop}_Perturbed ../SCRIPT/PRESTIMULUS_PREPARATION_DATA.sbatch 2_0 ${sub_pop} Perturbed $root $rootCLUSTER)
    jobs_one_group_pre="${jobs_one_group_pre}${jid16##* },${jid17##* },"
    for pop in $POPS; do 
        jid18=$(sbatch --dependency=aftercorr:${jid4##* } -J 1_1_${sub_pop}_${pop} ../SCRIPT/PRESTIMULUS_PREPARATION_DATA.sbatch 1_1 ${sub_pop} ${pop} $root $rootCLUSTER)
        jobs_two_group_pre="${jobs_two_group_pre}${jid18##* },"
        for g_ampa in $(seq 2 ${G_AMPA_STEP} ${G_AMPA_END}); do 
            for beta in $(seq 1 ${BETA_STEP} ${BETA_END}); do
                jid19=$(sbatch --dependency=aftercorr:${JOB_G_AMPA_BETA[${g_ampa}_${beta}]} -J ${g_ampa}_${beta}_${sub_pop}_${pop} ../SCRIPT/PRESTIMULUS_PREPARATION_DATA.sbatch ${g_ampa}_${beta} ${sub_pop} ${pop} $root $rootCLUSTER)
                jobs_two_group_pre="${jobs_two_group_pre}${jid19##* },"
                sleep 2
            done
        done
    done
done
# 
# 
# 
# 
sbatch --dependency=aftercorr:${jobs_one_group_pre%%,} -J ONE_CORTICAL_COLUMN ../SCRIPT/PRESTIMULUS_DATA.sbatch ONE_CORTICAL_COLUMN ${G_AMPA_STEP} ${G_AMPA_END} ${BETA_STEP} ${BETA_END} $root $rootCLUSTER
sbatch --dependency=aftercorr:${jobs_two_group_pre%%,} -J TWO_CORTICAL_COLUMN ../SCRIPT/PRESTIMULUS_DATA.sbatch TWO_CORTICAL_COLUMN ${G_AMPA_STEP} ${G_AMPA_END} ${BETA_STEP} ${BETA_END} $root $rootCLUSTER
sbatch --dependency=aftercorr:${jobs_one_group%%,} -J ONE_CORTICAL_COLUMN ../SCRIPT/POSTSTIMULUS_DATA.sbatch ONE_CORTICAL_COLUMN ${G_AMPA_STEP} ${G_AMPA_END} ${BETA_STEP} ${BETA_END} $root $rootCLUSTER
sbatch --dependency=aftercorr:${jobs_two_group%%,} -J TWO_CORTICAL_COLUMN ../SCRIPT/POSTSTIMULUS_DATA.sbatch TWO_CORTICAL_COLUMN ${G_AMPA_STEP} ${G_AMPA_END} ${BETA_STEP} ${BETA_END} $root $rootCLUSTER
# 
# 
# 
squeue -u $USER
