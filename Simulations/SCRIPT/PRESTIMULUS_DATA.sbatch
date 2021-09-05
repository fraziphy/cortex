#!/bin/bash
##############################################################
#                         Data Curation                      #
#                                                            #
#                   -Prestimulus intervals                   #
##############################################################
#SBATCH --time=0-00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=farhad.razi01@estudiant.upf.edu
#SBATCH --output=output_prestimulus_data_%x.out
#SBATCH --error=error_prestimulus_data_%x.err
# 
# 
# 
# 
ml matplotlib/2.1.2-foss-2018a-Python-3.6.4
# 
# 
# 
# 
export COND=$1 G_AMPA_STEP=$2 G_AMPA_END=$3 BETA_STEP=$4 BETA_END=$5 ROOT=$6 ROOT_C=$7
python <<-EOF
import numpy as np
import pickle
import os
import time
# 
# 
# 
# 
cond = os.environ.get("COND")
g_ampa_step = int(os.environ.get("G_AMPA_STEP"))
g_ampa_end = int(os.environ.get("G_AMPA_END"))
beta_step = int(os.environ.get("BETA_STEP"))
beta_end = int(os.environ.get("BETA_END"))
root_suff = os.environ.get("ROOT")
root_c = os.environ.get("ROOT_C")
root = "/" + root_c +"/scratch/" + root_suff
# 
# 
# 
# 
# 
# 
G_AMPA = [2]
BETA = [0]
pops=["Perturbed"]
beta_NREM = 0
if cond == "TWO_CORTICAL_COLUMN":
    G_AMPA = [i for i in range(2,g_ampa_end+1,g_ampa_step)]
    BETA = [i for i in range(1,beta_end+1,beta_step)]
    pops=["Perturbed","Unperturbed"]
    beta_NREM = 1
# 
# 
signals=["lfp","firing","v"]
sub_population = ["pyr","inh"]
# 
# 
# 
# 
# 
# 
# 
# 
DATA_PLOT = {}
DATA_PLOT["BETA"] = BETA
DATA_PLOT["G_AMPA"] = G_AMPA
# 
# 
# 
# 
# 
# 
for pop in pops:
    DATA_PLOT[pop] = {}
    for sub_pop in sub_population:
        DATA_PLOT[pop][sub_pop] = {}
        for gampa in G_AMPA:
            for beta in BETA:
                vig = "{}_{}".format(gampa,beta)
                DATA_PLOT[pop][sub_pop][vig] = {}
                a = None
                with open(root+'PRESTIMULUS_PREPARATION_DATA_{}_{}_{}.pickle'.format(vig,pop,sub_pop), 'rb') as handle:
                    a = pickle.load(handle)
                for signal in signals:
                    DATA_PLOT[pop][sub_pop][vig][signal] = {}
                    DATA_PLOT[pop][sub_pop][vig][signal]["SAMPEL"] = a[signal]["SAMPEL"]
                    if signal == "lfp":
                        DATA_PLOT[pop][sub_pop][vig][signal]["BINS"] = a[signal]["BINS"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["HIST"] = a[signal]["HIST"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["FREQUENCY"] = a[signal]["FREQUENCY"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["FFT_M"] = a[signal]["FFT_M"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["FFT_SD"] = a[signal]["FFT_SD"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["FFT_RATIO_BINS"] = a[signal]["FFT_RATIO_BINS"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["FFT_RATIO"] = a[signal]["FFT_RATIO"]
                    else:
                        DATA_PLOT[pop][sub_pop][vig][signal]["BINS"] = a[signal]["BINS"]
                        DATA_PLOT[pop][sub_pop][vig][signal]["HIST"] = a[signal]["HIST"]
        vig = "1_{}".format(beta_NREM)
        DATA_PLOT[pop][sub_pop][vig] = {}
        a = None
        with open(root+'PRESTIMULUS_PREPARATION_DATA_{}_{}_{}.pickle'.format(vig,pop,sub_pop), 'rb') as handle:
            a = pickle.load(handle)
        for signal in signals:
            DATA_PLOT[pop][sub_pop][vig][signal] = {}
            DATA_PLOT[pop][sub_pop][vig][signal]["SAMPEL"] = a[signal]["SAMPEL"]
            if signal == "lfp":
                DATA_PLOT[pop][sub_pop][vig][signal]["BINS"] = a[signal]["BINS"]
                DATA_PLOT[pop][sub_pop][vig][signal]["HIST"] = a[signal]["HIST"]
                DATA_PLOT[pop][sub_pop][vig][signal]["FREQUENCY"] = a[signal]["FREQUENCY"]
                DATA_PLOT[pop][sub_pop][vig][signal]["FFT_M"] = a[signal]["FFT_M"]
                DATA_PLOT[pop][sub_pop][vig][signal]["FFT_SD"] = a[signal]["FFT_SD"]
                DATA_PLOT[pop][sub_pop][vig][signal]["FFT_RATIO_BINS"] = a[signal]["FFT_RATIO_BINS"]
                DATA_PLOT[pop][sub_pop][vig][signal]["FFT_RATIO"] = a[signal]["FFT_RATIO"]
            else:
                DATA_PLOT[pop][sub_pop][vig][signal]["UP_BINS"] = a[signal]["UP_BINS"]
                DATA_PLOT[pop][sub_pop][vig][signal]["UP_HIST"] = a[signal]["UP_HIST"]
                DATA_PLOT[pop][sub_pop][vig][signal]["DOWN_BINS"] = a[signal]["DOWN_BINS"]
                DATA_PLOT[pop][sub_pop][vig][signal]["DOWN_HIST"] = a[signal]["DOWN_HIST"]
# 
# 
# 
# 
# 
# 
with open(root+'PRESTIMULUS_DATA_{}.pickle'.format(cond), 'wb') as handle:
    pickle.dump(DATA_PLOT, handle, protocol=pickle.HIGHEST_PROTOCOL)
# 
# 
# 
# 
# 
# 
for pop in pops:
    for sub_pop in sub_population:
        vig = "1_{}".format(beta_NREM)
        os.remove(root+'PRESTIMULUS_PREPARATION_DATA_{}_{}_{}.pickle'.format(vig,pop,sub_pop))
        time.sleep(5)
        for gampa in G_AMPA:
            for beta in BETA:
                vig = "{}_{}".format(gampa,beta)
                os.remove(root+'PRESTIMULUS_PREPARATION_DATA_{}_{}_{}.pickle'.format(vig,pop,sub_pop))
                time.sleep(5)
# 
# 
# 
# 
# 
# 
EOF
