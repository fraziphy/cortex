#!/bin/bash
##############################################################
#                         Data Curation                      #
#                                                            #
#                   -Poststimulus intervals                  #
##############################################################
#SBATCH --time=0-00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=farhad.razi01@estudiant.upf.edu
#SBATCH --output=output_poststimulus_data_%x.out
#SBATCH --error=error_poststimulus_data_%x.err
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
DATA_PLOT["WITHIN_POP"] = {}
for pop in pops:
    DATA_PLOT["WITHIN_POP"][pop] = {}
    for sub_pop in sub_population:
        DATA_PLOT["WITHIN_POP"][pop][sub_pop] = {}
        for signal in signals:
            DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal] = {}
            vig = "1_{}".format(beta_NREM)
            DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig] = {}
            data_cluster_test = None
            with open(root+'Non_Parametric_Test_Within_{}_Population_{}_{}_{}.pickle'.format(pop,signal,vig,sub_pop), 'rb') as handle:
                data_cluster_test = pickle.load(handle)
            if signal == "firing":
                DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Difference"] = data_cluster_test["Difference"] * 1000
            else:
                DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Difference"] = data_cluster_test["Difference"]
            DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Significant_Index"] = data_cluster_test["Significant_Index"]
            DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["P_Value"] = data_cluster_test["P_Value"]
            DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Initial_T_Cluster_AND_Index_Info"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info"]
            DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Initial_T_Cluster_AND_Index_Info_BON"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info_BON"]
            for gampa in G_AMPA:
                for i,beta in enumerate(BETA):
                    vig = "{}_{}".format(gampa,beta)
                    DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig] = {}
                    data_cluster_test = None
                    with open(root+'Non_Parametric_Test_Within_{}_Population_{}_{}_{}.pickle'.format(pop,signal,vig,sub_pop), 'rb') as handle:
                        data_cluster_test = pickle.load(handle)
                    if signal == "firing":
                        DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Difference"] = data_cluster_test["Difference"] * 1000
                    else:
                        DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Difference"] = data_cluster_test["Difference"]
                    DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Significant_Index"] = data_cluster_test["Significant_Index"]
                    DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["P_Value"] = data_cluster_test["P_Value"]
                    DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Initial_T_Cluster_AND_Index_Info"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info"]
                    DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]["Initial_T_Cluster_AND_Index_Info_BON"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info_BON"]
# 
# 
# 
# 
if cond == "ONE_CORTICAL_COLUMN":
    DATA_PLOT["ACROSS_STATE"] = {}
    for pop in pops:
        DATA_PLOT["ACROSS_STATE"][pop] = {}
        for sub_pop in sub_population:
            DATA_PLOT["ACROSS_STATE"][pop][sub_pop] = {}
            for signal in signals:
                DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal] = {}
                for gampa in G_AMPA:
                    for i,beta in enumerate(BETA):
                        vig = "{}_{}".format(gampa,beta)
                        DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig] = {}
                        data_cluster_test = None
                        with open(root+'Non_Parametric_Test_Across_State_{}_{}_{}_{}.pickle'.format(pop,signal,vig,sub_pop), 'rb') as handle:
                            data_cluster_test = pickle.load(handle)
                        if signal == "firing":
                            DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig]["Difference"] = data_cluster_test["Difference"] * 1000
                        else:
                            DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig]["Difference"] = data_cluster_test["Difference"]
                        DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig]["Significant_Index"] = data_cluster_test["Significant_Index"]
                        DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig]["P_Value"] = data_cluster_test["P_Value"]
                        DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig]["Initial_T_Cluster_AND_Index_Info"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info"]
                        DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig]["Initial_T_Cluster_AND_Index_Info_BON"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info_BON"]
# 
# 
# 
# 
if cond == "TWO_CORTICAL_COLUMN":
    # 
    # 
    # 
    #   
    # 
    # 
    # 
    # 
    # 
    pop,signal = "Unperturbed", "firing"
    vig = "2_5"
    DATA_PLOT["EXAMPLE_T_CLUSTER"] = {}
    DATA_PLOT["EXAMPLE_T_CLUSTER"]["info"] = "{} -- {} -- {}".format(pop,signal,vig)
    for sub_pop in sub_population:
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop] = {}
        data_cluster_test = None
        with open(root+'Non_Parametric_Test_Within_{}_Population_{}_{}_{}.pickle'.format(pop,signal,vig,sub_pop), 'rb') as handle:
            data_cluster_test = pickle.load(handle)
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["t_stat"] = data_cluster_test["t_stat"]
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["T_Critical"] = data_cluster_test["T_Critical"]
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["Initial_T_Cluster_AND_Index_Info"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info"]
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["P_Value"] = data_cluster_test["P_Value"]
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["T_Critical_BON"] = data_cluster_test["T_Critical_BON"]
        DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["Initial_T_Cluster_AND_Index_Info_BON"] = data_cluster_test["Initial_T_Cluster_AND_Index_Info_BON"]
# 
# 
# 
with open(root+'POSTSTIMULUS_DATA_{}.pickle'.format(cond), 'wb') as handle:
    pickle.dump(DATA_PLOT, handle, protocol=pickle.HIGHEST_PROTOCOL)
EOF
