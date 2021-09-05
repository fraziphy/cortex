#!/bin/bash
##############################################################
#         Non parametric data_test within populations        #
#                                                            #
##############################################################
#SBATCH --time=0-02:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3500M
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=farhad.razi01@estudiant.upf.edu
#SBATCH --output=output_non_parametric_within_pop_%x.out
#SBATCH --error=error_non_parametric_within_pop_%x.err
# 
# 
# 
# 
ml Python
# 
# 
# 
# 
export VIG_1=$1 POP=$2 SUB_POP=$3 SIGNAL=$4 ROOT=$5 ROOT_C=$6
python <<-EOF
import sys
sys.path.append("../MODULE/")
import functions_non_parametric_test as FNPT
import numpy as np
from scipy import stats
import pickle
import os
# 
# 
# 
# 
vig_state = os.environ.get("VIG_1")
population = os.environ.get("POP")
sub_population = os.environ.get("SUB_POP")
signal = os.environ.get("SIGNAL")
root_suff = os.environ.get("ROOT")
root_c = os.environ.get("ROOT_C")
root = "/" + root_c +"/scratch/" + root_suff
g_ampa = int(vig_state.split("_")[0])
beta = int(vig_state.split("_")[1])
# 
# 
# 
n_pop = 2
if beta == 0:
    n_pop = 1
# 
# 
key_list = ["Perturbed_Population_lfp_pyr" ,"Unperturbed_Population_lfp_pyr", "Perturbed_Population_lfp_inh", "Unperturbed_Population_lfp_inh", "Perturbed_Population_firing_pyr", "Unperturbed_Population_firing_pyr", "Perturbed_Population_firing_inh", "Unperturbed_Population_firing_inh","Perturbed_Population_v_pyr" ,"Unperturbed_Population_v_pyr", "Perturbed_Population_v_inh", "Unperturbed_Population_v_inh"]
index_data = key_list.index("{}_Population_{}_{}".format(population,signal,sub_population))
# 
# 
# 
data = np.load(root+"{}_Population_{}_{}_{}.npy".format(population,signal,vig_state,sub_population))
n = data.shape[1]//2
n_trial = data.shape[0]
baseline_data = (data[:,:n].reshape(-1,n//2)).copy()
data_test = (np.vstack((data[:,n:],data[:,:n]))).copy()
# 
# 
# 
# 
# 
# 
entropy_seed = 12345
RNG_init = np.random.default_rng(np.random.SeedSequence(entropy=entropy_seed,spawn_key=(n_pop-1,0,g_ampa,beta,1,index_data),))
RNG = np.random.default_rng(np.random.SeedSequence(entropy=entropy_seed,spawn_key=(n_pop-1,1,g_ampa,beta,1,index_data),))
# 
# 
# 
# 
df = n_trial * 2 - 2
t_crit = stats.t.ppf(1-FNPT.significance_level/2,df)
t_stat = stats.ttest_ind(data_test[:n_trial],data_test[n_trial:])[0]
t_cluster = FNPT.T_CLUSTER(t_stat,t_crit)
# 
# 
# 
# 
bon_t_crit = stats.t.ppf(1-FNPT.significance_level/data_test.shape[1]/2,df)
bon_t_cluster= FNPT.T_CLUSTER(t_stat,bon_t_crit)
# 
# 
# 
# 
a = {}
# a["Difference"] = data_test[:n_trial].mean(axis=0)-data_test[n_trial:].mean(axis=0)
# n_initial = 0
a["Difference"] = data.mean(axis=0)
n_initial = data.shape[1]//2
a["t_stat"] = t_stat
a["Degree of freedom"] = df
a["T_Critical"] = t_crit
a["Initial_T_Cluster_AND_Index_Info"] = t_cluster
# 
# 
# 
# 
a["T_Critical_BON"] = bon_t_crit
a["Initial_T_Cluster_AND_Index_Info_BON"] = bon_t_cluster
# 
# 
# 
# 
# 
#
if t_cluster is not None:
    t_cluster_NOT = FNPT.NON_SIGNIFICANT_CLUSTER(baseline_data,t_crit,n_trial,RNG_init)
    index = t_cluster[1][t_cluster[0]>t_cluster_NOT]
    t_cluster = t_cluster[0][t_cluster[0]>t_cluster_NOT]
    t_cluster_ordered = (np.sort(t_cluster)[::-1]).reshape(1,-1)
    index_sort = np.argsort(t_cluster)[::-1]
    t_n = len(t_cluster)
    if t_n > 0:
        t_distribution = FNPT.CLUSTER_DISTRIBUTION(data_test,t_crit,n_trial,t_cluster_ordered,t_n,RNG)
        N_t_dist = (~np.isnan(t_distribution)).sum(axis=0)
        t_distribution[np.isnan(t_distribution)] = -1.
        critical_level = (t_distribution>=t_cluster).sum(axis=0)
        p_value = np.round(critical_level/N_t_dist,4)
        p_value = p_value[index_sort]
        # 
        # 
        # 
        # 
        a["T_Cluster_Info"] = t_cluster,index + n_initial
        a["P_Value"] = p_value
        a["Distribution"] = t_distribution.T[index_sort].T
        a["Significant_Index"] = index[p_value<=FNPT.significance_level] + n_initial
        a["Non_Significant_Cluster"] = t_cluster_NOT
        a["Significant_Info"] = "There is significanct cluster"
    else:
        a["T_Cluster_Info"] = None
        a["P_Value"] = None
        a["Distribution"] = None
        a["Significant_Index"] = np.nan
        a["Non_Significant_Cluster"] = t_cluster_NOT
        a["Significant_Info"] = "There is cluster but NOT significant"
else:
    a["T_Cluster_Info"] = None
    a["P_Value"] = None
    a["Distribution"] = None
    a["Significant_Index"] = np.nan
    a["Non_Significant_Cluster"] = None
    a["Significant_Info"] = "There is NOT cluster"

a["General_info"] = "Within population",population,vig_state, signal, sub_population
# 
# 
# 
# 
with open(root+'Non_Parametric_Test_Within_{}_Population_{}_{}_{}.pickle'.format(population,signal,vig_state,sub_population), 'wb') as handle:
    pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
EOF
