import numpy as np
from scipy import stats
# 
# 
# 
# 
#
# 
#
permut = 1000
# 
# 
# 
# 
significance_level = 0.01
# 
# 
# 
# 
# 
def FIND_IND(tstat,tcri):
    f_ind = next((i for i,x in enumerate(tstat) if x>tcri),None)
    all_lis = []
    while f_ind is not None:
        all_lis.append(f_ind)
        l_ind = next((i for i,x in enumerate(tstat) if x<=tcri and i>=f_ind),None)
        if l_ind is not None:
            all_lis.append(l_ind)
            f_ind = next((i for i,x in enumerate(tstat) if x>tcri and i>l_ind),None)
        else:
            all_lis.append(len(tstat)-1)
            f_ind = None
    return all_lis
# 
# 
# 
# 
# 
def FIND_PERIOUD(tstat,tcri):
    p = FIND_IND(tstat,tcri)
    n = FIND_IND(-tstat,tcri)
    all_ind = p+n
    if not all_ind:
        return None
    all_ind.sort()
    return np.array(all_ind).reshape(-1,2)
# 
# 
# 
# 
# 
def T_CLUSTER(tstat,tcri):
    index = FIND_PERIOUD(tstat,tcri)
    if index is None:
        return None
    t_clusters = np.zeros(index.shape[0],dtype=float)
    for i in range(index.shape[0]):
        t_clusters[i] = np.abs(tstat[index[i,0]:index[i,1]+1].sum())
    return t_clusters,index
# 
# 
# 
# 
def NON_SIGNIFICANT_CLUSTER(data,t_crit,n_trial,rng_init):
    data_aux = data.copy()
    t_cluster_init = np.zeros(permut,dtype=float)
    k=0
    while k<permut:
        t_stat_init = stats.ttest_ind(data_aux[:n_trial],data_aux[n_trial:])[0]
        t_cluster = T_CLUSTER(t_stat_init,t_crit)
        if t_cluster is None:
            rng_init.shuffle(data_aux)
            continue
        else:
            t_cluster_init[k] = t_cluster[0].max()
            rng_init.shuffle(data_aux)
            k += 1
    return t_cluster_init.max()
# 
# 
# 
# 
# 
def CLUSTER_DISTRIBUTION(data,t_crit,n_trial,t_cluster_ordered,t_n,rng):
    data_aux = data.copy()
    t_distribution = t_cluster_ordered.copy()
    collect = np.full(t_n,np.nan)
    while (~np.isnan(t_distribution[:,-1])).sum() < permut:
        collect[:] = np.nan
        rng.shuffle(data_aux)
        t_stat_new = stats.ttest_ind(data_aux[:n_trial],data_aux[n_trial:])[0]
        t_cluster_aux = T_CLUSTER(t_stat_new,t_crit)
        if t_cluster_aux is not None:
            t_cluster_aux = np.sort(t_cluster_aux[0])[::-1][:t_n]
            collect[:len(t_cluster_aux)] = t_cluster_aux
            t_distribution = np.vstack((t_distribution,collect))
        else:
            continue
    return t_distribution
# 
# 
# 
# 
# 
