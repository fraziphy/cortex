#!/bin/bash
##############################################################
#                       G_GABA calculation                   #
#                                                            #
##############################################################
#SBATCH --time=0-00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=farhad.razi@upf.edu
#SBATCH --output=output_new_g_gaba.out
#SBATCH --error=error_new_g_gaba.err
# 
# 
# 
# 
ml Python
# 
# 
# 
# 
export G_AMPA_STEP=$1 G_AMPA_END=$2 BETA_STEP=$3 BETA_END=$4 ROOT=$5 ROOT_C=$6
python <<-EOF
import sys
sys.path.append("../MODULE/")
import functions_mass_model as FMM
import numpy as np
from scipy.signal import argrelextrema
import os
import time
# 
# 
# 
# 
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
T = 10
dt = 0.1
n_total = int(T*1000/dt)
n = int(n_total/2)
dur = 100
dur = int(dur/dt)
# 
# 
# 
# 
n_realizations = 500
# 
# 
# 
# 
def SMOTH(y,subpop):
    z=np.linspace(-60,60,121)
    zz=np.exp(-1/500*z**2)
    if subpop == "inh":
        zz=np.exp(-1/100*z**2)
    zz /= zz.sum()
    zzz=np.convolve(y,zz,mode="same")
    return zzz
# 
# 
# 
# 
def SMOTH_V(y):
    z=np.linspace(-10,10,21)
    zz=np.exp(-1/20*z**2)
    zz /= zz.sum()
    zzz=np.convolve(y,zz,mode="same")
    return zzz
# 
# 
# 
# 
def FIND_PEAK_V(hist):
    z = SMOTH_V(hist)
    return np.argmax(z)
# 
# 
# 
# 
def FIND_PEAK(hist,subpop):
    z = SMOTH(hist,subpop)
    return argrelextrema(z, np.greater)[0]
# 
# 
# 
# 
def FIND_G_GABA(g_ampa,beta,vp,vi):
    qp,qi = FMM.Qp(vp),FMM.Qi(vi)

    A = FMM.alpha_Na*qp/FMM.R_pump+(FMM.Na_eq**3)/(FMM.Na_eq**3+3375)
    NA = np.cbrt(A*3375 / (1-A))
    W = 0.37 / (1+(38.7/NA)**3.5)
    
    g_gaba_exc_wake = -(FMM.g_l*(vp - FMM.E_p_l)+ g_ampa* FMM.N_pp*qp*(vp-FMM.E_ampa)+ g_ampa * beta* FMM.N_pP*qp*(vp-FMM.E_ampa)+ FMM.tau_p/FMM.c_m*FMM.g_k_Na*W*(vp-FMM.E_k)) / (FMM.N_pi*qi*(vp-FMM.E_gaba) )
    g_gaba_inh_wake = -(FMM.g_l*(vi - FMM.E_i_l)+ g_ampa* FMM.N_ip*qp*(vi-FMM.E_ampa) + g_ampa * beta* FMM.N_iP*qp*(vi-FMM.E_ampa)) / ( FMM.N_ii*qi*(vi-FMM.E_gaba))
    
    return g_gaba_exc_wake,g_gaba_inh_wake
# 
# 
# 
# 
data_vp = np.load(root+"Perturbed_Population_lfp_1_0_pyr.npy")
data_vi = np.load(root+"Perturbed_Population_lfp_1_0_inh.npy")
v_p = np.load(root+"Perturbed_Population_v_1_0_pyr.npy")
v_i = np.load(root+"Perturbed_Population_v_1_0_inh.npy")
# 
# 
# 
# 
bins = np.arange(data_vp[:,:n].min(),data_vp[:,:n].max()+0.5,0.5)
HIST = np.histogram(data_vp[:,:n].ravel(),bins=bins)[0]
ind1,ind2  = FIND_PEAK(HIST,"pyr")
threshold = (bins[ind1]+bins[ind2])/2
mask_v_p = data_vp[:,:n]>=threshold
v_p_up = v_p[:,:n][mask_v_p]
bins_1 = np.arange(v_p_up.min(),v_p_up.max()+0.5,0.5)
v_p_up_hist = np.histogram(v_p_up,bins=bins_1)[0]
ind1  = FIND_PEAK_V(v_p_up_hist)
v_p_up_mean = bins_1[ind1]
# 
# 
# 
# 
bins = np.arange(data_vi[:,:n].min(),data_vi[:,:n].max()+0.5,0.5)
HIST = np.histogram(data_vi[:,:n].ravel(),bins=bins)[0]
ind1,ind2  = FIND_PEAK(HIST,"inh")
threshold = (bins[ind1]+bins[ind2])/2
mask_v_i = data_vi[:,:n]>=threshold
v_i_up = v_i[:,:n][mask_v_i]
bins_1 = np.arange(v_i_up.min(),v_i_up.max()+0.5,0.5)
v_i_up_hist = np.histogram(v_i_up,bins=bins_1)[0]
ind1  = FIND_PEAK_V(v_i_up_hist)
v_i_up_mean = bins_1[ind1]
# 
# 
# 
# 
new_g_gaba = np.round(FIND_G_GABA(2,0,v_p_up_mean,v_i_up_mean),3)
np.save(root+"G_GABA_2_0.npy",new_g_gaba)
# 
# 
# 
# 
for i in range(2,g_ampa_end+1,g_ampa_step):
    for j in range(1,beta_end+1,beta_step):
        new_g_gaba = np.round(FIND_G_GABA(i,j,v_p_up_mean,v_i_up_mean),3)
        np.save(root+"G_GABA_{}_{}.npy".format(i,j),new_g_gaba)
        time.sleep(5)
# 
# 
# 
# 
# 
# 
# 
# 
FMM.n_pop = 1
# 
# 
# 
# 
FMM.G_AMPA = np.array([1,1]).reshape(2,-1)
FMM.G_AMPA_INPUT = np.array([0,0]).reshape(2,-1)
FMM.G_AMPA_NOISE = np.array([0,0]).reshape(2,-1)
FMM.G_GABA = np.array([1,1]).reshape(2,-1)
FMM.beta = 0
# 
# 
# 
# 
T = 40 # simulation time in seconds
dt = 0.1 # timestep in ms
n = int(T*1000/dt)
# 
# 
# 
# 
#     
input_no = np.zeros((2,FMM.n_pop))
# 
# 
# 
# 
# building variables
data_cor = np.zeros((23, FMM.n_pop,2), dtype=float)  #second 2 is the number of populations
l1 = np.zeros((23,FMM.n_pop),dtype=float)           #stocastic term to be added for the SDE solution in the cortical network
# 
# 
# 
# initial value for cortical network
entropy_seed = 12345
RNG_init = np.random.default_rng(np.random.SeedSequence(entropy=entropy_seed,spawn_key=(0,),))
data_cor[0,:,0] = -10*RNG_init.random(FMM.n_pop) + FMM.theta_p
data_cor[1,:,0] = -10*RNG_init.random(FMM.n_pop) + FMM.theta_i
data_cor[2:,:,0] = 0.01*RNG_init.random((21,FMM.n_pop))
data_cor[[15,16,17,18],:,0] = 0
#  
# 
# 
# 
# simulation letting system to reach steady state solution
for i in range(n - 1):
    # implimenting sensory noise
    # 
    # 
    # 
    # 
    data_cor[:,:, 1] = FMM.RK2order_Cor(dt,data_cor[:,:, 0],l1,input_no)
    # preparing for the next iteration
    data_cor[:,:, 0] = data_cor[:,:, 1]
#  
# 
# 
# 
new_g_gaba = np.round(FIND_G_GABA(1,1,data_cor[0,0,-1],data_cor[1,0,-1]),3)
np.save(root+"G_GABA_1_1.npy",new_g_gaba)
#  
# 
# 
# 
EOF
