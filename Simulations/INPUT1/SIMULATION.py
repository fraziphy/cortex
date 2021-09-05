import sys
sys.path.append("../MODULE/")
import functions_mass_model as FMM
import os
import numpy as np
from scipy import stats
from mpi4py import MPI
# 
# 
# 
# 
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()
# 
# 
# 
# 
n_trial = 500
sub_n_trial = n_trial//p
# 
# 
# 
# 
vig_state = os.environ.get("VIG_STATE")
root_suff = os.environ.get("ROOT")
root_c = os.environ.get("ROOT_C")
root = "/" + root_c +"/scratch/" + root_suff
g_ampa = int(vig_state.split("_")[0])
beta = int(vig_state.split("_")[1])
# 
# 
# 
# 
# 
# 
# 
# 
FMM.n_pop = 2
if beta == 0:
    FMM.n_pop = 1
# 
# 
# 
# 
# 
# 
# 
RNG = None
RNG_init = None
g_gaba = None
if my_rank == 0:
    # 
    # 
    # 
    # 
    entropy_seed = 12345
    RNG = [np.random.default_rng(np.random.SeedSequence(entropy=entropy_seed,spawn_key=(FMM.n_pop-1,0,g_ampa,beta,0,k),)) for k in range(n_trial * 2 * FMM.n_pop)]
    RNG = [[RNG[i*FMM.n_pop*2*sub_n_trial + j*FMM.n_pop*2:i*FMM.n_pop*2*sub_n_trial + (j+1)*FMM.n_pop*2] for j in range(sub_n_trial)] for i in range(p)]
    # 
    # 
    RNG_init = [np.random.default_rng(np.random.SeedSequence(entropy=entropy_seed,spawn_key=(FMM.n_pop-1,1,g_ampa,beta,0,k),)) for k in range(p)]
    # 
    # 
    # 
    # 
    g_gaba = 1, 1
    if beta != 0 or g_ampa == 2:
        g_gaba = np.load(root+"G_GABA_{}.npy".format(vig_state))
    g_gaba = g_gaba[0],g_gaba[1]
RNG = comm.scatter(RNG, root=0)
RNG_init = comm.scatter(RNG_init, root=0)
# 
# 
# 
# 
g_gaba = comm.bcast(g_gaba, root=0)
# 
# 
# 
# 
FMM.G_AMPA = np.array([g_ampa,g_ampa]).reshape(2,-1)
FMM.G_AMPA_INPUT = np.array([g_ampa,g_ampa]).reshape(2,-1)
FMM.G_AMPA_NOISE = np.array([g_ampa,g_ampa]).reshape(2,-1)
FMM.G_GABA = np.array(g_gaba).reshape(2,-1)
FMM.beta = beta
# 
# 
# 
# 
phi_n_sd = 1.8
if g_ampa != 1:
    phi_n_sd = 1.
FMM.phi_n_sd = phi_n_sd
# 
# 
# 
# 
T = 10 # simulation time in seconds
dt = 0.1 # timestep in ms
n = int(T*1000/dt)
sens_amp = 1.
sens_duration = 100
sens_input_iter = int(sens_duration/dt)
n_1 = int(n/2)          #sensory implimentation time
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#     
input_no = np.zeros((2,FMM.n_pop))
input_1 = np.zeros((2,FMM.n_pop))
input_1[0,0] = 1
# 
# 
# 
# 
# building variables
data = np.zeros((sub_n_trial,23,FMM.n_pop,n),dtype=float)
data_initi = np.zeros((23, FMM.n_pop,2), dtype=float)  #second 2 is the number of populations
cortical_white_noise = np.zeros((2,FMM.n_pop),dtype=float)          #noise to cortical population
l1 = np.zeros((23,FMM.n_pop),dtype=float)           #stocastic term to be added for the SDE solution in the cortical network
# 
# 
# 
#  
# 
# 
# 
n_initial = int(100000 / dt)
# simulation letting system to reach steady state solution
for j in range(sub_n_trial):
    # initial value for cortical network
    data_initi[0,:,0] = -10*RNG_init.random(FMM.n_pop) + FMM.theta_p
    data_initi[1,:,0] = -10*RNG_init.random(FMM.n_pop) + FMM.theta_i
    data_initi[2:,:,0] = 0.01*RNG_init.random((21,FMM.n_pop))
    data_initi[11:19,:,0] = 0
    for i in range(n_initial - 1):
        # Runge–Kutta method 2nd order for cortical network
        cortical_white_noise[:] = RNG_init.normal(0,FMM.phi_n_sd,size=(2,FMM.n_pop))
        l1[[12,14]] = FMM.y_e * FMM.y_e * np.sqrt(dt) * cortical_white_noise
        data_initi[:,:, 1] = FMM.RK2order_Cor(dt,data_initi[:,:, 0],l1,input_no)
        # preparing for the next iteration
        data_initi[:, :,0] = data_initi[:,:, 1]
    # 
    # 
    # 
    # 
    data[j,:,:,0] = data_initi[:,:,0]
    # 
    #
    # 
    # 
    for i in range(n - 1):
        # implimenting sensory noise
        # 
        # 
        # 
        #    
        cortical_white_noise[:] = [[RNG[j][k+kk*FMM.n_pop].normal(0,FMM.phi_n_sd) for k in range(FMM.n_pop)]for kk in range(2)]
        l1[[12, 14]] = FMM.y_e * FMM.y_e * np.sqrt(dt) * cortical_white_noise
        input_yes = int(i >= n_1 and i <= n_1+sens_input_iter)
        data[j,:,:, i+1] = FMM.RK2order_Cor(dt,data[j,:,:, i],l1,input_yes*sens_amp*input_1)
# 
# 
# 
# 
recvbuf = None
if my_rank == 0:
    recvbuf =np.empty((n_trial,23,FMM.n_pop,n),dtype=float)
# 
# 
# 
# 
comm.Gather(data, recvbuf, root=0)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
if my_rank == 0:
    new_arang = np.transpose(recvbuf,(1,2,0,3))
    np.save(root+"Simulation_{}.npy".format(vig_state),new_arang)
    # 
    # 
    # 
    # 
    i_ampa_pyr = FMM.I_AMPA(FMM.G_AMPA[0]*new_arang[2] + FMM.G_AMPA_NOISE[0]*new_arang[11] + FMM.G_AMPA_INPUT[0]*new_arang[15] + FMM.beta*FMM.G_AMPA[0]*new_arang[19], new_arang[0])
    i_ampa_inh = FMM.I_AMPA(FMM.G_AMPA[1]*new_arang[6] + FMM.G_AMPA_NOISE[1]*new_arang[13] + FMM.G_AMPA_INPUT[1]*new_arang[17] + FMM.beta*FMM.G_AMPA[1]*new_arang[21], new_arang[1])
    i_gaba_pyr = FMM.I_GABA(FMM.G_GABA[0]*new_arang[4], new_arang[0])
    i_gaba_inh = FMM.I_GABA(FMM.G_GABA[1]*new_arang[8], new_arang[1])
    # 
    LFP_pyr = np.abs(i_ampa_pyr) + np.abs(i_gaba_pyr)
    # 
    LFP_inh = np.abs(i_ampa_inh) + np.abs(i_gaba_inh)
    # 
    # 
    # 
    # 
    Firing_pyr = FMM.Qp(new_arang[0])
    # 
    Firing_inh = FMM.Qi(new_arang[1])
    # 
    # 
    # 
    # 
    np.save(root+"Perturbed_Population_lfp_{}_pyr.npy".format(vig_state),LFP_pyr[0]-LFP_pyr[0][:,:n_1].mean(axis=1).reshape(n_trial,-1))
    np.save(root+"Perturbed_Population_v_{}_pyr.npy".format(vig_state),new_arang[0,0])
    np.save(root+"Perturbed_Population_firing_{}_pyr.npy".format(vig_state),Firing_pyr[0])
    np.save(root+"Perturbed_Population_lfp_{}_inh.npy".format(vig_state),LFP_inh[0]-LFP_inh[0][:,:n_1].mean(axis=1).reshape(n_trial,-1))
    np.save(root+"Perturbed_Population_v_{}_inh.npy".format(vig_state),new_arang[1,0])
    np.save(root+"Perturbed_Population_firing_{}_inh.npy".format(vig_state),Firing_inh[0])
    if beta != 0:
        np.save(root+"Unperturbed_Population_lfp_{}_pyr.npy".format(vig_state),LFP_pyr[1]-LFP_pyr[1][:,:n_1].mean(axis=1).reshape(n_trial,-1))
        np.save(root+"Unperturbed_Population_v_{}_pyr.npy".format(vig_state),new_arang[0,1])
        np.save(root+"Unperturbed_Population_firing_{}_pyr.npy".format(vig_state),Firing_pyr[1])
        np.save(root+"Unperturbed_Population_lfp_{}_inh.npy".format(vig_state),LFP_inh[1]-LFP_inh[1][:,:n_1].mean(axis=1).reshape(n_trial,-1))
        np.save(root+"Unperturbed_Population_v_{}_inh.npy".format(vig_state),new_arang[1,1])
        np.save(root+"Unperturbed_Population_firing_{}_inh.npy".format(vig_state),Firing_inh[1])       
# 
# 
# 
# 
# 
# 
# 
MPI.Finalize
