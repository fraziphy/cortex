import os
import numpy as np
import pickle
import sympy
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time
# 
# 
# 
# 
root_c = os.environ.get("rootMACHINE")
root = root_c +"/DATA/" 
root_fig = root_c+"/" 
# 
# 
# 
# 
def LAMDA(g_ampa_pyr,g_ampa_inh,g_gaba_pyr,g_gaba_inh,vp_steady,vi_steady):


    var = sympy.symbols('v_p v_i sep sep_dot sgp sgp_dot sei sei_dot sgi sgi_dot na')
    qp = 0.5*q_max_p *(1 + sympy.tanh((var[0]-theta_p)*0.5 * np.pi/sigma_p/np.sqrt(3)))
    qi = 0.5*q_max_i *(1 + sympy.tanh((var[1]-theta_i)*0.5 * np.pi/sigma_i/np.sqrt(3)))

    ipl = g_l*(var[0]-E_p_l)
    iil = g_l*(var[1]-E_i_l)

    iep = g_ampa_pyr*var[2]*(var[0]-E_ampa)
    igp = g_gaba_pyr*var[4]*(var[0]-E_gaba)

    iei = g_ampa_inh*var[6]*(var[1]-E_ampa)
    igi = g_gaba_inh*var[8]*(var[1]-E_gaba)

    wkna = 0.37/(1+(38.7/var[10])**3.5)
    ikna = g_k_Na*wkna*(var[0]-E_k)

    Na_pump = R_pump * (var[10]**3 / (var[10]**3 + 3375)) - Na_pump_0

    f0 = 1/tau_p*(-ipl-iep-igp-tau_p/c_m*ikna)
    f1 = 1/tau_i*(-iil-iei-igi)
    f2 = var[3]
    f3 = y_e**2 * (N_pp*qp - var[2]) - 2*y_e*var[3]
    f4 = var[5]
    f5 = y_g**2 * (N_pi*qi - var[4]) - 2*y_g*var[5]
    f6 = var[7]
    f7 = y_e**2 * (N_pp*qp - var[6]) - 2*y_e*var[7]
    f8 = var[9]
    f9 = y_g**2 * (N_pi*qi - var[8]) - 2*y_g*var[9]
    f10 = 1/tau_Na*(alpha_Na*qp - Na_pump)


    f = sympy.sympify([f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10])
    J = sympy.zeros(len(f),len(var))
    for i, fi in enumerate(f):
        for j, s in enumerate(var):
            J[i,j] = sympy.diff(fi, s)


    vp,vi = vp_steady,vi_steady

    qp = 0.5*q_max_p *(1 + np.tanh((vp-theta_p)*0.5 * np.pi/sigma_p/np.sqrt(3)))
    qi = 0.5*q_max_i *(1 + np.tanh((vi-theta_i)*0.5 * np.pi/sigma_i/np.sqrt(3)))

    sep = N_pp*qp
    sep_dot = 0
    sgp = N_pi*qi
    sgp_dot = 0
    sei = N_ip*qp
    sei_dot = 0
    sgi = N_ii*qi
    sgi_dot = 0

    na_coff = alpha_Na/R_pump*0.5*q_max_p *(1 + np.tanh((vp-theta_p)*0.5 * np.pi/sigma_p/np.sqrt(3)))+(Na_eq**3/(Na_eq**3+3375))
    na = np.cbrt(na_coff*3375/(1-na_coff))


    jacob = J.evalf(subs={var[0]: vp,
                          var[1]: vi,
                          var[2]: sep,
                          var[3]: sep_dot,
                          var[4]: sgp,
                          var[5]: sgp_dot,
                          var[6]: sei,
                          var[7]: sei_dot,
                          var[8]: sgi,
                          var[9]: sgi_dot,
                          var[10]: na})
    jacob = np.array(jacob).astype(np.float64)
    lamda=linalg.eigvals(jacob)
    
    return lamda
# 
# 
# 
# 
def I_AMPA(gs,v):
    return gs * (v - E_ampa)
# 
# 
# 
# 
def I_GABA(gs,v):
    return gs * (v - E_gaba)
# 
# 
# 
# 
def Qp(v):
    return 0.5*q_max_p *(1 + np.tanh((v-theta_p)*i_sigma_p))
# 
# 
# 
# 
def Qi(v):
    return 0.5*q_max_i *(1 + np.tanh((v-theta_i)*i_sigma_i))
# 
# 
# 
# 
# 
def Cortex_Field(yi,input_1):
    # array for storing the field value in "Cortical_Field" function for cortical network
    y = np.empty(yi.shape,dtype=float)
    q_p,q_i = Qp(yi[0]),Qi(yi[1])
    # 
    # 
    # 
    # 
    na_aux = yi[10] * yi[10] * yi[10]
    # 
    # 
    # 
    #  
    i_ampa = I_AMPA(G_AMPA*yi[[2,6]] + G_AMPA_NOISE*yi[[11,13]] + G_AMPA_INPUT*yi[[15,17]] + beta*G_AMPA*yi[[19,21]], yi[:2])
    i_gaba = I_GABA(G_GABA*yi[[4,8]], yi[:2])
    # 
    # 
    # 
    #  
    y[0] = (- g_l * (yi[0] - E_p_l) - i_ampa[0] - i_gaba[0]) * i_tau_p - g_k_Na * 0.37 / (1 + (38.7 / yi[10]) ** 3.5) * (yi[0] - E_k) * i_c_m
    y[1] = (- g_l * (yi[1] - E_i_l) - i_ampa[1] - i_gaba[1]) * i_tau_i
    y[2:9:2] = yi[3:10:2]
    y[[3,7]] = y_e * (y_e * (INTRA_CONN_EXC * q_p - yi[[2,6]]) - 2 * yi[[3,7]])
    y[[5,9]] = y_g * (y_g * (INTRA_CONN_INH * q_i - yi[[4,8]]) - 2 * yi[[5,9]])
    y[10] = (alpha_Na * q_p - (R_pump * (na_aux / (na_aux + 3375)) - Na_pump_0)) * i_tau_Na
    y[[11,13]] = yi[[12,14]]
    y[[12,14]] = y_e * (y_e * (- yi[[11,13]]) - 2 * yi[[12,14]])
    y[[15,17]] = yi[[16,18]]
    y[[16,18]] = y_e * (y_e * (input_1 - yi[[15,17]]) - 2 * yi[[16,18]])
    y[[19,21]] = yi[[20,22]]
    y[[20,22]] = y_e * (y_e * (INTER_CONN_EXC * q_p[::-1] - yi[[19,21]]) - 2 * yi[[20,22]])
    return y
# 
# 
# 
# 
def RK2order_Cor(dt,data_cor,l,input_1):
    k1_cor = dt * Cortex_Field(data_cor,input_1)
    k2_cor = dt * Cortex_Field(data_cor + k1_cor + l,input_1)
    return data_cor + 0.5 * (k1_cor + k2_cor) + l
# 
#
# 
# 
# 
#
# 
# 
n_pop = 1
# 
# 
# 
# 
# 
sigma_p = 6.7
g_k_Na = 1.9
# 
# 
# 
# 
# Parameter Space of Cortical Network
# 
# 
# 
# 
c_m = 1.
tau_p,tau_i = 30.,30.
q_max_p,q_max_i = 30e-3,60e-3
theta_p,theta_i = -58.5,-58.5
sigma_i = 6.
y_e,y_g = 70e-3,58.6e-3
g_l = 1.
E_p_l,E_i_l = -66.,-64.
E_k = -100.
E_ampa,E_gaba = 0., -70.
alpha_Na = 2.
tau_Na = 1.7
R_pump = 0.09
Na_eq = 9.5
Na_pump_0 = R_pump * (Na_eq**3 / (Na_eq**3 + 3375))
# 
# 
# 
# 
# 
N_pp = 160
N_ii = 40
N_ip = 40
N_pi = 160
N_pP = 8 
N_iP = 2 
# 
# 
# 
# 
g_ampa_pyr,g_ampa_inh = 1, 1
g_ampa_noise = 1
g_ampa_input = 1
g_gaba_pyr,g_gaba_inh = 1, 1
phi_n_sd = 2.5
# 
# 
# 
# 
i_c_m = 1/c_m
i_tau_Na = 1/tau_Na
i_tau_p,i_tau_i = 1/tau_p,1/tau_i
INTRA_CONN_EXC = np.array([N_pp,N_ip]).reshape(2,-1)
INTRA_CONN_INH = np.array([N_pi,N_ii]).reshape(2,-1)
INTER_CONN_EXC = np.array([N_pP,N_iP]).reshape(2,-1)
# 
# 
# 
# 
G_AMPA = np.array([g_ampa_pyr,g_ampa_inh]).reshape(2,-1)
G_AMPA_NOISE = np.array([g_ampa_noise,g_ampa_noise]).reshape(2,-1)
G_AMPA_INPUT = np.array([g_ampa_input,g_ampa_input]).reshape(2,-1)
G_GABA = np.array([g_gaba_pyr,g_gaba_inh]).reshape(2,-1)
# 
# 
# 
# 
beta = 0
i_sigma_p = 0.5 * np.pi/sigma_p/np.sqrt(3)
i_sigma_i = 0.5 * np.pi/sigma_i/np.sqrt(3)

n_pop = 1
# 
# 
# 
# 
G_AMPA_INPUT = np.array([0,0]).reshape(2,-1)
G_AMPA_NOISE = np.array([0,0]).reshape(2,-1)
beta = 0
# 
# 
# 
# 
T = 10 # simulation time in seconds
dt = 0.1 # timestep in ms
n = int(T*1000/dt)
# 
# 
# 
# 
#     
input_no = np.zeros((2,n_pop))
# 
# 
# 
# 
# 
# 
# 
# 
G_AMPA = [1,2]
Dynamic = {}
Dynamic["STEADY_STATE"] = {}
Dynamic["LAMDA"] = {}




entropy_seed = 12345
RNG_init = np.random.default_rng(np.random.SeedSequence(entropy=entropy_seed,spawn_key=(0,),))



for g_ampa in G_AMPA:
    Dynamic["STEADY_STATE"]["{}".format(g_ampa)] = {}
    if g_ampa==1:
        g_ampa_pyr,g_ampa_inh=g_ampa,g_ampa
        g_gaba_pyr,g_gaba_inh=1,1
    else:
        g_ampa_pyr,g_ampa_inh=g_ampa,g_ampa
        g_gaba_pyr,g_gaba_inh=2.294,2.313

    G_AMPA = np.array([g_ampa_pyr,g_ampa_inh]).reshape(2,-1)
    G_GABA = np.array([g_gaba_pyr,g_gaba_inh]).reshape(2,-1) 
    l1 = np.zeros((23,n_pop),dtype=float)           #stocastic term to be added for the SDE solution in the cortical network
    
    for XX in [0,10,50,100,200]:
        data_cor = np.zeros((23, n_pop,n), dtype=float)
        data_cor[0,:,0] = -10*RNG_init.random(n_pop) + theta_p
        data_cor[1,:,0] = -10*RNG_init.random(n_pop) + theta_i
        data_cor[2:,:,0] = 0.01*RNG_init.random((21,n_pop))
        data_cor[[15,16,17,18],:,0] = 0
        for i in range(n - 1):
            # implimenting sensory noise
            # 
            # 
            # 
            # 
            if i==50000:
                data_cor[3,:, i] += y_e**2*XX
                data_cor[7,:, i] += y_e**2*XX
            data_cor[:,:, i+1] = RK2order_Cor(dt,data_cor[:,:, i],l1,input_no)
        # 
        Dynamic["STEADY_STATE"]["{}".format(g_ampa)]["XX_{}".format(XX)] = data_cor
        if XX==0:
            Dynamic["LAMDA"]["{}".format(g_ampa)] = LAMDA(g_ampa_pyr,g_ampa_inh,g_gaba_pyr,g_gaba_inh,Dynamic["STEADY_STATE"]["{}".format(g_ampa)]["XX_0"][0,0,-1],Dynamic["STEADY_STATE"]["{}".format(g_ampa)]["XX_0"][1,0,-1])
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
# 
# 
# 
# 
scale=0.45
labelsize_font = 16.5 *scale
panel_font = 36 *scale
label_font = 21 *scale


fix_x = 7 #unit in inch
ws = 0.2 #of the axis x lenght
distx = 0.4 #of the axis x lenght
disty = 0.8
ax_x = fix_x/(2+1*ws)
ax_y = 0.525 * ax_x

fix_y = (3+ws+ disty) * ax_y

gs1 = GridSpec(1,2, bottom=(2+ws+disty)*ax_y/fix_y, top=1, left=0., right=1, wspace=ws)
gs2 = GridSpec(2,2, bottom=0, top=(2+ws)*ax_y/fix_y, left=0, right=1,wspace=ws)

fig = plt.figure(figsize=(fix_x, fix_y))


ax1 = []
for i in range(1):
    for j in range(2):
        ax1.append(fig.add_subplot(gs1[i,j]))
ax1 = np.array(ax1).reshape(1,2)
ax2 = []
for i in range(2):
    for j in range(2):
        ax2.append(fig.add_subplot(gs2[i,j]))
ax2 = np.array(ax2).reshape(2,2)

TITLE=["NREM", "Wakefulness"]

COL = ["lightblue","cyan","deepskyblue","b"]
COL = ["cyan","deepskyblue","b","darkblue"]

for j in range(2):
    ax1[0,j].plot([-1,0.2],[0,0],"k",lw=0.1)
    ax1[0,j].plot([0,0],[-0.07,0.07],"k",lw=0.1)
    ax1[0,j].plot(Dynamic["LAMDA"]["{}".format(j+1)].real,Dynamic["LAMDA"]["{}".format(j+1)].imag,"or",markersize=3)
    ax1[0,j].tick_params(axis='both', labelsize=labelsize_font)
    ax1[0,j].set_xticks([-0.8,-0.4,0])
    ax1[0,j].set_yticks([-0.05,0,0.05])
    ax1[0,j].set_xlim([-1,0.2])
    ax1[0,j].set_ylim([-0.07,0.07])
    ax1[0,j].set_title(TITLE[j],fontsize=label_font)
    ax1[0,j].set_xlabel("Real part",fontsize=label_font)
ax1[0,0].set_ylabel("Imaginary part",fontsize=label_font)
   
for i in range(2):
    for j in range(2):  
        for k,XX in enumerate([10,50,100,200]):
            ax2[i,j].plot(Dynamic["STEADY_STATE"]["{}".format(j+1)]["XX_{}".format(XX)][i,0,:],color=COL[k])
        ax2[i,j].tick_params(axis='both', labelsize=labelsize_font)
        ax2[i,j].set_xticks([50000,55000,60000,65000])
        ax2[i,j].set_xticklabels([0,0.5,1,1.5])
        ax2[i,j].set_xlim(48000-100,68000)
        
        ax2[i,j].legend([10,50,100,200],fontsize=labelsize_font,loc='upper right')
        
    if i==0:   
        ax2[i,0].set_ylabel(r"$V_{p}$"+" (mV)",fontsize=label_font)
        
    else:
        ax2[i,0].set_ylabel(r"$V_{i}$"+" (mV)",fontsize=label_font)
        
for j in range(2):
    ax2[0,j].set_ylim(-62,-39)
    ax2[1,j].set_ylim(-59,-20)
    ax2[0,j].set_title(TITLE[j],fontsize=label_font)
    ax2[-1,j].set_xlabel("Time (s)",fontsize=label_font)



fig.align_ylabels([ax1[0,0],ax2[0,0],ax2[1,0]])


plt.text(-0.091,1.03,"A",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(-0.091,0.61,"B",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.savefig("./DYNAMICAL_ANALYSIS.pdf",bbox_inches = 'tight', pad_inches = 0)
