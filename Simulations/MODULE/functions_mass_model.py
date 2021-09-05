import numpy as np
# 
# 
# 
# 
#
# 
# 
# 
# 
n_pop = 2
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
beta = 1
# 
# 
# 
# 
i_sigma_p = 0.5 * np.pi/sigma_p/np.sqrt(3)
i_sigma_i = 0.5 * np.pi/sigma_i/np.sqrt(3)
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
