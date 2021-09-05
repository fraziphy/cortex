import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
import pickle
import time
import math
import os
# 
# 
# 
# 
sub_pop = os.environ.get("SUB_POP")
root_c = os.environ.get("rootMACHINE")
root = root_c +"/DATA/" 
root_fig = root_c+"/" 
# 
# 
# 
# 
# 
# 
# 
# 
signals=["lfp","firing","v"]
SIGNAL=["LFP","Firing Rate","V"+r"$_{p}$"]
if sub_pop == "inh":
    SIGNAL=["LFP","Firing Rate","V"+r"$_{i}$"]
unit=["mV", "Hz", "mV"]
TITLE=["NREM", "Wakefulness"]
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

significance_level = 0.01
# 
# 
# 
# 
with open(root+"POSTSTIMULUS_DATA_TWO_CORTICAL_COLUMN.pickle", 'rb') as handle:
    DATA_PLOT = pickle.load(handle)


t_crit = DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["T_Critical_BON"]
t_stat = DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["t_stat"]
index=DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["Initial_T_Cluster_AND_Index_Info_BON"][1]
t_n=len(DATA_PLOT["EXAMPLE_T_CLUSTER"][sub_pop]["Initial_T_Cluster_AND_Index_Info_BON"][0])



scale=0.45
labelsize_font = 16.5 *scale
panel_font = 36 *scale
label_font = 21 *scale


fix_x = 7 #unit in inch
ws = 0.3 #of the axis x lenght
distx = 0.4 #of the axis x lenght
disty = 0.5
ax_x = fix_x/(2+ws+distx)
ax_y = 0.625 * ax_x

fix_y = (2+disty) * ax_y

gs1 = GridSpec(1,1, bottom=0, top=1, left=0., right=(1)*ax_x/fix_x)
gs2 = GridSpec(2,1, bottom=0, top=1, left=(1+distx)*ax_x/fix_x, right=1,hspace=disty)

fig = plt.figure(figsize=(fix_x, fix_y))


ax1 = []
ax1.append(fig.add_subplot(gs1[0]))
ax2 = []
for i in range(2):
    for j in range(1):
        ax2.append(fig.add_subplot(gs2[i,j]))
ax2 = np.array(ax2).reshape(2,1)


tt=np.arange(50000)

ax1[0].plot(t_stat,"k")
ax1[0].plot([0,50000],[t_crit,t_crit],"--",color="y")
ax1[0].plot([0,50000],[-t_crit,-t_crit],"--",color="y")
ax1[0].plot([0,t_stat.shape[0]-1],[0 ,0],zorder=1,color="gray",linewidth=0.5)
ax1[0].fill_between(tt[index[0][0]:index[0][1]],t_stat[index[0][0]:index[0][1]],[t_crit]*(index[0][1]-index[0][0]),color="gray",alpha=0.5)
ylim= ax1[0].get_ylim()
ax1[0].fill_between([0,1000],[ylim[0]]*2,[ylim[0]-0.01*(ylim[1]-ylim[0])]*2,color="k")
ax1[0].tick_params(axis='both', labelsize=labelsize_font)
ax1[0].set_xlim(0,50000)
ax1[0].set_xticks(np.linspace(0,50000,6))
ax1[0].set_xticklabels([0,1,2,3,4,5])
ax1[0].set_yticks([-t_crit,0,t_crit])
ax1[0].set_ylabel("t-statistics",fontsize=label_font)
ax1[0].set_xlabel("Time (s)",fontsize=label_font)
divider = make_axes_locatable(ax1[0])
ax_1 = divider.append_axes("bottom", size="20%", pad=0.4)
ax_1.plot([0,t_stat.shape[0]-1],[0 ,0],zorder=1,color="gray",linewidth=0.5)


# 
# 
# 
#       
dist=0.4



for i in range(t_n):
    ax_1.fill_between(index[i],[0.2 ,0.2],[-0.2 ,-0.2],zorder=2,color="r")
ax_1.tick_params(axis='both', labelsize=labelsize_font)
ax_1.set_xlim(0,50000)
ax_1.set_ylim(-2,2)
ax_1.text(25000,-2.5,"Bonferroni correction cluster",ha='center',fontsize=label_font)
ax_1.axis('off')



signal ="firing"
pop = "Unperturbed"
g_ampa = DATA_PLOT["G_AMPA"]
beta = DATA_PLOT["BETA"]
data = np.zeros((len(g_ampa),len(beta)),dtype=float)
data1 = np.zeros((len(g_ampa),len(beta)),dtype=float)
color = [[None]*5]*3
for i,g in enumerate(g_ampa):
    col = [None]*5
    data_aux = [None]*5
    data_aux1 = [None]*5
    for j,b in enumerate(beta):
        vig = "{}_{}".format(g,b)
        a = None
        a = DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig]
        BON = a['Initial_T_Cluster_AND_Index_Info_BON']
        data_aux[j] = np.diff(a['Initial_T_Cluster_AND_Index_Info'][1])[0]/10
        data_aux1[j] = a['Initial_T_Cluster_AND_Index_Info'][0][0]/10000
        if BON is not None:
            data_aux[j] = np.diff(BON[1])[0]/10
            data_aux1[j] = BON[0][0]/10000
            col[j] = "r"
        else:
            data_aux[j] = 0
            data_aux1[j] = 0
            col[j] = "k"
    data[i] = data_aux
    data1[i] = data_aux1
    color[i]=col


g_color=["violet","darkviolet","indigo"]

for i in range(3):
    ax2[0,0].plot(beta,data1[i],color=g_color[i],zorder=1,alpha=0.8)
    ax2[0,0].legend([r"$\bar{g}_{AMPA} = 2$",r"$\bar{g}_{AMPA} = 6$",r"$\bar{g}_{AMPA} = 10$"], fontsize=labelsize_font)
    ax2[0,0].scatter(beta,data1[i],marker="s",s=[7]*5,zorder=2,color=color[i],alpha=0.8)
    
    ax2[0,0].set_xticks([1,2,3,4,5])
    #ax2[0,0].set_yticks([0,0.4,0.8,1.2])
    ax2[0,0].tick_params(axis='both', labelsize=labelsize_font)
ax2[0,0].set_xlabel(r"$\beta$",fontsize=label_font)
ax2[0,0].set_ylabel("Area (s)",fontsize=label_font)

for i in range(3):
    ax2[1,0].plot(beta,data[i],color=g_color[i],zorder=1,alpha=0.8)
    ax2[1,0].legend([r"$\bar{g}_{AMPA} = 2$",r"$\bar{g}_{AMPA} = 6$",r"$\bar{g}_{AMPA} = 10$"], fontsize=labelsize_font)
    ax2[1,0].scatter(beta,data[i],marker="s",s=[7]*5,zorder=2,color=color[i],alpha=0.8)
    
    ax2[1,0].set_xticks([1,2,3,4,5])
    #ax2[1,0].set_yticks([25,75,125,175])
    ax2[1,0].tick_params(axis='both', labelsize=labelsize_font)
ax2[1,0].set_xlabel(r"$\beta$",fontsize=label_font)
ax2[1,0].set_ylabel("Duration (ms)",fontsize=label_font)
fig.align_ylabels(ax2[:,0])

titlefig = {"pyr":"Pyramidal", "inh":"Inhibitory"}
plt.text(0.5,1.15,"Unperturbed Column --- {} Population".format(titlefig[sub_pop]),horizontalalignment='center', transform=fig.transFigure, fontsize=label_font*1.5)
plt.text(-0.08,1.025,"A",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.44,1.025,"B",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.44,0.455,"C",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.savefig(root_fig+"POSTSTIMULUS_PLOT_UNPERTURBED_CORTICAL_COLUMN_SCREENING_Bonferroni_{}.pdf".format(sub_pop),bbox_inches = 'tight', pad_inches = 0)
