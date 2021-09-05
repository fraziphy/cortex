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
# 
# 
# 
# 

scale=0.4
labelsize_font = 16.5 *scale
panel_font = 36 *scale
label_font = 21 *scale

t = np.arange(n_total)
xticks = np.linspace(0,n_total,11)


lfp_ticks=[0,100,200]
firing_ticks=[15,20,24,28]
v_ticks = [-58,-54,-50]
yticks_ax=[lfp_ticks,firing_ticks,v_ticks]






fix_x = 4 #unit in inch
ws = 0.3 #of the axis x lenght
distx = 0.7 #of the axis x lenght

ax_x = fix_x/(2+ws)
ax_y = 0.625 * ax_x

fix_y = ((2 + 1*ws)) * ax_y


gs1 = GridSpec(2,2,bottom=0 , top=1, left=0., right=1, wspace=ws,hspace = ws)

fig = plt.figure(figsize=(fix_x, fix_y))

ax = []
for i in range(2):
    for j in range(2):
        ax.append(fig.add_subplot(gs1[i,j]))
ax = np.array(ax).reshape(2,2)

        
G_AMPA=[2,6]
BETA=[1,2]

ylim_sync = [None] * 3
cc=["g","k"]
pop="Unperturbed"
signal="firing"
f = open(root_fig+'p_value_POSTSTIMULUS_PLOT_UNPERTURBED_CORTICAL_COLUMN_{}.txt'.format(sub_pop), 'w')
f.write(" "*20+"{:<10}".format("beta=1")+" "*20+"{:<10}\n\n".format("beta=2"))
ylim_min,ylim_max = np.inf,-np.inf
for i,g_ampa in enumerate(G_AMPA):
    for j,beta in enumerate(BETA):
        vig = "{}_{}".format(g_ampa,beta)
        ax[i,j].plot(DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig]["Difference"],color="k")
        ax[i,j].set_xticks(xticks)
        ax[i,j].set_xticklabels((xticks/1000*dt).astype(int) - 5)
        ax[i,j].tick_params(axis='both', labelsize=labelsize_font)
        ax[i,j].set_xlim(0,xticks[-1])
        ylim= ax[i,j].get_ylim()
        ylim_min,ylim_max = min(ylim_min,ylim[0]),max(ylim_max,ylim[1])
ylim = ylim_min,ylim_max
for i,g_ampa in enumerate(G_AMPA):
    f.write("{:^20}".format("g_ampa = {}".format(g_ampa)))
    for j,beta in enumerate(BETA):
        vig = "{}_{}".format(g_ampa,beta)
        ax[i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[0]-0.02*(ylim[1]-ylim[0])]*2,color="k")
        ax[i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[1]]*2,color="gray",alpha=0.2)
        ax[i,j].plot([n]*2,[ylim[0],ylim[1]],"--k",alpha=0.9)
        index = DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig]["Significant_Index"]
        if ~np.isnan(index).all():         
            for k in range(index.shape[0]):
                f.write("{:<10}".format(DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig]["P_Value"][k]))
                ax[i,j].fill_between(index[k],[ylim[1]]*2,[ylim[1]+0.02*(ylim[1]-ylim[0])]*2,color="red",label="Significance")

        else:
            f.write(" "*20)
        ax[i,j].set_ylim(ylim[0]-0.04*(ylim[1]-ylim[0]),ylim[1]+0.04*(ylim[1]-ylim[0]))
        ax[i,j].set_xlim(xticks[4]+8000,xticks[-4])
        ax[i,j].annotate(r'$\bar{g}_{AMPA} =$'+"{}".format(g_ampa), (0.60, 0.850),xycoords='axes fraction', ha='left',fontsize=label_font,color="k")
        ax[0,j].set_title(r"$\beta=$"+"{}".format(beta),loc="center",fontsize=label_font)
        f.write(" "*10)
    f.write("\n\n")
    ax[i,0].set_ylabel("Average\n Firing rate (Hz)",fontsize=label_font)
ax[1,0].set_xlabel("Time (s)",fontsize=label_font)
ax[1,1].set_xlabel("Time (s)",fontsize=label_font)
f.close()


fig.align_ylabels([ax[:,0]])

titlefig = {"pyr":"Pyramidal", "inh":"Inhibitory"}
plt.text(0.5,1.15,"Unperturbed Column --- {} Population".format(titlefig[sub_pop]),horizontalalignment='center', transform=fig.transFigure, fontsize=label_font*1.5)
plt.savefig(root_fig+"POSTSTIMULUS_PLOT_UNPERTURBED_CORTICAL_COLUMN_{}.pdf".format(sub_pop),bbox_inches = 'tight', pad_inches = 0)
