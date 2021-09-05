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
pop = "Perturbed"
vig_state = {"NREM":"1_0", "WAKE":"2_0"}
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
with open(root+"POSTSTIMULUS_DATA_ONE_CORTICAL_COLUMN.pickle", 'rb') as handle:
    DATA_PLOT = pickle.load(handle)
# 
# 
# 
# 
scale=0.45
labelsize_font = 16.5 *scale
panel_font = 36 *scale
label_font = 21 *scale

t = np.arange(n_total)
xticks = np.linspace(0,n_total,11)


fix_x = 7 #unit in inch
ws = 0.25 #of the axis x lenght
dist = 0.7 #of the axis x lenght
ax_x = fix_x/(3+ws+dist)
ax_y = 0.625 * ax_x

fix_y = (3 + 2*ws) * ax_y


gs1 = GridSpec(3,2,bottom=0 , top=1, left=0., right=(2+ws)*ax_x/fix_x, wspace=ws,hspace = ws)
gs2 = GridSpec(3,1,bottom=0 , top=1, left=(2+ws+dist)*ax_x/fix_x, right=1,hspace = ws)

fig = plt.figure(figsize=(fix_x, fix_y))

ax = []
for i in range(3):
    for j in range(2):
        ax.append(fig.add_subplot(gs1[i,j]))
ax = np.array(ax).reshape(3,2)
ax_1 = []
for i in range(3):
    ax_1.append(fig.add_subplot(gs2[i]))

cc=["g","k"]
for i,signal in enumerate(signals):
    ylim_min,ylim_max = np.inf,-np.inf
    for j,state in enumerate(vig_state):
        ax[i,j].plot(DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig_state[state]]["Difference"],color=cc[j])
        ax[i,j].set_xticks(xticks)
        ax[i,j].set_xticklabels((xticks/1000*dt).astype(int) - 5)
        ax[i,j].tick_params(axis='both', labelsize=labelsize_font)
        ylim= ax[i,j].get_ylim()
        ylim_min,ylim_max = min(ylim[0],ylim_min),max(ylim[1],ylim_max)
    ax[i,0].set_ylim(ylim_min,ylim_max)
    ax[i,1].set_ylim(ylim_min,ylim_max)
    ax[i,0].set_ylabel("Average\n {} ({})".format(SIGNAL[i],unit[i]),fontsize=label_font)
ax[2,0].set_xlabel("Time (s)",fontsize=label_font)
ax[2,1].set_xlabel("Time (s)",fontsize=label_font)
ax[0,0].set_title("{}".format("NREM"),loc="center",fontsize=label_font)
ax[0,1].set_title("{}".format("Wakefulness"),loc="center",fontsize=label_font)
for i in range(3):
    for j in range(2):
        ax[i,j].set_xlim(xticks[4]+8000,xticks[-4])
f = open('p_value_POSTSTIMULUS_PLOT_ONE_CORTICAL_COLUMN_{}.txt'.format(sub_pop), 'w')
f.write("Panel A\n\n")
f.write(" "*20+"{:<10}".format("NREM")+" "*20+"{:<10}\n\n".format("WAKE"))
for i,signal in enumerate(signals):
    f.write("{:^20}".format(signal))
    for j,state in enumerate(vig_state):
        ylim= ax[i,j].get_ylim()
        ax[i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[0]-0.02*(ylim[1]-ylim[0])]*2,color="k")
        ax[i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[1]]*2,color="gray",alpha=0.2)
        ax[i,j].plot([n]*2,[ylim[0],ylim[1]],"--k",alpha=0.9)
        index = DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig_state[state]]["Significant_Index"]
        ylim_min,ylim_max = ax[i,j].get_ylim()
        if ~np.isnan(index).all():         
            for k in range(index.shape[0]):
                f.write("{:<10}".format(DATA_PLOT["WITHIN_POP"][pop][sub_pop][signal][vig_state[state]]["P_Value"][k]))
                ax[i,j].fill_between(index[k],[ylim_max]*2,[ylim_max+0.02*(ylim_max-ylim_min)]*2,color="red",label="Significance")
        ylim= ax[i,j].get_ylim()
        ylim_min,ylim_max = min(ylim[0],ylim_min),max(ylim[1],ylim_max)
        f.write(" "*10)
    f.write("\n\n")
    ax[i,0].set_ylim(ylim[0]-0.04*(ylim[1]-ylim[0]),ylim_max+0.04*(ylim_max-ylim_min))
    ax[i,1].set_ylim(ylim[0]-0.04*(ylim[1]-ylim[0]),ylim_max+0.04*(ylim_max-ylim_min))

f.write("-"*80)
f.write("\n\n")
f.write("Panel B\n\n")
for i,signal in enumerate(signals):
    ax_1[i].plot(DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig_state[state]]["Difference"],color="orange")
    ax_1[i].set_xticks(xticks)
    ax_1[i].set_xticklabels((xticks/1000*dt).astype(int) - 5)
    ax_1[i].tick_params(axis='both', labelsize=labelsize_font)
    ax_1[i].set_ylabel("Average\n"+r"$\Delta$"+ "{} ({})".format(SIGNAL[i],unit[i]),fontsize=label_font)
    ylim= ax_1[i].get_ylim()
    ylim_min,ylim_max = ax_1[i].get_ylim()
    ax_1[i].fill_between([n,n + dur],[ylim[0]]*2,[ylim[0]-0.02*(ylim[1]-ylim[0])]*2,color="k")
    ax_1[i].fill_between([n,n + dur],[ylim[0]]*2,[ylim[1]]*2,color="gray",alpha=0.2)
    ax_1[i].plot([n]*2,[ylim[0],ylim[1]],"--k",alpha=0.9)
    index = DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig_state[state]]["Significant_Index"]
    f.write("{:^20}".format(signal))
    if ~np.isnan(index).all():         
        for k in range(index.shape[0]):
            f.write("{:<10}".format(DATA_PLOT["ACROSS_STATE"][pop][sub_pop][signal][vig_state[state]]["P_Value"][k]))
            ax_1[i].fill_between(index[k],[ylim_max]*2,[ylim_max+0.02*(ylim_max-ylim_min)]*2,color="red",label="Significance")
    ax_1[i].set_ylim(ylim[0]-0.04*(ylim[1]-ylim[0]),ylim_max+0.04*(ylim_max-ylim_min))
    ax_1[i].set_xlim(xticks[4]+8000,xticks[-4])
    f.write("\n\n")
ax_1[2].set_xlabel("Time (s)",fontsize=label_font)
ax_1[0].set_title("{}".format("NREM vs. wakefulness"),loc="center",fontsize=label_font)
f.close()


fig.align_ylabels(ax[:,0])
fig.align_ylabels(ax_1[:])

titlefig = {"pyr":"Pyramidal", "inh":"Inhibitory"}
plt.text(0.5,1.15,"One Cortical Column Model --- {} Population".format(titlefig[sub_pop]),horizontalalignment='center', transform=fig.transFigure, fontsize=label_font*1.3)
plt.text(-0.09,1.047,"A",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.64,1.047,"B",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)

plt.savefig(root_fig+"POSTSTIMULUS_PLOT_ONE_CORTICAL_COLUMN_{}.pdf".format(sub_pop),bbox_inches = 'tight', pad_inches = 0)
