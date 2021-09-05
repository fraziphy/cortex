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
vig_state = {"NREM":"1_1", "WAKE":"2_1"}
# 
# 
# 
# 
pops=["Perturbed","Unperturbed"]
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
scale=0.5
labelsize_font = 16.5 *scale
panel_font = 36 *scale
label_font = 21 *scale

t = np.arange(n_total)
xticks = np.linspace(0,n_total,11)



fix_x = 7 #unit in inch
ws = 0.3 #of the axis x lenght
distx = 0.7 #of the axis x lenght
disty = 1
ax_x = fix_x/(3+ws+distx)
ax_y = 0.625 * ax_x

fix_y = (2*(3 + 2*ws) + disty) * ax_y


gs1 = GridSpec(3,2,bottom=(3+2*ws+disty)*ax_y/fix_y , top=1, left=0., right=(2+ws)*ax_x/fix_x, wspace=ws,hspace = ws)

gs2 = GridSpec(3,2,bottom=(3+2*ws+disty)*ax_y/fix_y , top=1 , left=(2+ws+distx)*ax_x/fix_x, right=(2+ws+distx)*ax_x/fix_x+(2+ws)*ax_x/fix_x, wspace=ws,hspace = ws)

fig = plt.figure(figsize=(fix_x, fix_y))

ax1 = []
for i in range(3):
    for j in range(2):
        ax1.append(fig.add_subplot(gs1[i,j]))
ax1 = np.array(ax1).reshape(3,2)
ax2 = []
for i in range(3):
    for j in range(2):
        ax2.append(fig.add_subplot(gs2[i,j]))
ax2 = np.array(ax2).reshape(3,2)

ax = [ax1,ax2]

ylim_sync = [None] * 3
cc=["g","k"]
f = open(root_fig+'p_value_POSTSTIMULUS_PLOT_TWO_CORTICAL_COLUMN_{}.txt'.format(sub_pop), 'w')
f.write("Panel A\n\n")
f.write(" "*20+"{:<10}".format("NREM")+" "*20+"{:<10}\n\n".format("WAKE"))
for ii,pop in enumerate(pops):
    for i,signal in enumerate(signals):
        for j,state in enumerate(vig_state):
            ax[ii][i,j].plot(DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig_state[state]]["Difference"],color=cc[j])
            ax[ii][i,j].set_xticks(xticks)
            ax[ii][i,j].set_xticklabels((xticks/1000*dt).astype(int) - 5)
            ax[ii][i,j].tick_params(axis='both', labelsize=labelsize_font)
            ax[ii][i,j].set_xlim(0,xticks[-1])
        ax[ii][i,0].set_ylabel("Average\n {} ({})".format(SIGNAL[i],unit[i]),fontsize=label_font)
    ax[ii][2,0].set_xlabel("Time (s)",fontsize=label_font)
    ax[ii][2,1].set_xlabel("Time (s)",fontsize=label_font)
    ax[ii][0,0].set_title("{}".format("NREM"),loc="center",fontsize=label_font)
    ax[ii][0,1].set_title("{}".format("Wakefulness"),loc="center",fontsize=label_font)
    for i in range(3):
        for j in range(2):
            ax[ii][i,j].set_xlim(xticks[4]+8000,xticks[-4])
    if ii==0:
        for i,signal in enumerate(signals):
            f.write("{:^20}".format(signal))
            ylim_min,ylim_max = np.inf,-np.inf
            for j,state in enumerate(vig_state):
                ylim= ax[ii][i,j].get_ylim()
                ylim_min,ylim_max = min(ylim_min,ylim[0]),max(ylim_max,ylim[1])
            ylim = ylim_min,ylim_max
            for j,state in enumerate(vig_state):
                ax[ii][i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[0]-0.02*(ylim[1]-ylim[0])]*2,color="k")
                ax[ii][i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[1]]*2,color="gray",alpha=0.2)
                ax[ii][i,j].plot([n]*2,[ylim[0],ylim[1]],"--k",alpha=0.9)
                index = DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig_state[state]]["Significant_Index"]
                if ~np.isnan(index).all():         
                    for k in range(index.shape[0]):
                        f.write("{:<10}".format(DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig_state[state]]["P_Value"][k]))
                        ax[ii][i,j].fill_between(index[k],[ylim[1]]*2,[ylim[1]+0.02*(ylim[1]-ylim[0])]*2,color="red",label="Significance")

                else:
                    f.write(" "*20)
                ax[ii][i,j].set_ylim(ylim[0]-0.04*(ylim[1]-ylim[0]),ylim[1]+0.04*(ylim[1]-ylim[0]))
                f.write(" "*10)
            f.write("\n\n")
        f.write("-"*80)
        f.write("\n\n")
        f.write("Panel B\n\n")
        f.write(" "*20+"{:<10}".format("NREM")+" "*20+"{:<10}\n\n".format("WAKE"))
    else:
        for i,signal in enumerate(signals):
            f.write("{:^20}".format(signal))
            ylim_min,ylim_max = np.inf,-np.inf
            for j,state in enumerate(vig_state):
                ylim= ax[ii][i,j].get_ylim()
                ylim_min,ylim_max = min(ylim_min,ylim[0]),max(ylim_max,ylim[1])
            ylim = ylim_min,ylim_max
            for j,state in enumerate(vig_state):
                ax[ii][i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[0]-0.02*(ylim[1]-ylim[0])]*2,color="k")
                ax[ii][i,j].fill_between([n,n + dur],[ylim[0]]*2,[ylim[1]]*2,color="gray",alpha=0.2)
                ax[ii][i,j].plot([n]*2,[ylim[0],ylim[1]],"--k",alpha=0.9)
                index = DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig_state[state]]["Significant_Index"]
                if ~np.isnan(index).all():         
                    for k in range(index.shape[0]):
                        f.write("{:<10}".format(DATA_PLOT['WITHIN_POP'][pop][sub_pop][signal][vig_state[state]]["P_Value"][k]))
                        ax[ii][i,j].fill_between(index[k],[ylim[1]]*2,[ylim[1]+0.02*(ylim[1]-ylim[0])]*2,color="red",label="Significance")

                else:
                    f.write(" "*20)
                ax[ii][i,j].set_ylim(ylim[0]-0.04*(ylim[1]-ylim[0]),ylim[1]+0.04*(ylim[1]-ylim[0]))
                f.write(" "*10)
            f.write("\n\n")
f.close()

fig.align_ylabels([ax[0][:,0]])
fig.align_ylabels([ax[1][:,0]])

plt.text(0.6,1.12,"Two Cortical Column Model",horizontalalignment='center', transform=fig.transFigure, fontsize=label_font*1.3)

titlefig = {"pyr":"Pyramidal", "inh":"Inhibitory"}
plt.text(-0.09,1.06,"A",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.009,1.06,"Perturbed Column --- {} Population".format(titlefig[sub_pop]),horizontalalignment='left', transform=fig.transFigure, fontsize=label_font*1.2)
plt.text(0.65,1.06,"B",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.74,1.06,"Unperturbed Column --- {} Population".format(titlefig[sub_pop]),horizontalalignment='left', transform=fig.transFigure, fontsize=label_font*1.2)
plt.savefig(root_fig+"POSTSTIMULUS_PLOT_TWO_CORTICAL_COLUMN_{}.pdf".format(sub_pop),bbox_inches = 'tight', pad_inches = 0)
