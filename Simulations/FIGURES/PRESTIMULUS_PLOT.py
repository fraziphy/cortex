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
cond = os.environ.get("COND")
sub_pop = os.environ.get("SUB_POP")
pop = os.environ.get("POP")
root_c = os.environ.get("rootMACHINE")
root = root_c +"/DATA/" 
root_fig = root_c+"/" 
# 
# 
# 
# 
vig_state = {"NREM":"1_0", "WAKE":"2_0"}
NAME = "ONE_CORTICAL_COLUMN"
if cond == "TWO_CORTICAL_COLUMN":
    vig_state = {"NREM":"1_1", "WAKE":"2_1"}
    NAME = "TWO_CORTICAL_COLUMN"
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
with open(root+"PRESTIMULUS_DATA_{}.pickle".format(NAME), 'rb') as handle:
    DATA_PLOT = pickle.load(handle)
# 
# 
# 
# 

def YTICKS(yticks):
    c = yticks[-1]
    if c//10 !=0:
        yticks = np.round(yticks)
    else:
        k = 0
        while c//10 ==0 and c%10<=10:
            c *= 10
            k += 1
        yticks = [math.trunc(i*10**k)/10**k for i in yticks]
    return yticks
# 
# 
# 
# 
scale=0.45
labelsize_font = 16.5 *scale
panel_font = 36 *scale
label_font = 21 *scale


fix_x = 7 #unit in inch
ws = 0.3 #of the axis x lenght
distx = 0.7 #of the axis x lenght
disty = 1
ax_x = fix_x/(4+2*ws+distx)
ax_y = 0.625 * ax_x

fix_y = (3+1.6 + 2*ws + disty) * ax_y

gs1 = GridSpec(3,2, bottom=(1.6+disty)*ax_y/fix_y, top=1, left=0., right=(2+ws)*ax_x/fix_x, wspace=ws,hspace = ws)
gs3 = GridSpec(2,2, bottom=(1.6+disty)*ax_y/fix_y, top=1, left=(2+ws+distx)*ax_x/fix_x, right=1,wspace=ws, hspace=2.5*ws)
gs2 = GridSpec(1,2, bottom=0, top=(1.6)*ax_y/fix_y, left=0, right=(2+ws)*ax_x/fix_x,wspace=ws)
gs4 = GridSpec(1,2, bottom=0, top=(1.6)*ax_y/fix_y, left=(2+ws+distx)*ax_x/fix_x, right=1, wspace=2*ws)

fig = plt.figure(figsize=(fix_x, fix_y))


ax1 = []
for i in range(3):
    for j in range(2):
        ax1.append(fig.add_subplot(gs1[i,j]))
ax1 = np.array(ax1).reshape(3,2)
ax2 = []
for i in range(2):
    ax2.append(fig.add_subplot(gs2[i]))
ax3 = []
for i in range(2):
    for j in range(2):
        ax3.append(fig.add_subplot(gs3[i,j]))
ax3 = np.array(ax3).reshape(2,2)
ax4 = []
for i in range(2):
    ax4.append(fig.add_subplot(gs4[i]))
# 
# 
# 
#
cc=["g","k"]
xticks = np.linspace(0,n,6)
for i,signal in enumerate(signals): 
    ylim_min,ylim_max = np.inf,-np.inf
    for j,state in enumerate(vig_state):
        ax1[i,j].plot(DATA_PLOT[pop][sub_pop][vig_state[state]][signal]["SAMPEL"],cc[j])
        ax1[i,j].tick_params(axis='both', labelsize=labelsize_font)
        ylim = ax1[i,j].get_ylim()
        ylim_min,ylim_max = min(ylim_min,ylim[0]),max(ylim_max,ylim[1])
        ax1[i,j].set_xticks([0,25000,50000])
        ax1[i,j].set_xticklabels([-5,-2.5,0])
        ax1[i,j].set_xlim(0,n)
    ax1[i,0].set_ylabel("{} ({})".format(SIGNAL[i],unit[i]),fontsize=label_font)
    if signal == "lfp":
        yticks_ax1 = ylim_min + 0.1*(ylim_max - ylim_min),0 ,ylim_max - 0.1*(ylim_max - ylim_min)
        yticks_ax1 = np.round(yticks_ax1).astype(int)
    elif signal == "v":
        yticks_ax1 = ylim_min + 0.1*(ylim_max - ylim_min), DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["SAMPEL"].mean() ,ylim_max - 0.1*(ylim_max - ylim_min)
        yticks_ax1 = np.round(yticks_ax1).astype(int)
    else:
        yticks_ax1 = 0,DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["SAMPEL"].mean(),30
        if sub_pop == "inh":
            yticks_ax1 = 0,DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["SAMPEL"].mean(),60
        yticks_ax1 = np.round(yticks_ax1).astype(int)
    #for j,state in enumerate(vig_state):
            #ax1[i,j].set_yticks(yticks_ax1)
            #ax1[i,j].set_yticklabels(yticks_ax1)
    for j in range(2):
        ax1[i,j].set_ylim(ylim_min,ylim_max)
for j in range(2):
    ax1[0,j].set_title(TITLE[j],fontsize=label_font)
    ax1[-1,j].set_xlabel("Time (s)",fontsize=label_font)
    ax1[1,j].set_ylim(-1,31)
    if sub_pop == "inh":
        ax1[1,j].set_ylim(-1,61)
# 
# 
# 
#
for i,state in enumerate(vig_state):
    ylim_min,ylim_max = np.inf,-np.inf
    xlim_min,xlim_max = np.inf,-np.inf
    ax2[i].plot(DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["BINS"],DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["HIST"],color=cc[i])
    ax2[i].tick_params(axis='both', labelsize=labelsize_font)
    ylim = ax2[i].get_ylim()
    ylim_min,ylim_max = min(ylim_min,ylim[0]),max(ylim_max,ylim[1])
    xlim = ax2[i].get_xlim()
    xlim_min,xlim_max = min(xlim_min,xlim[0]),max(xlim_max,xlim[1])
for i in range(2):
    yticks_ax2 = YTICKS(np.linspace(0,ylim_max,4))
    #ax2[i].set_yticks(yticks_ax2)
    #ax2[i].set_yticklabels(yticks_ax2)
    ax2[i].set_ylim(0,ylim_max)
    xticks_ax2 = xlim_min + 0.1*(xlim_max - xlim_min),0 ,xlim_max - 0.1*(xlim_max - xlim_min)
    xticks_ax2 = np.round(xticks_ax2).astype(int)
    #ax2[i].set_xticks(xticks_ax2)
    ax2[i].set_xlim(xlim_min,xlim_max)
    ax2[i].set_title(TITLE[i],fontsize=label_font)
    ax2[i].set_xlabel("LFP (mV)",fontsize=label_font)
ax2[0].set_ylabel("Occurrence %",fontsize=label_font)
# 
# 
# 
#
for i,signal in enumerate(["firing","v"]):
    ylim_min,ylim_max = np.inf,-np.inf
    xlim_min,xlim_max = np.inf,-np.inf
    for j,state in enumerate(vig_state):
        if state=="WAKE":
            ax3[i,1].plot(DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["BINS"],DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["HIST"],"k")
        else:
            ax3[i,0].plot(DATA_PLOT[pop][sub_pop][vig_state["NREM"]][signal]["UP_BINS"],DATA_PLOT[pop][sub_pop][vig_state["NREM"]][signal]["UP_HIST"],"darkgreen")
            ax3[i,0].plot(DATA_PLOT[pop][sub_pop][vig_state["NREM"]][signal]["DOWN_BINS"],DATA_PLOT[pop][sub_pop][vig_state["NREM"]][signal]["DOWN_HIST"],"limegreen")
            ax3[i,0].annotate('Up', (0.70, 0.50),xycoords='axes fraction', ha='left',fontsize=label_font,color="darkgreen")
            ax3[i,0].annotate('Down', (0.20, 0.50),xycoords='axes fraction', ha='left',fontsize=label_font,color="limegreen")
        ax3[i,j].tick_params(axis='both', labelsize=labelsize_font)
        ylim = ax3[i,j].get_ylim()
        ylim_min,ylim_max = min(ylim_min,ylim[0]),max(ylim_max,ylim[1])
        xlim = ax3[i,j].get_xlim()
        xlim_min,xlim_max = min(xlim_min,xlim[0]),max(xlim_max,xlim[1])
    arg = np.argmax(DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["HIST"])
    arg = DATA_PLOT[pop][sub_pop][vig_state["WAKE"]][signal]["BINS"][arg]
    if signal == "v":
        xticks_ax3 = xlim_min + 0.01*(xlim_max - xlim_min), arg ,xlim_max - 0.01*(xlim_max - xlim_min)
        xticks_ax3 = np.round(xticks_ax3).astype(int)
    else:
        xticks_ax3 = 0, arg ,30
        if sub_pop == "inh":
            xticks_ax3 = 0, arg ,60
        xticks_ax3 = np.round(xticks_ax3).astype(int)
    for j,state in enumerate(vig_state):
        #ax3[i,j].set_xticks(xticks_ax3)
        #ax3[i,j].set_xticklabels(xticks_ax3)
        yticks_ax3 = YTICKS(np.linspace(0,ylim_max,4))
        ax3[i,j].set_yticks(yticks_ax3)
        ax3[i,j].set_yticklabels(yticks_ax3)
        if signal == "v":
            ax3[i,j].set_xlim(xlim_min,xlim_max)
        else:
            ax3[i,j].set_xlim(-1,31)
            if sub_pop == "inh":
                ax3[i,j].set_xlim(-1,61)
            
    ylim0 = ax3[i,0].get_ylim()
    ylim1 = ax3[i,1].get_ylim()
    ax3[i,0].set_ylim(0,max(ylim0[1],ylim1[1]))
    ax3[i,1].set_ylim(0,max(ylim0[1],ylim1[1]))

for i in range(2):
    for j in range(2):
        ax3[i,j].set_xlabel("{} ({})".format(SIGNAL[i+1],unit[i+1]),fontsize=label_font)
    ax3[0,i].set_title(TITLE[i],fontsize=label_font)
    ax3[i,0].set_ylabel("Occurrence %",fontsize=label_font)  
# 
# 
# 
# 
colors = ["g","k"]
for i,state in enumerate(vig_state):
    ax4[0].plot(DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FREQUENCY"],DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_M"],color=colors[i])
    ax4[0].fill_between(DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FREQUENCY"],DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_M"]+DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_SD"],DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_M"]-DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_SD"],color=colors[i],alpha=0.2)
ylim = ax4[0].get_ylim()
yticks_ax4 = np.linspace(ylim[0]+0.1*(ylim[1] - ylim[0]),ylim[1]-0.1*(ylim[1] - ylim[0]),4)
yticks_ax4 = np.round(yticks_ax4,1)
ax4[0].set_yticks(yticks_ax4)
ax4[0].set_xlim(-5,100)
ax4[0].set_ylabel("LFP Log Power",fontsize=label_font)
ax4[0].set_xlabel("Frequency (Hz)",fontsize=label_font)
ax4[0].tick_params(axis='both', labelsize=labelsize_font)
ax4[0].legend(["NREM","W"],fontsize=labelsize_font,loc='upper right')
for i,state in enumerate(vig_state):
    ax4[1].plot(DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_RATIO_BINS"],DATA_PLOT[pop][sub_pop][vig_state[state]]["lfp"]["FFT_RATIO"],color=colors[i])
ax4[1].tick_params(axis='both', labelsize=labelsize_font)
ax4[1].set_ylabel("Occurrence %",fontsize=label_font)
ax4[1].set_xlabel("High-/low-frequency \npower ratio (log scale)",fontsize=label_font)
ylim = ax4[1].get_ylim()
yticks_ax4 = YTICKS(np.linspace(0,ylim[1],4))
ax4[1].set_yticks(yticks_ax4)
ax4[1].set_ylim(0,ylim[1])
xlim = ax4[1].get_xlim()
xticks_ax4 = np.linspace(xlim[0]+0.1*(xlim[1] - xlim[0]),xlim[1]-0.1*(xlim[1] - xlim[0]),3)
xticks_ax4 = np.round(xticks_ax4,1)
ax4[1].set_xticks(xticks_ax4)


fig.align_ylabels([ax1[0,0],ax1[1,0],ax1[2,0],ax2[0]])
fig.align_ylabels([ax3[0,0],ax3[1,0],ax4[0]])

ylbl = ax1[2,0].yaxis.get_label()
ylbl = ylbl.get_position()

titlefig = {"pyr":"Pyramidal", "inh":"Inhibitory"}
if NAME=="ONE_CORTICAL_COLUMN":
    plt.text(0.5,1.12,"One Cortical Column Model --- {} Population".format(titlefig[sub_pop]),horizontalalignment='center', transform=fig.transFigure, fontsize=label_font*1.3)
else:
    plt.text(0.5,1.12,"Two Cortical Column Model --- {} Population".format(titlefig[sub_pop]),horizontalalignment='center', transform=fig.transFigure, fontsize=label_font*1.3)
plt.text(-0.075,1.035,"A",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.485,1.035,"C",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(-0.075,0.3,"B",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.text(0.485,0.3,"D",horizontalalignment='left', transform=fig.transFigure, fontsize=panel_font)
plt.savefig(root_fig+"PRESTIMULUS_PLOT_{}_{}_{}.pdf".format(NAME,pop,sub_pop),bbox_inches = 'tight', pad_inches = 0)
