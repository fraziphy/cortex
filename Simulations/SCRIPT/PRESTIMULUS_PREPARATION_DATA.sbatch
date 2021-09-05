#!/bin/bash
##############################################################
#                         Data Curation                      #
#                                                            #
#              -Prestimulus intervals preparation            #
##############################################################
#SBATCH --time=0-00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=farhad.razi01@estudiant.upf.edu
#SBATCH --output=output_prestimulus_preparation_%x.out
#SBATCH --error=error_prestimulus_preparation_%x.err
# 
# 
# 
# 
ml matplotlib/2.1.2-foss-2018a-Python-3.6.4
# 
# 
# 
# 
export G_AMPA_BETA=$1 SUB_POP=$2 POP=$3 ROOT=$4 ROOT_C=$5
python <<-EOF
import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
import pickle
import math
import time
import os
# 
# 
# 
# 
g_ampa_beta = os.environ.get("G_AMPA_BETA")
sub_population = os.environ.get("SUB_POP")
pop = os.environ.get("POP")
root_suff = os.environ.get("ROOT")
root_c = os.environ.get("ROOT_C")
root = "/" + root_c +"/scratch/" + root_suff
root_fig = "/" + root_c +"/FIGURE/" + root_suff
g_ampa = int(g_ampa_beta.split("_")[0])
beta = int(g_ampa_beta.split("_")[1])
# 
# 
# 
# 
# 
# 
# 
signals=["lfp","firing","v"]
bin_size = {"lfp":0.5,"firing":0.5,"v":0.5}
# 
# 
# 
# 
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
n_trial = 500
# 
# 
# 
# 
def Spectrogram(signal,taper,shift,sampling_rate):
    #overlap in percent
    N = taper.shape[0]
    M = 1+ (signal.shape[0] - N)//shift
    lfp = np.empty((M,N),dtype=float)
    for i in range(M):
        lfp[i,:] = signal[i*shift:i*shift+N]
    lfp = lfp - np.mean(lfp,axis=1).reshape(M,-1)
    lfp = np.fft.rfft(lfp*taper)
    spec = np.abs(lfp)**2/N/sampling_rate
    spec[:,1:] *= 2
    return spec.T

time_window = 2000 #size in ms
overlap  = 90
sampling_rate = 1000/dt   #sampling rate in Hz

taper_size = int(time_window/dt) #size in number of data_lfppoint
shift = (100 - overlap)*taper_size//100
hann = (np.sin(np.pi*np.arange(taper_size)/taper_size))**2
frequency = (sampling_rate*np.arange(taper_size)/taper_size)
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
data={}
for signal in signals:
    if signal == "firing":
        data[signal] = 1000 * np.load(root+"{}_Population_{}_{}_{}.npy".format(pop,signal,g_ampa_beta,sub_population))
    else:
        data[signal] = np.load(root+"{}_Population_{}_{}_{}.npy".format(pop,signal,g_ampa_beta,sub_population))
# 
# 
# 
# 
DATA_PLOT = {}
spec = np.empty((n_trial,len(frequency[:201])))
np.random.seed(10000*beta+1000*g_ampa+100*beta)
random_trial = np.random.randint(n_trial)
for signal in signals:
    DATA_PLOT[signal] = {}
    DATA_PLOT[signal]["SAMPEL"] = data[signal][random_trial,:n]
    if signal == "lfp":
        bins = np.arange(data[signal][:,:n].min(),data[signal][:,:n].max()+bin_size[signal],bin_size[signal])
        DATA_PLOT[signal]["BINS"] = bins[:-1]
        for i in range(n_trial):
            a=Spectrogram(data[signal][i,:n],hann,shift,sampling_rate)
            ab=a.sum(axis=1)
            spec[i] = ab[0:201]
        DATA_PLOT[signal]["HIST"] = np.histogram(data[signal][:,:n].ravel(),bins=bins)[0] * 100/data[signal][:,:n].size
        DATA_PLOT[signal]["FREQUENCY"] = frequency[:201]
        DATA_PLOT[signal]["FFT_M"] = np.log10(spec).mean(axis=0)
        DATA_PLOT[signal]["FFT_SD"] = np.log10(spec).std(axis=0)
        ah=np.log10(spec[:,np.argwhere(frequency>=30)[0][0]:].sum(axis=1)/spec[:,:np.argwhere(frequency>4)[0][0]].sum(axis=1))
        fft_ratio_bins = np.arange(ah.min(),ah.max()+0.05,0.05)
        aah=np.histogram(ah,bins=fft_ratio_bins)[0]
        DATA_PLOT[signal]["FFT_RATIO_BINS"] = fft_ratio_bins[:-1]
        DATA_PLOT[signal]["FFT_RATIO"] = aah/aah.sum()
    else:
        if g_ampa !=1:
            bins = np.arange(data[signal][:,:n].min(),data[signal][:,:n].max()+bin_size[signal],bin_size[signal])
            DATA_PLOT[signal]["BINS"] = bins[:-1]
            DATA_PLOT[signal]["HIST"] = np.histogram(data[signal][:,:n].ravel(),bins=bins)[0] * 100/data[signal][:,:n].size
        else:
            bins = np.arange(data["lfp"][:,:n].min(),data["lfp"][:,:n].max()+bin_size["lfp"],bin_size["lfp"])
            HIST = np.histogram(data["lfp"][:,:n].ravel(),bins=bins)[0]
            ind1,ind2  = FIND_PEAK(HIST,sub_population)
            threshold = (bins[ind1]+bins[ind2])/2
            mask = data["lfp"][:,:n]>=threshold
            signal_up = (data[signal][:,:n][mask]).ravel()
            signal_down = (data[signal][:,:n][~mask]).ravel()
            
            bins_1 = np.arange(signal_up.min(),signal_up.max(),bin_size[signal])
            DATA_PLOT[signal]["UP_BINS"] = bins_1[:-1]
            DATA_PLOT[signal]["UP_HIST"] = np.histogram(signal_up,bins=bins_1)[0] * 100/data["lfp"][:,:n].size
            bins_1 = np.arange(signal_down.min(),signal_down.max(),bin_size[signal])
            DATA_PLOT[signal]["DOWN_BINS"] = bins_1[:-1]
            DATA_PLOT[signal]["DOWN_HIST"] = np.histogram(signal_down,bins=bins_1)[0] * 100/data["lfp"][:,:n].size
# 
# 
# 
# 
# 
# 
# 
with open(root+'PRESTIMULUS_PREPARATION_DATA_{}_{}_{}.pickle'.format(g_ampa_beta,pop,sub_population), 'wb') as handle:
    pickle.dump(DATA_PLOT, handle, protocol=pickle.HIGHEST_PROTOCOL)
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
EOF
