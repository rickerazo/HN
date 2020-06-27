####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps
# BURST analysis -> discover if there is a bifurcation of some kind in the data

############# experiment by experiment - manual clustering
#   Using mean spike F, discriminate between experiments in which HN neurons went into high-spike frequency bursting mode
#   and which experiments did not
#   IMPORTANT: The reason to use mean spike frequency is because spiking frequency is well-defined and easy to understand
#   Also important: Vm excursion and spiking frequency are strongly positively correlated.
#   The goals of this script are:
#       1. analyze data prep-by-prep
#          1.1 identifiy what characteristics of spiking frequency separate between hi/lo frequency bursting modes
#          1.2 use these characteristics to cluster data
#       2. apply clustering
#           2.1 all experiments
#           2.2 experiment-by-experiment
#
#       output saved to data_analyses/ebe/*png
#
#       EXPERIMENT-BY-EXPERIMENT
#
#
####  Imported data: only if 3 successful repeated experiments
#
###    - INPUT: curated experiment list
###     - Novel: running mean data filter (Savitzky-Golay) to clean Nai and detect peaks
############################################################################

#### libraries necessary for code
from scipy.signal import savgol_filter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import OrderedDict
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import linregress
from scipy.stats.stats import pearsonr
from mpl_toolkits import mplot3d

import pandas as pd
import os
import sys

cmaps = OrderedDict()

##############################################################################
####	graphix stuff
mpl.rcParams['axes.linewidth']=10


markers = ['1', '2', 'x0', '1','4','3','p','+','h','|','.', '2', 'x0', '1','4','3','p','+','h','|']
# colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia']
colores = ['purple','lawngreen','gold','steelblue','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','firebrick','fuchsia']


font = {'weight' : 'bold',
        'size'   : 50}
plt.rc('font', **font)
plt.rcParams['agg.path.chunksize'] = 10000

plt.close('all')

#### data analysis parameters
## coefficient of variance standard: reduce alpha
coef_std1 = 0.25
## coefficient of variane absolute standard: do not plot
covar = 0.3
length_standard = 3

minimum_spikes = 2

##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
experiment_list = np.array([
    18902004,
    18914000,
    18929002,
    19522001,
    19529000,
    19603000,
    19607000,       #DCC = 0.0
    19607001,       #DCC = -0.05
    19614000,       #DCC = -0.1
    19625000,       #DCC = 0         ######### only I and Vm were saved by pClamp protocol
    19626000,       #DCC = -0.1     ######### only I and Vm were saved by pClamp protocol
    19805000,
    19809001,     ######### 
    19826000,
    19918001,
    19918002,
    19919000,
    19930000,
    19930001,
    '19n05001',
    '19n06000',
    '19d12000',
    '19d12001',
    20204006,
    20205003,
    20210002,
    20225002,
    20226001,
    20303001,
    20310000,
    20311002,     ################
    20317002
	])
# np.save('experiment_list',experiment_list)

def coefficient(coefficient_of_variation, covar,coef_std1):
	if coefficient_of_variation < covar:
		if coefficient_of_variation>coef_std1:
			alpha2 = 0.6
			fst1 = 'none'
		else:
			alpha2 = 1
			fst1 = 'full'
        	# elif coefficient_of_variation >= covar and coefficient_of_variation < 0.5:
        	# 	alpha2 = 0.2
        	# 	fst1 = 'none'
	else:
		alpha2 = 0.001
		fst1 = 'none'
	return alpha2, fst1

def normalize_Nai(list10):
    control = list10[-1]
    protocol = list10[-4]
    coupling = list10[-3]
    list1 = list10[0:-5]

    if control=='I':
        IpumpMax = list10[-5]
        IpumpMax = float(IpumpMax)

        control_param = IpumpMax
        for j in range(0,len(list1)):
            
            burst_mean_v = np.array([])
            ibi_min_v = np.array([])

            gp = float(list1[j])
            x_all = int(gp)
            time = np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')


    if control=='G':
        gp = list10[-5]
        gp = float(gp)
        control_param=gp

        for j in range(0,len(list1)):
            burst_mean_v = np.array([])
            ibi_min_v = np.array([])

            IpumpMax = list1[j]
            IpumpMax = float(IpumpMax)
            x_all = int(IpumpMax*10)
            time = np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')


    # norm_cytNa = savgol_filter(cytNa,20001,2)
    norm_cytNa = savgol_filter(cytNa,51,2)
    peaks, locs = find_peaks(norm_cytNa, prominence=np.std(norm_cytNa))
    min_peaks, min_locs = find_peaks(-norm_cytNa, prominence=np.std(norm_cytNa))

    cytNa_excursion = np.mean(cytNa[peaks]-cytNa[min_peaks])

    return time, cytNa,norm_cytNa, V, cytNa[peaks], cytNa[min_peaks], time[peaks], time[min_peaks], cytNa_excursion
    # return time, cytNa,norm_cytNa, V, norm_cytNa[peaks], norm_cytNa[min_peaks], time[peaks], time[min_peaks]
    # #######################################  cytNa filter
    # fig2 = plt.figure(2,figsize=(20,20))
    # ax2= fig2.add_subplot(111)
    # time, cytNa,norm_cytNa, V, cytNa_peaks, cytNa_min_peaks, time_peaks, time_min_peaks, cytNa_excursion = normalize_Nai(list10)
    # # time, cytNa,norm_cytNa, V, cytNa_peaks, cytNa_min_peaks, time_peaks, time_min_peaks = normalize_Nai(list10)


    # ax2.plot(time,cytNa)
    # ax2.plot(time,norm_cytNa)

    # # ax2.plot(time[peaks], norm_cytNa[peaks],'*',markersize=15)
    # # ax2.plot(time[inv_peaks], norm_cytNa[inv_peaks],'*',markersize=15)

    # ax2.plot(time_peaks, cytNa_peaks,'*',markersize=25)
    # ax2.plot(time_min_peaks, cytNa_min_peaks,'*',markersize=25)

def amplitude_oscillations_sodium(list10):
    control = list10[-1]
    protocol = list10[-4]
    coupling = list10[-3]
    list1 = list10[0:-5]

    cytNa_mean = np.array([])
    cytNa_std = np.array([])
    cytNa_max = np.array([])
    cytNa_min = np.array([])

    gp_vec = np.array([])
    IpumpMax_vec = np.array([])

    bd_vec_std = np.array([])
    bd_vec = np.array([])
    ibi_vec_std = np.array([])
    ibi_vec = np.array([])
    hz_vec = np.array([])
    period_vec = np.array([])

    ct_vmean = np.array([])
    ct_vmin = np.array([])
    burst_vm_excursion = np.array([])
    cytNa_excursion = np.array([])

    Nai_peaks = np.array([])
    Nai_troughs = np.array([])

    fID = np.array([])

    if control=='I':
        IpumpMax = list10[-5]
        IpumpMax = float(IpumpMax)

        control_param = IpumpMax
        for j in range(0,len(list1)):
            
            burst_mean_v = np.array([])
            ibi_min_v = np.array([])

            gp = float(list1[j])
            x_all = int(gp)
            time = np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')            
            Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')

            cytNa_mean = np.append(cytNa_mean, np.mean(cytNa))
            cytNa_std = np.append(cytNa_std, np.std(cytNa))
            cytNa_max = np.append(cytNa_max, np.max(cytNa))
            cytNa_min = np.append(cytNa_min, np.min(cytNa))
            
            gp_vec = np.append(gp_vec, gp)
            IpumpMax_vec = np.append(IpumpMax_vec, IpumpMax)

            bd_vec_std = np.append(bd_vec_std, np.std(burst_duration))
            bd_vec = np.append(bd_vec, np.mean(burst_duration))
            ibi_vec_std = np.append(ibi_vec_std, np.std(burst_duration))
            ibi_vec = np.append(ibi_vec, np.mean(interburst_interval))
            hz_vec = np.append(hz_vec, np.mean(Hz))
            period_vec = np.append(period_vec, np.mean(cycle_period))


                        # 1. spike detection
            Vthreshold = -35
            spikes, vpks = find_peaks(V, height=Vthreshold, prominence=5)
            spike_times = time[spikes]
            spike_volts = V[spikes]
            spike_lag = np.diff(spike_times)

            # last spike in bursts -> mechanism: identify spiking frequency, then when the interspike interval is greater than the mean spiking frequency + 4 std deviations
            # interspike interval tolerance for burst discrimination:

            p1 = np.nonzero(spike_lag>interspike_tolerance)
            p1 = p1[0]
            last_spike = np.append(p1,np.size(spike_times)-1)

            # first spike in bursts
            first_spike = 0
            first_spike = np.append(first_spike,p1+1)

            events1 = np.nonzero((last_spike-first_spike>=minimum_spikes)) #minimum five spikes to be considered a burst
            first_spike = first_spike[events1]
            last_spike = last_spike[events1]
            isi = np.mean(np.diff(spike_times))
            # Burst_volt 
            for j in range(0,len(first_spike)): # a burst is the time interval between first spike and last spike
                t_ini = spike_times[first_spike[j]]
                t_end = spike_times[last_spike[j]]
                
                q0=np.nonzero(time<=t_end)
                q1=np.nonzero(time[q0]>=t_ini)

                # isolate t and v within one burst
                v1 = V[q1]
                t1 = time[q1]

                burst_mean_v = np.append(burst_mean_v,np.mean(v1))

            for j in range(0,len(first_spike)-1):
                t_ini = spike_times[last_spike[j]]
                t_end= spike_times[first_spike[j+1]]

                q0=np.nonzero(time<=t_end)
                q1=np.nonzero(time[q0]>=t_ini)

               # isolate t and v within one burst
                v1 = V[q1]
                t1 = time[q1]
                ibi_min_v = np.append(ibi_min_v, np.min(v1))

            ct_vmean = np.append(ct_vmean, np.mean(burst_mean_v))
            ct_vmin = np.append(ct_vmin, np.mean(ibi_min_v))

            fID = np.append(fID,fh[0:-10])
            burst_vm_excursion = np.append(burst_vm_excursion, np.mean(ct_vmean - ct_vmin))

            norm_cytNa = savgol_filter(cytNa,51,2)
            peaks, locs = find_peaks(norm_cytNa, prominence=np.std(norm_cytNa))
            troughs, trough_locs = find_peaks(-norm_cytNa, prominence=np.std(norm_cytNa))

            p1 = len(cytNa[peaks])
            p2 = len(cytNa[troughs])
            p3 = np.array([p1,p2])
            p4 = np.min(p3)
            
            Nai_p = cytNa[peaks]
            Nai_t = cytNa[troughs]
            # cytNa_excursion = 0
            if p4 == 0:
                cytNa_excursion = np.append(cytNa_excursion, 2*np.std(cytNa))

                Nai_peaks = np.append(Nai_peaks, np.mean(cytNa)+np.std(cytNa))
                Nai_troughs = np.append(Nai_troughs, np.mean(cytNa)-np.std(cytNa))

            else:
                cytNa_excursion = np.append(cytNa_excursion, np.mean(Nai_p[0:p4] - Nai_t[0:p4]))

                Nai_peaks = np.append(Nai_peaks, np.mean(Nai_p[0:p4]))
                Nai_troughs = np.append(Nai_troughs, np.mean(Nai_t[0:p4]))

    if control=='G':
        gp = list10[-5]
        gp = float(gp)
        control_param=gp

        for j in range(0,len(list1)):
            burst_mean_v = np.array([])
            ibi_min_v = np.array([])

            IpumpMax = list1[j]
            IpumpMax = float(IpumpMax)
            x_all = int(IpumpMax*10)
            time = np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')

            cytNa_mean = np.append(cytNa_mean, np.mean(cytNa))
            cytNa_std = np.append(cytNa_std, np.std(cytNa))
            cytNa_max = np.append(cytNa_max, np.max(cytNa))
            cytNa_min = np.append(cytNa_min, np.min(cytNa))

            gp_vec = np.append(gp_vec, gp)
            IpumpMax_vec = np.append(IpumpMax_vec, IpumpMax)

            bd_vec_std = np.append(bd_vec_std, np.std(burst_duration))
            bd_vec = np.append(bd_vec, np.mean(burst_duration))
            ibi_vec_std = np.append(ibi_vec_std, np.std(burst_duration))
            ibi_vec = np.append(ibi_vec, np.mean(interburst_interval))
            hz_vec = np.append(hz_vec, np.mean(Hz))
            period_vec = np.append(period_vec, np.mean(cycle_period))


                        # 1. spike detection
            Vthreshold = -35
            spikes, vpks = find_peaks(V, height=Vthreshold, prominence=5)
            spike_times = time[spikes]
            spike_volts = V[spikes]
            spike_lag = np.diff(spike_times)

            # last spike in bursts -> mechanism: identify spiking frequency, then when the interspike interval is greater than the mean spiking frequency + 4 std deviations
            # interspike interval tolerance for burst discrimination:

            p1 = np.nonzero(spike_lag>interspike_tolerance)
            p1 = p1[0]
            last_spike = np.append(p1,np.size(spike_times)-1)

            # first spike in bursts
            first_spike = 0
            first_spike = np.append(first_spike,p1+1)

            events1 = np.nonzero((last_spike-first_spike>=minimum_spikes)) #minimum five spikes to be considered a burst
            first_spike = first_spike[events1]
            last_spike = last_spike[events1]
            isi = np.mean(np.diff(spike_times))
            # Burst_volt 
            for j in range(0,len(first_spike)): # a burst is the time interval between first spike and last spike
                t_ini = spike_times[first_spike[j]]
                t_end = spike_times[last_spike[j]]
                
                q0=np.nonzero(time<=t_end)
                q1=np.nonzero(time[q0]>=t_ini)

                # isolate t and v within one burst
                v1 = V[q1]
                t1 = time[q1]

                burst_mean_v = np.append(burst_mean_v,np.mean(v1))

            for j in range(0,len(first_spike)-1):
                t_ini = spike_times[last_spike[j]]
                t_end= spike_times[first_spike[j+1]]

                q0=np.nonzero(time<=t_end)
                q1=np.nonzero(time[q0]>=t_ini)

               # isolate t and v within one burst
                v1 = V[q1]
                t1 = time[q1]
                ibi_min_v = np.append(ibi_min_v, np.min(v1))

            ct_vmean = np.append(ct_vmean, np.mean(burst_mean_v))
            ct_vmin = np.append(ct_vmin, np.mean(ibi_min_v))

            fID = np.append(fID,fh[0:-10])
            burst_vm_excursion = np.append(burst_vm_excursion, np.mean(ct_vmean - ct_vmin))



            norm_cytNa = savgol_filter(cytNa,51,2)
            peaks, locs = find_peaks(norm_cytNa, prominence=np.std(norm_cytNa))
            troughs, trough_locs = find_peaks(-norm_cytNa, prominence=np.std(norm_cytNa))

            p1 = len(cytNa[peaks])
            p2 = len(cytNa[troughs])
            p3 = np.array([p1,p2])
            p4 = np.min(p3)
            
            Nai_p = cytNa[peaks]
            Nai_t = cytNa[troughs]
            # cytNa_excursion = 0
            # if fh[0:-10]=='19826000':
            #     cytNa_excursion = np.append(cytNa_excursion, float('NaN'))
            # elif p4 == 0:
            if p4 == 0:
                cytNa_excursion = np.append(cytNa_excursion, 2*np.std(cytNa))

                Nai_peaks = np.append(Nai_peaks, np.mean(cytNa)+np.std(cytNa))
                Nai_troughs = np.append(Nai_troughs, np.mean(cytNa)-np.std(cytNa))

            else:
                cytNa_excursion = np.append(cytNa_excursion, np.mean(Nai_p[0:p4] - Nai_t[0:p4]))

                Nai_peaks = np.append(Nai_peaks, np.mean(Nai_p[0:p4]))
                Nai_troughs = np.append(Nai_troughs, np.mean(Nai_t[0:p4]))


    return gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion, Nai_peaks,Nai_troughs,fID
    # return gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, ct_vmean, ct_vmin, fID

def labelfun(list10):   #function to determine legend labels
    if list10[-1]=='I':
        lab1 = r'g$_P$'
        title1 = r'I$_{pump}^{max}$'
        constant = float(list10[-5])
    elif list10[-1]=='G':
        lab1 = r'I$_{pump}^{max}$'
        title1 = r'g_$P$'
        constant = float(list10[-5])
    return lab1, title1, constant

def cluster_colr(p1):
    if p1>=4:
        colr = 'blue'
    else:
        colr = 'red'

    return colr


############## IMPORTANT : minimum 3 successful experiment repeats to consider the data for subsequent analysis

# compute_raw_data= True
compute_raw_data= False
if compute_raw_data==True:
    print('Parameters from experiments analyzed:')
    id0 = np.array([])

    x0=np.array([])
    y0=np.array([])
    a0=np.array([])
    b0=np.array([])
    c0=np.array([])

    z0=np.array([])
    w0=np.array([])
    u0=np.array([])
    t0=np.array([])


    d0=np.array([])
    # e0=np.array([])
    # f0=np.array([])
    # g0=np.array([])
    # h0=np.array([])

    bd_z = np.array([])
    ibi_z = np.array([])
    exp_ID = np.array([])

    ######################################### Survey the experiments, figure out which experiments have enough successful repeats: appropriate 'n' sample size for subsequent analyses
    list_n = np.array([])
    for i in range(0,len(experiment_list)):

        fileID= str(experiment_list[i])
        fh = fileID+'/'+fileID+'_'
        params = list(np.load(fh+'param.npy'))
        
        for k in range(0,len(params)):# import file
            list20 = params[k]
            han2 = fh+str(list20)+'.npy'
            list10 =  list(np.load(han2))
            # if list10[-1]=='G' and list10[-4]=='nm':
            if list10[-4]=='nm':
                constant_param = float(list10[-5])
                list_n = np.append(list_n, constant_param)


    ###########################################
    for i in range(0,len(experiment_list)):

        fileID= str(experiment_list[i])
        fh = fileID+'/'+fileID+'_'
        params = list(np.load(fh+'param.npy'))

        standard_bd = np.array([])
        standard_bd_std = np.array([])

        standard_ibi = np.array([])
        standard_ibi_std = np.array([])

        for k in range(0,len(params)):# import file
            list20 = params[k]
            han2 = fh+str(list20)+'.npy'
            list10 =  list(np.load(han2))

            n = float(list10[-5])
            nn = np.nonzero(n==list_n)
            nn = nn[0]

            if len(nn)>=1:          # minimum successful repeats (minimum n) to accept the data and analyze it.
            # print(len(nn))
                # if list10[-1]=='G' and list10[-4]=='nm':# and list10[-5]=='2.0':
                if list10[-4]=='nm':# and list10[-5]=='2.0':
                    # if list10[-4]=='nm':
                    # if float(list10[-5])==1.0 or float(list10[-5])==2.0 or float(list10[-5])==5.0 or float(list10[-5])==6.0:

                    print(fileID,list10)
                    lab1, title1, constant = labelfun(list10)

                    ######################## V
                    # gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, ct_vmean, ct_vmin, fID = amplitude_oscillations_sodium(list10)
                    # gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion,fID = amplitude_oscillations_sodium(list10)
                    gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion, Nai_peaks,Nai_troughs,fID = amplitude_oscillations_sodium(list10)
                    ################# 
                    x0 = np.append(x0, gp_vec)
                    y0 = np.append(y0, IpumpMax_vec)

                    a0 = np.append(a0, bd_vec)
                    b0 = np.append(b0, ibi_vec)
                    c0 = np.append(c0, hz_vec)
                    d0 = np.append(d0, period_vec)

                    z0 = np.append(z0, burst_vm_excursion)
                    w0 = np.append(w0, cytNa_excursion)

                    u0 = np.append(u0,Nai_peaks)
                    t0 = np.append(t0, Nai_troughs)

                    exp_ID = np.append(exp_ID, fID)

    ###########################################################################################
    ######### save the variables to make subsequent computations faster
    dir0 = 'data_analyses/ebe/vars/'
    np.save(dir0+'x0.npy',x0)
    np.save(dir0+'y0.npy',y0)

    np.save(dir0+'a0.npy',a0)
    np.save(dir0+'b0.npy',b0)
    np.save(dir0+'c0.npy',c0)
    np.save(dir0+'d0.npy',d0)

    np.save(dir0+'z0.npy',z0)
    np.save(dir0+'w0.npy',w0)

    np.save(dir0+'u0.npy',u0)
    np.save(dir0+'t0.npy',t0)

    np.save(dir0+'exp0.npy',exp_ID)

############################################################
######### load the variables 
############################################################
dir0 = 'data_analyses/ebe/vars/'

x0 = np.load(dir0+'x0.npy')             ### gP
y0 = np.load(dir0+'y0.npy')             ### IpumpMax

a0 = np.load(dir0+'a0.npy')             ### burst duration
b0 = np.load(dir0+'b0.npy')             ### interburst interval
c0 = np.load(dir0+'c0.npy')             ### spike f 
d0 = np.load(dir0+'d0.npy')             ### burst period

z0 = np.load(dir0+'z0.npy')             ### Vm excursion

u0 = np.load(dir0+'u0.npy')             ### Nai peaks
t0 = np.load(dir0+'t0.npy')             ### Nai troughs

w0 = np.load(dir0+'w0.npy')             ### Nai excursion
# w0 = u0-t0                              ### Nai excursion

#############   there is no difference in the data when we compare
#############   (Nai peaks - Nai troughs) : averaged    
#############           vs
#############   average Nai peaks - average Nai troughs

############    therefore it's safe to assume there isnt a difference in Vm when either approach is used

exp_ID = np.load(dir0+'exp0.npy')       ### file ID

###########################################################################################
### normalize the variables : z-scores
raw_c0 = c0
raw_z0 = z0

a0 = (a0-np.mean(a0))/np.std(a0)
b0 = (b0-np.mean(b0))/np.std(b0)

c0 = (c0-np.mean(c0))/np.std(c0)
z0 = (z0-np.mean(z0))/np.std(z0)
w0 = (w0-np.mean(w0))/np.std(w0)

ncluj=2

spikef_threshold = 8
Vm_excursion_threshold = 10


np.save(dir0+'spikef_threshold.npy', spikef_threshold)
np.save(dir0+'Vm_excursion_threshold.npy', Vm_excursion_threshold)
#############################################################################################
###########     EXPERIMENT - BY- EXPERIMENT CLUSTERING      #################################
##################################################################  SPIKE FREQUENCY     #####
cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

for i in range(0,len(x0)):

    if raw_c0[i]> spikef_threshold:
        cluster_ID[i]=1
    else:
        cluster_ID[i]=0

# plt.ion()
# plt.show()

h1 = plt.figure(figsize=(20,20))
ax1 = h1.add_subplot(111)

p0 = np.nonzero(cluster_ID==0)
p1 = np.nonzero(cluster_ID==1)

# ax1.scatter(x0[p0], y0[p0], c=colores[0], s=150, alpha = 0.5, marker='*')
# ax1.scatter(x0[p1], y0[p1], c=colores[1], s=150, alpha = 0.5, marker='1')

ax1.scatter(x0[p0], y0[p0], c='r', s=1500, alpha = 0.5, marker='2')
ax1.scatter(x0[p1], y0[p1], c='b', s=1500, alpha = 0.5, marker='1')

x1 = x0[p1]
y1 = y0[p1]
c1 = raw_c0[p1]
ex1 = exp_ID[p1]

p2 = np.nonzero(x1<4)

x2 = x1[p2]
y2 = y1[p2]
c2 = c1[p2]
ex2 = ex1[p2]
print('High spike frequency with low gP. Spike f discrimination')
print('x2',x2)
print('y2',y2)
print('c2',c2)
print('ex2',ex2)

h1.savefig('data_analyses/ebe/spikef/all_exp.png')
##  discriminate when a prep always remains in low spike f mode
##      when mean of spike f is below 5 and 
##  discriminate when a prep always remains in high spike f mode
##      mean is above 10 and standard dev is greater than 2

# ############################################################################
# plt.ion()
# plt.show()

# plt.figure()
# ax=plt.axes(projection='3d')
# ax.plot_trisurf(x0,y0,u0-t0, cmap='Spectral')
# ax.scatter(x0,y0,u0-t0)



# ######################################
# plt.figure()
# ax1=plt.axes(projection='3d')

# ax1.plot_trisurf(x0,y0,c0, cmap='Spectral')
# ax1.scatter(x0,y0,c0)
# ax1.set_title('c0')
# # ###############################
# plt.figure()
# ax2=plt.axes(projection='3d')

# ax2.plot_trisurf(x0,y0,z0, cmap='Spectral')
# ax2.scatter(x0,y0,z0)
# ax2.set_title('z0')

# # ###################################
# plt.figure()
# ax3=plt.axes(projection='3d')

# ax3.scatter(x0,y0,w0)
# ax3.plot_trisurf(x0,y0,w0, cmap='Spectral')
# ax3.set_title('w0')


############################################################################


cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

p10 = np.unique(exp_ID)

exp_cluster = np.zeros((1,len(x0)))
exp_cluster=exp_cluster[0]

for i in range(0,len(p10)):
    p2 = np.nonzero(exp_ID==p10[i])
    x2 = x0[p2]
    y2 = y0[p2]
    c2 = raw_c0[p2]

    temp1 = np.zeros((1,len(c2)))
    temp1 = temp1[0]
    for j in range(0,len(c2)):
        # if c2[j]>9:
        if c2[j]>spikef_threshold:
            temp1[j] = 1
        else:
            temp1[j] = 0

    exp_cluster[p2] = temp1



    p0 = np.nonzero(exp_cluster[p2]==0)
    p1 = np.nonzero(exp_cluster[p2]==1)

    ###

    h2 = plt.figure(figsize=(20,20))
    ax2 = h2.add_subplot(111)

    ax2.scatter(x2[p0],y2[p0],c='r',marker='1',s=5000)
    ax2.scatter(x2[p1],y2[p1],c='b',marker='2',s=5000)

    ax2.set_xlim(-1,10.5)
    ax2.set_ylim(-0.1,1.1)

    ax2.tick_params(axis='both', length=20, width=30)

    ax2.set_ylabel(r'I$_{pump}^{max}$ (nA)')
    ax2.set_xlabel(r'$\bar{g}_P$ (nS)')
    ax2.set_title('Experiment ID ='+str(p10[i]))

    h2.savefig('data_analyses/ebe/spikef/'+str(p10[i])+'.png')
    
    plt.close()


#############################################################################################
###########     EXPERIMENT - BY- EXPERIMENT CLUSTERING      #################################
##################################################################  Vm EXCURSION     ########

cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

for i in range(0,len(x0)):

    if raw_z0[i]> Vm_excursion_threshold:
        cluster_ID[i]=1
    else:
        cluster_ID[i]=0

# plt.ion()
# plt.show()

h1 = plt.figure(figsize=(20,20))
ax1 = h1.add_subplot(111)

p0 = np.nonzero(cluster_ID==0)
p1 = np.nonzero(cluster_ID==1)

# ax1.scatter(x0[p0], y0[p0], c=colores[0], s=150, alpha = 0.5, marker='*')
# ax1.scatter(x0[p1], y0[p1], c=colores[1], s=150, alpha = 0.5, marker='1')

ax1.scatter(x0[p0], y0[p0], c='r', s=1750, alpha = 0.5, marker='2')
ax1.scatter(x0[p1], y0[p1], c='b', s=1750, alpha = 0.5, marker='1')

x1 = x0[p1]
y1 = y0[p1]
c1 = raw_c0[p1]
z1 = raw_z0[p1]
ex1 = exp_ID[p1]

p2 = np.nonzero(x1<4)

x2 = x1[p2]
y2 = y1[p2]
c2 = c1[p2]
z2 = z1[p2]
ex2 = ex1[p2]
print('High spike frequency: low gP. Vm excursion discrimination')
print('x2',x2)
print('y2',y2)
print('z2',z2)
print('ex2',ex2)

h1.savefig('data_analyses/ebe/Vm/all_exp.png')
##  discriminate when a prep always remains in low spike f mode
##      when mean of spike f is below 5 and 
##  discriminate when a prep always remains in high spike f mode
##      mean is above 10 and standard dev is greater than 2

# ############################################################################
# plt.ion()
# plt.show()

# plt.figure()
# ax=plt.axes(projection='3d')
# ax.plot_trisurf(x0,y0,u0-t0, cmap='Spectral')
# ax.scatter(x0,y0,u0-t0)



# ######################################
# plt.figure()
# ax1=plt.axes(projection='3d')

# ax1.plot_trisurf(x0,y0,c0, cmap='Spectral')
# ax1.scatter(x0,y0,c0)
# ax1.set_title('c0')
# # ###############################
# plt.figure()
# ax2=plt.axes(projection='3d')

# ax2.plot_trisurf(x0,y0,z0, cmap='Spectral')
# ax2.scatter(x0,y0,z0)
# ax2.set_title('z0')

# # ###################################
# plt.figure()
# ax3=plt.axes(projection='3d')

# ax3.scatter(x0,y0,w0)
# ax3.plot_trisurf(x0,y0,w0, cmap='Spectral')
# ax3.set_title('w0')


############################################################################


cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

p10 = np.unique(exp_ID)

exp_cluster = np.zeros((1,len(x0)))
exp_cluster=exp_cluster[0]

for i in range(0,len(p10)):
    p2 = np.nonzero(exp_ID==p10[i])
    x2 = x0[p2]
    y2 = y0[p2]
    z2 = raw_z0[p2]


    temp1 = np.zeros((1,len(z2)))
    temp1 = temp1[0]
    for j in range(0,len(z2)):
        # if c2[j]>9:
        if z2[j]>Vm_excursion_threshold:
            temp1[j] = 1
        else:
            temp1[j] = 0

    exp_cluster[p2] = temp1



    p0 = np.nonzero(exp_cluster[p2]==0)
    p1 = np.nonzero(exp_cluster[p2]==1)

    ###

    h2 = plt.figure(figsize=(20,20))
    ax2 = h2.add_subplot(111)

    ax2.scatter(x2[p0],y2[p0],c='r',marker='1',s=5000)
    ax2.scatter(x2[p1],y2[p1],c='b',marker='2',s=5000)

    ax2.set_xlim(-1,10.5)
    ax2.set_ylim(-0.1,1.1)

    ax2.tick_params(axis='both', length=20, width=30)

    ax2.set_ylabel(r'I$_{pump}^{max}$ (nA)')
    ax2.set_xlabel(r'$\bar{g}_P$ (nS)')
    ax2.set_title('Experiment ID ='+str(p10[i]))

    h2.savefig('data_analyses/ebe/Vm/'+str(p10[i])+'.png')
    
    plt.close()    


######################################################################################################################################
################################ probably the most robust way to identify when a neuron goes into high spike frequency bursting mode
################################ is to use both spiking frequency and Vm excursion
################################ That's next:
########################################################################### Vm excursion + Spike F 
######################################################################################################################################
# Vm_excursion_threshold = 12
# spikef_threshold = 9

cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

for i in range(0,len(x0)):

    if raw_z0[i]> Vm_excursion_threshold and raw_c0[i] > spikef_threshold:
        cluster_ID[i]=1
    else:
        cluster_ID[i]=0

h1 = plt.figure(figsize=(20,20))
ax1 = h1.add_subplot(111)

p0 = np.nonzero(cluster_ID==0)
p1 = np.nonzero(cluster_ID==1)

ax1.scatter(x0[p0], y0[p0], c='r', s=1750, alpha = 0.5, marker='2')
ax1.scatter(x0[p1], y0[p1], c='b', s=1750, alpha = 0.5, marker='1')

ax1.set_ylabel(r'I$_{pump}^{max}$ (nA)')
ax1.set_xlabel(r'$\bar{g}_P$ (nS)')
ax1.set_title('All experiments')


x1 = x0[p1]
y1 = y0[p1]
c1 = raw_c0[p1]
z1 = raw_z0[p1]
ex1 = exp_ID[p1]

p2 = np.nonzero(x1<4)

x2 = x1[p2]
y2 = y1[p2]
c2 = c1[p2]
z2 = z1[p2]
ex2 = ex1[p2]

print('High spike frequency. low gP. Vm and spike f discrimination')
print('x2',x2)
print('y2',y2)
print('c2',c2)
print('z2',z2)
print('ex2',ex2)

h1.savefig('data_analyses/ebe/Vm_spikef/all_exp.png')

############################################################################


cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

p10 = np.unique(exp_ID)

exp_cluster = np.zeros((1,len(x0)))
exp_cluster=exp_cluster[0]

for i in range(0,len(p10)):
    p2 = np.nonzero(exp_ID==p10[i])
    x2 = x0[p2]
    y2 = y0[p2]
    z2 = raw_z0[p2]
    c2 = raw_c0[p2]


    temp1 = np.zeros((1,len(z2)))
    temp1 = temp1[0]
    for j in range(0,len(z2)):
        # if c2[j]>9:
        if z2[j]>Vm_excursion_threshold and c2[j]> spikef_threshold:
            temp1[j] = 1
        else:
            temp1[j] = 0

    exp_cluster[p2] = temp1



    p0 = np.nonzero(exp_cluster[p2]==0)
    p1 = np.nonzero(exp_cluster[p2]==1)

    ###

    h2 = plt.figure(figsize=(20,20))
    ax2 = h2.add_subplot(111)

    ax2.scatter(x2[p0],y2[p0],c='r',marker='1',s=5000)
    ax2.scatter(x2[p1],y2[p1],c='b',marker='2',s=5000)

    ax2.set_xlim(-1,10.5)
    ax2.set_ylim(-0.1,1.1)

    ax2.tick_params(axis='both', length=20, width=30)

    ax2.set_ylabel(r'I$_{pump}^{max}$ (nA)')
    ax2.set_xlabel(r'$\bar{g}_P$ (nS)')
    ax2.set_title('Experiment ID ='+str(p10[i]))

    h2.savefig('data_analyses/ebe/Vm_spikef/'+str(p10[i])+'.png')
    
    plt.close()    

