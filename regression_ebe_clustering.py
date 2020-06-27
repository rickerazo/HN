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
    # 19625000,       #DCC = 0         ######### only I and Vm were saved by pClamp protocol
    # 19626000,       #DCC = -0.1     ######### only I and Vm were saved by pClamp protocol
    19805000,
    19809001,     #########     outlier
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
    20311002,     ################  outlier
    20317002
	])
# np.save('experiment_list',experiment_list)

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
                    # lab1, title1, constant = labelfun(list10)

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

# spikef_threshold = 8
# Vm_excursion_threshold = 10
spikef_threshold = np.load(dir0+'spikef_threshold.npy')
Vm_excursion_threshold = np.load(dir0+'Vm_excursion_threshold.npy')

################################    Simple      -->"threshold-based"<--    clustering algorithm
cluster_ID = np.zeros((1,len(x0)))
cluster_ID = cluster_ID[0]

for i in range(0,len(x0)):

    if raw_z0[i]> Vm_excursion_threshold and raw_c0[i] > spikef_threshold:
    # if raw_c0[i] > spikef_threshold:
        cluster_ID[i]=1
    else:
        cluster_ID[i]=0

## quality of the clustering:
dv = {'x':raw_c0,'y':raw_z0}
dv_dfs = pd.DataFrame(dv,columns=['x','y'])
Vm_spikef_score = silhouette_score(dv_dfs, cluster_ID, metric='euclidean', sample_size=None, random_state=1)
print('Score='+str(Vm_spikef_score))

### score = 0.5 -> decent!

# silhouette_score(c0.reshape(-1,1), kmeans_spikef.labels_, metric='euclidean', sample_size=None, random_state=i)
######################################################################################################################################
######################################################################################################################################
################################    Based on the previous work from ebe_clustering.py, here we use both Vm excursion and mean spike F
################################    to cluster the experiments
################################    Next is to begin curve fitting:
################################    1. BD vs IBI
################################    2. BD vs spike f
################################    3. BD vs Vm excursion
################################    4. BD vs Nai
################################    5. BD vs T
######################################################################################################################################

### variables of interest are: a0 (BD), b0(IBI), c0(spike F), z0(Vm excursion), w0(Nai excursion), d0(T)
p0 = np.nonzero(cluster_ID==0)
p1 = np.nonzero(cluster_ID==1)

bd_ibi_r0 = linregress(a0[p0],b0[p0])
bd_ibi_r1 = linregress(a0[p1],b0[p1])

bd_spikef_r0 = linregress(a0[p0],c0[p0])
bd_spikef_r1 = linregress(a0[p1],c0[p1])

bd_vm_r0 = linregress(a0[p0],z0[p0])
bd_vm_r1 = linregress(a0[p1],z0[p1])

bd_nai_r0 = linregress(a0[p0],w0[p0])
bd_nai_r1 = linregress(a0[p1],w0[p1])

bd_t_r0 = linregress(a0[p0],d0[p0])
bd_t_r1 = linregress(a0[p1],d0[p1])


print('\n')
### print results
#############################
print('BD + IBI regression')
print('cluster 0 ')
print('r2 = ',bd_ibi_r0.rvalue, ', p= ', bd_ibi_r0.pvalue)
print('cluster 1')
print('r2 = ',bd_ibi_r1.rvalue, ', p= ', bd_ibi_r1.pvalue)
print('\n')

#############################
print('BD + spike f regression')
print('cluster 0 ')
print('r2 = ', bd_spikef_r0.rvalue, ', p= ', bd_spikef_r0.pvalue)
print('cluster 1')
print('r2 = ', bd_spikef_r1.rvalue, ', p= ', bd_spikef_r1.pvalue)
print('\n')

#############################
print('BD + Vm excursion regression')
print('cluster 0 ')
print('r2 = ', bd_vm_r0.rvalue, ', p= ', bd_vm_r0.pvalue)
print('cluster 1')
print('r2 = ', bd_vm_r1.rvalue, ', p= ', bd_vm_r1.pvalue)
print('\n')

#############################
print('BD + Nai excursion regression')
print('cluster 0 ')
print('r2 = ', bd_nai_r0.rvalue, ', p=', bd_nai_r0.pvalue)
print('cluster 1')
print('r2 = ', bd_nai_r1.rvalue, ', p=', bd_nai_r1.pvalue)
print('\n')

#############################
print('BD + T regression')
print('cluster 0 ')
print('r2 = ', bd_t_r0.rvalue, ', p= ', bd_t_r0.pvalue)
print('cluster 1')
print('r2 = ', bd_t_r1.rvalue, ', p= ', bd_t_r1.pvalue)
print('\n')



################# break in data

x6 = x0[p0]
y6 = y0[p0]
a6 = a0[p0]
c6 = c0[p0]
z6 = raw_z0[p0]
ex6 = exp_ID[p0]

p6 = np.nonzero(x6==6)
p7 = np.nonzero(x6==7)

p10 = np.nonzero(x6==10)
###a0

x1 = x0[p1]
y1 = y0[p1]
a1 = a0[p1]
c1 = c0[p1]
z1 = z0[p1]
ex1 = exp_ID[p1]

p4 = np.nonzero(x6==4)
p5 = np.nonzero(x6==5)

# p6 = np.nonzero(x1==6)
p3 = np.nonzero(x1==3)


plt.ion()
plt.show()
