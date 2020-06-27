####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps
# BURST analysis -> discover if there is a bifurcation of some kind in the data

############# gP constant, sweep IpumpMax
# clustered_regression.py
#   The goals of this script are:
#       1. run linear regression analysis in DVs: DB, IBI, T, spike f
#       2. detect qualitative data changes


###    - INPUT: curated experiment list

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
import pandas as pd
import os
import sys

cmaps = OrderedDict()

##############################################################################
####	graphix stuff
mpl.rcParams['axes.linewidth']=10


markers = ['1', '2', 'x0', '1','4','3','p','+','h','|','.', '2', 'x0', '1','4','3','p','+','h','|']
# colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia']
colores = ['forestgreen','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia']


font = {'weight' : 'bold',
        'size'   : 50}
plt.rc('font', **font)
plt.rcParams['agg.path.chunksize'] = 10000

mk1 = 20
mk2 = 15
mk3 = 3
mk10 = 25
lw1= 6
alpha0 = 1
plt.close('all')


# coefficient of variance standard: reduce alpha
coef_std1 = 0.25
# coefficient of variane absolute standard: do not plot
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
    19625000,       #DCC = 0
    19626000,       #DCC = -0.1
    19805000,
    # 19809001,
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
    # 20311002,
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
            if fh[0:-10]=='19826000':
                cytNa_excursion = np.append(cytNa_excursion, float('NaN'))
            elif p4 == 0:
                cytNa_excursion = np.append(cytNa_excursion, np.mean(cytNa))
            else:
                cytNa_excursion = np.append(cytNa_excursion, np.mean(Nai_p[0:p4] - Nai_t[0:p4]))

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
            if fh[0:-10]=='19826000':
                cytNa_excursion = np.append(cytNa_excursion, float('NaN'))
            elif p4 == 0:
                cytNa_excursion = np.append(cytNa_excursion, np.mean(cytNa))
            else:
                cytNa_excursion = np.append(cytNa_excursion, np.mean(Nai_p[0:p4] - Nai_t[0:p4]))

    return gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion,fID
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

# time, cytNa,norm_cytNa, V, cytNa_peaks, cytNa_min_peaks, time_peaks, time_min_peaks, cytNa_excursion = normalize_Nai(list10)
# # time, cytNa,norm_cytNa, V, cytNa_peaks, cytNa_min_peaks, time_peaks, time_min_peaks = normalize_Nai(list10)


# ax2.plot(time,cytNa)
# ax2.plot(time,norm_cytNa)

# # ax2.plot(time[peaks], norm_cytNa[peaks],'*',markersize=15)
# # ax2.plot(time[inv_peaks], norm_cytNa[inv_peaks],'*',markersize=15)

# ax2.plot(time_peaks, cytNa_peaks,'*',markersize=25)
# ax2.plot(time_min_peaks, cytNa_min_peaks,'*',markersize=25)



print('Parameters from experiments analyzed:')
id0 = np.array([])

x0=np.array([])
y0=np.array([])
a0=np.array([])
b0=np.array([])
c0=np.array([])

z0=np.array([])
w0=np.array([])

d0=np.array([])
e0=np.array([])
f0=np.array([])
g0=np.array([])
h0=np.array([])

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

        if len(nn)>=0:          # minimum successful repeats (minimum n) to accept the data and analyze it.
        # print(len(nn))
            # if list10[-1]=='G' and list10[-4]=='nm':# and list10[-5]=='2.0':
            # if list10[-4]=='nm':# and list10[-5]=='2.0':
            if list10[-5]=='1.0':# and list10[-5]=='2.0':
                # if list10[-4]=='nm':
                # if float(list10[-5])==1.0 or float(list10[-5])==2.0 or float(list10[-5])==5.0 or float(list10[-5])==6.0:

                print(fileID,list10[-5],list10[-1],list10[0:-5])
                lab1, title1, constant = labelfun(list10)

                ######################## V
                # gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, ct_vmean, ct_vmin, fID = amplitude_oscillations_sodium(list10)
                gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion,fID = amplitude_oscillations_sodium(list10)
                ################# 
                x0 = np.append(x0, gp_vec)
                y0 = np.append(y0, IpumpMax_vec)

                a0 = np.append(a0, bd_vec)
                b0 = np.append(b0, ibi_vec)
                c0 = np.append(c0, hz_vec)
                d0 = np.append(d0, period_vec)

                ################ 
                # z0 = np.append(z0, ct_vmean)
                # w0 = np.append(w0, ct_vmin)
                z0 = np.append(z0, burst_vm_excursion)

                cytNa_min[np.nonzero(cytNa_max>=25)]=np.mean(cytNa_min)
                cytNa_max[np.nonzero(cytNa_max>=25)]=np.mean(cytNa_max)

                e0 = np.append(e0, cytNa_max)
                f0 = np.append(f0, cytNa_std)
                g0 = np.append(g0, cytNa_mean)
                h0 = np.append(h0, cytNa_min)

                exp_ID = np.append(exp_ID, fID)


#################################### How many experiments are there for each parameter combination?
print('\n')

x1 = np.unique(x0)
y1 = np.unique(y0)
for i in range(0,len(x1)):
    x2= x0[np.nonzero(x0==x1[i])]

    for j in range(0,len(y1)):
        y2 = y0[np.nonzero(y0==y1[j])]


        print(x1[i],'   ',len(x2),'   ',y1[j],'   ',len(y2))

print('\n')

for j in range(0,len(y1)):
    y2 = y0[np.nonzero(y0==y1[j])]
    for i in range(0,len(x1)):
        x2= x0[np.nonzero(x0==x1[i])]

        print(y1[j],'   ',len(y2),'   ',x1[i],'   ',len(x2))

# plt.ion()
# plt.show()

# ####
# h1 = plt.figure(1,figsize=(30,30))
# ax = h1.add_subplot(111)
# ax.scatter(y0,a0,s=100)