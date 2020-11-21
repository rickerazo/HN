#get_data.py

#This is a repository to store functions that capture data from the database

from scipy.signal import savgol_filter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
# from collections import OrderedDict
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
# from sklearn.metrics import silhouette_score
# from scipy.stats import linregress
# from scipy.stats.stats import pearsonr
# from mpl_toolkits import mplot3d
# from scipy.stats import ttest_ind
# from statsmodels.sandbox.stats.multicomp import multipletests 


def amplitude_oscillations_sodium(list10, ax10, ax20,plist, fh, minimum_spikes, colores, countr):
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

    CoV = np.array([])

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
            ##################
            if len(plist)>=4:
                ax10.plot([IpumpMax, IpumpMax], [np.mean(burst_duration)-np.std(burst_duration), np.mean(burst_duration)+np.std(burst_duration)], color = colores[countr-1], alpha=0.5, linewidth=3)
                ax20.plot([IpumpMax, IpumpMax], [np.mean(interburst_interval)-np.std(interburst_interval), np.mean(interburst_interval)+np.std(interburst_interval)], color = colores[countr-1], alpha=0.5, linewidth=3)

            ##############
            
            gp_vec = np.append(gp_vec, gp)
            IpumpMax_vec = np.append(IpumpMax_vec, IpumpMax)

            bd_vec_std = np.append(bd_vec_std, np.std(burst_duration))
            bd_vec = np.append(bd_vec, np.mean(burst_duration))
            ibi_vec_std = np.append(ibi_vec_std, np.std(interburst_interval))
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

            coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
            CoV = np.append(CoV, coefficient_of_variation)

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
            ##################
            if len(plist)>=4:
                ax10.plot([IpumpMax, IpumpMax], [np.mean(burst_duration)-np.std(burst_duration), np.mean(burst_duration)+np.std(burst_duration)], color = colores[countr-1], alpha=0.5, linewidth=3)
                ax20.plot([IpumpMax, IpumpMax], [np.mean(interburst_interval)-np.std(interburst_interval), np.mean(interburst_interval)+np.std(interburst_interval)], color = colores[countr-1], alpha=0.5, linewidth=3)

            ##############

            gp_vec = np.append(gp_vec, gp)
            IpumpMax_vec = np.append(IpumpMax_vec, IpumpMax)

            bd_vec_std = np.append(bd_vec_std, np.std(burst_duration))
            bd_vec = np.append(bd_vec, np.mean(burst_duration))
            ibi_vec_std = np.append(ibi_vec_std, np.std(interburst_interval))
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

            coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
            CoV = np.append(CoV, coefficient_of_variation)

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


    return gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion, Nai_peaks,Nai_troughs,fID, CoV
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

def experimental_variables(control, protocol, coupling, list10, list1, j, fh):
    # print(list1, j)
    if control=='I':
        IpumpMax = list10[-5]
        IpumpMax = float(IpumpMax)
        control_param = IpumpMax

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
    
    if control=='G':
        gp = list10[-5]
        gp = float(gp)
        control_param=gp
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


    return IpumpMax, gp, x_all, time,V, cytNa, interspike_tolerance, burst_duration, interburst_interval, cycle_period, interspike_tolerance

def experimental_variables_classic(control, protocol, coupling, list10, list1, j, fh):

    if control=='I':
        IpumpMax = list10[-5]
        IpumpMax = float(IpumpMax)
        control_param = IpumpMax

        gp = float(list1[j])
        x_all = int(gp)
        time = np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')            
        # Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        mean_spike_f = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
    
    if control=='G':
        gp = list10[-5]
        gp = float(gp)
        control_param=gp

        burst_mean_v = np.array([])
        ibi_min_v = np.array([])

        IpumpMax = list1[j]
        IpumpMax = float(IpumpMax)
        x_all = int(IpumpMax*10)
        time = np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        # Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
        mean_spike_f= np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')

    return IpumpMax, gp, x_all, time,V, cytNa, interspike_tolerance, burst_duration, interburst_interval, cycle_period, mean_spike_f    