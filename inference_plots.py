# inference_plots.py

#This script will store all functions associated with plotting

from scipy.signal import savgol_filter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
# from collections import OrderedDict
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import linregress
from scipy.stats.stats import pearsonr
from mpl_toolkits import mplot3d
from scipy.stats import ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests 

import custom_analysis
import get_data

def robust_stats(list10, ax11, ax21,plist, fh, minimum_spikes, colores, countr):
    control = list10[-1]
    protocol = list10[-4]
    coupling = list10[-3]
    list1 = list10[0:-5]

    aggregate_bd = np.array([])
    aggregate_ibi = np.array([])

    gp_vec = np.array([])
    IpumpMax_vec = np.array([])

    bd_mean = np.array([])
    bd_median = np.array([])

    ibi_mean = np.array([])
    ibi_median = np.array([])

    period_mean = np.array([])
    period_median = np.array([])

    Nai_ex = np.array([])
    # cytNa_std = np.array([])
    # cytNa_max = np.array([])
    # cytNa_min = np.array([])

    bd_first_quartile = np.array([])
    bd_third_quartile = np.array([])
    bd_median = np.array([])

    ibi_first_quartile = np.array([])
    ibi_third_quartile = np.array([])
    ibi_median = np.array([])

    mean_f_vec = np.array([])
    median_f_vec = np.array([])

    Vm_excursion_vec = np.array([])

    fID = np.array([])

    cluster_vec = np.array([])


    if len(plist)>=4:    

        for j in range(0,len(list1)):
            IpumpMax, gp, x_all, time,V, cytNa, interspike_tolerance, burst_duration, interburst_interval, cycle_period, interspike_tolerance = get_data.experimental_variables(control, protocol, coupling, list10, list1, j, fh)
            aggregate_bd = np.append(aggregate_bd, burst_duration)
            aggregate_ibi = np.append(aggregate_ibi, interburst_interval)


            ###### Nai
            cytNa_excursion = custom_analysis.normalize_Nai(list10, fh)
            Nai_ex = np.append(Nai_ex, cytNa_excursion)
            # cytNa_mean = np.append(cytNa_mean, np.mean(cytNa))
            # cytNa_std = np.append(cytNa_std, np.std(cytNa))
            # cytNa_max = np.append(cytNa_max, np.max(cytNa))
            # cytNa_min = np.append(cytNa_min, np.min(cytNa))
            ########### instant spiking frequency
            mean_v = V
            time1 = time - time[0]	

            mean_f, median_f = custom_analysis.instant_spiking_frequency(mean_v, time1, cytNa, interspike_tolerance, x_all)

            ################## Sweep parameters

            gp_vec = np.append(gp_vec, gp)
            IpumpMax_vec = np.append(IpumpMax_vec, IpumpMax)
            fID = np.append(fID, fh)

            ################## Rhythmic characteristics
            bd_first_quartile = np.append(bd_first_quartile, np.percentile(burst_duration, 25, interpolation='midpoint'))
            bd_third_quartile = np.append(bd_third_quartile, np.percentile(burst_duration, 75, interpolation='midpoint'))

            bd_median = np.append(bd_median, np.median(burst_duration))

            ibi_first_quartile = np.append(ibi_first_quartile, np.percentile(interburst_interval, 25, interpolation='midpoint'))
            ibi_third_quartile = np.append(ibi_third_quartile, np.percentile(interburst_interval, 75, interpolation='midpoint'))

            ibi_median = np.append(ibi_median, np.median(interburst_interval))

            mean_f_vec = np.append(mean_f_vec, mean_f)
            median_f_vec = np.append(median_f_vec, median_f)

            ################# Vm excursion
            # print(IpumpMax, gp)
            # x1, y1, y3, y4 = custom_analysis.flat_spikes(V, time, -30, interspike_tolerance, burst_duration, cytNa)
            mean_v, time, y1 = custom_analysis.average_voltage(V, time, -30, interspike_tolerance, burst_duration, cytNa)
            quiescent_potential = np.min(mean_v)
            bursting_potential = np.max(mean_v)
            Vm_excursion = np.abs(bursting_potential - quiescent_potential)

            # plt.figure()
            # plt.plot(time, y1)
            # plt.plot(time, mean_v)

            # print(IpumpMax, bursting_potential, quiescent_potential,Vm_excursion)
            # quiescent_potential, bursting_potential = custom_analysis.Voltage_excursion(V, time, -25, interspike_tolerance, fh)
            # Vm_excursion = np.abs(bursting_potential[0:-1] - quiescent_potential)
            
            Vm_excursion_vec = np.append(Vm_excursion_vec, np.median(Vm_excursion))

            # # # Clustering
            cluster_ID = custom_analysis.two_threshold_clustering(median_f, Vm_excursion, plist, bd_median, ibi_median)
            cluster_vec = np.append(cluster_vec, cluster_ID)

            ########## Robust plots
            ax11, ax21 = median_interquartile_ranges(cluster_ID, IpumpMax, burst_duration, interburst_interval,colores, countr, plist, ax11, ax21)
    else:
        mean_f = 0
        median_f = 0
        Vm_excursion = 0
        cluster_vec = 0

        cluster_vec = np.append(cluster_vec, 0)
        Vm_excursion_vec = np.append(Vm_excursion_vec, Vm_excursion)
        median_f_vec = np.append(median_f_vec, median_f)
        mean_f_vec = np.append(mean_f_vec, mean_f)
        aggregate_bd= np.append(aggregate_bd, 0) 
        aggregate_ibi = np.append(aggregate_ibi, 0)

    # return gp_vec, IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_median, bd_first_quartile, bd_third_quartile, ibi_median, ibi_first_quartile, ibi_third_quartile, mean_f_vec, median_f_vec, Vm_excursion_vec, fID, ax11, ax21
    # print(aggregate_bd)
    if aggregate_bd.all()!=0:
        normalization_constant_BD = np.max(aggregate_bd[0:7])
        normalization_constant_IBI = np.max(aggregate_ibi[0:7])
    else:
        normalization_constant_BD = np.max(aggregate_bd)
        normalization_constant_IBI = np.max(aggregate_ibi)

    # return gp_vec, IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_median, bd_first_quartile, bd_third_quartile, ibi_median, ibi_first_quartile, ibi_third_quartile, mean_f, median_f, Vm_excursion_vec, cluster_vec ,fID, ax11, ax21, normalization_constant_BD, normalization_constant_IBI
    return gp_vec, IpumpMax_vec, Nai_ex, bd_median, bd_first_quartile, bd_third_quartile, ibi_median, ibi_first_quartile, ibi_third_quartile, mean_f_vec, median_f_vec, Vm_excursion_vec, cluster_vec ,fID, ax11, ax21, normalization_constant_BD, normalization_constant_IBI

def classic_stats(list10, ax10, ax20,plist, fh, minimum_spikes, colores, countr):
    control = list10[-1]
    protocol = list10[-4]
    coupling = list10[-3]
    list1 = list10[0:-5]

    gp_vec = np.array([])
    IpumpMax_vec = np.array([])

    bd_mean = np.array([])
    bd_median = np.array([])

    ibi_mean = np.array([])
    ibi_median = np.array([])

    period_mean = np.array([])
    period_median = np.array([])

    cytNa_mean = np.array([])
    cytNa_std = np.array([])
    cytNa_max = np.array([])
    cytNa_min = np.array([])

    bd_vec = np.array([])
    bd_vec_std = np.array([])
    ibi_vec = np.array([])
    ibi_vec_std = np.array([])
    hz_vec = np.array([])
    period_vec = np.array([])

    burst_vm_excursion = np.array([])
    burst_mean_v = np.array([])

    ibi_min_v = np.array([])

    ct_vmean = np.array([])
    ct_vmin = np.array([])


    cytNa_excursion = np.array([])
    Nai_peaks = np.array([])
    Nai_troughs = np.array([])

    CoV = np.array([])

    fID = np.array([])

    if len(plist)>=4:    

        for j in range(0,len(list1)):
            IpumpMax, gp, x_all, time,V, cytNa, interspike_tolerance, burst_duration, interburst_interval, cycle_period, Hz = get_data.experimental_variables_classic(control, protocol, coupling, list10, list1, j, fh)


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


def classic_trend(IpumpMax_vec, bd_vec, bd_vec_std,ibi_vec, ibi_vec_std, colores, countr, km, ax10, ax20, plist):
	if len(plist)>=4:
		# print('\n',plist)
		# print(len(IpumpMax_vec), len(bd_vec), '\n')
		ax10.plot(IpumpMax_vec, bd_vec, 'd',color = colores[countr-1], alpha=0.5, markersize=km)
		ax10.plot(IpumpMax_vec, bd_vec+bd_vec_std, '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
		ax10.plot(IpumpMax_vec, bd_vec-bd_vec_std, '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
		ax10.plot(IpumpMax_vec, bd_vec ,color = colores[countr-1], linewidth=3)

		ax20.plot(IpumpMax_vec, ibi_vec, 'd',color = colores[countr-1], alpha=0.5, markersize=km)
		ax20.plot(IpumpMax_vec, ibi_vec+ibi_vec_std, '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
		ax20.plot(IpumpMax_vec, ibi_vec-ibi_vec_std, '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
		ax20.plot(IpumpMax_vec, ibi_vec ,color = colores[countr-1], linewidth=3)

		countr = countr+1


	return ax10, ax20, countr 

def robust_trend(cluster_ID, IpumpMax_vec, bd_vec, bd_first_quartile, bd_third_quartile ,ibi_vec, ibi_first_quartile, ibi_third_quartile, colores, countr, km, ax11, ax21, plist):
    p1 = (cluster_ID == 1)
    if len(plist)>=4:

        ax11.plot(IpumpMax_vec[p1], bd_vec[p1], 'd',color = colores[countr-1], alpha=0.5, markersize=km)
        ax11.plot(IpumpMax_vec[p1], bd_first_quartile[p1], '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
        ax11.plot(IpumpMax_vec[p1], bd_third_quartile[p1], '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
        ax11.plot(IpumpMax_vec[p1], bd_vec[p1] ,color = colores[countr-1], linewidth=3)

        ax21.plot(IpumpMax_vec[p1], ibi_vec[p1], 'd',color = colores[countr-1], alpha=0.5, markersize=km)
        ax21.plot(IpumpMax_vec[p1], ibi_first_quartile[p1], '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
        ax21.plot(IpumpMax_vec[p1], ibi_third_quartile[p1], '_',color = colores[countr-1], alpha=0.5, markersize=km, markeredgewidth=4)
        ax21.plot(IpumpMax_vec[p1], ibi_vec[p1] ,color = colores[countr-1], linewidth=3)

        countr = countr+1


    return ax11, ax21, countr 	


def median_interquartile_ranges(cluster_ID, IpumpMax, burst_duration, interburst_interval,colores, countr, plist, ax11, ax21):
    if cluster_ID == 1:
        if len(plist)>=4:
            # ax11.plot([IpumpMax, IpumpMax], [np.percentile(burst_duration, 25, interpolation='midpoint'), np.percentile(burst_duration, 75, interpolation='midpoint')], color = colores[countr-1], alpha=0.5, linewidth=3)
            # ax21.plot([IpumpMax, IpumpMax], [np.percentile(interburst_interval, 25, interpolation='midpoint'), np.percentile(interburst_interval, 75, interpolation='midpoint')], color = colores[countr-1], alpha=0.5, linewidth=3)

            ax11.plot([IpumpMax, IpumpMax], [np.percentile(burst_duration, 25, interpolation='midpoint'), np.percentile(burst_duration, 75, interpolation='midpoint')], color = colores[countr-1], alpha=cluster_ID, linewidth=3)
            ax21.plot([IpumpMax, IpumpMax], [np.percentile(interburst_interval, 25, interpolation='midpoint'), np.percentile(interburst_interval, 75, interpolation='midpoint')], color = colores[countr-1], alpha=cluster_ID, linewidth=3)

    return ax11, ax21




# interburst_interval right after the burst
# def bd_ibi_stats(list10, ax13, ax23,plist, fh, minimum_spikes, colores, countr):
# if len(plist)>=4:    

def bd_ibi_stats(list10, cluster_ID, ax12, ax22, plist, fh, minimum_spikes, colores, countr0):
    control = list10[-1]
    protocol = list10[-4]
    coupling = list10[-3]
    list1 = list10[0:-5]

    gp_vec = np.array([])
    IpumpMax_vec = np.array([])

    bd_mean = np.array([])
    bd_median = np.array([])

    ibi_mean = np.array([])
    ibi_median = np.array([])
    
    period_mean = np.array([])
    period_median = np.array([])

    if len(plist)>=4:    

        for j in range(0,len(list1)):
            cluster = cluster_ID[j]

            if cluster!=0:
                IpumpMax, gp, x_all, time,V, cytNa, interspike_tolerance, burst_duration, interburst_interval, cycle_period, interspike_tolerance = get_data.experimental_variables(control, protocol, coupling, list10, list1, j, fh)


                gp_vec = np.append(gp_vec, gp)
                IpumpMax_vec = np.append(IpumpMax_vec, IpumpMax)

                bd_mean = np.append(bd_mean, np.mean(burst_duration))
                ibi_mean = np.append(ibi_mean, np.mean(interburst_interval))
                period_mean = np.append(period_mean, np.mean(cycle_period))


                bd_median = np.append(bd_median, np.median(burst_duration))
                ibi_median = np.append(ibi_median, np.median(interburst_interval))
                period_median = np.append(period_median, np.median(cycle_period))


                p1 = len(interburst_interval)
                # print(len(burst_duration), len(interburst_interval))
                kolor = int(IpumpMax*10)

                ax12.plot(burst_duration[0:p1], interburst_interval[0:p1], 'd', color = colores[countr0])
                ax12.plot(burst_duration[0:p1], interburst_interval[0:p1], color = colores[countr0])

                ax22.plot(np.median(burst_duration), np.median(interburst_interval), 'd', color = colores[countr0])
                # ax22.plot(np.median(burst_duration), np.median(interburst_interval), color = colores[countr0])


                # ax12.plot(burst_duration[0:p1], interburst_interval[0:p1], 'd', color = colores[kolor])
                # ax12.plot(burst_duration[0:p1], interburst_interval[0:p1], color = colores[kolor])

                # ax22.plot(np.median(burst_duration), np.median(interburst_interval), 'd', color = colores[kolor])
                # ax22.plot(np.median(burst_duration), np.median(interburst_interval), color = colores[kolor])
        ax22.plot(bd_median, ibi_median, color = colores[countr0])

        countr0 = countr0+1
    # return ax12, ax22

    return ax12, ax22, countr0


def fancy_trends(ax1, ax2, f1, dir2, gP, protokol, suf1, ID, sweep_length):
    if sweep_length>=4:
        ax1.tick_params(axis='both', which='major', length=10, width=10)
        ax2.tick_params(axis='both', which='major', length=10, width=10)

        ax1.set_yticks(np.arange(0,7,2))
        ax2.set_yticks(np.arange(0,9,2))
        ax2.set_xticks(np.arange(0,1.1,0.2))
        ax1.set_xticks(np.arange(0,1.1,0.2))

        ax1.set_ylabel('BD(s)')
        ax2.set_ylabel('IBI(s)')

        ax1.set_xlim(0.05, 0.95)
        ax2.set_xlim(0.05, 0.95)

        ax1.set_ylim(0,8)
        ax2.set_ylim(0,8)


        f1.suptitle(r'$\bar{g}_P$= '+str(gP)+' (nS)')
        ax1.set_xlabel(r'I$^{pump}_{max}$ (nA)')
        ax2.set_xlabel(r'I$^{pump}_{max}$ (nA)')

        f1.savefig(dir2+str(int(gP))+protokol+'_'+str(ID)+'_'+suf1+'.png')    

def fancy_characteristics(ax1, ax2, f1, dir2, gP, protokol, xlabl, ylabl):
    ax1.tick_params(axis='both', which='major', length=10, width=10)
    ax2.tick_params(axis='both', which='major', length=10, width=10)

    # ax1.set_yticks(np.arange(0.5,3.6,1.0))
    # ax2.set_yticks(np.arange(0.5,3.6,1.0))
    # ax2.set_xticks(np.arange(0.5,3.6,1.0))
    # ax1.set_xticks(np.arange(0.5,3.6,1.0))

    tick_marks(0, 9, 0, 9, 1, 2, ax1)
    tick_marks(0, 9, 0, 9, 1, 2, ax2)

    ax1.set_xlim(0, 10)
    ax2.set_xlim(0, 10)

    ax1.set_ylim(0, 10)
    ax2.set_ylim(0, 10)


    ax1.set_xlabel(xlabl)
    ax2.set_xlabel(xlabl)
    ax1.set_ylabel(ylabl)

    ax1.set_title('Cycle-by-cycle')
    ax2.set_title('Trial Median')

    f1.suptitle(r'$\bar{g}_P$= '+str(gP)+' (nS)')

    f1.savefig(dir2+str(int(gP))+protokol+'_'+xlabl[0:-3]+'_'+ylabl[0:-3]+'.png')

def tick_marks(min_x, max_x, min_y, max_y, h_x, h_y, ax):
    ax.set_xticks(np.arange(min_x, max_x , h_x))
    ax.set_yticks(np.arange(min_y, max_y , h_y))
    return ax


def regress_bd_ibi_g(protocol, gP, predictor_x, predicted_y, ax, f, colores, kolor, sweep_length, sample_label, xlabl, ylabl, cluster_ID,dir3, choice):
    p1 = cluster_ID== choice
    if sweep_length >=4:
        r0 = linregress(predictor_x[p1], predicted_y[p1])
        print('Linear regression output. Params: gP= ', str(gP), ', IpumpMax Dynamic Clamp sweeps. ', sample_label)
        print('r2 = ', str('%.3f' % r0.rvalue))
        print('pval = ',str('%.3f' % r0.pvalue))
        print('predictor x=',str(xlabl))
        print('predicted y=',str(ylabl))
        print('y = ', str(r0.intercept), '+ x*',str(r0.slope), '\n')

        ax.scatter(predictor_x[p1], predicted_y[p1], color =colores[kolor])
        ax.plot(predictor_x[p1], r0.intercept+predictor_x[p1]* r0.slope, color=colores[kolor])
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)

        ax.tick_params(axis='both', which='major', length=10, width=10)
        ax.tick_params(axis='both', which='major', length=10, width=10)

        # ax1.set_yticks(np.arange(0.5,3.6,1.0))
        # ax2.set_yticks(np.arange(0.5,3.6,1.0))
        # ax2.set_xticks(np.arange(0.5,3.6,1.0))
        # ax1.set_xticks(np.arange(0.5,3.6,1.0))


        # tick_marks(0, 9, 1, 10, 2, 2, ax)
        # ax.set_xlim(0, 10)
        # ax.set_ylim(0, 10)

        ### ##raw
        # px = [0, 4]
        # py = [1, 4]
        # tick_marks(px[0], px[1], py[0]+1, py[1], 1, 1, ax)

        ## ##normalized
        px = [0, 1]
        py = [0, 1]
        tick_marks(px[0], px[1], py[0]+0.2, py[1], 0.2, 0.2, ax)

        ax.set_xlim(px[0], px[1])
        ax.set_ylim(py[0], py[1])

        ax.set_xlabel(xlabl)
        ax.set_ylabel(ylabl)

        f.suptitle('r2 = '+str('%.3f' % r0.rvalue)+', p = '+str('%.3f' % r0.pvalue)+', gP= '+str(gP)+' (nS)')


        f.savefig(dir3+str(int(gP))+str(protocol)+'_'+sample_label+'.png')

    else:
        r0 = 0
        # ax = 0

    return r0, ax



def polyfit_bd_ibi(x0,y0, degree, ax, countr, colores, cluster_ID, choice):
    p1 = cluster_ID ==choice
    x0 = x0[p1]
    y0 = y0[p1]
    if len(x0!=0):
        z = np.polyfit(x0,y0, degree)

        p = np.poly1d(z)
        # print(p)

        print('polynomial regression. degree= ',str(degree))
        print(p, '\n')


        ax.scatter(x0, p(x0), color='cyan', s=100)
        ax.scatter(x0,y0, color=colores[countr])

    else: 
        p=0

    return ax, p