# custom_analysis.py

#this script will store all scripts that do some sort of custom analysis: For example intrabust instant spiking frequency

from scipy.signal import find_peaks
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.signal import find_peaks

from scipy.signal import savgol_filter

import matplotlib.pyplot as plt
import pandas as pd


import get_data

minimum_spikes = 3

spike_lag_standard_deviation_adjusted = 1.0

def instant_spiking_frequency(mean_v, time, cytNa, interspike_tolerance, x_all):
	max_spike = np.max(mean_v)

	# 1. spike detection
	Vthreshold = -30
	spikes, vpks = find_peaks(mean_v, height=Vthreshold, prominence = 4)
	spike_times = time[spikes]
	spike_volts = mean_v[spikes]
	spike_lag = np.diff(spike_times)	
	#######
	ap = np.zeros((1,np.size(spike_times)))+Vthreshold
	ap=ap[0]
	# ax1.plot(spike_times,ap,'*')
	spike_lag = np.diff(spike_times)
	interspike_tolerance = np.mean(spike_lag)+ np.std(spike_lag)*spike_lag_standard_deviation_adjusted

	p1 = np.nonzero(spike_lag>interspike_tolerance)
	p1 = p1[0]
	last_spike = np.append(p1,np.size(spike_times)-1)
	first_spike = 0
	first_spike = np.append(first_spike,p1+1)

	events1 = np.nonzero((last_spike-first_spike>=minimum_spikes)) #minimum five spikes to be considered a burst
	first_spike = first_spike[events1]
	last_spike = last_spike[events1]
	isi = np.mean(np.diff(spike_times))

	ap = np.zeros((1,np.size(last_spike)))-13 #plotting arctifact
	ap=ap[0]
	# ax1.plot(spike_times[last_spike],ap,'*', label='last spike', markersize=15)
	# ax1.plot(spike_times[first_spike],ap,'*', label='first spike', markersize=15)

	middle_spike = np.zeros((1,np.size(first_spike)))
	middle_spike = middle_spike[0]
	Hz = np.zeros((1,np.size(first_spike)))
	mean_f = Hz[0]
	median_f = Hz[0]

	counter = 0

	for i in range(0,np.size(first_spike)):
		burst = spike_times[first_spike[i]:last_spike[i]+1]
		middle_spike[i] = np.median(burst)
		# mean_Hz[i] = np.mean(1/np.diff(burst))
		f_instant = 1/np.diff(burst)
		
		mean_f[i] = np.mean(f_instant)

		median_f[i] = np.median(f_instant)

	return mean_f, median_f




def clustering_KMeans(x1, x2, y1, K):
	data = {'x':x1,'y':x2,'z':y1}
	df_spikef = pd.DataFrame(data,columns=['x','y','z'])
	# kmeans_spikef = KMeans(n_clusters=ncluj, random_state=1).fit(df_spikef)
	# kmeans_data = KMeans(n_clusters=K , random_state=1).fit(y1.reshape(-1,1))
	DV = pd.Series(df_spikef.z)
	kmeans_data = KMeans(n_clusters=K , random_state=1).fit(df_spikef)
	centroid_data = kmeans_data.cluster_centers_

	## print silhoutte score
	spikef_score = silhouette_score(y1.reshape(-1,1), kmeans_data.labels_, metric='euclidean', sample_size=None, random_state=1)
	print('Score of spike f clustering: '+str(spikef_score))

	return kmeans_data.labels_

# (median_f, mean_f, Vm_excursion, bd_median, ibi_median, plist)
def two_threshold_clustering(spike_frequency_characteristic, Vm_excursion, plist, bd_median, ibi_median):	
	if len(plist) >= 4:
		spikef_threshold = 9
		# spikef_threshold = 10
		Vm_excursion_threshold = 19

		# spikef_threshold = 7
		# Vm_excursion_threshold = 15

		# print( '\n', np.median(spike_frequency_characteristic), np.median(Vm_excursion), '\n')
		if np.median(Vm_excursion)> Vm_excursion_threshold and np.median(spike_frequency_characteristic) > spikef_threshold:
			cluster_ID = 1
		else:
			cluster_ID = 0
	else:
		cluster_ID= 0



		# the next line termines whether or not to discriminate between 
	return cluster_ID
	# return 1



def Voltage_excursion(mean_v, time, Vthreshold, interspike_tolerance, fh):

	max_spike = np.max(mean_v)
	# mean_v = mean_v
	time = time - time[0]
	# 1. spike detection
	spikes, vpks = find_peaks(mean_v, threshold=Vthreshold, prominence=5)
	spike_times = time[spikes]
	spike_volts = mean_v[spikes]
	spike_lag = np.diff(spike_times)	
	#######
	spike_lag = np.diff(spike_times)
	interspike_tolerance = np.mean(spike_lag)+ np.std(spike_lag)*spike_lag_standard_deviation_adjusted

	p1 = np.nonzero(spike_lag>interspike_tolerance)
	p1 = p1[0]
	last_spike = np.append(p1,np.size(spike_times)-1)
	first_spike = 0
	first_spike = np.append(first_spike,p1+1)

	events1 = np.nonzero((last_spike-first_spike>=minimum_spikes)) #minimum five spikes to be considered a burst
	first_spike = first_spike[events1]
	last_spike = last_spike[events1]
	isi = np.mean(np.diff(spike_times))

	bursting_potential = np.array([])
	quiescent_potential = np.array([])

	# Burst_volt 
	for k in range(0,len(first_spike)): # burst by burst fashion
		t_ini = spike_times[first_spike[k]]
		t_end = spike_times[last_spike[k]]

		q0=np.nonzero(time<=t_end- isi )
		q1=np.nonzero(time[q0]>=t_ini+isi)

		v1 = mean_v[q1]
		t1 = time[q1]
		burst_peaks, burst_v = find_peaks(v1, threshold=Vthreshold, prominence=4)

		if np.size(burst_peaks)>=2:
			for nm in range(0,len(burst_peaks)):	#spike by spike inside burst
				thisspike = burst_peaks[nm]

				if nm < len(burst_peaks)-1:
					nextspike = burst_peaks[nm+1]
					v2 = v1[thisspike:nextspike]
					t2 = t1[thisspike:nextspike]
					local_minima = np.min(v2)
					p0 = np.nonzero(v2==local_minima)
					t_minima = t2[p0]

					if len(t_minima)>1:
						t_minima = t_minima[-1]

					lag1 = t_minima - t1[burst_peaks[nm]]
					t_up = t1[burst_peaks[nm]] - lag1
					t_down = t1[burst_peaks[nm]]+ lag1
					p0=np.nonzero(time<=t_down)
					p1=np.nonzero(time[p0]>=t_up)


				else:
					v2 = v1[burst_peaks[nm-1]:thisspike]
					t2 = t1[burst_peaks[nm-1]:thisspike]

					local_minima = np.min(v2)
					p0 = np.nonzero(v2==local_minima)
					t_minima = t2[p0]

					if len(t_minima)>1:
						t_minima = t_minima[-1]

					lag1 = t_minima - t1[burst_peaks[nm]]
					
					t_up = t1[burst_peaks[nm]] + lag1
					t_down = t1[burst_peaks[nm]]- lag1
					p0=np.nonzero(time<=t_down)
					p1=np.nonzero(time[p0]>=t_up)

				spike = mean_v[p1]
				av_spike = np.mean(spike)
				mean_v[p1] = av_spike
				
				mean_v[p1] = av_spike
				ts = time[p1]

		bursting_base_potential = mean_v[q1]

		bursting_potential = np.append(bursting_potential, np.median(bursting_base_potential))


	for k in range(0,len(first_spike)-1):
		t_ini = spike_times[last_spike[k]]
		t_end = spike_times[first_spike[k+1]]

		q0=np.nonzero(time<=t_end- isi )
		q1=np.nonzero(time[q0]>=t_ini+isi)

		a0, b0 = np.shape(q1)

		if b0!=0:
			trough_potential = np.min(mean_v[q1])
		else:
			trough_potential = np.min(mean_v)
		t1 = time[q1]

		quiescent_potential = np.append(quiescent_potential, trough_potential)

	return quiescent_potential, bursting_potential



def flat_spikes(mean_v, t, Vthreshold, interspike_tolerance, burst_duration, Na_in):
	y2 = mean_v.copy()

	# if np.mean(burst_duration) >=3: 
	# 	b4 = 25
	# if np.mean(burst_duration) <3: 
	# 	b4 = 16
	
	b4 = 16
	q00=np.nonzero(t<b4)

	max_spike = np.max(mean_v)
	# mean_v = mean_v
	time = t - t[0]
	# 1. spike detection
	# Vthreshold = -20
	spikes, vpks = find_peaks(mean_v, height=Vthreshold, prominence=4)
	spike_times = time[spikes]
	spike_volts = mean_v[spikes]
	spike_lag = np.diff(spike_times)	
	#######
	ap = np.zeros((1,np.size(spike_times)))+Vthreshold
	ap=ap[0]
	# ax1.plot(spike_times,ap,'*')
	spike_lag = np.diff(spike_times)
	interspike_tolerance = np.mean(spike_lag)+ np.std(spike_lag)*spike_lag_standard_deviation_adjusted

	p1 = np.nonzero(spike_lag>interspike_tolerance)
	p1 = p1[0]
	last_spike = np.append(p1,np.size(spike_times)-1)
	first_spike = 0
	first_spike = np.append(first_spike,p1+1)

	events1 = np.nonzero((last_spike-first_spike>=0)) #minimum five spikes to be considered a burst
	first_spike = first_spike[events1]
	last_spike = last_spike[events1]
	isi = np.mean(np.diff(spike_times))

	ap = np.zeros((1,np.size(last_spike)))-13 #plotting arctifact
	ap=ap[0]
	# ax1.plot(spike_times[last_spike],ap,'*', label='last spike', markersize=15)
	# ax1.plot(spike_times[first_spike],ap,'*', label='first spike', markersize=15)

	middle_spike = np.zeros((1,np.size(first_spike)))
	middle_spike = middle_spike[0]
	mean_Hz = np.zeros((1,np.size(first_spike)))
	mean_Hz = mean_Hz[0]

	counter = 0
	for i in range(0,np.size(first_spike)):
		burst = spike_times[first_spike[i]:last_spike[i]+1]
		middle_spike[i] = np.median(burst)
		mean_Hz[i] = np.mean(1/np.diff(burst))
		# f_instant = 1/np.diff(burst)

		# ax2.plot(burst[0:-1], f_instant, c='k', linewidth=3)

	burst_duration = spike_times[last_spike] - spike_times[first_spike]
	interburst_interval = spike_times[first_spike[1:np.size(first_spike)]] - spike_times[last_spike[0:-1]]
	cycle_period = np.diff(middle_spike)

	burst_duration = burst_duration[~np.isnan(burst_duration)]
	interburst_interval = interburst_interval[~np.isnan(interburst_interval)]
	cycle_period = cycle_period[~np.isnan(cycle_period)]
	mean_Hz = mean_Hz[~np.isnan(mean_Hz)]

	guide_T = [0]
	guide_BD = [2]
	guide_IBI = [1]


	# ax1.plot([middle_spike[guide_T], middle_spike[guide_T]+cycle_period[guide_T]], [max_spike+6, max_spike+6],label='cycle period', linewidth=4)

	# ax1.plot([spike_times[first_spike[guide_BD]],spike_times[first_spike[guide_BD]]+burst_duration[guide_BD]],[max_spike+6, max_spike+6],color='blue',label='burst duration', linewidth=4)
	# ax1.plot([spike_times[last_spike[guide_IBI]],spike_times[last_spike[guide_IBI]]+interburst_interval[guide_IBI]], [max_spike+4 , max_spike+4],color='gray',label='interburst interval', linewidth=4)

	# ax1.plot(middle_spike[0], max_spike+3, 'D', markersize=5, c='k', fillstyle='none', markeredgewidth = 2)
	# ax1.plot(middle_spike[1], max_spike+3, 'D', markersize=5, c='k', fillstyle='none', markeredgewidth = 2)
	# ax1.plot(middle_spike[2], max_spike+3, 'D', markersize=5, c='k', fillstyle='none', markeredgewidth = 2)



	# Burst_volt 
	for k in range(0,len(first_spike)): # burst by burst fashion
		t_ini = spike_times[first_spike[k]]
		t_end = spike_times[last_spike[k]]

		q0=np.nonzero(time<=t_end+isi)
		q1=np.nonzero(time[q0]>=t_ini-isi)

		v1 = mean_v[q1]
		t1 = time[q1]
		burst_peaks, burst_v = find_peaks(v1, height=Vthreshold)

		if np.size(burst_peaks)>=2:
			for nm in range(0,len(burst_peaks)):	#spike by spike inside burst
				thisspike = burst_peaks[nm]

				if nm < len(burst_peaks)-1:
					nextspike = burst_peaks[nm+1]
					v2 = v1[thisspike:nextspike]
					t2 = t1[thisspike:nextspike]
					local_minima = np.min(v2)
					p0 = np.nonzero(v2==local_minima)
					t_minima = t2[p0]

					if len(t_minima)>1:
						t_minima = t_minima[-1]

					lag1 = t_minima - t1[burst_peaks[nm]]
					t_up = t1[burst_peaks[nm]] - lag1
					t_down = t1[burst_peaks[nm]]+ lag1
					p0=np.nonzero(time<=t_down)
					p1=np.nonzero(time[p0]>=t_up)


				else:
					v2 = v1[burst_peaks[nm-1]:thisspike]
					t2 = t1[burst_peaks[nm-1]:thisspike]

					local_minima = np.min(v2)
					p0 = np.nonzero(v2==local_minima)
					t_minima = t2[p0]

					if len(t_minima)>1:
						t_minima = t_minima[-1]

					lag1 = t_minima - t1[burst_peaks[nm]]
					
					t_up = t1[burst_peaks[nm]] + lag1
					t_down = t1[burst_peaks[nm]]- lag1
					p0=np.nonzero(time<=t_down)
					p1=np.nonzero(time[p0]>=t_up)

				spike = mean_v[p1]
				av_spike = np.mean(spike)
				mean_v[p1] = av_spike
				
				# mean_v[p1] = av_spike
				ts = time[p1]


	################# xlim
	if np.mean(burst_duration) >=3: 
		# b4 = 25

		b1 = b4+4	#x_lim
		b2 = 14	#bar
		b3 = 14.5	#label

		
	if np.mean(burst_duration) <3: 
		# b4 = 14

		b1 = b4+4	#x_lim
		b2 = 4	#bar
		b3 = 4.5	#label


	# q00=np.nonzero(time<20)
	q00=np.nonzero(time<b4)
	# q10=np.nonzero(t[q00]>t_ini)
			
	x1 = time[q00]
	y1 = mean_v[q00]
	y3 = y2[q00]
	y4 = Na_in[q00]
	# y1 = mean_v[q00]
	# ax1.legend()
	# ax1.plot(x1,y1,'.-',markersize=5,color=colores[x_all],alpha=0.65)
	# font1 = 15
	# ax1.plot(x1,y1,linewidth=5,color=colores[x_all],alpha=0.65)
	# # ax1.plot([20.5, 20.5], [np.min(y1), np.max(y1)], linewidth=10, color=colores[x_all])

	# ax1.plot([-1.5, -1.5],[np.min(y1), np.min(y1)+20], color=colores[x_all], linewidth=6)
	# ax2.plot([-1.5, -1.5],[5,15], color=colores[x_all], linewidth=6)

	# ax1.text(-3, np.min(y1)+1, '20 mV', rotation = 90, fontsize= font1)
	# ax2.text(-3, 7, '10 Hz', rotation = 90, fontsize= font1)

	# ax1.text(b4+0.25, -47, '-45 mV', fontsize=font1)
	# ax2.text(b4+0.25, 10, '10 Hz', fontsize=font1)

	# ax1.set_xlim(-5,b1)
	# ax2.set_xlim(-5,b1)

	# ax1.set_ylim(-70, 20)
	# ax2.set_ylim(0, 20)

	# ############ TEN second
	# ax1.plot([b2,b2+10],[19,19],linewidth=6,color='black')
	# ax1.text(b3, 23,'10 seconds' , fontsize = font1)

	# ########## 45 second dashed
	# ## 50 mV dotted line
	# ax1.plot([time[0], b4],[-45, -45],'--',color='black',linewidth=3)
	# ## 10 Hz dash
	# ax2.plot([time[0], b4],[10, 10],'--',color='black',linewidth=3)


	return x1, y1, y3, y4





def normalize_Nai(list10, fh):
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

    p1 = len(peaks)
    p2 = len(min_peaks)
    p3 = np.min([p1,p2])
    cytNa_excursion = np.mean(cytNa[peaks[0:p3]]-cytNa[min_peaks[0:p3]])

    # return norm_cytNa, cytNa[peaks], cytNa[min_peaks], time[peaks], time[min_peaks], cytNa_excursion
    return cytNa_excursion




def average_voltage(mean_v, t, Vthreshold, interspike_tolerance, burst_duration, Na_in):
	y2 = mean_v.copy()

	b4 = t[-1]
	q00=np.nonzero(t<b4)

	time = t - t[0]

	# 1. spike detection
	# Vthreshold = -20
	spikes, vpks = find_peaks(mean_v, height=Vthreshold, prominence=4)
	spike_times = time[spikes]
	spike_volts = mean_v[spikes]
	spike_lag = np.diff(spike_times)	
	#######
	# print(np.shape(spike_times))
	ap = np.zeros((1,np.size(spike_times)))+Vthreshold
	ap=ap[0]


	# interspike_tolerance = np.mean(spike_lag)+ np.std(spike_lag)*spike_lag_standard_deviation_adjusted
	interspike_tolerance = 1e-1

	for i in range(0,len(ap)):
		thisspike = spike_times[i]
		previous_point = spike_times[i]- interspike_tolerance
		next_point = spike_times[i]+ interspike_tolerance

		# print(spike_times[i])
		# print(next_point)
		q0=np.nonzero(time<=next_point)
		q1=np.nonzero(time[q0]>=previous_point)

		v1 = mean_v[q1]

		mean_v[q1] = np.mean(v1)

	return mean_v, time, y2