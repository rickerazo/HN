# Average_voltage vs Cyt Na .py
####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps
# GOAL: plot Vm vs CytNa, in these ecperiments one parameter stays constant, and the other varies in a sweep fashion.
# Vm was averaged during the burst
# plots show that as CytNa increases, the average burst potential and the hyperpolarized phase are more depolarized and hyperpolarized, respectively. These demonstrate an intersting dynamic -> depolarizing and hyperpolarzing currents gain strength as gp increases

#### libraries necessary for code
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import find_peaks
from scipy.signal import peak_prominences
from matplotlib import pyplot as plt
# from mpl_toolkits import mplot3d
from scipy.stats import mode
import matplotlib
import numpy as np
# import itertools
import os
import sys
####	graphix stuff
matplotlib.rcParams['axes.linewidth']=10
font = {'weight' : 'bold',
        'size'   : 100}
plt.rc('font', **font)
plt.rcParams['agg.path.chunksize'] = 10000


mk1 = 16
mk2 = 15
mk3 = 3
mk10 = 35
lw1= 6
alpha0 = 1
markers = ['1', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia']

# coefficient of variance standard: reduce alpha
coef_std1 = 0.2
# coefficient of variane absolute standard: do not plot
covar = 0.55
minimum_spikes = 3
######################################################### file ID stuff
cwd = os.getcwd()

def compute(list10):
	constant1 = list10[-1]
	list0 = list10[0:-5]
	protocol = list10[-4]
	coupling = list10[-3]

	output_params = np.array([])
	if constant1=='I':
		IpumpMax = list10[-5]
		IpumpMax = float(IpumpMax)

		control_param = IpumpMax
		list1 = list0
		# list1 = list0[::-1]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			output_params = np.append(output_params, int(gp))

			x_all = int(gp)
			# print(list1)
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			time= np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			# 1. spike detection
			Vthreshold = -20
			spikes, vpks = find_peaks(V, height=Vthreshold)
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
			for j in range(0,len(first_spike)):
				t_ini = spike_times[first_spike[j]]
				t_end = spike_times[last_spike[j]]
				
				q0=np.nonzero(time<=t_end+isi)
				q1=np.nonzero(time[q0]>=t_ini-isi)

				v1 = V[q1]
				t1 = time[q1]
				# plt.plot(t1,v1)
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

						spike = V[p1]
						av_spike = np.mean(spike)
						V[p1] = av_spike
						ts = time[p1]
			x1 = cytNa
			y1 = V
			AVpmax.plot(x1,y1,'.-',markersize=20,color=colores[x_all],marker=markers[x_all],alpha=0.75)

	if constant1=='G':
		gp = list10[-5]
		gp = float(gp)
	
		control_param=gp
		x_all = int(control_param)

		list1 = list0[::-1]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			output_params = np.append(output_params, int(IpumpMax*10))

			x_all = int(IpumpMax*10)
			# print(str(gp)+','+str(IpumpMax))
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			time= np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			# 1. spike detection
			Vthreshold = -20
			spikes, vpks = find_peaks(V, height=Vthreshold)
			spike_times = time[spikes]
			spike_volts = V[spikes]
			spike_lag = np.diff(spike_times)	

			# last spike in bursts -> mechanism: identify spiking frequency, then when the interspike interval is greater than the mean spiking frequency + 4 std deviations
			# interspike interval tolerance for burst discrimination:
			p1 = np.nonzero(spike_lag>interspike_tolerance)
			p1 = p1[0]
			last_spike = np.append(p1,np.size(spike_times)-1)
			first_spike = 0
			first_spike = np.append(first_spike,p1+1)

			events1 = np.nonzero((last_spike-first_spike>=minimum_spikes)) #minimum five spikes to be considered a burst
			first_spike = first_spike[events1]
			last_spike = last_spike[events1]
			isi = np.mean(np.diff(spike_times))
					# Burst_volt 
			for k in range(0,len(first_spike)): # burst by burst fashion
				t_ini = spike_times[first_spike[k]]
				t_end = spike_times[last_spike[k]]

				q0=np.nonzero(time<=t_end+isi)
				q1=np.nonzero(time[q0]>=t_ini-isi)

				v1 = V[q1]
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

						spike = V[p1]
						av_spike = np.mean(spike)
						V[p1] = av_spike
						ts = time[p1]
					
			x1 = cytNa
			y1 = V

			AVgpax.plot(x1,y1,'.-',markersize=20,color=colores[x_all],alpha=0.65)

	if constant1=='C':
		gp = list10[-5]
		gp = float(gp)
	
		control_param=gp
		x_all = int(control_param)

		list1 = list0[::-1]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			output_params = np.append(output_params, int(IpumpMax*10))

			x_all = int(IpumpMax*10)
			# print(str(gp)+','+str(IpumpMax))
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			time= np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			# 1. spike detection
			Vthreshold = -20
			spikes, vpks = find_peaks(V, height=Vthreshold)
			spike_times = time[spikes]
			spike_volts = V[spikes]
			spike_lag = np.diff(spike_times)	

			# last spike in bursts -> mechanism: identify spiking frequency, then when the interspike interval is greater than the mean spiking frequency + 4 std deviations
			# interspike interval tolerance for burst discrimination:
			p1 = np.nonzero(spike_lag>interspike_tolerance)
			p1 = p1[0]
			last_spike = np.append(p1,np.size(spike_times)-1)
			first_spike = 0
			first_spike = np.append(first_spike,p1+1)

			events1 = np.nonzero((last_spike-first_spike>=minimum_spikes)) #minimum five spikes to be considered a burst
			first_spike = first_spike[events1]
			last_spike = last_spike[events1]
			isi = np.mean(np.diff(spike_times))
					# Burst_volt 
			for k in range(0,len(first_spike)): # burst by burst fashion
				t_ini = spike_times[first_spike[k]]
				t_end = spike_times[last_spike[k]]

				q0=np.nonzero(time<=t_end+isi)
				q1=np.nonzero(time[q0]>=t_ini-isi)

				v1 = V[q1]
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

						spike = V[p1]
						av_spike = np.mean(spike)
						V[p1] = av_spike
						ts = time[p1]
					
			x1 = cytNa
			y1 = V

			AVbaseline.plot(x1,y1,'.-',markersize=20,color=colores[x_all],alpha=0.65)

			
	return V,cytNa,time,x_all, float(list10[-5]), output_params, constant1, protocol, coupling

riff_gp = []
riff_pump = []
riff_ctrl = []
	
fileID= cwd[-8:len(cwd)]
fh = fileID+'_'
params = list(np.load(fh+'param.npy'))

for k in range(0,len(params)):# import file
	list20 = params[k]
	han2 = fh+str(list20)+'.npy'
	list10 =  list(np.load(han2))
	Ibias = float(list10[-2])
	if Ibias <=0.1:
		# print(str(fileID)+'	'+str(list10))

		if list10[-1] =='G':
			AV_gp = plt.figure(figsize=(35,30))
			AVgpax = plt.axes()
			AVgpax.yaxis.tick_right()
		elif list10[-1] =='I':
			AV_pm = plt.figure(figsize=(35,30))
			AVpmax = plt.axes()
			AVpmax.yaxis.tick_right()
		elif list10[-1] == 'C':
			AV_baseline = plt.figure(figsize=(35,30))
			AVbaseline = plt.axes()
			AVbaseline.yaxis.tick_right()
		# ax.yaxis.tick_right()
		V,cytNa,time,x_all, control_param, output_params, constant1, protocol, coupling = compute(list10)	
		
		print(constant1)
		print(control_param)
		print(output_params)
		print(protocol)
		print(coupling)
		# print('\n')
		for i in range(0,len(output_params)):
			index0 = int(output_params[i])
			if constant1 =='G':
				param_label = format(output_params[i]*0.1, '.1f')
				label1 = r'I$_{pump}^{max}$='+str(param_label)+' nA'
				label2 = str(int(control_param))
			elif constant1 == 'I':
				label1 = r'G$_P$='+str(output_params[i])+' nS'
				label2 = str(int(control_param*10))
			else:
				label1 = 'Baseline rhythm'
				label2 = 'baseline'

			plt.plot(cytNa[-1],V[-1], 'o',markersize = 1e-1, color = colores[index0], label = label1)

		# print(np.min(V),np.max(V), np.min(cytNa), np.max(cytNa))
		# xmiN = np.round(np.min(cytNa))
		# xmaX = np.round(np.max(cytNa))
		# ymiN = np.round(np.min(V))
		# ymaX = np.round(np.max(V))

		xmiN = 5
		xmaX = 40
		ymiN = -70
		ymaX = 10

		plt.tick_params(axis='both', which='major', length=35, width=15,labelsize=60)
		
		plt.xticks(np.arange(xmiN,xmaX,5),fontsize=80)
		plt.yticks(np.arange(ymiN,ymaX,15),fontsize=80)
		plt.xlim(5,np.round(np.max(cytNa))+2)

		plt.ylabel(r'V$_m$ (mV)',labelpad=10, fontweight = 'bold')
		plt.xlabel(r'[Na]$_i$ (mM)',labelpad=10, fontweight = 'bold')


		if protocol =='nm':
			label3 = 'CTRL'
		elif protocol == 'cs':
			label3 = r'Cs$^+$'
		else:
			label3 = 'error'

		if constant1 == 'I':
			plt.legend(ncol = 4, markerscale = 3e2, loc = [-0.05,1.005], fontsize = 52, fancybox=True, shadow=False)

			plt.text(5.5,-65,r'I$_{pump}^{max}$ = '+str(control_param)+ 'nA 	'+label3)

		elif constant1 == 'G':
			plt.legend(ncol = 4, markerscale = 3e2, loc = [-0.1,1.005], fontsize = 48, fancybox=True, shadow=False)

			plt.text(5.5,-65,r'G$_P$ = '+str(control_param)+ ' nS 	'+label3)

		else:
			plt.legend(ncol = 4, markerscale = 3e2, loc = [-0.1,1.005], fontsize = 48, fancybox=True, shadow=False)			

			plt.text(5.5,-65,'BL 	'+label3)

		plt.savefig('averageV_'+str(constant1)+label2+'_'+str(protocol)+'_'+str(coupling)+'.png')
		plt.close()
		print('\n')
