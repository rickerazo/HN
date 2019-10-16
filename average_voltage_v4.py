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
from mpl_toolkits import mplot3d
from scipy.stats import mode
import matplotlib
import numpy as np
import itertools
import os
import sys
####	graphix stuff
matplotlib.rcParams['axes.linewidth']=10



font = {'weight' : 'bold',
        'size'   : 250}
plt.rc('font', **font)
mk1 = 16
mk2 = 15
mk3 = 3
mk10 = 35
lw1= 6
alpha0 = 1
markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = [' ', 'brown', 'r', 'c','m','k','y','b', 'orange', 'brown','b', 'g', 'r', 'c','m','k','y','b', 'g', 'b']

# coefficient of variance standard: reduce alpha
coef_std1 = 0.2
# coefficient of variane absolute standard: do not plot
covar = 0.55
# AV_gp = plt.figure(1,figsize=(50,50))
# AVgpax = plt.axes()

# AV_pm = plt.figure(6,figsize=(75,50))
# AVpmax = plt.axes()

# gra1= plt.figure(1,figsize=(50,50))
# ggax = plt.axes()
# gra2= plt.figure(1,figsize=(50,50))
# gpax = plt.axes()
######################################################### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
experiment_list = np.load('experiment_list.npy')
# 
def compute(list10):
	constant1 = list10[-1]
	list0 = list10[0:-2]
	if constant1=='I':
		IpumpMax = list10[-2]
		IpumpMax = float(IpumpMax)

		control_param = IpumpMax
		# x_all = int(control_param)*10
		# print(x_all)
		list1 = list0[::-1]

		control_param = IpumpMax

		for j in range(0,len(list1)):
			gp = float(list1[j])
			x_all = int(gp)
			# print(list1)
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			time= np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
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

			events1 = np.nonzero((last_spike-first_spike>=1)) #minimum five spikes to be considered a burst
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
							# plt.figure(10)
							# plt.plot(t2,v2,'k')
							local_minima = np.min(v2)
							p0 = np.nonzero(v2==local_minima)
							t_minima = t2[p0]

							if len(t_minima)>1:
								t_minima = t_minima[-1]

							lag1 = t_minima - t1[burst_peaks[nm]]
							# plt.figure(10)
							# plt.plot(t2[p0],v2[p0],'*k',markersize=10)
							# plt.figure(13)
							# plt.plot(t2,v2)
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
		gp = list10[-2]
		gp = float(gp)
	
		control_param=gp
		x_all = int(control_param)

		list1 = list0[::-1]

		control_param = gp

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			x_all = int(IpumpMax*10)
			# print(str(gp)+','+str(IpumpMax))
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			time= np.load(fh+'time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			interspike_tolerance = np.load(fh+'interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
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

			events1 = np.nonzero((last_spike-first_spike>=1)) #minimum five spikes to be considered a burst
			first_spike = first_spike[events1]
			last_spike = last_spike[events1]
			# nr_spike = np.zeros((1,len(first_spike)))
			# nr_spike = nr_spike[0]
			isi = np.mean(np.diff(spike_times))
# time - raw data, t1 - data to be analyzed, t2 = burst time.t3 spike time
					# Burst_volt 
			for k in range(0,len(first_spike)): # burst by burst fashion
				t_ini = spike_times[first_spike[k]]
				t_end = spike_times[last_spike[k]]
				# plt.figure(11)
				# plt.plot([t_ini, t_end],[0,0])

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
							# plt.figure(10)
							# plt.plot(t2,v2,'k')
							local_minima = np.min(v2)
							p0 = np.nonzero(v2==local_minima)
							t_minima = t2[p0]

							if len(t_minima)>1:
								t_minima = t_minima[-1]

							lag1 = t_minima - t1[burst_peaks[nm]]
							# plt.figure(10)
							# plt.plot(t2[p0],v2[p0],'*k',markersize=10)
							# plt.figure(13)
							# plt.plot(t2,v2)
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
			# if np.mean(cytNa>0):
				
			# 	t1 = float(list1[j])
			# 	t1 = np.array(t1)*10
			# 	t2 = int(t1)

			AVgpax.plot(x1,y1,'.-',markersize=20,color=colores[x_all],alpha=0.65)
	return V,cytNa,time,x_all



# IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
# experiment_list = ['19809001']
experiment_list = ['19826000']
ctr_prm=[]
str_prm=[]
for i in range(0,len(experiment_list)):
	
	fileID= str(experiment_list[i])
	fh = fileID+'/'+fileID+'_'
	params = list(np.load(fh+'param.npy'))

	for k in range(0,len(params)):# import file
		# list20 = params[0]
		list20 = params[k]
		han2 = fh+str(list20)+'.npy'
		list10 =  list(np.load(han2))
		print(str(fileID)+'	'+str(list10))
		# print(list10)
		ctr_prm.append(list10[-2])
		str_prm.append(list10[-1])
		if list10[-1] =='G':
			AV_gp = plt.figure(k,figsize=(100,75))
			AVgpax = plt.axes()
		if list10[-1] =='I':
			AV_pm = plt.figure(k,figsize=(100,75))
			AVpmax = plt.axes()			

		V,cytNa,time,x_all=compute(list10)
		# spike,ts,V,time = compute(list10)


riff1 = [1,2,3,4]#,5,6,7]#,8]
color_riff1 = [colores[1],colores[2], colores[3], colores[4]]#, colores[5], colores[6], colores[7]]#, colores[8]]

riff2 = [2,3,4,5,6,7,8,9]
color_riff2 = [colores[2], colores[3], colores[4], colores[5], colores[6], colores[7], colores[8],colores[9]]

riff3 = [2,3,4,5,6,7,8,9]
color_riff3 = [colores[2], colores[3], colores[4], colores[5], colores[6], colores[7], colores[8],colores[9]]

riff4 = [2,3,4]
color_riff4 = [colores[2],colores[3],colores[4]]

riff5 = [4,5,6,7,]
color_riff5 = [colores[4],colores[5],colores[6],colores[7]]

titles = [' ','BD','IBI','T','Hz','Ipump_cytNa','BD','IBI','T','Hz','Ipump_cytNa']
phrase= [' ','_gp','_pm']

plt.figure(0)
for j in range(0,len(riff1)):
	p1 = float(riff1[j])
	p1 = format(format(p1,'.1f'))
	str1 = str(p1)
	plt.plot(0,0,'o',markersize=1e-3,color=color_riff1[j],label=r'$\bar{g}_P$= '+str1+' nS')
# plt.axis([7,15,-65,15])
plt.title(r'I$_{pump}^{max}$='+str(ctr_prm[0])+' nA')
plt.legend(loc=2,bbox_to_anchor=(-0.1,1.2),markerscale=1e+5,ncol=3,fontsize=200)
# plt.savefig(fileID+'_'+str(str_prm[0]+str(ctr_prm[0]))+'V_Na_average'+'.png')#+phrase[i])

plt.figure(1)
for j in range(0,len(riff2)):
	p1 = float(riff2[j])*0.1
	p1 = format(format(p1,'.1f'))

	str1 = str(p1)
	plt.plot(0,0,'o',markersize=1e-3,color=color_riff2[j],label=r'I$_{pump}^{max}$= '+str1+' nA')
# plt.axis([7,15,-75,15])	
plt.title(r'$\bar{g}_P$='+str(ctr_prm[1])+' nS')
plt.legend(loc=2,bbox_to_anchor=(-0.1,1.2),markerscale=1e+5,ncol=3,fontsize=200)
# plt.savefig(fileID+'_'+str(str_prm[1]+str(ctr_prm[1]))+'V_Na_average'+'.png')#+phrase[i])


plt.figure(2)
for j in range(0,len(riff3)):
	p1 = float(riff3[j])*0.1
	p1 = format(format(p1,'.1f'))

	str1 = str(p1)
	plt.plot(0,0,'o',markersize=1e-3,color=color_riff3[j],label=r'I$_{pump}^{max}$= '+str1+' nA')
# plt.axis([7,15,-75,15])
plt.title(r'\bar{g}$_P$='+str(ctr_prm[2])+' nS')
plt.legend(loc=2,bbox_to_anchor=(-0.1,1.2),markerscale=1e+5,ncol=3,fontsize=200)
# plt.savefig(fileID+'_'+str(str_prm[2]+str(ctr_prm[2]))+'V_Na_average'+'.png')#+phrase[i])

plt.figure(3)
for j in range(0,len(riff4)):
	p1 = float(riff4[j])*0.1
	p1 = format(format(p1,'.1f'))

	str1 = str(p1)
	plt.plot(0,0,'o',markersize=1e-3,color=color_riff4[j],label=r'$I_{pump}^{max}$= '+str1+' nA')
# plt.axis([7,15,-75,15])
plt.title(r'$\bar{g}_P$='+str(ctr_prm[3])+' nS')
plt.legend(loc=2,bbox_to_anchor=(-0.1,1.2),markerscale=1e+5,ncol=3,fontsize=200)
# plt.savefig(fileID+'_'+str(str_prm[2]+str(ctr_prm[2]))+'V_Na_average'+'.png')#+phrase[i])
plt.figure(4)
for j in range(0,len(riff5)):
	p1 = float(riff5[j])
	p1 = format(format(p1,'.1f'))

	str1 = str(p1)
	plt.plot(0,0,'o',markersize=1e-3,color=color_riff5[j],label=r'$\bar{g}_P$= '+str1+' nA')
# plt.axis([9,15,-75,15])
plt.title(r'$\bar{g}_P$='+str(ctr_prm[4])+' nS')
plt.legend(loc=2,bbox_to_anchor=(-0.1,1.2),markerscale=1e+5,ncol=3,fontsize=200)
# plt.savefig(fileID+'_'+str(str_prm[2]+str(ctr_prm[2]))+'V_Na_average'+'.png')#+phrase[i])

# for c1 in range(3,6):
# 	plt.figure(c1)
# 	for c2 in range(len(riff4)):
# 		p1 = float(riff4[c2])
# 		p1 = format(format(p1,'.1f'))

# 		str1 = str(p1)
# 		plt.plot(0,0,'o',markersize=1e-3,color=color_riff4[c2],label=r'G$_P$= '+str1+' nS')
# 	plt.axis([9,20,-75,15])
# 	plt.title(r'I$_{pump}^{max}$='+str(ctr_prm[c1])+' nA')
# 	plt.legend(loc=2,bbox_to_anchor=(-0.1,1.15),markerscale=1e+5,ncol=4,fontsize=100)

for i in range(0,len(str_prm)):
	plt.figure(i)
	plt.ylabel('Membrane potential V(mV)',weight='bold')
	plt.yticks(np.arange(-70,-10,20))
	plt.xticks(np.arange(10,30,5))
	# plt.axis([10,21,-85,20])
	plt.axis([9,15,-65,10])
	plt.xlabel('Cytosol Sodium concentration (mM)',weight='bold')
	plt.savefig(fileID+'_'+str(str_prm[i]+str(ctr_prm[i]))+'V_Na_average'+'.png')#+phrase[i])

# plt.figure(6)
# plt.savefig(fileID+'_'+str(str_prm[5]+str(ctr_prm[5]))+'V_Na_average'+'.png')#+phrase[i])
plt.ion()
plt.show()