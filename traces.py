####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to take text file *.atf from pCLamp and import to analysis with numpy. Output of this script are: saved state variables: Ipump, Intracellular Na, injected membrane current, and neuron membrane potential; in addition, the neuron activity is also processed: neuron burst, interburst interval, cycle period, and spiking frequency.

# Libraries necessary for data import and processing:
import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import find_peaks
import os
## graphs stuff settings
font = {'weight' : 'bold',
        'size'   : 30}
plt.rc('font', **font)


# Import text file atf from clampfit : code assumes input file has format: time, Ipump, I_membrane, Na_in,Vm in that order
# fileID is the name of the input text .atf file
cwd = os.getcwd()

fileID= cwd[-8:len(cwd)]
f=open(fileID+'.atf')
dat = f.readlines()

f.close()
t = np.zeros((1,np.size(dat)))
t=t[0]
Ipump = np.zeros((1,np.size(dat)))
Ipump=Ipump[0]
I_mem = np.zeros((1,np.size(dat)))
I_mem=I_mem[0]
Na_in = np.zeros((1,np.size(dat)))
Na_in=Na_in[0]
Vm = np.zeros((1,np.size(dat)))
Vm = Vm[0]

for i in range(10,np.size(dat)):
	p0 = dat[i]
	p1=p0.split()
	t[i] = float(p1[0]) # seconds
	Ipump[i] = float(p1[1]) # nanoAmps
	I_mem[i] = float(p1[2]) # nanoAmps
	Na_in[i] = float(p1[3]) # nanoMoles
	Vm[i] = float(p1[4])	# miliVolts

######################################	TIME BINS 	########################################################################################################
# select a time range to plot (timebin), all the data can be overwhelming since it's a long trace
# carefully set the IpumpMax and gp from the readme file and/or timetags at the *.abf file:

IpumpMax = 0.8
gp= 6.0
# time bin: specifically analyze a time window after experimental manipulation
tm = [	3864.42, 3970.302	]
t_ini = tm[0]
t_end = tm[1]

q0=np.nonzero(t<t_end)
q1=np.nonzero(t[q0]>t_ini)

# assignment of state variables:
V = Vm[q1]
time = t[q1]
pump = Ipump[q1]
cyt_na = Na_in[q1]
I_tot = I_mem[q1]

######################################## VOLTAGE TRACES ######################################################################################
## plot the data

plt.figure(figsize=(20,10))
plt.plot(time,V)
plt.title(r'Membrane potential $HN_7$ $I_{pump}^{max}$='+str(IpumpMax)+' nA, $g_p$= '+str(gp)+' nS')
plt.ylabel(r'$V_m$(mV)')
plt.xlabel('Time(s)')


# data analysis
# 1. spike detection
Vthreshold = -20
spikes = find_peaks(V, height=Vthreshold, prominence=10)
spike_times = time[spikes[0]]

## plot spikes
# x-axis
ap = np.zeros((1,np.size(spike_times)))+Vthreshold
ap=ap[0]

plt.plot(spike_times,ap,'*')

spike_lag = np.diff(spike_times)

# last spike in bursts -> mechanism: identify spiking frequency, then when the interspike interval is greater than the mean spiking frequency + 4 std deviations
# interspike interval tolerance for burst discrimination:
interspike_tolerance = np.mean(spike_lag)+np.std(spike_lag)
p1 = np.nonzero(spike_lag>interspike_tolerance)
p1 = p1[0]
last_spike = np.append(p1,np.size(spike_times)-1)

# first spike in bursts
first_spike = 0
first_spike = np.append(first_spike,p1+1)

events1 = np.nonzero((last_spike-first_spike>=1)) #minimum five spikes to be considered a burst
first_spike = first_spike[events1]
last_spike = last_spike[events1]

ap = np.zeros((1,np.size(last_spike)))-13 #plotting arctifact
ap=ap[0]
plt.plot(spike_times[last_spike],ap,'*', label='last spike')
plt.plot(spike_times[first_spike],ap,'*', label='first spike')


# Id median spikes, burst by burst
middle_spike = np.zeros((1,np.size(first_spike)))
middle_spike = middle_spike[0]
mean_Hz = np.zeros((1,np.size(first_spike)))
mean_Hz = mean_Hz[0]
for i in range(0,np.size(first_spike)):
	burst = spike_times[first_spike[i]:last_spike[i]]
	middle_spike[i] = np.median(burst)
	mean_Hz[i] = np.mean(1/np.diff(burst))

plt.plot(middle_spike,ap,'*',label='median spike')

burst_duration = spike_times[last_spike] - spike_times[first_spike]
interburst_interval = spike_times[first_spike[1:np.size(first_spike)]] - spike_times[last_spike[0:-1]]
cycle_period = np.diff(middle_spike)

burst_duration = burst_duration[~np.isnan(burst_duration)]
interburst_interval = interburst_interval[~np.isnan(interburst_interval)]
cycle_period = cycle_period[~np.isnan(cycle_period)]
mean_Hz = mean_Hz[~np.isnan(mean_Hz)]


plt.plot([spike_times[first_spike],spike_times[first_spike]+burst_duration],[-10,-10],color='black',label='burst duration')
plt.plot([spike_times[last_spike[0:-1]],spike_times[last_spike[0:-1]]+interburst_interval],	[-8,-8],color='gray',label='interburst interval')
plt.plot([middle_spike[0:-1], middle_spike[0:-1]+cycle_period], [-9,-9],label='cycle period')


np.save(fileID+'_BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),burst_duration)
np.save(fileID+'_IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),interburst_interval)
np.save(fileID+'_T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),cycle_period)
np.save(fileID+'_Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),mean_Hz)

plt.savefig(fileID+'_Vm_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.png')
######################################## PUMP CURRENT ######################################################################################
## descriptive stats on the data
pump_mean = np.mean(pump)
pump_std = np.std(pump)

# detect peaks in the oscillations
pump_peaks = find_peaks(pump,height=[pump_mean+pump_std/2, np.max(pump)],distance=20000) #current issue with code: the "distance" must be adjusted to each data set to capture peaks correctly: need to debug this problem!
pump_peaks = pump_peaks[0]

#plot the data
plt.figure(figsize=(20,10))
plt.plot(time,pump)															# plot all data 

# plt.plot(time[pump_peaks],pump[pump_peaks],'*',markersize=20)				# plot peaks

plt.plot([time[0],time[-1]],[pump_mean, pump_mean])							# plot mean
plt.plot([time[0],time[-1]],[pump_mean+pump_std, pump_mean+pump_std])		# plot mean+std
plt.plot([time[0],time[-1]],[pump_mean-pump_std, pump_mean-pump_std])		# plot mean-std
plt.title(r'Pump current $HN_7$ $I_{pump}^{max}$='+str(IpumpMax)+' nA, $g_p$= '+str(gp)+' nS')
plt.ylabel(r'$I_{pump}$(nA)')
plt.xlabel('Time(s)')

#save the plot
plt.savefig(fileID+'_Ipump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.png')

######################################## CYTOSOL NA CONCENTRATION ######################################################################################
## descriptive stats on data
cyt_na_mean = np.mean(cyt_na)
cyt_na_std = np.std(cyt_na)
#filter the data
cyt_na_filtered = gaussian_filter(cyt_na,sigma=cyt_na_std*5)
cyt_na_mean = np.mean(cyt_na_filtered)
cyt_na_std = np.std(cyt_na_filtered)
# detect peaks in the oscillations
cyt_na_peaks = find_peaks(cyt_na_filtered,height=cyt_na_mean+cyt_na_std/2,distance=10000) #current issue with code: the "distance" must be adjusted to each data set to capture peaks correctly: need to debug this problem!
cyt_na_peaks = cyt_na_peaks[0]

## plot the data
plt.figure(figsize=(20,10))
plt.plot(time,cyt_na)																# plot raw data (noisy)
plt.plot(time,cyt_na_filtered)														# plot filtered data

# plt.plot(time[cyt_na_peaks],cyt_na[cyt_na_peaks],'*',markersize=20)					# plot pleaks

plt.plot([time[0],time[-1]],[cyt_na_mean, cyt_na_mean])								# plot mean
plt.plot([time[0],time[-1]],[cyt_na_mean+cyt_na_std, cyt_na_mean+cyt_na_std])		# plot mean + std
plt.plot([time[0],time[-1]],[cyt_na_mean-cyt_na_std, cyt_na_mean-cyt_na_std])		# plot mean - std

plt.title(r'Intracellular Na concentration $HN_7$ $I_{pump}^{max}$='+str(IpumpMax)+' nA, $g_p$= '+str(gp)+' nS')
plt.ylabel(r'Cytosol $Na_{in}$(nM)')
plt.xlabel('Time(s)')
#save the plot
plt.savefig(fileID+'_Na_in_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.png')


##########################################################################################################################################

np.save(fileID+'_V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),V)
np.save(fileID+'_time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),time)
np.save(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),pump)
np.save(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),cyt_na)
np.save(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),I_tot)
np.save(fileID+'_interspikeTol_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),interspike_tolerance)
##########################################################################################################################################
plt.ion()
plt.show()

# subsequent script is burst.py
