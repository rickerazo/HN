####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps
# 

#### libraries necessary for code
from matplotlib import pyplot as plt
import numpy as np
import itertools
import os
####	graphix stuff
font = {'weight' : 'bold',
        'size'   : 60}
plt.rc('font', **font)
mk = 20
mk2 = 8
lw= 15
alpha0 = 1
markers = ['.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['b', 'g', 'r', 'c','m','k','y','b', 'g', 'b']

# The next section analyzes the data files specified by the max and min parameters above. There is one section for gp sweeps and another for Ipumpmax sweeps.
def compute(list10):
	control = list10[-1]
	list1 = list10[0:-2]
	####################################################################################	
	############################# Sodium sweep ###########################################	
	if control == 'I':
		IpumpMax = list10[-2]	
		list0=list1[::-1] #reverse the order of the list, for plotting purposes
		##
		# File ID stuff. It must be after parameters are set
		str1 = 'BD_IpumpMax='+str(IpumpMax)+'_gp='
		str2 = 'IBI_IpumpMax='+str(IpumpMax)+'_gp='
		str3 = 'T_IpumpMax='+str(IpumpMax)+'_gp='
		str4 = 'Hz_IpumpMax='+str(IpumpMax)+'_gp='
		strend = '.npy'

		##
		p1 = len(list0)
		BD_mean = np.zeros((1,p1))
		BD_mean = BD_mean[0]
		BD_std =  np.zeros((1,p1))
		BD_std = BD_std[0]

		IBI_mean = np.zeros((1,p1))
		IBI_mean = IBI_mean[0]
		IBI_std = np.zeros((1,p1))
		IBI_std = IBI_std[0]

		T_mean = np.zeros((1,p1))
		T_mean = T_mean[0]
		T_std =  np.zeros((1,p1))
		T_std = T_std[0]

		Hz_mean = np.zeros((1,p1))
		Hz_mean = Hz_mean[0]
		Hz_std = np.zeros((1,p1))
		Hz_std = Hz_std[0]

		for i in range(0,np.size(list0)):
			j=list1[i]
			gp = j
			#BD
			filename1 = fileID+'_'+str1+str(gp)+strend
			burst_duration = np.load(filename1)

			x = np.zeros((1,np.size(burst_duration)))
			x=x[0]+gp

			plt.figure(1)
			# plt.plot([gp,gp],[np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration)+np.std(burst_duration)])
			plt.plot(gp,np.max(burst_duration),'^',markersize=mk,color='green')
			plt.plot(gp,np.min(burst_duration),'v',markersize=mk,color='red')
			plt.plot(x,burst_duration,'.',markersize=mk-10)
			plt.plot(gp,np.mean(burst_duration),'*',markersize=mk,color='black')
			plt.title(r'Burst duration $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('time(s)')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			BD_mean[i] = np.mean(burst_duration)
			BD_std[i] = np.std(burst_duration)

			#IBI
			filename2 = fileID+'_'+str2+str(gp)+strend
			interburst_interval = np.load(filename2)

			x = np.zeros((1,np.size(interburst_interval)))
			x=x[0]+gp
			plt.figure(2)
			plt.plot(gp,np.max(interburst_interval),'^',markersize=mk,color='green')
			plt.plot(gp,np.min(interburst_interval),'v',markersize=mk,color='red')
			plt.plot(x,interburst_interval,'.',markersize=mk-10)
			plt.plot(gp,np.mean(interburst_interval),'*',markersize=mk,color='black')
			plt.title(r'Interburst interval $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('time(s)')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			IBI_mean[i] = np.mean(interburst_interval)
			IBI_std[i] = np.std(interburst_interval)

			#T
			filename3 = fileID+'_'+str3+str(gp)+strend
			cycle_period = np.load(filename3)

			x = np.zeros((1,np.size(cycle_period)))
			x=x[0]+gp

			plt.figure(3)
			plt.plot(gp,np.max(cycle_period),'^',markersize=mk,color='green')
			plt.plot(gp,np.min(cycle_period),'v',markersize=mk,color='red')
			plt.plot(x,cycle_period,'.',markersize=mk-10)
			plt.plot(gp,np.mean(cycle_period),'*',markersize=mk,color='black')
			plt.title(r'Cycle period $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('time(s)')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			T_mean[i] = np.mean(cycle_period)
			T_std[i] = np.std(cycle_period)

			#Hz
			filename4 = fileID+'_'+str4+str(gp)+strend
			mean_hz = np.load(filename4)

			x = np.zeros((1,np.size(mean_hz)))
			x=x[0]+gp

			plt.figure(4)
			plt.plot(gp,np.max(mean_hz),'^',markersize=mk,color='green')
			plt.plot(gp,np.min(mean_hz),'v',markersize=mk,color='red')
			plt.plot(x,mean_hz,'.',markersize=mk-10)
			plt.plot(gp,np.mean(mean_hz),'*',markersize=mk,color='black')

			plt.title(r'$\mu$ intrabust spiking frequency $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('Hz')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			Hz_mean[i] = np.mean(mean_hz)
			Hz_std[i] = np.std(mean_hz)

			# Pump Na
			shift = 0.01
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Ipump = Ipump/IpumpMax
			CytNa = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			plt.figure(5)
			plt.plot(CytNa,Ipump, marker=markers[i], color=colores[i],linewidth=lw,linestyle='', markersize=mk2 ,label=r'$g_p$='+str(gp))
			plt.plot(np.mean(CytNa),np.min(Ipump), '*',markersize=mk, color=colores[i])
			plt.plot(np.min(CytNa),np.min(Ipump), marker=markers[i], markersize=mk, color=colores[i])
			plt.plot(np.max(CytNa),np.min(Ipump), marker=markers[i], markersize=mk, color=colores[i])

		plt.figure(1)
		# plt.plot(gp,np.mean(burst_duration)+np.std(burst_duration),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		# plt.plot(gp,np.mean(burst_duration)-np.std(burst_duration),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		# plt.plot(gp,np.mean(burst_duration),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_gp,max_gp,1),BD_mean,linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),BD_mean+BD_std,':',linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),BD_mean-BD_std,':',linewidth=lw)
		plt.legend()
		plt.savefig(fileID+'_burst_duration.png')

		plt.figure(2)
		# plt.plot(gp,np.mean(interburst_interval)+np.std(interburst_interval),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		# plt.plot(gp,np.mean(interburst_interval)-np.std(interburst_interval),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		# plt.plot(gp,np.mean(interburst_interval),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_gp,max_gp,1),IBI_mean,linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),IBI_mean+IBI_std,':',linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),IBI_mean-IBI_std,':',linewidth=lw)
		plt.legend()
		plt.savefig(fileID+'_interburst_interval.png')

		plt.figure(3)
		# plt.plot(gp,np.mean(cycle_period)+np.std(cycle_period),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		# plt.plot(gp,np.mean(cycle_period)-np.std(cycle_period),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		# plt.plot(gp,np.mean(cycle_period),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_gp,max_gp,1),T_mean,linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),T_mean+T_std,':',linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),T_mean-T_std,':',linewidth=lw)
		plt.legend()
		plt.savefig(fileID+'_cycle_period.png')

		plt.figure(4)
		# plt.plot(gp,np.mean(mean_hz)+np.std(mean_hz),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		# plt.plot(gp,np.mean(mean_hz)-np.std(mean_hz),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		# plt.plot(gp,np.mean(mean_hz),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_gp,max_gp,1),Hz_mean,linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),Hz_mean+Hz_std,':',linewidth=lw)
		# plt.plot(np.arange(min_gp,max_gp,1),Hz_mean-Hz_std,':',linewidth=lw)
		plt.legend()
		plt.savefig(fileID+'_mean_hz.png')

		plt.figure(5)
		plt.axis([np.min(CytNa)-2.5,np.max(CytNa)+5,0,1.1])
		plt.title(r'Na/K pump current and Cytosol Sodium $I_{pump}^{max}$='+str(IpumpMax))
		plt.xlabel(r'$[Na^+]_{in}$ nM')
		plt.ylabel(r'Activation of $I_{pump}$ nA')
		plt.legend()
		plt.savefig(fileID+'_CytNa_Ipump.png')

		plt.ion()
		plt.show()

	####################################################################################
	################################ pump sweep ###########################################
	if control == 'G':
		gp = list10[-2]
		list1=list1[::-1] #reverse the order of the list, for plotting purposes
		p1 = len(list1)

		BD_mean = np.zeros((1,p1))
		BD_mean = BD_mean[0]
		BD_std =  np.zeros((1,p1))
		BD_std = BD_std[0]

		IBI_mean = np.zeros((1,p1))
		IBI_mean = IBI_mean[0]
		IBI_std = np.zeros((1,p1))
		IBI_std = IBI_std[0]

		T_mean = np.zeros((1,p1))
		T_mean = T_mean[0]
		T_std =  np.zeros((1,p1))
		T_std = T_std[0]

		Hz_mean = np.zeros((1,p1))
		Hz_mean = Hz_mean[0]
		Hz_std = np.zeros((1,p1))
		Hz_std = Hz_std[0]

		for i in range(0,np.size(list1)):
			j=list1[i]
			IpumpMax = j
			# p2 = int(np.round(10*(IpumpMax-min_pump)))
			p2 = i

			str1 = 'BD_IpumpMax='+str(IpumpMax)+'_gp='
			str2 = 'IBI_IpumpMax='+str(IpumpMax)+'_gp='
			str3 = 'T_IpumpMax='+str(IpumpMax)+'_gp='
			str4 = 'Hz_IpumpMax='+str(IpumpMax)+'_gp='
			strend = '.npy'

			#BD
			filename1 = fileID+'_'+str1+str(gp)+strend
			burst_duration = np.load(filename1)

			x = np.zeros((1,np.size(burst_duration)))
			x=x[0]+IpumpMax

			plt.figure(1)
			plt.plot(IpumpMax,np.mean(burst_duration)+np.std(burst_duration),'^',markersize=mk,color='green')
			plt.plot(IpumpMax,np.mean(burst_duration)-np.std(burst_duration),'v',markersize=mk,color='red')
			plt.plot(x,burst_duration,'.',markersize=mk-10)
			plt.plot(IpumpMax,np.mean(burst_duration),'*',markersize=mk,color='black')
			plt.title(r'Burst duration $g_p$='+str(gp))
			plt.ylabel('time(s)')
			plt.xlabel(r'Na,K Pump Strength $I_{pump}^{max}$ (nA)')

			BD_mean[i] = np.mean(burst_duration)
			BD_std[i] = np.std(burst_duration)

			#IBI
			filename2 = fileID+'_'+str2+str(gp)+strend
			interburst_interval = np.load(filename2)

			x = np.zeros((1,np.size(interburst_interval)))
			x=x[0]+IpumpMax

			plt.figure(2)
			plt.plot(IpumpMax,np.mean(interburst_interval)+np.std(interburst_interval),'^',markersize=mk,color='green')
			plt.plot(IpumpMax,np.mean(interburst_interval)-np.std(interburst_interval),'v',markersize=mk,color='red')
			plt.plot(x,interburst_interval,'.',markersize=mk-10)
			plt.plot(IpumpMax,np.mean(interburst_interval),'*',markersize=mk,color='black')
			plt.title(r'Interburst interval $g_p$='+str(gp))
			plt.ylabel('time(s)')
			plt.xlabel(r'Na,K Pump Strength $I_{pump}^{max}$ (nA)')

			IBI_mean[i] = np.mean(interburst_interval)
			IBI_std[i] = np.std(interburst_interval)

			#T
			filename3 = fileID+'_'+str3+str(gp)+strend
			cycle_period = np.load(filename3)

			x = np.zeros((1,np.size(cycle_period)))
			x=x[0]+IpumpMax

			plt.figure(3)
			plt.plot(IpumpMax,np.mean(cycle_period)+np.std(cycle_period),'^',markersize=mk,color='green')
			plt.plot(IpumpMax,np.mean(cycle_period)-np.std(cycle_period),'v',markersize=mk,color='red')
			plt.plot(x,cycle_period,'.',markersize=mk-10)
			plt.plot(IpumpMax,np.mean(cycle_period),'*',markersize=mk,color='black')
			plt.title(r'Cycle period $g_p$='+str(gp))
			plt.ylabel('time(s)')
			plt.xlabel(r'Na,K Pump Strength $I_{pump}^{max}$ (nA)')

			T_mean[i] = np.mean(cycle_period)
			T_std[i] = np.std(cycle_period)

			#Hz
			filename4 = fileID+'_'+str4+str(gp)+strend
			mean_hz = np.load(filename4)

			x = np.zeros((1,np.size(mean_hz)))
			x=x[0]+IpumpMax

			plt.figure(4)
			plt.plot(IpumpMax,np.mean(mean_hz)+np.std(mean_hz),'^',markersize=mk,color='green')
			plt.plot(IpumpMax,np.mean(mean_hz)-np.std(mean_hz),'v',markersize=mk,color='red')
			plt.plot(x,mean_hz,'.',markersize=mk-10)
			plt.plot(IpumpMax,np.mean(mean_hz),'*',markersize=mk,color='black')

			plt.title(r'$\mu$ intrabust spiking frequency $g_p$='+str(gp))
			plt.ylabel('Hz')
			plt.xlabel(r'Na,K Pump Strength $I_{pump}^{max}$ (nA)')

			Hz_mean[i] = np.mean(mean_hz)
			Hz_std[i] = np.std(mean_hz)


			# Pump Na
			shift = 0.0
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Ipump = Ipump/IpumpMax
			CytNa = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			plt.figure(5)
			# plt.plot(CytNa,Ipump, linewidth=lw/2,label='gp='+str(gp))
			plt.plot(CytNa,Ipump, marker=markers[i], color=colores[i],linewidth=lw,linestyle='', markersize=mk2 ,label=r'$I_{pump}^{max}$='+str(IpumpMax))
			plt.plot(np.mean(CytNa),np.min(Ipump), '*',markersize=mk, color=colores[i])
			plt.plot(np.min(CytNa),np.min(Ipump), marker=markers[i], markersize=mk, color=colores[i])
			plt.plot(np.max(CytNa),np.min(Ipump), marker=markers[i], markersize=mk, color=colores[i])
			# plt.plot(np.mean(CytNa)+np.std(CytNa),np.mean(Ipump)+np.std(Ipump) +shift, '^',markersize=mk, color='green')
			# plt.plot(np.mean(CytNa)-np.std(CytNa),np.mean(Ipump)-np.std(Ipump) +shift, 'v',markersize=mk, color='red')

		plt.figure(1)
		plt.plot(IpumpMax,np.mean(burst_duration)+np.std(burst_duration),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		plt.plot(IpumpMax,np.mean(burst_duration)-np.std(burst_duration),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		plt.plot(IpumpMax,np.mean(burst_duration),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_pump,max_pump,0.1),BD_mean,linewidth=lw)
		# plt.plot(np.arange(min_pump,max_pump,0.1),BD_mean+BD_std,':',linewidth=lw)
		# plt.plot(np.arange(min_pump,max_pump,0.1),BD_mean-BD_std,':',linewidth=lw)
		plt.plot(list1,BD_mean,linewidth=lw)
		plt.plot(list1,BD_mean+BD_std,':',linewidth=lw)
		plt.plot(list1,BD_mean-BD_std,':',linewidth=lw)

		plt.legend()
		plt.savefig(fileID+'_burst_duration.png')

		plt.figure(2)
		plt.plot(IpumpMax,np.mean(interburst_interval)+np.std(interburst_interval),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		plt.plot(IpumpMax,np.mean(interburst_interval)-np.std(interburst_interval),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		plt.plot(IpumpMax,np.mean(interburst_interval),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_pump,max_pump,0.1),IBI_mean,linewidth=lw)
		# plt.plot(np.arange(min_pump,max_pump,0.1),IBI_mean+IBI_std,':',linewidth=lw)
		# plt.plot(np.arange(min_pump,max_pump,0.1),IBI_mean-IBI_std,':',linewidth=lw)
		plt.plot(list1,IBI_mean,linewidth=lw)
		plt.plot(list1,IBI_mean+IBI_std,':',linewidth=lw)
		plt.plot(list1,IBI_mean-IBI_std,':',linewidth=lw)

		plt.legend()
		plt.savefig(fileID+'_interburst_interval.png')

		plt.figure(3)
		plt.plot(IpumpMax,np.mean(cycle_period)+np.std(cycle_period),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		plt.plot(IpumpMax,np.mean(cycle_period)-np.std(cycle_period),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		plt.plot(IpumpMax,np.mean(cycle_period),'*',markersize=mk,color='black',label=r'$\mu$')
		# plt.plot(np.arange(min_pump,max_pump,0.1),T_mean,linewidth=lw)
		# plt.plot(np.arange(min_pump,max_pump,0.1),T_mean+T_std,':',linewidth=lw)
		# plt.plot(np.arange(min_pump,max_pump,0.1),T_mean-T_std,':',linewidth=lw)
		plt.plot(list1,T_mean,linewidth=lw)
		plt.plot(list1,T_mean+T_std,':',linewidth=lw)
		plt.plot(list1,T_mean-T_std,':',linewidth=lw)

		plt.legend()
		plt.savefig(fileID+'_cycle_period.png')

		plt.figure(4)
		plt.plot(IpumpMax,np.mean(mean_hz)+np.std(mean_hz),'^',markersize=mk,color='green',label=r'$\mu + \sigma$')
		plt.plot(IpumpMax,np.mean(mean_hz)-np.std(mean_hz),'v',markersize=mk,color='red',label=r'$\mu - \sigma$')
		plt.plot(IpumpMax,np.mean(mean_hz),'*',markersize=mk,color='black',label=r'$\mu$')
		plt.plot(list1,Hz_mean,linewidth=lw)
		plt.plot(list1,Hz_mean+Hz_std,':',linewidth=lw)
		plt.plot(list1,Hz_mean-Hz_std,':',linewidth=lw)
		plt.legend()
		plt.savefig(fileID+'_mean_hz.png')

		plt.figure(5)
		plt.title(r'Na/K pump current and Cytosol Sodium $g_p$='+str(gp))
		plt.xlabel(r'$[Na^+]_{in}$ nM')
		plt.ylabel(r'Activation of $I_{pump}$ nA')
		plt.legend()
		plt.savefig(fileID+'_CytNa_Ipump.png')

		plt.ion()
		plt.show()



##### file ID stuff
cwd = os.getcwd()
fileID= cwd[-8:len(cwd)]

############## MANIPULATE to modify behavior of script, later the code should be more interactive and receive these parameters as user-input in the command-line
list1 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7]
list2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7]
list3 = [0.4,0.5,0.7,0.8,0.9]

list4 = [4.0,5.0,6.0]
list5 = [4.0,5.0,6.0]
list6 = [4.0,5.0,6.0]

I_bias = 0.
#list1
gp = 4.0			#used exclusively for Gp sweep
# IpumpMax = 0.3
control = 'G'

if control=='I':
	list1.append(IpumpMax)
	cntrl = IpumpMax
elif control =='G':
	list1.append(gp)
	cntrl = gp
list1.append(control)

#list2
gp = 5.0
# IpumpMax = 0.75
control = 'G'

if control=='I':
	list2.append(IpumpMax)
	cntrl = IpumpMax
elif control =='G':
	list2.append(gp)
	cntrl = gp
list2.append(control)

#list3
# IpumpMax = 0.8
gp = 6.0
control = 'G'

if control=='I':
	list3.append(IpumpMax)
	cntrl = IpumpMax
elif control =='G':
	list3.append(gp)
	cntrl = gp
list3.append(control)

#list4
IpumpMax = 0.4
# gp = 6.0
control = 'I'

if control=='I':
	list4.append(IpumpMax)
	cntrl = IpumpMax
elif control =='G':
	list4.append(gp)
	cntrl = gp
list4.append(control)

#list5
IpumpMax = 0.5
# gp = 6.0
control = 'I'

if control=='I':
	list5.append(IpumpMax)
	cntrl = IpumpMax
elif control =='G':
	list5.append(gp)
	cntrl = gp
list5.append(control)

#list6
IpumpMax = 0.7
# gp = 6.0
control = 'I'

if control=='I':
	list6.append(IpumpMax)
	cntrl = IpumpMax
elif control =='G':
	list6.append(gp)
	cntrl = gp
list6.append(control)


############			save data parameters for subsequent use
np.save(fileID+'_list1',list1)
np.save(fileID+'_list2',list2)
np.save(fileID+'_list3',list3)

np.save(fileID+'_list4',list4)
np.save(fileID+'_list5',list5)
np.save(fileID+'_list6',list6)

np.save(fileID+'_Ibias',I_bias)

param0 = ['list1','list2','list3','list4','list5','list6']
np.save(fileID+'_param',param0)
################################################### 1
plt.figure(1,figsize=(45,30))
# IBI
plt.figure(2,figsize=(45,30))
# T
plt.figure(3,figsize=(45,30))
# Hz
plt.figure(4,figsize=(45,30))
# Pump vs Na_in
plt.figure(5,figsize=(45,30))


compute(list1)

# plt.close('all')
compute(list2)
compute(list3)
# ################################################### 2
plt.figure(6,figsize=(45,30))
# IBI
plt.figure(7,figsize=(45,30))
# T
plt.figure(8,figsize=(45,30))
# Hz
plt.figure(9,figsize=(45,30))
# Pump vs Na_in
plt.figure(10,figsize=(45,30))
compute(list4)
compute(list5)
compute(list6)

# plt.close('all')


# ################################################### 3
# plt.figure(1,figsize=(45,30))
# # IBI
# plt.figure(2,figsize=(45,30))
# # T
# plt.figure(3,figsize=(45,30))
# # Hz
# plt.figure(4,figsize=(45,30))
# # Pump vs Na_in
# plt.figure(5,figsize=(45,30))
# compute(list3)

# plt.close('all')