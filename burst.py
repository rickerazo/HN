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
        'size'   : 70}
plt.rc('font', **font)
mk = 20
mk2 = 8
lw= 15
alpha0 = 1
markers = ['1', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia']

# The next section analyzes the data files specified by the max and min parameters above. There is one section for gp sweeps and another for Ipumpmax sweeps.
def compute(list10):
	control = list10[-1]
	list1 = list10[0:-5]
	protocol = list10[-4]
	coupling = list10[-3]
	# print(protocol,coupling,list1)
	################pl####################################################################	
	############################# Sodium sweep ###########################################	
	if control == 'I':
		IpumpMax = list10[-5]	
		list0=list1[::-1] #reverse the order of the list, for plotting purposes

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
			x_all = int(gp)

			# File ID stuff. It must be after parameters are set
			str1 = 'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			str2 = 'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			str3 = 'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			str4 = 'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			strend = '.npy'


			#BD
			filename1 = fileID+'_'+str1+strend
			burst_duration = np.load(filename1)

			x = np.zeros((1,np.size(burst_duration)))
			x=x[0]+gp

			plt.figure(1)
			# plt.plot([gp,gp],[np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration)+np.std(burst_duration)])
			plt.plot(gp,np.mean(burst_duration)+np.std(burst_duration),'^',markersize=mk,color='green')
			plt.plot(gp,np.mean(burst_duration)-np.std(burst_duration),'v',markersize=mk,color='red')
			plt.plot(x,burst_duration,'.',markersize=mk-10)
			plt.plot(gp,np.mean(burst_duration),'*',markersize=mk,color='black')
			plt.title(r'Burst duration $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('time(s)')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			BD_mean[i] = np.mean(burst_duration)
			BD_std[i] = np.std(burst_duration)

			#IBI
			filename2 = fileID+'_'+str2+strend
			interburst_interval = np.load(filename2)

			x = np.zeros((1,np.size(interburst_interval)))
			x=x[0]+gp
			plt.figure(2)
			plt.plot(gp,np.mean(interburst_interval)+np.std(interburst_interval),'^',markersize=mk,color='green')
			plt.plot(gp,np.mean(interburst_interval)-np.std(interburst_interval),'v',markersize=mk,color='red')
			plt.plot(x,interburst_interval,'.',markersize=mk-10)
			plt.plot(gp,np.mean(interburst_interval),'*',markersize=mk,color='black')
			plt.title(r'Interburst interval $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('time(s)')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			IBI_mean[i] = np.mean(interburst_interval)
			IBI_std[i] = np.std(interburst_interval)

			#T
			filename3 = fileID+'_'+str3+strend
			cycle_period = np.load(filename3)

			x = np.zeros((1,np.size(cycle_period)))
			x=x[0]+gp

			plt.figure(3)
			plt.plot(gp,np.mean(cycle_period)+np.std(cycle_period),'^',markersize=mk,color='green')
			plt.plot(gp,np.mean(cycle_period)-np.std(cycle_period),'v',markersize=mk,color='red')
			plt.plot(x,cycle_period,'.',markersize=mk-10)
			plt.plot(gp,np.mean(cycle_period),'*',markersize=mk,color='black')
			plt.title(r'Cycle period $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('time(s)')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			T_mean[i] = np.mean(cycle_period)
			T_std[i] = np.std(cycle_period)

			#Hz
			filename4 = fileID+'_'+str4+strend
			mean_hz = np.load(filename4)

			x = np.zeros((1,np.size(mean_hz)))
			x=x[0]+gp

			plt.figure(4)
			plt.plot(gp,np.mean(mean_hz)+np.std(mean_hz),'^',markersize=mk,color='green')
			plt.plot(gp,np.mean(mean_hz)-np.std(mean_hz),'v',markersize=mk,color='red')
			plt.plot(x,mean_hz,'.',markersize=mk-10)
			plt.plot(gp,np.mean(mean_hz),'*',markersize=mk,color='black')

			plt.title(r'$\mu$ intrabust spiking frequency $I_{pump}^{max}$='+str(IpumpMax))
			plt.ylabel('Hz')
			plt.xlabel(r'Maximal conductance of $I_p$ $g_p$')

			Hz_mean[i] = np.mean(mean_hz)
			Hz_std[i] = np.std(mean_hz)

			# Pump Na
			shift = 0.01
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			Ipump = Ipump/IpumpMax
			CytNa = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			plt.figure(5)
			plt.plot(CytNa,Ipump, marker='o', color=colores[x_all],linewidth=lw,linestyle='', markersize=mk2 ,label=r'$g_p$='+str(gp))
			plt.plot([np.min(CytNa),np.max(CytNa)],[np.min(Ipump),np.min(Ipump)],color=colores[x_all],linewidth=3)
			plt.plot([np.min(CytNa),np.max(CytNa)],[np.max(Ipump),np.max(Ipump)],color=colores[x_all],linewidth=3)
			# plt.plot(np.mean(CytNa),np.min(Ipump), '*',markersize=mk, color=colores[i])
			# plt.plot(np.min(CytNa),np.min(Ipump), marker='o', markersize=mk, color=colores[i])
			# plt.plot(np.max(CytNa),np.min(Ipump), marker='o', markersize=mk, color=colores[i])

		plt.figure(1)
		plt.plot(list1,BD_mean,linewidth=lw, label = protocol)
		plt.plot(list1,BD_mean+BD_std,':',linewidth=lw)
		plt.plot(list1,BD_mean-BD_std,':',linewidth=lw)
		plt.legend()
		plt.savefig(fileID+'_burst_duration.png')

		plt.figure(2)
		plt.plot(list1,IBI_mean,linewidth=lw, label = protocol)
		plt.plot(list1,IBI_mean+IBI_std,':',linewidth=lw)
		plt.plot(list1,IBI_mean-IBI_std,':',linewidth=lw)		
		plt.legend()
		plt.savefig(fileID+'_interburst_interval.png')

		plt.figure(3)
		plt.plot(list1,T_mean,linewidth=lw, label = protocol)
		plt.plot(list1,T_mean+T_std,':',linewidth=lw)
		plt.plot(list1,T_mean-T_std,':',linewidth=lw)		
		plt.legend()
		plt.savefig(fileID+'_cycle_period.png')

		plt.figure(4)
		plt.plot(list1,Hz_mean,linewidth=lw, label = protocol)
		plt.plot(list1,Hz_mean+Hz_std,':',linewidth=lw)
		plt.plot(list1,Hz_mean-Hz_std,':',linewidth=lw)		
		plt.legend()
		plt.savefig(fileID+'_mean_hz.png')

		plt.figure(5)
		# plt.axis([np.min(CytNa)-2.5,np.max(CytNa)+5,0,1.1])
		plt.xlim(10,26)
		plt.ylim(0.,1.1)
		plt.title(r'$Na^+/K^+$ and $[Na^+]_{cyt}$$I_{pump}^{max}$='+str(IpumpMax))
		plt.xlabel(r'$[Na^+]_{in}$ nM')
		plt.ylabel(r'Activation of $I_{pump}$ nA')
		plt.legend(markerscale=5,loc=(1.2,1))
		plt.savefig(fileID+'_CytNa_Ipump.png')

		plt.ion()
		plt.show()

	####################################################################################
	################################ pump sweep ###########################################
	if control == 'G':
		gp = list10[-5]
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
			# x_all = int(IpumpMax*10)
			# p2 = int(np.round(10*(IpumpMax-min_pump)))
			p2 = i
			x_all = int(IpumpMax*10)

			#file handles
			str1 = 'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			str2 = 'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			str3 = 'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			str4 = 'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)
			strend = '.npy'

			#BD
			filename1 = fileID+'_'+str1+strend
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
			filename2 = fileID+'_'+str2+strend
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
			filename3 = fileID+'_'+str3+strend
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
			filename4 = fileID+'_'+str4+strend
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
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			Ipump = Ipump/IpumpMax
			CytNa = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			plt.figure(5)
			# plt.plot(CytNa,Ipump, linewidth=lw/2,label='gp='+str(gp))
			plt.plot(CytNa,Ipump, marker='o', color=colores[x_all],linewidth=lw,linestyle='', markersize=mk2 ,label=r'$I_{pump}^{max}$='+str(IpumpMax))
			# plt.plot(np.mean(CytNa),np.min(Ipump), '*',markersize=mk, color=colores[x_all])
			plt.plot([np.min(CytNa),np.max(CytNa)],[np.min(Ipump),np.min(Ipump)],color=colores[x_all],linewidth=3)
			plt.plot([np.min(CytNa),np.max(CytNa)],[np.max(Ipump),np.max(Ipump)],color=colores[x_all],linewidth=3)

			# plt.plot(np.min(CytNa),np.min(Ipump), marker=markers[i], markersize=mk, color=colores[x_all])
			# plt.plot(np.max(CytNa),np.min(Ipump), marker=markers[i], markersize=mk, color=colores[x_all])
			# plt.plot(np.mean(CytNa)+np.std(CytNa),np.mean(Ipump)+np.std(Ipump) +shift, '^',markersize=mk, color='green')
			# plt.plot(np.mean(CytNa)-np.std(CytNa),np.mean(Ipump)-np.std(Ipump) +shift, 'v',markersize=mk, color='red')

		plt.figure(1)
		plt.plot(IpumpMax,np.mean(burst_duration)+np.std(burst_duration),'^',markersize=mk,color='green')
		plt.plot(IpumpMax,np.mean(burst_duration)-np.std(burst_duration),'v',markersize=mk,color='red')
		plt.plot(IpumpMax,np.mean(burst_duration),'*',markersize=mk,color='black')

		plt.plot(list1,BD_mean,linewidth=lw, label = protocol)
		plt.plot(list1,BD_mean+BD_std,':',linewidth=lw)
		plt.plot(list1,BD_mean-BD_std,':',linewidth=lw)

		plt.legend()
		plt.savefig(fileID+'_burst_duration.png')

		plt.figure(2)
		plt.plot(IpumpMax,np.mean(interburst_interval)+np.std(interburst_interval),'^',markersize=mk,color='green')
		plt.plot(IpumpMax,np.mean(interburst_interval)-np.std(interburst_interval),'v',markersize=mk,color='red')
		plt.plot(IpumpMax,np.mean(interburst_interval),'*',markersize=mk,color='black')

		plt.plot(list1,IBI_mean,linewidth=lw, label = protocol)
		plt.plot(list1,IBI_mean+IBI_std,':',linewidth=lw)
		plt.plot(list1,IBI_mean-IBI_std,':',linewidth=lw)

		plt.legend()
		plt.savefig(fileID+'_interburst_interval.png')

		plt.figure(3)
		plt.plot(IpumpMax,np.mean(cycle_period)+np.std(cycle_period),'^',markersize=mk,color='green')
		plt.plot(IpumpMax,np.mean(cycle_period)-np.std(cycle_period),'v',markersize=mk,color='red')
		plt.plot(IpumpMax,np.mean(cycle_period),'*',markersize=mk,color='black')

		plt.plot(list1,T_mean,linewidth=lw, label = protocol)
		plt.plot(list1,T_mean+T_std,':',linewidth=lw)
		plt.plot(list1,T_mean-T_std,':',linewidth=lw)

		plt.legend()
		plt.savefig(fileID+'_cycle_period.png')

		plt.figure(4)
		plt.plot(IpumpMax,np.mean(mean_hz)+np.std(mean_hz),'^',markersize=mk,color='green')
		plt.plot(IpumpMax,np.mean(mean_hz)-np.std(mean_hz),'v',markersize=mk,color='red')
		plt.plot(IpumpMax,np.mean(mean_hz),'*',markersize=mk,color='black')
		plt.plot(list1,Hz_mean,linewidth=lw, label = protocol)
		plt.plot(list1,Hz_mean+Hz_std,':',linewidth=lw)
		plt.plot(list1,Hz_mean-Hz_std,':',linewidth=lw)
		plt.legend(loc=(1.2,1))
		plt.savefig(fileID+'_mean_hz.png')

		plt.figure(5)
		plt.title(r'$Na^+/K^+$ and $[Na^+]_{cyt}$ $\bar{g}_P$='+str(gp))
		plt.ylim(0,1.1)
		plt.xlabel(r'$[Na^+]_{in}$ nM')
		plt.ylabel(r'Activation of $I_{pump}$ nA')
		plt.legend(markerscale=5, loc=(1.2,1))
		plt.savefig(fileID+'_CytNa_Ipump.png')

		plt.ion()
		plt.show()


def append_list(control,list1,x1,x2,x3,x4):
	list1.append(x1)
	list1.append(x2)
	list1.append(x3)
	list1.append(x4)
	list1.append(control)
	return list1

##### file ID stuff
cwd = os.getcwd()
fileID= cwd[-8:len(cwd)]

############## MANIPULATE to modify behavior of script, later the code should be more interactive and receive these parameters as user-input in the command-line
list1 = [0.0] #gp = 7

I_bias = 0.1

gp = 0.0
control = 'C'

#### protocol settings : virtually any protocol that can be described in a two letter code word
protocol = 'nm'       # normal saline
# protocol = 'cs'         # Cs+ saline
# protocol = 'mm'       # myomodulin


#### coupling settings: virtually any protocol that can be described in a two letter code word
# coupling = 'hco'          # hybrid HCO neural network configuration
# coupling = 'bio'          # biological HCO : injected pump current to one neuron
coupling = 'one'            # single, isolated neuron


list10 = append_list(control, list1, gp, protocol, coupling, I_bias)
print(list10)
#######################################################
list2 = [0.1, 0.2]

I_bias = 0.1

gp = 0.0
control = 'G'

#### protocol settings : virtually any protocol that can be described in a two letter code word
protocol = 'nm'       # normal saline
# protocol = 'cs'         # Cs+ saline
# protocol = 'mm'       # myomodulin


#### coupling settings: virtually any protocol that can be described in a two letter code word
# coupling = 'hco'          # hybrid HCO neural network configuration
# coupling = 'bio'          # biological HCO : injected pump current to one neuron
coupling = 'one'            # single, isolated neuron


list20 = append_list(control, list2, gp, protocol, coupling, I_bias)
print(list20)

#######################################################
list3 = [0.2,0.4]

I_bias = 0.1

gp = 3.0
control = 'G'

#### protocol settings : virtually any protocol that can be described in a two letter code word
protocol = 'nm'       # normal saline
# protocol = 'cs'         # Cs+ saline
# protocol = 'mm'       # myomodulin


#### coupling settings: virtually any protocol that can be described in a two letter code word
# coupling = 'hco'          # hybrid HCO neural network configuration
# coupling = 'bio'          # biological HCO : injected pump current to one neuron
coupling = 'one'            # single, isolated neuron


list30 = append_list(control, list3, gp, protocol, coupling, I_bias)
print(list30)

#######################################################
list4 = [0.1, 0.4, 0.6, 0.7]

I_bias = 0.0

gp = 9.0
control = 'G'

#### protocol settings : virtually any protocol that can be described in a two letter code word
protocol = 'nm'       # normal saline
# protocol = 'cs'         # Cs+ saline
# protocol = 'mm'       # myomodulin


#### coupling settings: virtually any protocol that can be described in a two letter code word
# coupling = 'hco'          # hybrid HCO neural network configuration
# coupling = 'bio'          # biological HCO : injected pump current to one neuron
coupling = 'one'            # single, isolated neuron


list40 = append_list(control, list4, gp, protocol, coupling, I_bias)
print(list40)

#######################################################
list5 = [0.1,0.2,0.3, 0.4, 0.6, 0.7, 0.9]

I_bias = 0.0

gp = 9.0
control = 'G'

#### protocol settings : virtually any protocol that can be described in a two letter code word
# protocol = 'nm'       # normal saline
protocol = 'cs'         # Cs+ saline
# protocol = 'mm'       # myomodulin


#### coupling settings: virtually any protocol that can be described in a two letter code word
# coupling = 'hco'          # hybrid HCO neural network configuration
# coupling = 'bio'          # biological HCO : injected pump current to one neuron
coupling = 'one'            # single, isolated neuron


list50 = append_list(control, list5, gp, protocol, coupling, I_bias)
print(list50)


############			save data parameters for subsequent use
np.save(fileID+'_list1',list10)
np.save(fileID+'_list2',list20)
np.save(fileID+'_list3',list30)
np.save(fileID+'_list4',list40)
np.save(fileID+'_list5',list50)

# np.save(fileID+'_Ibias',I_bias)

param0 = ['list1','list2', 'list3', 'list4', 'list5']
np.save(fileID+'_param',param0)
################################################### 1
plt.figure(1,figsize=(30,25))
# IBI
plt.figure(2,figsize=(30,25))
# T
plt.figure(3,figsize=(30,25))
# Hz
plt.figure(4,figsize=(30,25))
# Pump vs Na_in
plt.figure(5,figsize=(30,25))

compute(list10)
compute(list20)
compute(list30)
compute(list40)
compute(list50)