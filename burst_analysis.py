####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps
# BURST analysis ->BD, IBI, Hz, T, Ipump_CytNa

# Selective data -> coefficient of variation = 25%

#### libraries necessary for code
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import itertools
import os
import sys
####	graphix stuff
matplotlib.rcParams['axes.linewidth']=10


markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['k','brown','b','g','c','gray','orange','r','violet','']
# colores = ['','#98FB98','#50C878','']
colores = colores[::-1]
# list1 = list0[::-1]

font = {'weight' : 'bold',
        'size'   : 70}
plt.rc('font', **font)
mk1 = 16
mk2 = 15
mk3 = 3
mk10 = 25
lw1= 6
alpha0 = 1


# coefficient of variance standard: reduce alpha
coef_std1 = 0.25
# coefficient of variane absolute standard: do not plot
covar = 0.3
length_standard = 3
for i in range(1,9):
	plt.figure(i,figsize=(30,20))

##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
experiment_list = np.array([
	# 18831003,
	# 18902004,
	# 18914000,
	# 18929002,
	19522001,
	19529000,
	# # 19603000,
	# # 19607000,
	# # 19607001,
	# # 19612001,
	19612002,
	19614000,
	19625000,
	19626000,
	19626002,
	# 19731000,
	19731001,
	19805000,
	19809001,
	19826000,
	19918002,
	19919000
	])
np.save('experiment_list',experiment_list)

def coefficient(coefficient_of_variation, covar,coef_std1):
	if coefficient_of_variation < covar:
		if coefficient_of_variation>coef_std1:
			alpha2 = 0.6
			fst1 = 'none'
		else:
			alpha2 = 1
			fst1 = 'full'
	elif coefficient_of_variation >= covar and coefficient_of_variation < 0.5:
		alpha2 = 0.2
		fst1 = 'none'
	else:
		alpha2 = 0.001
		fst1 = 'none'
	return alpha2, fst1

def compute(list10):		
	constant1 = list10[-1]
	list0 = list10[0:-2]
	pump_list1 = []
	pump_list2 = []
	pump_list3 = []
	pump_list4 = []
	# pump
	gp_list1 = []
	gp_list2 =[]
	gp_list3 =[]
	gp_list4 =[]

	BD_mean_list = []
	IBI_mean_list = []
	T_mean_list = []
	Hz_mean_list=[]	
	Ibias = np.load(fh+'Ibias.npy')
	if Ibias>0.1:
		constant1='Z'

	if constant1=='I':
		IpumpMax = list10[-2]
		IpumpMax = float(IpumpMax)

		list1 = list0[::-1]

		control_param = IpumpMax
		for j in range(0,len(list1)):
			gp = float(list1[j])
			x_all = int(control_param*10)

			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)
			
			plt.figure(1)
			y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
			x = np.zeros((1,len(y1)))+gp
			plt.plot(x[0],y1,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(gp,np.mean(burst_duration),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(burst_duration)+np.std(burst_duration),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(burst_duration)-np.std(burst_duration),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			gp_list1.append(gp)
			BD_mean_list.append(np.mean(burst_duration))#/np.max(burst_duration))

			#IBI
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)
			plt.figure(2)
			y2 = np.array((np.mean(interburst_interval)-np.std(interburst_interval),np.mean(interburst_interval),np.mean(interburst_interval)+np.std(interburst_interval)))
			x = np.zeros((1,len(y2)))+gp
			plt.plot(x[0],y2,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(gp,np.mean(interburst_interval),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(interburst_interval)+np.std(interburst_interval),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(interburst_interval)-np.std(interburst_interval),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			gp_list2.append(gp)
			IBI_mean_list.append(np.mean(interburst_interval))#/np.max(cycle_period))


			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)
			plt.figure(3)
			y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
			x = np.zeros((1,len(y3)))+gp
			plt.plot(x[0],y3,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(gp,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(cycle_period)+np.std(cycle_period),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(cycle_period)-np.std(cycle_period),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			gp_list3.append(gp)
			T_mean_list.append(np.mean(cycle_period))#/np.max(cycle_period))

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)
			plt.figure(4)
			y4 = np.array((np.mean(Hz)-np.std(Hz), np.mean(Hz) ,np.mean(Hz)+np.std(Hz) ))
			x = np.zeros((1,len(y4)))+gp
			plt.plot(x[0],y4,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(gp,np.mean(Hz),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(Hz)+np.std(Hz),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(gp,np.mean(Hz)-np.std(Hz),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			gp_list4.append(gp)
			Hz_mean_list.append(np.mean(Hz))



		plt.figure(1)
		plt.plot(gp_list1, BD_mean_list,color=colores[x_all])
		plt.figure(2)
		plt.plot(gp_list2, IBI_mean_list,color=colores[x_all])
		plt.figure(3)
		plt.plot(gp_list3, T_mean_list,color=colores[x_all])
		plt.figure(4)
		plt.plot(gp_list4, Hz_mean_list,color=colores[x_all])



	if constant1=='G':
		gp = list10[-2]
		gp = float(gp)
		control_param=gp

		list1 = list0[::-1]

		for j in range(0,len(list1)):
			IpumpMax = list1[j]
			IpumpMax = float(IpumpMax)
			x_all = int(control_param)

			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			plt.figure(5)
			y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
			x = np.zeros((1,len(y1)))+IpumpMax
			plt.plot(x[0],y1,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(IpumpMax,np.mean(burst_duration),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(burst_duration)+np.std(burst_duration),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(burst_duration)-np.std(burst_duration),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			pump_list1.append(IpumpMax)
			BD_mean_list.append(np.mean(burst_duration))

			#IBI
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			plt.figure(6)
			y2 = np.array((np.mean(interburst_interval)-np.std(interburst_interval),np.mean(interburst_interval),np.mean(interburst_interval)+np.std(interburst_interval)))
			x = np.zeros((1,len(y2)))+IpumpMax
			plt.plot(x[0],y2,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(IpumpMax,np.mean(interburst_interval),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(interburst_interval)+np.std(interburst_interval),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(interburst_interval)-np.std(interburst_interval),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			pump_list2.append(IpumpMax)
			IBI_mean_list.append(np.mean(interburst_interval))

			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			plt.figure(7)
			y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
			x = np.zeros((1,len(y3)))+IpumpMax
			plt.plot(x[0],y3,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(IpumpMax,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(cycle_period)+np.std(cycle_period),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(cycle_period)-np.std(cycle_period),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])

			pump_list3.append(IpumpMax)
			T_mean_list.append(np.mean(cycle_period))


			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			plt.figure(8)
			y4 = np.array((np.mean(Hz)-np.std(Hz),np.mean(Hz),np.mean(Hz)+np.std(Hz)))
			x = np.zeros((1,len(y4)))+IpumpMax
			plt.plot(x[0],y4,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plt.plot(IpumpMax,np.mean(Hz),marker='o',markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(Hz)+np.std(Hz),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plt.plot(IpumpMax,np.mean(Hz)-np.std(Hz),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			pump_list4.append(IpumpMax)
			Hz_mean_list.append(np.mean(Hz))

		plt.figure(5)
		plt.plot(pump_list1, BD_mean_list,color=colores[x_all])
		plt.figure(6)
		plt.plot(pump_list2, IBI_mean_list,color=colores[x_all])
		plt.figure(7)
		plt.plot(pump_list3, T_mean_list,color=colores[x_all])
		plt.figure(8)
		plt.plot(pump_list4, Hz_mean_list,color=colores[x_all])
	return float(list10[-2])
riff_gp = []
riff_pump = []
print('Parameters from experiments analyzed:')
for i in range(0,len(experiment_list)):
	
	fileID= str(experiment_list[i])
	# fh1 = cwd+fileID+'/'
	fh = fileID+'/'+fileID+'_'
	# print(cwd+str(experiment_list[i]))
	# sys.path.insert(0,fh1)
	params = list(np.load(fh+'param.npy'))
	# params = ['list1']
	for k in range(0,len(params)):# import file
		list20 = params[k]
		han2 = fh+str(list20)+'.npy'
		list10 =  list(np.load(han2))
		Ibias = np.load(fh+'Ibias.npy')
		# print(fileID,list10,Ibias)
		if Ibias<=0.1:
			print(fileID,list10,Ibias)
			control_param = compute(list10)
			if list10[-1] == 'I':
				riff_pump=np.append(riff_pump,control_param)
			if list10[-1] == 'G':
				riff_gp = np.append(riff_gp,int(control_param))


###############################################################################################################
###############################################################################################################

plt.figure(1)
plt.ylabel('Burst duration (s)',weight='bold',fontsize=80)
plt.xlabel(r'$\bar{g}_P$ (nS)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(1,12,3))
plt.ylim(0,11)
plt.xticks(ticks=np.arange(1,11,3))
plt.xlim(0.5,10.5)


plt.figure(2)
plt.ylabel('IBI (s)',weight='bold',fontsize=80)
plt.xlabel(r'$\bar{g}_P$ (nS)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(1,12,3))
plt.ylim(0,11)
plt.xticks(ticks=np.arange(1,11,3))
plt.xlim(0.5,10.5)


plt.figure(3)
plt.ylabel('Period (s)',weight='bold',fontsize=80)
plt.xlabel(r'$\bar{g}_P$ (nS)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(1,21,4))
plt.ylim(0,14)
plt.xticks(ticks=np.arange(1,11,3))
plt.xlim(0.5,10.5)


plt.figure(4)
plt.ylabel('Spike frequency (Hz)',weight='bold',fontsize=80)
plt.xlabel(r'$\bar{g}_P$ (nS)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(0,60,10))
plt.xticks(ticks=np.arange(1,11,3))
plt.xlim(0.5,10.5)
plt.ylim(0,45)


plt.figure(5)
plt.ylabel('Burst duration (s)',weight='bold',fontsize=80)
plt.xlabel(r'$I_{pump}^{max}$ (nA)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(1,12,3))
plt.ylim(0,12.5)
plt.xticks(ticks=np.arange(0.2,1.2,0.2))
plt.xlim(0.05,1)


plt.figure(6)
plt.ylabel('IBI (s)',weight='bold',fontsize=80)
plt.xlabel(r'$I_{pump}^{max}$ (nA)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(1,12,3))
plt.ylim(0,12.5)
plt.xticks(ticks=np.arange(0.2,1.2,0.2))
plt.xlim(0.05,1)


plt.figure(7)
plt.ylabel('Period (s)',weight='bold',fontsize=80)
plt.xlabel(r'$I_{pump}^{max}$ (nA)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(1,21,4))
plt.ylim(0,16)
plt.xticks(ticks=np.arange(0.2,1.2,0.2))
plt.xlim(0.05,1)


plt.figure(8)
plt.ylabel('Spike frequency (Hz)',weight='bold',fontsize=80)
plt.xlabel(r'$I_{pump}^{max}$ (nA)',weight='bold',fontsize=80)
plt.yticks(ticks=np.arange(0,60,10))
plt.ylim(0,50)
plt.xticks(ticks=np.arange(0.2,1.2,0.2))
plt.xlim(0.05,1)
#############################################################################################################

markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
# colores = ['k','b','brown','c','g','gray','m','orange','pink','violet','r','yellow']
# colores = colores[::-1]

riff2 = np.unique(riff_pump)

color_riff2= [colores[1],colores[2],colores[3],colores[4],colores[5],colores[6],colores[7],colores[8],colores[9]]#,colores[6],colores[8],colores[9]]
titles = [' ','BD','IBI','T','Hz','Ipump_cytNa','BD','IBI','T','Hz','Ipump_cytNa']
#gp 
for i in range(1,5):
	plt.figure(i)
	fig = plt.figure(i)
	plot = fig.add_subplot(111)
	for j in range(0,len(riff2)):
		p1 = riff2[j]
		str1 = str(p1)
		plt.plot(0,0,'o',markersize=1e-3,color=color_riff2[j],label=r'$I_{pump}^{max}$='+str1+' nA')
		plt.legend(loc=2,bbox_to_anchor=(-0.014,1.15),markerscale=2e+4,ncol=3,fontsize=45)
		# plt.legend(framealpha=0.25)
		axes = plt.axes()
		axes.set_frame_on(False)
		# plot.tick_params(axis='both',which='major',labelsize=20)
		plt.savefig('gp_'+titles[i])

#		4,5,5.5,6,6.5,7
riff1 = np.unique(riff_gp)
color_riff1 = [colores[1],colores[2],colores[3],colores[4],colores[5],colores[6], colores[7],colores[8],colores[9]]

#pump
for i in range(5,9):
	plt.figure(i)
	fig = plt.figure(i)
	plot = fig.add_subplot(111)

	for j in range(0,len(riff1)):
		p1 = riff1[j]
		str1 = str(p1)
		plt.plot(0,0,'o',markersize=1e-3,color=color_riff1[j],label=r'$\bar{g}_P$='+str1+' nS')
		plt.legend(loc=2,bbox_to_anchor=(-0.16,1.15),markerscale=2e+4,ncol=4,fontsize=45)
		# plt.legend(framealpha=0.25)
		axes = plt.axes()
		axes.set_frame_on(False)
		# plot.tick_params(axis='both',which='major',labelsize=20)
		plt.savefig('pump_'+titles[i-4])


# plt.ion()
# plt.show()
