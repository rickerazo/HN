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



font = {'weight' : 'bold',
        'size'   : 50}
plt.rc('font', **font)
mk1 = 16
mk2 = 15
mk3 = 3
mk10 = 35
lw1= 6
alpha0 = 1

markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = [' ', 'g', 'r', 'c','m','k','y','b', 'g', 'b','b', 'g', 'r', 'c','m','k','y','b', 'g', 'b']

# coefficient of variance standard: reduce alpha
coef_std1 = 0.25
# coefficient of variane absolute standard: do not plot
covar = 0.25
length_standard = 2
# for i in range(1,13):
# 	plt.figure(i,figsize=(45,30))

##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
experiment_list = np.array([
	18831003,
	18902004,
	18914000,
	18929002,
	19522001,
	19529000,
	# # 19603000,
	# # 19607000,
	# # 19607001,
	# # 19612001,
	# 19612002,
	19614000,
	19625000,
	# 19626000,
	19626002,
	# 19731000,
	# 19731001,
	# 19805000,
	19809001,
	])
np.save('experiment_list',experiment_list)
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

	I_pump_mean = [] 
	I_pump_plus = []
	I_pump_minus = []
	cytNa_mean = []
	cytNa_plus = []
	cytNa_minus = []

	counter1 = 0
	counter2=0
	counter3 = 0
	output1 = []
	if constant1=='I':
		IpumpMax = list10[-2]
		IpumpMax = float(IpumpMax)

		list1 = list0[::-1]

		control_param = IpumpMax
		for j in range(0,len(list1)):
			gp = float(list1[j])
			x_all = int(control_param*10)
			# print(x_all)

			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# burst_duration = burst_duration/np.min(burst_duration)
			# burst_duration = burst_duration/np.max(burst_duration)
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			
			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				# plt.figure(1)
				y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
				# y1 = y1/np.min(burst_duration)
				# y1 = y1/np.max(burst_duration)
				x = np.zeros((1,len(y1)))+gp
				x=x[0]
				gp_list1.append(gp)
				# BD_mean_list.append(np.mean(burst_duration)/np.min(burst_duration))
				BD_mean_list.append(np.mean(burst_duration))#/np.max(burst_duration))

				counter1 = counter1+1


			#IBI
			interbust_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# interbust_interval = interbust_interval/np.min(interbust_interval)
			# interbust_interval = interbust_interval/np.max(interbust_interval)
			coefficient_of_variation = np.std(interbust_interval)/np.mean(interbust_interval)

			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				y2 = np.array((np.mean(interbust_interval)-np.std(interbust_interval),np.mean(interbust_interval),np.mean(interbust_interval)+np.std(interbust_interval)))
				# y2 = y2/np.min(interbust_interval)
				# y2 = y2/np.max(interbust_interval)
				x = np.zeros((1,len(y2)))+gp
				x=x[0]
				gp_list2.append(gp)

				# IBI_mean_list.append(np.mean(interbust_interval)/np.min(interbust_interval))
				IBI_mean_list.append(np.mean(interbust_interval))#/np.max(interbust_interval))

				counter2=counter2+1

			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# cycle_period = cycle_period/np.min(cycle_period)
			# cycle_period = cycle_period/np.max(cycle_period)
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)

			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
				# y3 = y3/np.min(cycle_period)
				# y3 = y3/np.max(cycle_period)
				x = np.zeros((1,len(y3)))+gp
				x=x[0]

				gp_list3.append(gp)

				T_mean_list.append(np.mean(cycle_period))#/np.max(cycle_period))

				counter3 = counter3+1

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			y4 = np.array((np.mean(Hz)-np.std(Hz),np.mean(Hz),np.mean(Hz)+np.std(Hz)))
			# y4 = y4/np.min(Hz)
			# y4 = y4/np.max(Hz)
			x = np.zeros((1,len(y4)))+gp
			x=x[0]

			gp_list4.append(gp)
			# Hz_mean_list.append(np.mean(Hz)/np.min(Hz))
			Hz_mean_list.append(np.mean(Hz))#/np.max(Hz))

			

			if len(gp_list1)>=length_standard:
				param1 = float(IpumpMax)
				delta_gp = float(gp_list1[0]-gp_list1[-1])
				delta_bd = float((BD_mean_list[0]-BD_mean_list[-1])*100/BD_mean_list[0])
				temp1 = [constant1,param1,delta_gp,delta_bd]
				output1.append([temp1])

			# if len(gp_list2)>=length_standard:
			# 	# plt.figure(2)
			# 	# plt.plot(gp_list2,IBI_mean_list,color = colores[x_all],linewidth=lw1)

			# if len(gp_list3)>=length_standard:
			# 	# plt.figure(3)
			# 	# plt.plot(gp_list3,T_mean_list,color=colores[x_all],linewidth=lw1)

			# if len(gp_list4)>=length_standard:
			# 	# plt.figure(4)
				# plt.plot(gp_list4,Hz_mean_list,color=colores[x_all],linewidth=lw1)


	if constant1=='G':
		gp = list10[-2]
		gp = float(gp)
		control_param=gp

		list1 = list0[::-1]

		for j in range(0,len(list1)):
			IpumpMax = list1[j]
			IpumpMax = float(IpumpMax)
			x_all = int(control_param)

			# BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# burst_duration = burst_duration/np.min(burst_duration)
			# burst_duration = burst_duration/np.max(burst_duration)
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)

			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
				# y1 = y1/np.min(burst_duration)
				# y1 = y1/np.max(burst_duration)
				x = np.zeros((1,len(y1)))+IpumpMax
				x=x[0]

				pump_list1.append(IpumpMax)
				# BD_mean_list.append(np.mean(burst_duration)/np,min(burst_duration))
				BD_mean_list.append(np.mean(burst_duration))#/np.max(burst_duration))

				counter1 = counter1+1

			#IBI
			interbust_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# interbust_interval = interbust_interval/np.min(interbust_interval)
			# interbust_interval = interbust_interval/np.max(interbust_interval)
			coefficient_of_variation = np.std(interbust_interval)/np.mean(interbust_interval)

			if coefficient_of_variation < covar:			
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				y2 = np.array((np.mean(interbust_interval)-np.std(interbust_interval),np.mean(interbust_interval),np.mean(interbust_interval)+np.std(interbust_interval)))
				# y2 = y2/np.min(interbust_interval)
				# y2 = y2/np.max(interbust_interval)
				x = np.zeros((1,len(y2)))+IpumpMax
				x=x[0]

				pump_list2.append(IpumpMax)

				# IBI_mean_list.append(np.mean(interbust_interval)/np.min(interbust_interval))
				IBI_mean_list.append(np.mean(interbust_interval))#/np.max(interbust_interval))

				counter2=counter2+1
				

			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# cycle_period = cycle_period/np.min(cycle_period)
			# cycle_period = cycle_period/np.max(cycle_period)
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)


			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
				# y3 = y3/np.min(y3)
				# y3 = y3/np.max(y3)
				x = np.zeros((1,len(y3)))+IpumpMax
				x=x[0]

				pump_list3.append(IpumpMax)

				# T_mean_list.append(np.mean(cycle_period)/np.min(cycle_period))
				T_mean_list.append(np.mean(cycle_period))#/np.max(cycle_period))

				counter3 = counter3+1
				

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			y4 = np.array((np.mean(Hz)-np.std(Hz),np.mean(Hz),np.mean(Hz)+np.std(Hz)))
			# y4 = y4/np.min(Hz)
			# y4 = y4/np.max(Hz)

			x = np.zeros((1,len(y4)))+IpumpMax
			x=x[0]

			pump_list4.append(IpumpMax)

			# Hz_mean_list.append(np.mean(Hz)/np.min(Hz))
			Hz_mean_list.append(np.mean(Hz))#/np.max(Hz))

			if len(pump_list1)>=length_standard:
				param1 = float(gp)
				delta_IpumpMax = float(pump_list1[0]-pump_list1[-1])
				delta_bd = float((BD_mean_list[0]-BD_mean_list[-1])*100/BD_mean_list[0])
				temp1 = [constant1,param1,delta_IpumpMax,delta_bd]
				# output1 = np.append(output1,temp1,axis=0)
				output1.append([temp1])

				# param1 = float(IpumpMax)
				# delta_gp = float(gp_list1[0]-gp_list1[-1])
				# delta_bd = 
				# temp1 = [constant1,param1,delta_gp,delta_bd]
				# output1.append([temp1])

			# if len(pump_list2)>=length_standard:
			# 	# plt.figure(7)
			# 	# plt.plot(pump_list2,IBI_mean_list,color = colores[x_all],linewidth=lw1)
			# if len(pump_list3)>=length_standard:
			# 	# plt.figure(8)
			# 	# plt.plot(pump_list3,T_mean_list,color=colores[x_all],linewidth=lw1)
			# if len(pump_list4)>=length_standard:
			# 	# plt.figure(9)
				# plt.plot(pump_list4,Hz_mean_list,color=colores[x_all],linewidth=lw1)
	# print(output1)
	return output1
# color_list=0
#loop on files inside experiment list
out2 =[]
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
		print(fileID,list10)
		out1 = compute(list10)
		if len(out1)>1:
			# print(out1)
			out2.append(out1)
			# out2 = np.append([out2],[out1],axis=0)

for j in range(0,len(out2)):
	temp1 = out2[j]
	for k in range(0,len(temp1)):
		temp2 = temp1[k]
		temp2 = temp2[0]
		constant1 = temp2[0]
		control = temp2[1]

		delta_bd = temp2[3]
		if constant1=='I':
			col0 = int(control*10)
			delta_gp = temp2[2]
			plt.figure(1,figsize=(30,30))
			plt.plot(delta_gp,delta_bd,'*',markersize=15,color = colores[col0])
			plt.ylabel('% change in BD')
			plt.xlabel(r'$\Delta$ G$_P$')
		if constant1=='G':
			col0 = int(control)
			delta_IpumpMax = temp2[2]
			plt.figure(2,figsize=(30,30))
			plt.plot(delta_IpumpMax,delta_bd,'*',markersize=15,color=colores[col0])
			plt.ylabel('% change in BD')
			plt.xlabel(r'$\Delta$ I$_{pump}^{max}$')

###############################################################################################################
###############################################################################################################
riff1 = [2,3,4,5,7,8,9]
for i in range(0,len(riff1)):
	str1 = float(riff1[i])*0.1
	str1 = format(str1,'.1f')
	plt.figure(1)
	plt.plot(0,0,'*',markersize=0.01,color = colores[riff1[i]],label='IpumpMax='+str(str1)+' nA')
plt.legend(markerscale=3000,bbox_to_anchor=[1,1.12],ncol=4,fontsize=30)
plt.savefig('deltaBD_pump')

riff2= [2,5,6]
for i in range(0,len(riff2)):
	str2 = float(riff2[i])
	col0 = int(control)
	plt.figure(2)
	plt.plot(0,0,'*',markersize=0.01,color = colores[riff2[i]],label='G$_P$='+str(str2)+' nS')
plt.legend(markerscale=3000,bbox_to_anchor=[1,1.12],ncol=4,fontsize=30)
plt.savefig('deltaBD_gp')

# plt.ion()
# plt.show()