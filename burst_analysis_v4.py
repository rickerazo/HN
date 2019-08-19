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
length_standard = 3
for i in range(1,13):
	plt.figure(i,figsize=(45,30))

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
	output1 = np.array([])
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

				plt.figure(1)
				y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
				# y1 = y1/np.min(burst_duration)
				# y1 = y1/np.max(burst_duration)
				x = np.zeros((1,len(y1)))+gp
				x=x[0]
				if counter1==0:
					plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
					# plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
					plt.plot(gp,np.mean(burst_duration),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
				else:	
					plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
					plt.plot(gp,np.mean(burst_duration),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])

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

				plt.figure(2)
				y2 = np.array((np.mean(interbust_interval)-np.std(interbust_interval),np.mean(interbust_interval),np.mean(interbust_interval)+np.std(interbust_interval)))
				# y2 = y2/np.min(interbust_interval)
				# y2 = y2/np.max(interbust_interval)
				x = np.zeros((1,len(y2)))+gp
				x=x[0]

				if counter2==0:
					plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
					# plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
					plt.plot(gp,np.mean(interbust_interval),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
					# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				else:
					plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
					plt.plot(gp,np.mean(interbust_interval),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
					# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1)

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

				plt.figure(3)
				y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
				# y3 = y3/np.min(cycle_period)
				# y3 = y3/np.max(cycle_period)
				x = np.zeros((1,len(y3)))+gp
				x=x[0]

				if counter3==0:
					# plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
					plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
					plt.plot(gp,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
					# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				else:
					plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
					plt.plot(gp,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
					# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				plt.plot(x[0],np.mean(cycle_period),'*',markersize=mk1,color=colores[x_all])
				plt.plot(x[0],np.min(cycle_period),'-',markersize=mk1,color=colores[x_all])
				plt.plot(x[0],np.max(cycle_period),'-',markersize=mk1,color=colores[x_all])

				gp_list3.append(gp)

				# T_mean_list.append(np.mean(cycle_period)/np.min(cycle_period))
				T_mean_list.append(np.mean(cycle_period))#/np.max(cycle_period))

				counter3 = counter3+1

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# Hz = Hz/np.min(Hz)
			# Hz = Hz/np.max(Hz)
			y4 = np.array((np.mean(Hz)-np.std(Hz),np.mean(Hz),np.mean(Hz)+np.std(Hz)))
			# y4 = y4/np.min(Hz)
			# y4 = y4/np.max(Hz)
			x = np.zeros((1,len(y4)))+gp
			x=x[0]
			plt.figure(4)
			if j ==0:
				plt.plot(x,y4,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])#,label=r'$I_{pump}^{max}$='+str(control_param))
				# plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],label=r'$I_{pump}^{max}$='+str(control_param))
				# plt.plot(x,Hz,markersize=mk2,marker=markers[i],color=colores[x_all],label=str(fileID))
			else:
				plt.plot(x,y4,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])
			# plt.plot(x[0],np.mean(Hz),'*',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.min(Hz),'-',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.max(Hz),'-',markersize=mk1,color=colores[x_all])

			gp_list4.append(gp)
			# Hz_mean_list.append(np.mean(Hz)/np.min(Hz))
			Hz_mean_list.append(np.mean(Hz))#/np.max(Hz))

			
			#Ipump
			Ipump = np.load(fh+'pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Ipump = Ipump/IpumpMax
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			plt.figure(5)
			# format(i,'.1f')
			# alpha1 = (i+1)*0.1
			alpha1 = 0.5
			if np.mean(Ipump)>0:
				if j ==0:
					# plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],label=r'$I_{pump}^{max}$='+str(control_param))
					plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all])#,label=r'$I_{pump}^{max}$='+str(control_param))
					# plt.plot(cytNa, Ipump,markersize=mk1,marker=markers[i],color=colores[x_all],label=fileID)
				else:
					plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],alpha=alpha1)

					I_pump_minus.append(np.mean(Ipump)-np.std(Ipump))
					I_pump_plus.append(np.mean(Ipump)+np.std(Ipump))
					I_pump_mean.append(np.mean(Ipump))

					cytNa_minus.append(np.mean(cytNa)- np.std(cytNa))
					cytNa_plus.append(np.mean(cytNa)+ np.std(cytNa))
					cytNa_mean.append(np.mean(cytNa))

					# plt.plot([19.9,19.9],[np.min(Ipump),np.max(Ipump)],marker=markers[i],markersize=mk3,color=colores[x_all])
					# plt.plot(19.9,np.mean(Ipump),marker=markers[i],markersize=20,color=colores[x_all])	

					# plt.plot([np.min(cytNa),np.max(cytNa)],[0.05,0.05],marker=markers[i],markersize=mk3,color=colores[x_all])
					# plt.plot(np.mean(cytNa),0.05,marker=markers[i],markersize=20,color=colores[x_all])	
			if len(gp_list1)>=length_standard:
				plt.figure(1)
				plt.plot(gp_list1,BD_mean_list,color=colores[x_all],linewidth=lw1)
				# temp2 = BD_mean_list[0]-BD_mean_list[-1]*100/BD_mean_list[0]
				temp1 = np.array([IpumpMax,gp_list1[0]-gp_list1[-1],BD_mean_list[0]-BD_mean_list[-1]*100/BD_mean_list[0]])
				output1 = np.append(output1,temp1)
				# print(temp1)
			if len(gp_list2)>=length_standard:
				plt.figure(2)
				plt.plot(gp_list2,IBI_mean_list,color = colores[x_all],linewidth=lw1)
			if len(gp_list3)>=length_standard:
				plt.figure(3)
				plt.plot(gp_list3,T_mean_list,color=colores[x_all],linewidth=lw1)
			if len(gp_list4)>=length_standard:
				plt.figure(4)
				plt.plot(gp_list4,Hz_mean_list,color=colores[x_all],linewidth=lw1)

			plt.figure(11)
			I_mem = np.load(fh+'Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Na_in = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')			
			Ip = I_mem +Ipump+0.1
			Ip= -Ip
			V1 = V/1000 # V is in miliVolts,V1 is in Volts
			Na_o = 0.115*1000 #M is a constant
			ENa = 0.02526*np.log(Na_o/(Na_in))
			mp = Ip /(gp*(V1-ENa))
			plt.plot(Ip,mp, markersize=mk1,color=colores[x_all])


	if constant1=='G':
		gp = list10[-2]
		gp = float(gp)
		control_param=gp

		list1 = list0[::-1]

		for j in range(0,len(list1)):
			IpumpMax = list1[j]
			IpumpMax = float(IpumpMax)
			x_all = int(control_param)
			# print(x_all)
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

				plt.figure(6)
				y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
				# y1 = y1/np.min(burst_duration)
				# y1 = y1/np.max(burst_duration)
				x = np.zeros((1,len(y1)))+IpumpMax
				x=x[0]

				if counter1==0:
					plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$g_p$='+str(control_param))
					# plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$g_p$='+str(control_param))
					# plt.plot(x,burst_duration,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,label=r'$g_p$='+str(control_param))
					# plt.plot(x,burst_duration,markersize=mk2,marker=markers[i],color=colores[x_all],label=str(fileID))
				else:
					plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
					# plt.plot(x,burst_duration,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2)
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


				plt.figure(7)
				y2 = np.array((np.mean(interbust_interval)-np.std(interbust_interval),np.mean(interbust_interval),np.mean(interbust_interval)+np.std(interbust_interval)))
				# y2 = y2/np.min(interbust_interval)
				# y2 = y2/np.max(interbust_interval)
				x = np.zeros((1,len(y2)))+IpumpMax
				x=x[0]

				if counter2==0:
					plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$g_p$='+str(control_param))
					# plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$g_p$='+str(control_param))
					# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,label=r'$g_p$='+str(control_param))
					# plt.plot(x,interbust_interval,markersize=mk2,marker=markers[i],color=colores[x_all],label=fileID)
				else:
					plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
					# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2)

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

				plt.figure(8)
				y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
				# y3 = y3/np.min(y3)
				# y3 = y3/np.max(y3)
				x = np.zeros((1,len(y3)))+IpumpMax
				x=x[0]

				if counter3==0:
					plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$g_p$='+str(control_param))
					# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2,label=r'$g_p$='+str(control_param))
				else:
					plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
					# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2)

				pump_list3.append(IpumpMax)

				# T_mean_list.append(np.mean(cycle_period)/np.min(cycle_period))
				T_mean_list.append(np.mean(cycle_period))#/np.max(cycle_period))

				counter3 = counter3+1
				

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# Hz = Hz/np.min(Hz)
			# Hz = Hz/np.max(Hz)
			y4 = np.array((np.mean(Hz)-np.std(Hz),np.mean(Hz),np.mean(Hz)+np.std(Hz)))
			# y4 = y4/np.min(Hz)
			# y4 = y4/np.max(Hz)

			x = np.zeros((1,len(y4)))+IpumpMax
			x=x[0]

			plt.figure(9)
			if j ==0:
				# plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],label=r'$g_p$='+str(control_param))
				plt.plot(x,y4,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])#,label=r'$g_p$='+str(control_param))
				# plt.plot(x,Hz,markersize=mk2,marker=markers[i],color=colores[x_all],label=str(fileID))
			else:
				plt.plot(x,y4,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])
			# plt.plot(x[0],np.mean(Hz),'*',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.min(Hz),'-',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.max(Hz),'-',markersize=mk1,color=colores[x_all])

			pump_list4.append(IpumpMax)

			# Hz_mean_list.append(np.mean(Hz)/np.min(Hz))
			Hz_mean_list.append(np.mean(Hz))#/np.max(Hz))
			
			#Ipump
			Ipump = np.load(fh+'pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Ipump = Ipump/IpumpMax
			cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			plt.figure(10)
			# alpha1 = (i+1)*0.1
			alpha1 = 0.5
			if np.mean(Ipump)>0:
				if j ==0:
					plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all])#,label=r'$I_{pump}^{max}$='+str(control_param))
					# plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],label=r'$I_{pump}^{max}$='+str(control_param))
					# plt.plot(cytNa, Ipump,markersize=mk1,marker=markers[i],color=colores[x_all],label=fileID)
				else:
					plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],alpha=alpha1)

					I_pump_minus.append(np.mean(Ipump)-np.std(Ipump))
					I_pump_plus.append(np.mean(Ipump)+np.std(Ipump))
					I_pump_mean.append(np.mean(Ipump))

					cytNa_minus.append(np.mean(cytNa)- np.std(cytNa))
					cytNa_plus.append(np.mean(cytNa)+ np.std(cytNa))
					cytNa_mean.append(np.mean(cytNa))
			if len(pump_list1)>=length_standard:
				plt.figure(6)
				plt.plot(pump_list1,BD_mean_list,color=colores[x_all],linewidth=lw1)
				temp1 = np.array([IpumpMax,pump_list1[0]-pump_list1[-1],BD_mean_list[0]-BD_mean_list[-1]*100/BD_mean_list[0]])
				output1 = np.append(output1,temp1)

			if len(pump_list2)>=length_standard:
				plt.figure(7)
				plt.plot(pump_list2,IBI_mean_list,color = colores[x_all],linewidth=lw1)
			if len(pump_list3)>=length_standard:
				plt.figure(8)
				plt.plot(pump_list3,T_mean_list,color=colores[x_all],linewidth=lw1)
			if len(pump_list4)>=length_standard:
				plt.figure(9)
				plt.plot(pump_list4,Hz_mean_list,color=colores[x_all],linewidth=lw1)

			plt.figure(12)
			I_mem = np.load(fh+'Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')			
			Na_in = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')			
			Ip = I_mem +Ipump+0.1
			Ip= -Ip
			V1 = V/1000 # V is in miliVolts,V1 is in Volts
			Na_o = 0.115*1000 #M is a constant
			ENa = 0.02526*np.log(Na_o/(Na_in))
			mp = Ip /(gp*(V1-ENa))
			plt.plot(Ip,mp)

	return output1
# color_list=0
#loop on files inside experiment list

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

###############################################################################################################
###############################################################################################################

plt.figure(1)
plt.ylabel('Burst duration (s)',fontsize=65,weight='bold')
plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,10,0,9])
plt.figure(2)
plt.ylabel('Interbust interval (s)',fontsize=65,weight='bold')
plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,10,0,9])
plt.figure(3)
plt.ylabel('Period',fontsize=65,weight='bold')
plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(1,12,3))
# plt.axis([0,10,0,12])
plt.figure(4)
plt.ylabel('Spiking frequency',fontsize=65,weight='bold')
plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(0,100,15))
# plt.axis([0,10,0,100])
plt.figure(5)
plt.ylabel(r'Activation of I$_{pump}$',fontsize=65,weight='bold')
plt.xlabel('Intracellular Na concentration (mM)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(0.15,1,0.15))
# plt.axis([0,20,0,1])


plt.figure(6)
plt.ylabel('Burst duration (s)',fontsize=65,weight='bold')
plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,1,0,9])
plt.figure(7)
plt.ylabel('Interbust interval (s)',fontsize=65,weight='bold')
plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,1,0,9])
plt.figure(8)
plt.ylabel('Period',fontsize=65,weight='bold')
plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(1,12,3))
# plt.axis([0,1,0,12])
plt.figure(9)
plt.ylabel('Spiking frequency',fontsize=65,weight='bold')
plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(0,100,15))
# plt.axis([0, 1, 0, 50])
plt.figure(10)
plt.ylabel(r'Activation of I$_{pump}$',fontsize=65,weight='bold')
plt.xlabel('Intracellular Na concentration (mM)',fontsize=65,weight='bold')
plt.yticks(ticks=np.arange(0.15,1,0.15))
# plt.axis([0,20,0,1])

#############################################################################################################

markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = [' ', 'g', 'r', 'c','m','k','y','b', 'g', 'b','b', 'g', 'r', 'c','m','k','y','b', 'g', 'b']

riff2 =[2,4,10]
color_riff2= [colores[2],colores[4],colores[10]]#,colores[6],colores[8],colores[9]]
titles = [' ','BD','IBI','T','Hz','Ipump_cytNa','BD','IBI','T','Hz','Ipump_cytNa']
#gp 
for i in range(1,6):
	plt.figure(i)
	for j in range(0,len(riff2)):
		p1 = float(riff2[j])*0.1
		p1 = format(format(p1,'.1f'))
		str1 = str(p1)
		plt.plot(0,0,'o',markersize=1e-3,color=color_riff2[j],label=r'I$_{pump}^{max}$= '+str1+' nA')
		plt.legend(loc=2,bbox_to_anchor=(-0.1,1.1),markerscale=2e+4,ncol=10,fontsize=40)
		plt.savefig('gp_'+titles[i])


riff1 = [6,7,8]
color_riff1 = [colores[6],colores[7], colores[8]]#,colores[7],colores[9]]

#pump
for i in range(6,11):
	plt.figure(i)
	for j in range(0,len(riff1)):
		p1 = float(riff1[j])
		str1 = str(p1)
		plt.plot(0,0,'o',markersize=1e-3,color=color_riff1[j],label=r'G$_P$= '+str1+' nS')
		plt.legend(loc=2,bbox_to_anchor=(-0.1,1.1),markerscale=2e+4,ncol=10,fontsize=50)
		plt.savefig('pump_'+titles[i])
# plt.ion()
# plt.show()

plt.figure(11)
plt.ylabel('mp')
plt.xlabel('Ip')
plt.savefig('gp_Ip_mp')

plt.figure(12)
plt.ylabel('mp')
plt.xlabel('Ip')
plt.savefig('p_Ip_mp')