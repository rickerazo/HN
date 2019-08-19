# 3D PLOTS.PY
####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps


#### libraries necessary for code
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
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
coef_std1 = 0.2
# coefficient of variane absolute standard: do not plot
covar = 0.35

#13,14, 15 -> BD, IBI, T
# for i in range(1,16):
# 	plt.figure(i,figsize=(45,30))

BD3 = plt.figure(13,figsize=(75,75))
BDax = plt.axes(projection='3d')

IBI3 = plt.figure(14,figsize=(75,75))
IBIax = plt.axes(projection='3d')

T3 = plt.figure(15,figsize=(75,75))
Tax = plt.axes(projection='3d')

VNa = plt.figure(16, figsize=(75,75))
VNax = plt.axes()

##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
experiment_list = np.load('experiment_list.npy')


def compute(list10):
	alpha2 = 0.1		
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
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			
			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				# plt.figure(1)
				# y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
				# x = np.zeros((1,len(y1)))+gp
				# x=x[0]
				# # if counter1==0:
				# 	# plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
				# 	# plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				# 	# plt.plot(gp,np.mean(burst_duration),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
				# # else:	
				# 	# plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				# 	# plt.plot(gp,np.mean(burst_duration),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])

				# gp_list1.append(gp)
				# BD_mean_list.append(np.mean(burst_duration))

				# counter1 = counter1+1

			# # trend
			# liga10 = np.zeros((1,len(gp)))
			# liga0=list0[0]+IpumpMax
			# print(liga10)
			# BDax.plot3D(liga10,gp_list1,BD_mean_list1,linewidth=lw1,color=colores[x_all],alpha=alpha2)

			# mean and standard dev
			BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
			BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration)+np.std(burst_duration),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
			BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration)-np.std(burst_duration),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)
			# BD_trend = 
			#IBI
			interbust_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			
			coefficient_of_variation = np.std(interbust_interval)/np.mean(interbust_interval)

			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				# plt.figure(2)
				# y2 = np.array((np.mean(interbust_interval)-np.std(interbust_interval),np.mean(interbust_interval),np.mean(interbust_interval)+np.std(interbust_interval)))
				# x = np.zeros((1,len(y2)))+gp
				# x=x[0]

				# if counter2==0:
				# 	plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
				# 	# plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				# 	plt.plot(gp,np.mean(interbust_interval),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
				# 	# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				# else:
				# 	plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				# 	plt.plot(gp,np.mean(interbust_interval),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
				# 	# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1)

			gp_list2.append(gp)

			IBI_mean_list.append(np.mean(interbust_interval))

				# counter2=counter2+1

			IBIax.scatter3D(IpumpMax,gp,np.mean(interbust_interval),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
			IBIax.scatter3D(IpumpMax,gp,np.mean(interbust_interval)+np.std(interbust_interval),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
			IBIax.scatter3D(IpumpMax,gp,np.mean(interbust_interval)-np.std(interbust_interval),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)


			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)

			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				# plt.figure(3)
				# y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
				# x = np.zeros((1,len(y3)))+gp
				# x=x[0]

				# if counter3==0:
				# 	# plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				# 	plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
				# 	plt.plot(gp,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
				# 	# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$I_{pump}^{max}$='+str(control_param))
				# else:
				# 	plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				# 	plt.plot(gp,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])
				# 	# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				# plt.plot(x[0],np.mean(cycle_period),'*',markersize=mk1,color=colores[x_all])
				# plt.plot(x[0],np.min(cycle_period),'-',markersize=mk1,color=colores[x_all])
				# plt.plot(x[0],np.max(cycle_period),'-',markersize=mk1,color=colores[x_all])

			gp_list3.append(gp)

			T_mean_list.append(np.mean(cycle_period))

				# counter3 = counter3+1

			Tax.scatter3D(IpumpMax,gp,np.mean(cycle_period),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
			Tax.scatter3D(IpumpMax,gp,np.mean(cycle_period)+np.std(cycle_period),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
			Tax.scatter3D(IpumpMax,gp,np.mean(cycle_period)-np.std(cycle_period),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)


			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'
				# plt.figure(4)
				# y4 = np.array((np.mean(Hz)-np.std(Hz), np.mean(Hz), np.mean(Hz)+np.std(Hz)))
				# x = np.zeros((1,len(y3)))+gp
				# x=x[0]
				# # plt.plot()
				# plt.plot(x,y4,'o-',markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))
				# plt.plot(gp,np.mean(cycle_period),marker='o',markersize=mk10,fillstyle=fst1,color = colores[x_all])

			# x = np.zeros((1,len(Hz)))+gp
			# x=x[0]
			# plt.figure(4)
			# if j ==0:
			# 	plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])#,label=r'$I_{pump}^{max}$='+str(control_param))
			# 	# plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],label=r'$I_{pump}^{max}$='+str(control_param))
			# 	# plt.plot(x,Hz,markersize=mk2,marker=markers[i],color=colores[x_all],label=str(fileID))
			# else:
			# 	plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])
			# plt.plot(x[0],np.mean(Hz),'*',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.min(Hz),'-',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.max(Hz),'-',markersize=mk1,color=colores[x_all])

			gp_list4.append(gp)
			Hz_mean_list.append(np.mean(Hz))

			
			#Ipump
			# Ipump = np.load(fh+'pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# Ipump = Ipump/IpumpMax
			# cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# # plt.figure(5)
			# # format(i,'.1f')
			# alpha1 = (i+1)*0.1
			# if np.mean(Ipump)>0:
			# 	if j ==0:
			# 		# plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],label=r'$I_{pump}^{max}$='+str(control_param))
			# 		plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all])#,label=r'$I_{pump}^{max}$='+str(control_param))
			# 		# plt.plot(cytNa, Ipump,markersize=mk1,marker=markers[i],color=colores[x_all],label=fileID)
			# 	else:
			# 		plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],alpha=alpha1)

			# 		I_pump_minus.append(np.mean(Ipump)-np.std(Ipump))
			# 		I_pump_plus.append(np.mean(Ipump)+np.std(Ipump))
			# 		I_pump_mean.append(np.mean(Ipump))

			# 		cytNa_minus.append(np.mean(cytNa)- np.std(cytNa))
			# 		cytNa_plus.append(np.mean(cytNa)+ np.std(cytNa))
			# 		cytNa_mean.append(np.mean(cytNa))

					# plt.plot([19.9,19.9],[np.min(Ipump),np.max(Ipump)],marker=markers[i],markersize=mk3,color=colores[x_all])
					# plt.plot(19.9,np.mean(Ipump),marker=markers[i],markersize=20,color=colores[x_all])	

					# plt.plot([np.min(cytNa),np.max(cytNa)],[0.05,0.05],marker=markers[i],markersize=mk3,color=colores[x_all])
					# plt.plot(np.mean(cytNa),0.05,marker=markers[i],markersize=20,color=colores[x_all])	

			# plt.figure(1)
			# plt.plot(gp_list1,BD_mean_list,color=colores[x_all],linewidth=lw1)
			# plt.figure(2)
			# plt.plot(gp_list2,IBI_mean_list,color = colores[x_all],linewidth=lw1)
			# plt.figure(3)
			# plt.plot(gp_list3,T_mean_list,color=colores[x_all],linewidth=lw1)
			# plt.figure(4)
			# plt.plot(gp_list4,Hz_mean_list,color=colores[x_all],linewidth=lw1)

			# plt.figure(11)
			# I_mem = np.load(fh+'Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# Na_in = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')			
			# Ip = I_mem +Ipump+0.1
			# Ip= -Ip
			# V1 = V/1000 # V is in miliVolts,V1 is in Volts
			# Na_o = 0.115*1000 #M is a constant
			# ENa = 0.02526*np.log(Na_o/(Na_in))
			# mp = Ip /(gp*(V1-ENa))
			# plt.plot(Ip,mp, markersize=mk1,color=colores[x_all])


			# V Na
			# if np.mean(Na_in):
			# 	VNax.plot(Na_in,V,color=colores[x_all])
		

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

			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)

			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				# plt.figure(6)
				# y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
				# x = np.zeros((1,len(y1)))+IpumpMax
				# x=x[0]

				# if counter1==0:
				# 	plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$g_p$='+str(control_param))
				# 	# plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$g_p$='+str(control_param))
				# 	# plt.plot(x,burst_duration,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,label=r'$g_p$='+str(control_param))
				# 	# plt.plot(x,burst_duration,markersize=mk2,marker=markers[i],color=colores[x_all],label=str(fileID))
				# else:
				# 	plt.plot(x,y1,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				# 	# plt.plot(x,burst_duration,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2)
			pump_list1.append(IpumpMax)
			BD_mean_list.append(np.mean(burst_duration))

				# counter1 = counter1+1
			# trend
			# liga10 = np.zeros((1,len(gp)))
			# liga0=list0[0]+IpumpMax
			# print(liga10)
			# BDax.plot3D(liga10,gp_list1,BD_mean_list1,linewidth=lw1,color=colores[x_all],alpha=alpha2)

			# mean and standard dev
			BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
			BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration)+np.std(burst_duration),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
			BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration)-np.std(burst_duration),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)



			#IBI
			interbust_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')

			coefficient_of_variation = np.std(interbust_interval)/np.mean(interbust_interval)

			# if coefficient_of_variation < covar:			
			# 	if coefficient_of_variation>coef_std1:
			# 		alpha2 = 0.3
			# 		fst1 = 'none'
			# 	else:
			# 		alpha2 = 1
			# 		fst1 = 'full'


			# 	plt.figure(7)
			# 	y2 = np.array((np.mean(interbust_interval)-np.std(interbust_interval),np.mean(interbust_interval),np.mean(interbust_interval)+np.std(interbust_interval)))
			# 	x = np.zeros((1,len(y2)))+IpumpMax
			# 	x=x[0]

			# 	if counter2==0:
			# 		plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$g_p$='+str(control_param))
			# 		# plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1,label=r'$g_p$='+str(control_param))
			# 		# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2,label=r'$g_p$='+str(control_param))
			# 		# plt.plot(x,interbust_interval,markersize=mk2,marker=markers[i],color=colores[x_all],label=fileID)
			# 	else:
			# 		plt.plot(x,y2,markersize=mk2,linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			# 		# plt.plot(x,interbust_interval,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],alpha=alpha2)

			pump_list2.append(IpumpMax)

			IBI_mean_list.append(np.mean(interbust_interval))

				# counter2=counter2+1

			IBIax.scatter3D(IpumpMax,gp,np.mean(interbust_interval),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
			IBIax.scatter3D(IpumpMax,gp,np.mean(interbust_interval)+np.std(interbust_interval),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
			IBIax.scatter3D(IpumpMax,gp,np.mean(interbust_interval)-np.std(interbust_interval),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)

				

			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')

			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)


			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'

				# plt.figure(8)
				# y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
				# x = np.zeros((1,len(y3)))+IpumpMax
				# x=x[0]

				# if counter3==0:
				# 	plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$g_p$='+str(control_param))
				# 	# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2,label=r'$g_p$='+str(control_param))
				# else:
				# 	plt.plot(x,y3,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)
				# 	# plt.plot(x,cycle_period,markersize=mk2,linewidth=lw1,marker = markers[i],color=colores[x_all],alpha=alpha2)

			pump_list3.append(IpumpMax)

			T_mean_list.append(np.mean(cycle_period))

				# counter3 = counter3+1
				
			Tax.scatter3D(IpumpMax,gp,np.mean(cycle_period),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
			Tax.scatter3D(IpumpMax,gp,np.mean(cycle_period)+np.std(cycle_period),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
			Tax.scatter3D(IpumpMax,gp,np.mean(cycle_period)-np.std(cycle_period),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)			

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			alpha2 = 0.1
			fst1 = 'none'
			
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			if coefficient_of_variation < covar:
				if coefficient_of_variation>coef_std1:
					alpha2 = 0.3
					fst1 = 'none'
				else:
					alpha2 = 1
					fst1 = 'full'
				
				# plt.figure(9)
				# y4 = np.array((np.mean(Hz)-np.std(Hz), np.mean(Hz), np.mean(Hz)+np.std(Hz)))
				# x = np.zeros((1,len(y4)))+gp
				# x=x[0]
				# # plt.plot()
				# # print(x,y4)
				# plt.plot(x,y4,markersize=mk2,linewidth=lw1,marker = 'o',color=colores[x_all],alpha=alpha2,fillstyle=fst1)#,label=r'$I_{pump}^{max}$='+str(control_param))

			# Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# x = np.zeros((1,len(Hz)))+IpumpMax
			# x=x[0]
			# plt.figure(9)
			# if j ==0:
			# 	# plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all],label=r'$g_p$='+str(control_param))
			# 	plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])#,label=r'$g_p$='+str(control_param))
			# 	# plt.plot(x,Hz,markersize=mk2,marker=markers[i],color=colores[x_all],label=str(fileID))
			# else:
			# 	plt.plot(x,Hz,markersize=mk2,linewidth=lw1,marker=markers[i],color=colores[x_all])
			# plt.plot(x[0],np.mean(Hz),'*',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.min(Hz),'-',markersize=mk1,color=colores[x_all])
			# plt.plot(x[0],np.max(Hz),'-',markersize=mk1,color=colores[x_all])

			pump_list4.append(IpumpMax)

			Hz_mean_list.append(np.mean(Hz))
			
			#Ipump
			# Ipump = np.load(fh+'pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# Ipump = Ipump/IpumpMax
			# cytNa = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# plt.figure(10)
			# alpha1 = (i+1)*0.1
			# if np.mean(Ipump)>0:
			# 	if j ==0:
			# 		plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all])#,label=r'$I_{pump}^{max}$='+str(control_param))
			# 		# plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],label=r'$I_{pump}^{max}$='+str(control_param))
			# 		# plt.plot(cytNa, Ipump,markersize=mk1,marker=markers[i],color=colores[x_all],label=fileID)
			# 	else:
			# 		plt.plot(cytNa, Ipump,markersize=mk3,marker=markers[i],color=colores[x_all],alpha=alpha1)

			# 		I_pump_minus.append(np.mean(Ipump)-np.std(Ipump))
			# 		I_pump_plus.append(np.mean(Ipump)+np.std(Ipump))
			# 		I_pump_mean.append(np.mean(Ipump))

			# 		cytNa_minus.append(np.mean(cytNa)- np.std(cytNa))
			# 		cytNa_plus.append(np.mean(cytNa)+ np.std(cytNa))
			# 		cytNa_mean.append(np.mean(cytNa))

			# plt.figure(6)
			# plt.plot(pump_list1,BD_mean_list,color=colores[x_all],linewidth=lw1)
			# plt.figure(7)
			# plt.plot(pump_list2,IBI_mean_list,color = colores[x_all],linewidth=lw1)
			# plt.figure(8)
			# plt.plot(pump_list3,T_mean_list,color=colores[x_all],linewidth=lw1)
			# plt.figure(9)
			# plt.plot(pump_list4,Hz_mean_list,color=colores[x_all],linewidth=lw1)

			# plt.figure(12)
			# I_mem = np.load(fh+'Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# V = np.load(fh+'V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')			
			# Na_in = np.load(fh+'cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')			
			# Ip = I_mem +Ipump+0.1
			# Ip= -Ip
			# V1 = V/1000 # V is in miliVolts,V1 is in Volts
			# Na_o = 0.115*1000 #M is a constant
			# ENa = 0.02526*np.log(Na_o/(Na_in))
			# mp = Ip /(gp*(V1-ENa))
			# plt.plot(Ip,mp)


			# V Na
			# if np.mean(Na_in):
			# 	VNax.plot(Na_in,V)
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
		compute(list10)




# plt.figure(1)
# plt.ylabel('Burst duration (s)',fontsize=65,weight='bold')
# plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,10,0,9])
# plt.figure(2)
# plt.ylabel('Interbust interval (s)',fontsize=65,weight='bold')
# plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,10,0,9])
# plt.figure(3)
# plt.ylabel('Period',fontsize=65,weight='bold')
# plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(1,12,3))
# plt.axis([0,10,0,12])
# plt.figure(4)
# plt.ylabel('Spiking frequency',fontsize=65,weight='bold')
# plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(0,100,15))
# plt.axis([0,11,0,50])
# plt.figure(5)
# plt.ylabel(r'Activation of I$_{pump}$',fontsize=65,weight='bold')
# plt.xlabel('Intracellular Na concentration (mM)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(0.15,1,0.15))
# plt.axis([0,20,0,1])


# plt.figure(6)
# plt.ylabel('Burst duration (s)',fontsize=65,weight='bold')
# plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,1,0,9])
# plt.figure(7)
# plt.ylabel('Interbust interval (s)',fontsize=65,weight='bold')
# plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(1,8,2))
# plt.axis([0,1,0,9])
# plt.figure(8)
# plt.ylabel('Period',fontsize=65,weight='bold')
# plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(1,12,3))
# plt.axis([0,1,0,12])
# plt.figure(9)
# plt.ylabel('Spiking frequency',fontsize=65,weight='bold')
# plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(0,100,15))
# plt.axis([0, 1, 0, 50])
# plt.figure(10)
# plt.ylabel(r'Activation of I$_{pump}$',fontsize=65,weight='bold')
# plt.xlabel('Intracellular Na concentration (mM)',fontsize=65,weight='bold')
# plt.yticks(ticks=np.arange(0.15,1,0.15))
# plt.axis([0,20,0,1])

riff1 = [2,5.5,6,9.0]
color_riff1 = [colores[2], colores[5],colores[6],colores[9]]

# riff2 =[2,8,3,4,2,3,3]
riff2 =[2,3,5,8,9]
color_riff2= [colores[2],colores[3],colores[5],colores[8],colores[9]]
titles = [' ','BD','IBI','T','Hz','Ipump_cytNa','BD','IBI','T','Hz','Ipump_cytNa']
#pump
# for i in range(1,6):
# 	plt.figure(i)
# 	for j in range(0,len(riff2)):
# 		p1 = float(riff2[j])*0.1
# 		p1 = format(format(p1,'.1f'))
# 		str1 = str(p1)
# 		plt.plot(0,0,'o',markersize=1e-3,color=color_riff2[j],label=r'I$_{pump}^{max}$= '+str1+' nA')
# 		plt.legend(loc=2,bbox_to_anchor=(0,1.1),markerscale=2e+4,ncol=10,fontsize=40)
# 		plt.savefig('gp_'+titles[i])



# #pump
# for i in range(6,11):
# 	plt.figure(i)
# 	for j in range(0,len(riff1)):
# 		p1 = float(riff1[j])
# 		str1 = str(p1)
# 		plt.plot(0,0,'o',markersize=1e-3,color=color_riff1[j],label=r'G$_P$= '+str1+' nS')
# 		# plt.ylabel('mp')
# 		# plt.xlabel('Ip')
# 		plt.legend(loc=2,bbox_to_anchor=(0,1.1),markerscale=2e+4,ncol=10,fontsize=50)
# 		plt.savefig('pump_'+titles[i])
# plt.ion()
# plt.show()

# for i in range(11,15):
# 	plt.figure(i)
# 	for j in range(len)

# plt.figure(11)
# plt.ylabel('mp')
# plt.xlabel('Ip')
# plt.savefig('gp_Ip_mp')

# plt.figure(12)
# plt.ylabel('mp')
# plt.xlabel('Ip')
# plt.savefig('p_Ip_mp')

# BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration),linewidth=lw1,marker='o',color=colores[x_all],alpha=alpha2)
# BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration)+np.std(burst_duration),linewidth=lw1,marker='^',color=colores[x_all],alpha=alpha2)
# BDax.scatter3D(IpumpMax,gp,np.mean(burst_duration)-np.std(burst_duration),linewidth=lw1,marker='v',color=colores[x_all],alpha=alpha2)
I_e1 = np.arange(0.2,1.0,0.2)
I_e2 = np.arange(1,10,2)
I_e3 = np.arange(2,11,2)
I_e4 = np.arange(1,12,3)
label1 = 25

plt.figure(13)
BDax.set_xlabel('IpumpMax (nA)',labelpad=label1)
BDax.set_xticks(I_e1)
BDax.set_ylabel('Gp (nS)',labelpad=label1)
BDax.set_yticks(I_e2)
BDax.set_zlabel('BD (s)',labelpad=label1)
BDax.set_zticks(I_e3)

plt.savefig('3D_BD')

plt.figure(14)
IBIax.set_xlabel('IpumpMax (nA)',labelpad=label1)
IBIax.set_xticks(I_e1)
IBIax.set_ylabel('Gp (nS)',labelpad=label1)
IBIax.set_yticks(I_e2)
IBIax.set_zlabel('IBI (s)',labelpad=label1)
IBIax.set_zticks(I_e3)

plt.savefig('3D_IBI')

plt.figure(15)
Tax.set_xlabel('IpumpMax (nA)',labelpad=label1)
Tax.set_xticks(I_e1)
Tax.set_ylabel('Gp (nS)',labelpad=label1)
Tax.set_yticks(I_e2)
Tax.set_zlabel('Period (s)',labelpad=label1)
Tax.set_zticks(I_e4)

plt.savefig('3D_T')


plt.figure(16)
VNax.set_xlabel('Na')
VNax.set_ylabel('V')
# savefig(fname, dpi=None, facecolor='w', edgecolor='w',
#       orientation='portrait', papertype=None, format=None,
#       transparent=False, bbox_inches=None, pad_inches=0.1,
#       frameon=None, metadata=None)
plt.savefig('V_Na')


plt.ion()
plt.show()