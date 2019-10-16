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
from sklearn.metrics import r2_score
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import itertools
import os
import sys
####	graphix stuff
matplotlib.rcParams['axes.linewidth']=10


markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['purple', 'g', 'r', 'c','m','k','y','b', 'pink','orange' ]


font = {'weight' : 'bold',
        'size'   : 50}
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
	# 19918002, #b
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
	# elif coefficient_of_variation >= covar and coefficient_of_variation < 0.5:
	# 	alpha2 = 0.2
	# 	fst1 = 'none'
	else:
		alpha2 = 0.001
		fst1 = 'none'
	return alpha2, fst1

def import_bd(list10):
	control = list10[-1]
	list1 = list10[0:-2]
	bd= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control == 'I':
		IpumpMax = list10[-2]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			# x_all = int(control_param*10)

			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			bd = np.append(bd,np.mean(burst_duration))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))
	
	if control == 'G':
		gp = list10[-2]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			bd = np.append(bd,np.mean(burst_duration))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'BD'
	return gp_out, IpumpMax_out, bd, coefficient_of_variation, label

def import_t(list10):
	control = list10[-1]
	list1 = list10[0:-2]
	t= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control == 'I':
		IpumpMax = list10[-2]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			# x_all = int(control_param*10)

			#BD
			period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(period)/np.mean(period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			t = np.append(t,np.mean(period))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))
	
	if control == 'G':
		gp = list10[-2]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#BD
			period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(period)/np.mean(period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			t = np.append(t,np.mean(period))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'T'
	return gp_out, IpumpMax_out, t, coefficient_of_variation,label

def import_ibi(list10):
	control = list10[-1]
	list1 = list10[0:-2]
	ibi= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control == 'I':
		IpumpMax = list10[-2]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			# x_all = int(control_param*10)

			#IBI
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			ibi = np.append(ibi,np.mean(interburst_interval))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))
	
	if control == 'G':
		gp = list10[-2]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#BD
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			ibi = np.append(ibi,np.mean(interburst_interval))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'IBI'
	return gp_out, IpumpMax_out, ibi, coefficient_of_variation,label

def import_Hz(list10):
	control = list10[-1]
	list1 = list10[0:-2]
	mean_Hz= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control == 'I':
		IpumpMax = list10[-2]

		for j in range(0,len(list1)):
			gp = float(list1[j])

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			mean_Hz = np.append(mean_Hz,np.mean(Hz))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))
	
	if control == 'G':
		gp = list10[-2]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			mean_Hz = np.append(mean_Hz,np.mean(Hz))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'Hz'
	return gp_out, IpumpMax_out, mean_Hz, coefficient_of_variation,label


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

# for i in range(0,len(experiment_list)):
	
# 	fileID= str(experiment_list[i])
# 	# fh1 = cwd+fileID+'/'
# 	fh = fileID+'/'+fileID+'_'
# 	# print(cwd+str(experiment_list[i]))
# 	# sys.path.insert(0,fh1)
# 	params = list(np.load(fh+'param.npy'))
# 	# params = ['list1']
# 	for k in range(0,len(params)):# import file
# 		list20 = params[k]
# 		han2 = fh+str(list20)+'.npy'
# 		list10 =  list(np.load(han2))
# 		print(fileID,list10)
# 		# print(list10[-1])
# 		if list10[-2]=='0.2':
# 			control_param = compute(list10)
# 			if list10[-1] == 'I':
# 				riff_pump=np.append(riff_pump,control_param)
# 			if list10[-1] == 'G':
# 				riff_gp = np.append(riff_gp,control_param)

x1_I_vec = np.array([])
x2_I_vec = np.array([])
y_I_vec = np.array([])

x1_G_vec = np.array([])
x2_G_vec = np.array([])
y_G_vec = np.array([])
for i in range(0,len(experiment_list)):
	
	fileID= str(experiment_list[i])
	fh = fileID+'/'+fileID+'_'
	params = list(np.load(fh+'param.npy'))
	for k in range(0,len(params)):# import file
		list20 = params[k]
		han2 = fh+str(list20)+'.npy'
		list10 =  list(np.load(han2))
		print(fileID,list10)
		if list10[-2]=='0.2':
			x1, x2, y1, coefficient_of_variation,label = import_bd(list10)
			# x1, x2, y1, coefficient_of_variation,label = import_t(list10)
			# x1, x2, y1, coefficient_of_variation,label = import_ibi(list10)
			# x1, x2, y1, coefficient_of_variation,label = import_Hz(list10)
			if coefficient_of_variation<0.25:
				if list10[-1]=='I':
					x1_I_vec = np.append(x1_I_vec,x1)
					x2_I_vec = np.append(x2_I_vec,x2)
					y_I_vec = np.append(y_I_vec,y1)
				if list10[-1]=='G':
					x1_G_vec = np.append(x1_G_vec,x1)
					x2_G_vec = np.append(x2_G_vec,x2)
					y_G_vec = np.append(y_G_vec,y1)

###############################################################################################################
###############################################################################################################
plot0=False
if plot0 == 'True':
	for i in range(1,9):
		plt.figure(i,figsize=(45,30))

	plt.figure(1)
	plt.ylabel('Burst duration (s)',fontsize=65,weight='bold')
	plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(1,12,3))
	# plt.axis([0,10,0,9])

	# plt.axis([0,10,0,9])
	plt.figure(2)
	plt.ylabel('IBI',fontsize=65,weight='bold')
	plt.xlabel(r'G$_P$ (nS)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(1,12,3))


	plt.figure(3)
	plt.ylabel('Period (s)',fontsize=65,weight='bold')
	plt.xlabel(r'$G_P$ (nS)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(1,21,4))
	# plt.axis([0,1,0,9])

	plt.figure(4)
	plt.ylabel('Spike frequency (Hz)',fontsize=65,weight='bold')
	plt.xlabel(r'$G_P$ (nS)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(0,60,10))
	# plt.ylim(0,50)

	plt.figure(5)
	plt.ylabel('Burst duration (s)',fontsize=65,weight='bold')
	plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(1,12,3))

	plt.figure(6)
	plt.ylabel('IBI',fontsize=65,weight='bold')
	plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(1,12,3))

	plt.figure(7)
	plt.ylabel('Period (s)',fontsize=65,weight='bold')
	plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(1,21,4))

	plt.figure(8)
	plt.ylabel('Spike frequency (Hz)',fontsize=65,weight='bold')
	plt.xlabel(r'$I_{pump}^{max}$ (nA)',fontsize=65,weight='bold')
	plt.yticks(ticks=np.arange(0,60,10))

#############################################################################################################

	markers = [' ', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
	colores = ['purple', 'g', 'r', 'c','m','k','y','b', 'pink','orange' ]

	riff2 = np.unique(riff_pump)

	# color_riff2= [colores[1],colores[2],colores[3],colores[4],colores[5],colores[6],colores[7],colores[8],colores[9]]#,colores[6],colores[8],colores[9]]
	color_riff2= [colores[2]]
	titles = [' ','BD','IBI','T','Hz','Ipump_cytNa','BD','IBI','T','Hz','Ipump_cytNa']
	#gp 
	for i in range(1,5):
		plt.figure(i)
		for j in range(0,len(riff2)):
			p1 = riff2[j]
			str1 = str(p1)
			plt.plot(0,0,'o',markersize=1e-3,color=color_riff2[j],label=r'I$_{pump}^{max}$= '+str1+' nA')
			plt.legend(loc=2,bbox_to_anchor=(0.,1.1),markerscale=2e+4,ncol=4,fontsize=40)
			# plt.savefig('gp_'+str(0.5)+'_'+titles[i]+'.png')

	#		4,5,5.5,6,6.5,7
	riff1 = np.unique(riff_gp)
	color_riff1 = [colores[1],colores[4],colores[5],colores[5],colores[6], colores[6],colores[7],colores[8],colores[9]]

	#pump
	for i in range(5,9):
		plt.figure(i)
		for j in range(0,len(riff1)):
			p1 = riff1[j]
			str1 = str(p1)
			plt.plot(0,0,'o',markersize=1e-3,color=color_riff1[j],label=r'G$_P$= '+str1+' nS')
			plt.legend(loc=2,bbox_to_anchor=(0.,1.1),markerscale=2e+4,ncol=4,fontsize=50)
			# plt.savefig('pump_'+titles[i-4])



################## LINEAR REGRESSION 
plt.figure(figsize=(20,20))
plt.plot(x1_I_vec,y_I_vec,'.',markersize=mk1)

slope, intercept, rvalue, pvalue, stderr = linregress(x1_I_vec,y_I_vec)

plt.plot(x1_I_vec,x1_I_vec*slope+intercept)
plt.text(1,5.5,'r='+str(rvalue),fontsize=25)
plt.text(1,5,'p='+str(pvalue),fontsize=25)
plt.text(1,4.5,'std error='+str(stderr),fontsize=25)

# plt.text()

plt.title(r'Linear regression, $I_{pump}^{max}$= 0.2')
plt.ylabel('BD (s)')
plt.xlabel('GP(nS)')

plt.savefig('0.2_linearegression_BD.png')
# residuals
def plot_residuals(x_vec,y_vec,predicted):
	plt.figure(figsize=(20,20))
	plt.plot(x_vec,y_vec,'.',markersize=mk1)
	plt.plot(x_vec,predicted,'.',markersize=mk1)
	for i in range(0,len(x_vec)):
		x = x_vec[i]
		y_pred = predicted[i]
		y_data = y_vec[i]
		plt.plot([x,x],[y_data, y_pred],linewidth=2,color='red')
	plt.ylabel(label)
# plt.plot(x1_I_vec,(x1_I_vec*slope+intercept)-y_I_vec,'.')
plot_residuals(x1_I_vec,y_I_vec,x1_I_vec*slope+intercept)
plt.title('residuals: linear')
plt.ylabel('BD')
plt.xlabel('GP')


################## POLYNOMIAL REGRESSION 

#quadratic function
def quadratic_poly(gp,a,b,c):
	A = a*pow(gp,2)+b*gp+c
	return A

def exponential_poly(gp,a,b):
	A = np.exp(b)* np.exp(a*gp)
	return A

# def squareroot_poly(gp,)	

popt_q,pcov_q = curve_fit(quadratic_poly, x1_I_vec,y_I_vec)
popt_e,pcov_e = curve_fit(exponential_poly, x1_I_vec,y_I_vec)

fit_q = quadratic_poly(x1_I_vec,*popt_q)
fit_e = exponential_poly(x1_I_vec,*popt_e) #need shift

# plt.figure(figsize=(20,20))
# plt.plot(x1_I_vec,y_I_vec,'.',markersize=15,label='raw data')
# plt.plot(x1_I_vec,fit_q,'.',markersize=15,label='quadratic')
plt.plot(x1_I_vec,fit_e,'.',markersize=15,label='exponential')

# plot_residuals(x1_I_vec,y_I_vec,fit_q)
# plt.title('residuals: quadratic')
# plt.ylabel('BD')
# plt.xlabel('GP')

# plot_residuals(x1_I_vec,y_I_vec,fit_e)
# plt.title('residuals : exponential')
# plt.ylabel('BD')
# plt.xlabel('GP')
####

a = np.polyfit((x1_I_vec),y_I_vec,0.5)
plt.figure(10)
plt.plot(x1_I_vec,(a*pow(x1_I_vec,0.5)),'.',label='sqrt',markersize=mk10)
# plot_residuals(x1_I_vec,y_I_vec,1/(a*pow(x1_I_vec,0.5)))


a,b,c = np.polyfit(x1_I_vec,y_I_vec, 2)
plt.figure(10,figsize=(20,20))
plt.plot(x1_I_vec,y_I_vec,'.',markersize=mk10)
plt.plot(x1_I_vec,a*pow(x1_I_vec,2)+b*x1_I_vec+c,'.',label='quad',markersize=mk10)
plt.title('polyfit')
r2_quad = r2_score(y_I_vec,a*pow(x1_I_vec,2)+b*x1_I_vec+c)

a,b,c,d = np.polyfit(x1_I_vec,y_I_vec,3)
plt.figure(10)
plt.plot(x1_I_vec,a*pow(x1_I_vec,3)+b*pow(x1_I_vec,2)+c*pow(x1_I_vec,1)+d,'.',label='cube',markersize=mk10)
r2_cube = r2_score(y_I_vec,a*pow(x1_I_vec,3)+b*pow(x1_I_vec,2)+c*pow(x1_I_vec,1)+d)

a,b,c,d,e = np.polyfit(x1_I_vec,y_I_vec,4)
plt.plot(x1_I_vec,a*pow(x1_I_vec,4)+b*pow(x1_I_vec,3)+c*pow(x1_I_vec,2)+d*pow(x1_I_vec,1)+e,'.',label='4',markersize=mk10)
r2_cuatro = r2_score(y_I_vec,a*pow(x1_I_vec,4)+b*pow(x1_I_vec,3)+c*pow(x1_I_vec,2)+d*pow(x1_I_vec,1)+e)


# a,b,c,d,e,f = np.polyfit(x1_I_vec,y_I_vec,5)
# plt.plot(x1_I_vec,a*pow(x1_I_vec,5)+b*pow(x1_I_vec,4)+c*pow(x1_I_vec,3)+d*pow(x1_I_vec,2)+e*pow(x1_I_vec,1)+f,'.',label='5',markersize=mk10)
# r2_penta = r2_score(y_I_vec,a*pow(x1_I_vec,5)+b*pow(x1_I_vec,4)+c*pow(x1_I_vec,3)+d*pow(x1_I_vec,2)+e*pow(x1_I_vec,1)+f)



# a,b,c,d,e,f,g = np.polyfit(x1_I_vec,y_I_vec,6)
# plt.plot(x1_I_vec,a*pow(x1_I_vec,6)+b*pow(x1_I_vec,5)+c*pow(x1_I_vec,4)+d*pow(x1_I_vec,3)+e*pow(x1_I_vec,2)+f*pow(x1_I_vec,1)+g,'.',label='6',markersize=mk10)
# r2_6 = r2_score(y_I_vec,a*pow(x1_I_vec,6)+b*pow(x1_I_vec,5)+c*pow(x1_I_vec,4)+d*pow(x1_I_vec,3)+e*pow(x1_I_vec,2)+f*pow(x1_I_vec,1)+g)
# plt.ylim(0,6)
plt.ion()
plt.show()
