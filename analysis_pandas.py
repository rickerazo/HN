####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps

# analysis_pandas.py
# 	- scans constant and sweep parameters
# 	- returns a list of parameters of performed experiments
# 	- saves experiment list

#### libraries necessary for code
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import itertools
import os
import sys
####	graphix stuff
plt.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['axes.linewidth']=10
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

markers = ['1', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia']

font = {'weight' : 'bold',
        'size'   : 70}
plt.rc('font', **font)
mk1 = 16
mk2 = 15
mk3 = 3
mk10 = 25
lw1= 6
alpha0 = 1

save_eps = False

# save_plot = 'cs'
save_plot = 'nm'

# coefficient of variance standard: reduce alpha
coef_std1 = 0.25
# coefficient of variane absolute standard: do not plot
covar = 0.45
length_standard = 3

##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
experiment_list = np.array([
    18902004,
    18914000,
    18929002,
    19522001,
    19529000,
    19603000,
    19607000,       #DCC = 0.0
    19607001,       #DCC = -0.05
    19614000,       #DCC = -0.1
    19625000,       #DCC = 0
    19626000,       #DCC = -0.1
    19805000,
    19809001,
    19826000,
    19918001,
    19918002,
    19919000,
    19930000,
    19930001,
    '19n05001',
    '19n06000',
    '19d12000',
    '19d12001',
    20204006,
    20205003,
    20210002,
    20225002,
    20226001,
    20303001,
    20310000,
    20311002,
    20317002
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
		alpha2 = 0.002
		fst1 = 'none'
	else:
		alpha2 = 0.001
		fst1 = 'none'
	return alpha2, fst1

def compute(list10):
	constant1 = list10[-1]
	list0 = list10[0:-5]
	protocol = list10[-4]
	coupling = list10[-3]
	Ibias = float(list10[-2])
	####################################
	if protocol == 'nm':
		mark_id = 'o'
	elif protocol == 'cs':
		mark_id = 'D'

	####################################
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
	# Ibias = np.load(fh+'Ibias.npy')
	if Ibias>0.1:
		constant1='Z'

	if constant1=='I' and protocol == save_plot:
	       # if constant1=='I' and protocol == save_plot and float(list10[-5])<0.5:
		IpumpMax = list10[-5]
		IpumpMax = float(IpumpMax)

		list1 = list0[::-1]

		control_param = IpumpMax


		for j in range(0,len(list1)):
			gp = float(list1[j])
			x_all = int(control_param*10)
			# print(x_all,'\n')
			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			# plot = fig.add_subplot(111)
			y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
			x = np.zeros((1,len(y1)))+gp
			plot1.plot(x[0],y1,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot1.plot(gp,np.mean(burst_duration),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot1.plot(gp,np.mean(burst_duration)+np.std(burst_duration),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot1.plot(gp,np.mean(burst_duration)-np.std(burst_duration),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				gp_list1.append(gp)
				BD_mean_list.append(np.mean(burst_duration))#/np.max(burst_duration))

			#IBI
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y2 = np.array((np.mean(interburst_interval)-np.std(interburst_interval),np.mean(interburst_interval),np.mean(interburst_interval)+np.std(interburst_interval)))
			x = np.zeros((1,len(y2)))+gp
			plot2.plot(x[0],y2,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot2.plot(gp,np.mean(interburst_interval),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot2.plot(gp,np.mean(interburst_interval)+np.std(interburst_interval),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot2.plot(gp,np.mean(interburst_interval)-np.std(interburst_interval),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				gp_list2.append(gp)
				IBI_mean_list.append(np.mean(interburst_interval))#/np.max(cycle_period))


			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
			x = np.zeros((1,len(y3)))+gp
			plot3.plot(x[0],y3,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot3.plot(gp,np.mean(cycle_period),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot3.plot(gp,np.mean(cycle_period)+np.std(cycle_period),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot3.plot(gp,np.mean(cycle_period)-np.std(cycle_period),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				gp_list3.append(gp)
				T_mean_list.append(np.mean(cycle_period))#/np.max(cycle_period))

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y4 = np.array((np.mean(Hz)-np.std(Hz), np.mean(Hz) ,np.mean(Hz)+np.std(Hz) ))
			x = np.zeros((1,len(y4)))+gp
			plot4.plot(x[0],y4,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot4.plot(gp,np.mean(Hz),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot4.plot(gp,np.mean(Hz)+np.std(Hz),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot4.plot(gp,np.mean(Hz)-np.std(Hz),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				gp_list4.append(gp)
				Hz_mean_list.append(np.mean(Hz))



		plot1.plot(gp_list1, BD_mean_list,color=colores[x_all])
		plot2.plot(gp_list2, IBI_mean_list,color=colores[x_all])
		plot3.plot(gp_list3, T_mean_list,color=colores[x_all])
		plot4.plot(gp_list4, Hz_mean_list,color=colores[x_all])



	# if constant1=='G' and protocol == save_plot and float(list10[-5])>4.0:
	if constant1=='G' and protocol == save_plot:
		gp = list10[-5]
		gp = float(gp)

		control_param=gp

		list1 = list0[::-1]


		for j in range(0,len(list1)):
			IpumpMax = list1[j]
			IpumpMax = float(IpumpMax)
			x_all = int(control_param)
			# print(x_all,'\n')
			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y1 = np.array((np.mean(burst_duration)-np.std(burst_duration),np.mean(burst_duration),np.mean(burst_duration)+np.std(burst_duration)))
			x = np.zeros((1,len(y1)))+IpumpMax
			plot5.plot(x[0],y1,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot5.plot(IpumpMax,np.mean(burst_duration),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot5.plot(IpumpMax,np.mean(burst_duration)+np.std(burst_duration),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot5.plot(IpumpMax,np.mean(burst_duration)-np.std(burst_duration),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				pump_list1.append(IpumpMax)
				BD_mean_list.append(np.mean(burst_duration))

			#IBI
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y2 = np.array((np.mean(interburst_interval)-np.std(interburst_interval),np.mean(interburst_interval),np.mean(interburst_interval)+np.std(interburst_interval)))
			x = np.zeros((1,len(y2)))+IpumpMax
			plot6.plot(x[0],y2,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot6.plot(IpumpMax,np.mean(interburst_interval),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot6.plot(IpumpMax,np.mean(interburst_interval)+np.std(interburst_interval),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot6.plot(IpumpMax,np.mean(interburst_interval)-np.std(interburst_interval),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				pump_list2.append(IpumpMax)
				IBI_mean_list.append(np.mean(interburst_interval))

			#T
			cycle_period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y3 = np.array((np.mean(cycle_period)-np.std(cycle_period),np.mean(cycle_period),np.mean(cycle_period)+np.std(cycle_period)))
			x = np.zeros((1,len(y3)))+IpumpMax
			plot7.plot(x[0],y3,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot7.plot(IpumpMax,np.mean(cycle_period),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot7.plot(IpumpMax,np.mean(cycle_period)+np.std(cycle_period),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot7.plot(IpumpMax,np.mean(cycle_period)-np.std(cycle_period),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				pump_list3.append(IpumpMax)
				T_mean_list.append(np.mean(cycle_period))


			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			y4 = np.array((np.mean(Hz)-np.std(Hz),np.mean(Hz),np.mean(Hz)+np.std(Hz)))
			x = np.zeros((1,len(y4)))+IpumpMax
			plot8.plot(x[0],y4,markersize=mk2,linewidth=lw1,color=colores[x_all],alpha=alpha2,fillstyle=fst1)
			plot8.plot(IpumpMax,np.mean(Hz),marker=mark_id,markersize=mk10,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot8.plot(IpumpMax,np.mean(Hz)+np.std(Hz),marker='^',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			plot8.plot(IpumpMax,np.mean(Hz)-np.std(Hz),marker='v',markersize=mk10/2,fillstyle=fst1,alpha=alpha2,color = colores[x_all])
			if coefficient_of_variation < covar:
				pump_list4.append(IpumpMax)
				Hz_mean_list.append(np.mean(Hz))

		plot5.plot(pump_list1, BD_mean_list,color=colores[x_all])
		plot6.plot(pump_list2, IBI_mean_list,color=colores[x_all])
		plot7.plot(pump_list3, T_mean_list,color=colores[x_all])
		plot8.plot(pump_list4, Hz_mean_list,color=colores[x_all])
	return float(list10[-5])

print('Parameters from experiments analyzed:\n')
print('Constant Ipm, sweep gP')
print('nA 	: 	nS')

col_counter = 0
for i in range(0,len(experiment_list)):

    fileID= str(experiment_list[i])
    fh = fileID+'/'+fileID+'_'
    params = list(np.load(fh+'param.npy'))
    # out = {}
    for k in range(0,len(params)):# import file
        list20 = params[k]
        han2 = fh+str(list20)+'.npy'
        list10 =  list(np.load(han2))
        protocol = str(list10[-4])

        list20 = {list10[-5]:list10[0:-5]}

        if list10[-1] == 'I' and protocol == save_plot:
            # print(fileID,',',list10)
            list20[list10[-5]]:list10[0:-5]
            out=list20
            df = pd.DataFrame(list20)
            print(out)


print('\nConstant gP, sweep Ipm')
print('nS 	: 	nA')

for i in range(0,len(experiment_list)):

    fileID= str(experiment_list[i])
    fh = fileID+'/'+fileID+'_'
    params = list(np.load(fh+'param.npy'))
    # out = {}
    for k in range(0,len(params)):# import file
        list20 = params[k]
        han2 = fh+str(list20)+'.npy'
        list10 =  list(np.load(han2))
        protocol = str(list10[-4])

        list20 = {list10[-5]:list10[0:-5]}

        if list10[-1] == 'G' and protocol == save_plot:
            list20[list10[-5]]:list10[0:-5]
            out=list20
            df = pd.DataFrame(list20)
            print(out)
