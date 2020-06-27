####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps

# threedee.py
# 	- create 3D plots: gP, Ipm, and:
# 									- BD
# 									- T
# 									- IBI
# 									- mean Hz
# 									- Duty cycle

#### libraries necessary for code
import pandas as pd
from sklearn.metrics import r2_score
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

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
        'size'   : 20}
plt.rc('font', **font)
mk1 = 16
mk2 = 15
mk3 = 3
mk10 = 25
lw1= 6
alpha0 = 1

plot_protocol = 'nm'
# plot_protocol = 'cs'

# coefficient of variance standard: reduce alpha
coef_std1 = 0.25
# coefficient of variane absolute standard: do not plot
covar = 0.3
length_standard = 3


##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
# experiment_list = np.load('experiment_list.npy')
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
    19625000,       #DCC = 0         ######### only I and Vm were saved by pClamp protocol
    19626000,       #DCC = -0.1     ######### only I and Vm were saved by pClamp protocol
    19805000,
    # 19809001,     ######### 
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
    # 20311002,     ################
    20317002
	])
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

	protocol = list10[-4]
	coupling = list10[-3]

	list1 = list10[0:-5]
	bd= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control =='C':
		IpumpMax = list10[-5]
		gp = list10[-5]		
		gp_out = np.append(gp_out,float(gp))
		IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))
		burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
		coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
		bd = np.mean(burst_duration)

	if control == 'I':
		IpumpMax = list10[-5]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			# x_all = int(control_param*10)

			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			bd = np.append(bd,np.mean(burst_duration))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	if control == 'G':
		gp = list10[-5]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#BD
			burst_duration = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(burst_duration)/np.mean(burst_duration)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			bd = np.append(bd,np.mean(burst_duration))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'BD'
	return gp_out, IpumpMax_out, bd, coefficient_of_variation, label

def import_t(list10):
	control = list10[-1]

	protocol = list10[-4]
	coupling = list10[-3]

	list1 = list10[0:-5]
	t= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control =='C':
		IpumpMax = list10[-5]
		gp = list10[-5]		
		gp_out = np.append(gp_out,float(gp))
		IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

		cycle_period = np.load(fh+'BD_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
		coefficient_of_variation = np.std(cycle_period)/np.mean(cycle_period)
		t = np.mean(cycle_period)

	if control == 'I':
		IpumpMax = list10[-5]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			# x_all = int(control_param*10)

			#BD
			period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(period)/np.mean(period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			t = np.append(t,np.mean(period))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	if control == 'G':
		gp = list10[-5]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#BD
			period = np.load(fh+'T_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(period)/np.mean(period)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			t = np.append(t,np.mean(period))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'T'
	return gp_out, IpumpMax_out, t, coefficient_of_variation,label

def import_ibi(list10):
	control = list10[-1]
	list1 = list10[0:-5]

	protocol = list10[-4]
	coupling = list10[-3]

	ibi= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control == 'I':
		IpumpMax = list10[-5]

		for j in range(0,len(list1)):
			gp = float(list1[j])
			# x_all = int(control_param*10)

			#IBI
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			ibi = np.append(ibi,np.mean(interburst_interval))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	if control == 'G':
		gp = list10[-5]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#BD
			interburst_interval = np.load(fh+'IBI_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(interburst_interval)/np.mean(interburst_interval)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			ibi = np.append(ibi,np.mean(interburst_interval))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'IBI'
	return gp_out, IpumpMax_out, ibi, coefficient_of_variation,label

def import_Hz(list10):
	control = list10[-1]

	protocol = list10[-4]
	coupling = list10[-3]

	list1 = list10[0:-5]
	mean_Hz= np.array([])
	gp_out= np.array([])
	IpumpMax_out= np.array([])

	if control == 'I':
		IpumpMax = list10[-5]

		for j in range(0,len(list1)):
			gp = float(list1[j])

			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			mean_Hz = np.append(mean_Hz,np.mean(Hz))
			gp_out = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	if control == 'G':
		gp = list10[-5]

		for j in range(0,len(list1)):
			IpumpMax = float(list1[j])
			#Hz
			Hz = np.load(fh+'Hz_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			coefficient_of_variation = np.std(Hz)/np.mean(Hz)
			alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

			mean_Hz = np.append(mean_Hz,np.mean(Hz))
			gp_oxut = np.append(gp_out,float(gp))
			IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

	label = 'spike f'
	return gp_out, IpumpMax_out, mean_Hz, coefficient_of_variation,label

def import_dc(list10):
    control = list10[-1]
    protocol = list10[-4]
    coupling = list10[-3]
    list1 = list10[0:-5]
    dc= np.array([])
    gp_out= np.array([])
    IpumpMax_out= np.array([])

    if control == 'I':
        IpumpMax = list10[-5]

        for j in range(0,len(list1)):
            gp = float(list1[j])
            # x_all = int(control_param*10)

            #BD
            duty_cycle = np.load(fh+'DC_IpumpIpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            coefficient_of_variation = np.std(duty_cycle)/np.mean(duty_cycle)
            alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

            dc = np.append(dc,np.mean(duty_cycle))
            gp_out = np.append(gp_out,float(gp))
            IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

    if control == 'G':
        gp = list10[-5]

        for j in range(0,len(list1)):
            IpumpMax = float(list1[j])
            #BD
            duty_cycle = np.load(fh+'DC_IpumpIpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
            coefficient_of_variation = np.std(duty_cycle)/np.mean(duty_cycle)
            alpha2, fst1 = coefficient(coefficient_of_variation, covar,coef_std1)

            dc = np.append(dc,np.mean(duty_cycle))
            gp_out = np.append(gp_out,float(gp))
            IpumpMax_out = np.append(IpumpMax_out,float(IpumpMax))

    label = 'DutyCycle'
    return gp_out, IpumpMax_out, dc, coefficient_of_variation, label


def data_vectors(xv,yv,zv,param):
	for i in range(0,len(yv)):
		p1 = np.nonzero(yv==param)
	return xv[p1],yv[p1],zv[p1]

x_I_vec = np.array([])
y_I_vec = np.array([])
z_I_vec = np.array([])

y_G_vec = np.array([])
x_G_vec = np.array([])
z_G_vec = np.array([])


print('Parameters from experiments analyzed:')
for i in range(0,len(experiment_list)):
	
	fileID= str(experiment_list[i])
	fh = fileID+'/'+fileID+'_'
	params = list(np.load(fh+'param.npy'))
	for k in range(0,len(params)):# import file
		list20 = params[k]
		han2 = fh+str(list20)+'.npy'
		list10 =  list(np.load(han2))
		protocol = list10[-4]
		if protocol == plot_protocol:
			print(fileID,list10)

			x, y, z, coefficient_of_variation,label = import_bd(list10)
			# x, y, z, coefficient_of_variation,label = import_t(list10)
			# x, y, z, coefficient_of_variation,label = import_ibi(list10)
			# x, y, z, coefficient_of_variation,label = import_Hz(list10)
			# x, y, z, coefficient_of_variation,label = import_dc(list10)

			if coefficient_of_variation<0.35:
				if list10[-1]=='I':
					x_I_vec = np.append(x_I_vec,x)
					y_I_vec = np.append(y_I_vec,y)
					z_I_vec = np.append(z_I_vec,z)
				if list10[-1]=='G':
					y_G_vec = np.append(y_G_vec,y)
					x_G_vec = np.append(x_G_vec,x)
					z_G_vec = np.append(z_G_vec,z)

BD3 = plt.figure(1,figsize=(75,75))
BDax = plt.axes(projection='3d')

bd_fig = plt.figure(figsize=(100,100))
bd_ax = bd_fig.add_subplot(111,projection='3d')

p_iv = np.array([])
nak_iv = np.array([])
dv = np.array([])

prm = np.unique(y_I_vec)
for i in range(0,len(prm)):
	x,y,z = data_vectors(x_I_vec,y_I_vec,z_I_vec,prm[i])
	p_iv = np.append(p_iv,x)
	nak_iv = np.append(nak_iv,y)
	dv = np.append(dv,z)


qrm = np.unique(x_G_vec)
for i in range(0,len(qrm)):
	a,b,c = data_vectors(y_G_vec,x_G_vec,z_G_vec,qrm[i])
	p_iv = np.append(p_iv,b)
	nak_iv = np.append(nak_iv,a)
	dv = np.append(dv,c)


p_iv0 = np.unique(p_iv)
nak_iv0 = np.unique(nak_iv)

x=np.array([])
y=np.array([])
z=np.array([])

for i in range(0,len(p_iv0)):
	index1 = np.nonzero(p_iv==p_iv0[i])

	nak_iv1 =nak_iv[index1] 
	nak_iv2 = np.unique(nak_iv1)

	dv1 = dv[index1]
	
	for j in range(0,len(nak_iv2)):
		index2 = np.nonzero(nak_iv1 == nak_iv2[j])
		# print(p_iv[index1])
		# print(nak_iv1[index2])
		# print(np.mean(dv1[index2]))
		# print('\n')
		x = np.append(x,p_iv0[i])
		y = np.append(y,nak_iv2[j])
		z = np.append(z,np.mean(dv1[index2]))

###########3

xu = np.unique(x)
for i in range(0,len(xu)):
	sweep = np.nonzero(x==xu[i])
	# print(x[sweep])
	BDax.plot(x[sweep],y[sweep],z[sweep],linewidth=3,color='k')

yu = np.unique(y)
for i in range(0,len(yu)):
	sweep = np.nonzero(y==yu[i])
	# print(y[sweep])
	BDax.plot(x[sweep],y[sweep],z[sweep],linewidth=3,color='k')


ytcks = np.arange(0.2,0.9,0.2)
xtcks = np.arange(0,10,2)

# ztcks = np.arange(2,np.max(z),int(np.max(z)/4))		# BD,IBI, T, spike f
ztcks = np.arange(0.2,1.1,0.2)							# duty cycle

top_angle = 370
fts2 = 30
padl = 30

BDax.plot(x,y,z,'.',markersize=20)
BDax.set_ylabel('IpumpMax (nA)',fontsize=fts2,labelpad = padl)
BDax.set_yticks(ytcks)
BDax.set_xlabel('gP max (nS)',fontsize=fts2,labelpad = padl)
BDax.set_xticks(xtcks)
# BDax.set_zlabel(label+' (s)',fontsize=fts2,labelpad = padl)			# BD, IBI, T
BDax.set_zlabel(label,fontsize=fts2,labelpad = padl)					# duty cycle
# BDax.set_zlabel(label+' (Hz)',fontsize=fts2,labelpad = padl)			# spike frequency
BDax.set_zticks(ztcks)

bd_ax.plot_trisurf(x,y,z,cmap='Spectral',alpha=0.5)
bd_ax.set_ylabel('IpumpMax (nA)',fontsize=fts2,labelpad = padl)
bd_ax.set_yticks(ytcks)
bd_ax.set_xlabel('gP max (nS)',fontsize=fts2,labelpad = padl)
bd_ax.set_xticks(xtcks)
# bd_ax.set_zlabel(label+' (s)',fontsize=fts2,labelpad = padl)			# BD, IBI, T
bd_ax.set_zlabel(label,fontsize=fts2,labelpad = padl)					# duty cycle
# bd_ax.set_zlabel(label+' (Hz)',fontsize=fts2,labelpad = padl)			# spike frequency

bd_ax.set_zticks(ztcks)


# ##########3
for angle in np.arange(0,top_angle,10):
	plt.figure(1)
	# BDax.view_init(40,angle)
	BDax.view_init(45,angle)
	plt.draw()
	plt.pause(.0001)

	plt.savefig(label+'_mesh_'+str(angle))


for angle in np.arange(0,top_angle,10):
	plt.figure(2)
	# bd_ax.view_init(40,angle)
	bd_ax.view_init(45,angle)
	plt.draw()
	plt.pause(.0001)

	plt.savefig(label+'_surface_'+str(angle))


