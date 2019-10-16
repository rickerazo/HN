# ANOVA.py

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.graphics.factorplots import interaction_plot

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

####	graphix stuff
markadores = ['0', '2', 'x', '1','4','3','p','+','h','|','.', '2', 'x', '1','4','3','p','+','h','|']
colores = ['purple', 'g', 'r', 'c','m','k','y','b', 'g', 'gray','orange', 'g', 'r', 'c','m','k','y','b', 'g', 'b']
font = {'weight' : 'bold',
        'size'   : 50}
plt.rc('font', **font)
# coefficient of variance standard: reduce alpha
coef_std1 = 0.25
# coefficient of variane absolute standard: do not plot
covar = 0.95
length_standard = 3

def coefficient(coefficient_of_variation, covar,coef_std1):
	if coefficient_of_variation < covar:
		if coefficient_of_variation>coef_std1:
			alpha2 = 0.3
			fst1 = 'none'
		else:
			alpha2 = 1
			fst1 = 'full'
	else:
		alpha2 = 0.0001
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


def eta_squared(aov):
    aov['eta_sq'] = 'NaN'
    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    return aov

def omega_squared(aov):
    mse = aov['sum_sq'][-1]/aov['df'][-1]
    aov['omega_sq'] = 'NaN'
    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*mse))/(sum(aov['sum_sq'])+mse)
    return aov

def sum_of_squares(x_vec,grand_mean,y_vec):
	x2 = np.unique(x_vec)
	for i in range(0,len(x2)):
		p1 = np.nonzero(x_vec==x2[i])
		p2 = np.square([y_vec[p1] - grand_mean])
		ssq_x = sum(p2[0])
	return ssq_x

def sum_of_squares_w(x_vec,y_vec):
	x2 = np.unique(x_vec)
	for i in range(0,len(x2)):
		p1 = np.nonzero(x_vec==x2[i])
		p2 = np.square([y_vec[p1] - np.mean(y_vec[p1])])
		ssq_x = sum(p2[0])
	return ssq_x

##### file ID stuff
cwd = os.getcwd()
cwd = cwd+'/'
# experiment_list = np.array([
# 	# 18831003,
# 	# 18902004,
# 	# 18914000,
# 	# 18929002,
# 	19522001,
# 	19529000,
# 	# # 19603000,
# 	# # 19607000,
# 	# # 19607001,
# 	# # 19612001,
# 	19612002,
# 	19614000,
# 	19625000,
# 	19626000,
# 	19626002,
# 	19731000,
# 	19731001,
# 	19805000,
# 	19809001,
# 	19826000
# 	])
experiment_list = np.load('experiment_list.npy')
cwd = os.getcwd()
x1_I_vec = np.array([])
x2_I_vec = np.array([])
y_I_vec = np.array([])

x1_G_vec = np.array([])
x2_G_vec = np.array([])
y_G_vec = np.array([])
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


fig_I = interaction_plot(x1_I_vec, x2_I_vec, y_I_vec,xlabel ='GP' , ylabel =label,linewidth=5, legendtitle='IpumpMax',legendloc='best')# ,colors=['purple', 'g', 'r', 'c','m','k','y','orange'], markers=['1', '2', 'x', '^','4','3','p','*'], ms=10)

fig_G = interaction_plot(x2_G_vec, x1_G_vec, y_G_vec,xlabel ='IpumpMax' , ylabel =label,linewidth=5, legendtitle='GP',legendloc='best')# ,colors=['purple', 'g', 'r', 'c','m','k','y'], markers=['1', '2', 'x', '^','4','3','p'], ms=10)
# fig_G = interaction_plot(x2_G_vec, x1_G_vec, y_G_vec,xlabel ='IpumpMax' , ylabel =label)# ,colors=['purple', 'g', 'r', 'c','m','k','y','orange','gray'], markers=['1', '2', 'x', '^','4','3','p','*','v'], ms=10)

# plt.figure(figsize=(50,50))
# plt.plot(x1_I_vec, y_I_vec,'.')


# plt.figure(figsize=(50,50))
# plt.plot(x2_G_vec, y_G_vec,'.')

plt.ion()
plt.show()

#######################

def ANOVA(x1_vec,x2_vec,y_vec,levels_df):
	# degrees of freedom
	N = len(y_vec)
	df_x1 = len(np.unique(x1_vec)) - 1
	df_x2 = len(np.unique(x2_vec)) - 1
	df_interaction = df_x1*df_x2 
	# df_w = N - (len(np.unique(x1_I_vec))*len(np.unique(x2_I_vec)))
	df_w = N - (levels_df*len(np.unique(x2_vec))) #*


	#sum of squares
	grand_mean = np.mean(y_vec)
	ssq_x2 = sum_of_squares(x2_vec,grand_mean,y_vec)
	ssq_x1 = sum_of_squares(x1_vec,grand_mean,y_vec)

	ssq_total =sum( np.square(y_vec-grand_mean))

	ssq_x2_w = sum_of_squares_w(x2_vec,y_vec)
	ssq_x1_w = sum_of_squares_w(x1_vec,y_vec)
	ssq_w = ssq_x2_w+ssq_x1_w

	# sum of squares interaction
	ssq_interaction = ssq_total - ssq_x1- ssq_x2 - ssq_w
	# Mean squares
	ms_x1 = ssq_x1/df_x1
	ms_x2 = ssq_x2/df_x2
	ms_x1x2 = ssq_interaction/df_interaction
	ms_w = ssq_w/df_w

	## F-ratio
	f_x1 = ms_x1/ms_w
	f_x2 = ms_x2/ms_w
	f_x1x2 = ms_x1x2/ms_w

	p_x1 = stats.f.sf(f_x1,df_x1,df_w)
	p_x2 = stats.f.sf(f_x2,df_x2,df_w)
	p_x1x2 = stats.f.sf(f_x1x2,df_interaction,df_w)


	results = {'sum_sq':[ssq_x1, ssq_x2, ssq_interaction, ssq_w],
	           'df':[df_x1, df_x2, df_interaction, df_w],
	           'F':[f_x1, f_x2, f_x1x2, 'NaN'],
	            'PR(>F)':[p_x1, p_x2, p_x1x2, 'NaN']}
	columns=['sum_sq', 'df', 'F', 'PR(>F)']

	aov_table1 = pd.DataFrame(results, columns=columns,
	                          index=['gp', 'IpumpMax', 
	                          'gp/Ipump', 'Residual'])


	eta_squared(aov_table1)
	omega_squared(aov_table1)
	print(aov_table1)

print('\n'+'Constant IpumpMax')
ANOVA(x1_I_vec,x2_I_vec,y_I_vec,3)

print('\n'+'Constant Gp')
ANOVA(x2_G_vec,x1_G_vec,y_G_vec,2)