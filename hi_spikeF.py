####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################
# Code developed in python 3.6.8: objective is to import saved data from traces.py to compare all data from dynamic clamp sweeps
# BURST analysis -> discover if there is a bifurcation of some kind in the data

############# experiment by experiment - manual clustering
#   Using mean spike F, discriminate between experiments in which HN neurons went into high-spike frequency bursting mode
#   and which experiments did not
#   IMPORTANT: The reason to use mean spike frequency is because spiking frequency is well-defined and easy to understand
#   Also important: Vm excursion and spiking frequency are strongly positively correlated.
#   The goals of this script are:
#       1. analyze data prep-by-prep
#          1.1 identifiy what characteristics of spiking frequency separate between HI/LO frequency bursting modes
#          1.2 use these characteristics to cluster data
#       2. apply clustering
#           2.1 all experiments
#           2.2 experiment-by-experiment
#
#       output saved to data_analyses/ebe/*png
#
#       EXPERIMENT-BY-EXPERIMENT
#
#
####  Imported data: only if 3 successful repeated experiments
#
###    - INPUT: curated experiment list
###     - Novel: running mean data filter (Savitzky-Golay) to clean Nai and detect peaks
############################################################################

#### libraries necessary for code
# from scipy.signal import savgol_filter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import OrderedDict
from scipy.signal import find_peaks
# from sklearn.cluster import KMeans
# from sklearn.metrics import silhouette_score
# from scipy.stats import linregress
# from scipy.stats.stats import pearsonr
from mpl_toolkits import mplot3d
# from scipy.stats import ttest_ind
# from statsmodels.sandbox.stats.multicomp import multipletests 

import get_data
import inference_plots
import custom_analysis

import pandas as pd
import os
import sys

cmaps = OrderedDict()

##############################################################################
####	graphix stuff
mpl.rcParams['axes.linewidth']=5


markers = ['1', '2', 'x0', '1','4','3','p','+','h','|','.', '2', 'x0', '1','4','3','p','+','h','|']
colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k']
colores_binario = ['red','blue']
# colores = ['purple','lawngreen','gold','steelblue','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','firebrick','fuchsia']


font = {'weight' : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rcParams['agg.path.chunksize'] = 10000

plt.close('all')

#### data analysis parameters
## coefficient of variance standard: reduce alpha
coef_std1 = 0.25
## coefficient of variane absolute standard: do not plot
covar = 0.3
length_standard = 3

minimum_spikes = 2

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
    19625000,       #DCC = 0         ######### only I and Vm were saved by pClamp protocol
    19626000,       #DCC = -0.1     ######### only I and Vm were saved by pClamp protocol
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
    # 20226001,         ###### defintely not rhythmic at all
    20303001,
    20310000,
    20311002,     ################
    20317002, 
    20730000,
    20804002,
    20806000,
    20806003,
    # 20813000      ######### outlier
    20902002,
    20903000,
    20909000,
    '20o14000',
    '20o14001',
    '20o26000',
    '20o30000'
	])
np.save('experiment_list',experiment_list)

############## IMPORTANT : minimum 3 successful experiment repeats to consider the data for subsequent analysis

dir0 = 'data_analyses/dependence_characteristics/vars/'
dir1 = 'data_analyses/dependence_characteristics/tests/'
dir2 = 'data_analyses/dependence_characteristics/'
dir3 = 'data_analyses/dependence_characteristics/regressions0/'
dir4 = 'data_analyses/dependence_characteristics/regressions1/'

dir10 = 'data_analyses/normalized_vars/'
### for gP in np.arange(5.0, 7.0, 2.0):

# for gP in np.arange(5.0, 7.0, 2.0):
# for gP in np.arange(6.0, 7.0, 2.0):
for gP in np.arange(5.0, 7.0, 1.0):
    tab1 = pd.DataFrame()

    f10 = plt.figure(figsize=(20,10))
    ax10 = f10.add_subplot(121)
    ax20 = f10.add_subplot(122)

    f11 = plt.figure(figsize=(20,10))
    ax11 = f11.add_subplot(121)
    ax21 = f11.add_subplot(122)

    f12 = plt.figure(figsize=(20,10))
    ax12 = f12.add_subplot(121)
    ax22 = f12.add_subplot(122)

    f14 = plt.figure(figsize=(10,10))
    ax14 = f14.add_subplot(111)

    f15 = plt.figure(figsize=(10,10))
    ax15 = f15.add_subplot(111)

    f17 = plt.figure(figsize=(20,10))
    ax17 = f17.add_subplot(121)
    ax27 = f17.add_subplot(122)


    countr0= 1
    countr1= 1
    countr2= 1


    # protokol = 'nm'
    protokol = 'cs'

    compute_raw_data= True
    # compute_raw_data= False
    if compute_raw_data==True:
        print('\n', 'gP ='+str(gP))
        print('Parameters from experiments analyzed:')
        id0 = np.array([])

        x0=np.array([])
        y0=np.array([])
        a0=np.array([])
        b0=np.array([])
        c0=np.array([])

        z0=np.array([])
        w0=np.array([])
        u0=np.array([])
        t0=np.array([])


        id0=np.array([])
        # e0=np.array([])
        # f0=np.array([])
        # g0=np.array([])
        h0=np.array([])

        bd_z = np.array([])
        ibi_z = np.array([])
        exp_ID = np.array([])

        ######################################### Survey the experiments, figure out which experiments have enough successful repeats: appropriate 'n' sample size for subsequent analyses
        list_n = np.array([])
        for i in range(0,len(experiment_list)):

            fileID= str(experiment_list[i])
            fh = fileID+'/'+fileID+'_'
            params = list(np.load(fh+'param.npy'))
            
            for k in range(0,len(params)):# import file
                list20 = params[k]
                han2 = fh+str(list20)+'.npy'
                list10 =  list(np.load(han2))
                # if list10[-1]=='G' and list10[-4]=='nm':
                if list10[-4]=='nm':
                    constant_param = float(list10[-5])
                    list_n = np.append(list_n, constant_param)


        ###########################################
        for i in range(0,len(experiment_list)):

            fileID= str(experiment_list[i])
            fh = fileID+'/'+fileID+'_'
            params = list(np.load(fh+'param.npy'))

            standard_bd = np.array([])
            standard_bd_std = np.array([])

            standard_ibi = np.array([])
            standard_ibi_std = np.array([])

            for k in range(0,len(params)):# import file
                list20 = params[k]
                han2 = fh+str(list20)+'.npy'
                list10 =  list(np.load(han2))

                n = float(list10[-5])
                nn = np.nonzero(n==list_n)
                nn = nn[0]

                if len(nn)>=1:          # minimum successful repeats (minimum n) to accept the data and analyze it.
                # print(len(nn))
                    # if list10[-1]=='G' and list10[-4]=='nm':# and list10[-5]=='2.0':
                    # if list10[-1]=='I':
                    if list10[-4]==protokol:
                        if float(list10[-5])==gP:
                            # if list10[-4]=='nm':
                            # if float(list10[-5])==1.0 or float(list10[-5])==2.0 or float(list10[-5])==5.0 or float(list10[-5])==6.0:

                            print(fileID,list10, countr0, countr1)
                            lab1, title1, constant = get_data.labelfun(list10)

                            plist = list10[0:-5]
                            km=10

                            ######################## V
                            # f16 = plt.figure(figsize=(20,10))
                            # ax16 = f16.add_subplot(121)
                            # ax26 = f16.add_subplot(122)

                            
                            # # generate classic statistics plots of trends. Upregulation of IPM
                            # gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion, Nai_peaks,Nai_troughs,fID, CoV = inference_plots.classic_stats(list10, ax10, ax20,plist, fh, minimum_spikes, colores, countr0)
                            # ax10, ax20, countr0 = inference_plots.classic_trend(IpumpMax_vec, bd_vec, bd_vec_std,ibi_vec, ibi_vec_std, colores, countr0, km, ax10, ax20, plist)
                            
                            # # generate robust statistics plots of trends. Upregulation of IPM
                            gp_vec, IpumpMax_vec, Nai_excursion, bd_median, bd_first_quartile, bd_third_quartile, ibi_median, ibi_first_quartile, ibi_third_quartile, mean_f, median_f, Vm_excursion, cluster_ID ,fID, ax11, ax21, normalization_constant_BD, normalization_constant_IBI= inference_plots.robust_stats(list10, ax11, ax21,plist, fh, minimum_spikes, colores, countr1)
                            ax11, ax21, countr1 = inference_plots.robust_trend(cluster_ID, IpumpMax_vec, bd_median, bd_first_quartile, bd_third_quartile ,ibi_median, ibi_first_quartile, ibi_third_quartile, colores, countr1, km, ax11, ax21, plist)

                            # Nai_excursion = cytNa_max - cytNa_min
                            # print(gp_vec,len(IpumpMax_vec),len(bd_median),len(median_f),len(Vm_excursion),len(ibi_median))
                            for ii in range(0,len(gp_vec)):
                                nrow = {'fileID':fileID, 'gp':gp_vec[ii], 'ipm': IpumpMax_vec[ii], 'Nai':Nai_excursion[ii], 'bd_median':bd_median[ii],'bd_1st_qtl':bd_first_quartile[ii], 'bd_3rd_qtl':bd_third_quartile[ii],'f_median':median_f[ii] , 'f_mean':mean_f,'V_ex': Vm_excursion[ii], 'ibi_median': ibi_median[ii], 'ibi_1st_qtl':ibi_first_quartile[ii], 'ibi_3rd_qtl':ibi_third_quartile[ii], 'normalization_constant_BD':normalization_constant_BD,'normalization_constant_IBI':normalization_constant_IBI, 'cluster':cluster_ID[ii]}
                                tab1 = tab1.append(nrow, ignore_index=True)
    pd.DataFrame.to_excel(tab1 , dir0+'data_table'+str(gP)+'.xlsx')
    # plt.close('all')
############

# #     #                         # # polyfit regression on BD data: experiment-by-experiment

# #     #                         polydeg = 2     # degree of polynomial to fit

# #     #                         ax16, p0 = inference_plots.polyfit_bd_ibi(IpumpMax_vec, bd_median, polydeg, ax16, countr1, colores, cluster_ID, 1)
# #     #                         ax26, p1 = inference_plots.polyfit_bd_ibi(IpumpMax_vec, ibi_median, polydeg, ax26, countr1, colores, cluster_ID, 1)

# #     #                         inference_plots.fancy_trends(ax16, ax26, f16, dir4, gP, protokol, 'robust', fileID, len(plist))

# #     #                         # # generate BD vs IBI characteristics
# #     #                         ax12, ax22, countr2 = inference_plots.bd_ibi_stats(list10, cluster_ID, ax12, ax22, plist, fh, minimum_spikes, colores, countr2)


# #     #                         # # LINEAR regression
# #     #                         r0, ax14 = inference_plots.regress_bd_ibi_g(protokol, gP, bd_median, ibi_median, ax14, f14,colores, countr2, len(plist), str(countr2-1)+'ID_'+str(fileID), 'BD(s)','IBI(s)' , cluster_ID ,dir3, 1)
# #     #                                                                   #(gP, predictor_x, predicted_y, ax, f, colores, kolor, sweep_length, sample_label, xlabl, ylabl, dir3)


# #     #                         # # Select the experiments that showed a significant linear relationship
# #     #                         if r0!=0:
# #     #                             if r0.pvalue <0.01:
# #     #                                 # print(r0.pvalue)
# #     #                                 # print(r0.rvalue)
# #     #                                 p1 = cluster_ID==1
# #     #                                 x0 = np.append(x0, bd_median[p1])
# #     #                                 y0 = np.append(y0, ibi_median[p1])

# #     #                                 a0 = np.append(a0, IpumpMax_vec[p1])

# #     #                                 id0 = np.append(id0, cluster_ID[p1])
# #     #                                 # p0 = inference_plots.polyfit_bd_ibi(bd_median,ibi_median, 2)                            

# #     #                         plt.close('all')
# #     # ######################################################################################################
# #     # ######################################################################################################
    ## classic stats plots
    # inference_plots.fancy_trends(ax10, ax20, f10, dir2, gP, protokol, 'classic', 'Cumulative', len(bd_vec))

    ## robust stats plots
    inference_plots.fancy_trends(ax11, ax21, f11, dir2, gP, protokol, 'robust', 'Cumulative', len(bd_median))

    # plt.close('all')

    # ## characteristiscs plots
    # inference_plots.fancy_characteristics(ax12, ax22, f12, dir2, gP, protokol, 'BD(s)', 'IBI(s)')    

    # ## cumulative linear regression
    # inference_plots.regress_bd_ibi_g(protokol, gP, x0, y0, ax15, f15,colores, -1, len(x0), 'Cumulative', 'BD(s)','IBI(s)' , id0  , dir3, 1)

    # ## cumulative polynomial fit
    # inference_plots.polyfit_bd_ibi(a0, x0, polydeg, ax17, -1, colores, id0, 1)
    # inference_plots.polyfit_bd_ibi(a0, y0, polydeg, ax27, -1, colores, id0, 1)
    # inference_plots.fancy_trends(ax17, ax27, f17, dir4, gP, protokol, 'robust', 'Cumulative', len(a0))

    # plt.close('all')

# gP = 5.0
# tab6 = pd.read_excel(dir0+'data_table5.0.xlsx')

# gP = 6.0
tab6 = pd.read_excel(dir0+'data_table6.0.xlsx')


tab1 = pd.DataFrame()

# h1 = plt.figure(figsize=(20,10))
# ax11 = h1.add_subplot(121)
# ax12 = h1.add_subplot(122)


h2 = plt.figure(figsize=(20,10))
ax20 = h2.add_subplot(121)
ax22 = h2.add_subplot(122)

# ipm_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
ipm_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for i in range(0,len(ipm_list)):
    ipm_ref = ipm_list[i]

    vars_of_interest = tab6.ipm==ipm_ref


    averagd_bd = np.mean(tab6.bd_median[vars_of_interest])
    averagd_bd_1st_qtl = np.mean(tab6.bd_1st_qtl[vars_of_interest])
    averagd_bd_3rd_qtl = np.mean(tab6.bd_3rd_qtl[vars_of_interest])


    averagd_ibi = np.mean(tab6.ibi_median[vars_of_interest])
    averagd_ibi_1st_qtl = np.mean(tab6.ibi_1st_qtl[vars_of_interest])
    averagd_ibi_3rd_qtl = np.mean(tab6.ibi_3rd_qtl[vars_of_interest])


    normalizd_averagd_bd = np.mean(tab6.bd_median[vars_of_interest]/tab6.normalization_constant_BD[vars_of_interest])
    normalizd_averagd_bd_1st_qtl = np.mean(tab6.bd_1st_qtl[vars_of_interest]/tab6.normalization_constant_BD[vars_of_interest])
    normalizd_averagd_bd_3rd_qtl = np.mean(tab6.bd_3rd_qtl[vars_of_interest]/tab6.normalization_constant_BD[vars_of_interest])


    normalizd_averagd_ibi = np.mean(tab6.ibi_median[vars_of_interest]/tab6.normalization_constant_IBI[vars_of_interest])
    normalizd_averagd_ibi_1st_qtl = np.mean(tab6.ibi_1st_qtl[vars_of_interest]/tab6.normalization_constant_IBI[vars_of_interest])
    normalizd_averagd_ibi_3rd_qtl = np.mean(tab6.ibi_3rd_qtl[vars_of_interest]/tab6.normalization_constant_IBI[vars_of_interest])

    # nrow = {'fileID':fileID, 'gp':gp_vec[ii], 'ipm': IpumpMax_vec[ii], 'Nai':Nai_excursion[ii], 'bd_median':bd_median[ii], 'f_median':median_f[ii] , 'f_mean':mean_f,'V_ex': Vm_excursion[ii], 'ibi_median': ibi_median[ii]}
    # tab1 = tab1.append(nrow, ignore_index=True)

    ax11.scatter(ipm_ref, averagd_bd ,c='k',s= 100)
    ax11.scatter(ipm_ref, averagd_bd_1st_qtl, marker='_',c='k',s= 100)
    ax11.scatter(ipm_ref, averagd_bd_3rd_qtl, marker='_',c='k',s= 100)
    ax11.plot([ipm_ref, ipm_ref],[averagd_bd_3rd_qtl, averagd_bd_1st_qtl], linewidth=4, color='k')

    ax21.scatter(ipm_ref, averagd_ibi,c='k',s= 100)
    ax21.scatter(ipm_ref, averagd_ibi_1st_qtl, marker='_',c='k',s= 100)
    ax21.scatter(ipm_ref, averagd_ibi_3rd_qtl, marker='_',c='k',s= 100)
    ax21.plot([ipm_ref, ipm_ref], [averagd_ibi_3rd_qtl, averagd_ibi_1st_qtl], linewidth=4, color='k')

    ax20.scatter(ipm_ref, normalizd_averagd_bd,c='k',s= 70)
    ax20.scatter(ipm_ref, normalizd_averagd_bd_1st_qtl, marker='_',c='k',s= 70)
    ax20.scatter(ipm_ref, normalizd_averagd_bd_3rd_qtl, marker='_',c='k',s= 70)

    ax22.scatter(ipm_ref, normalizd_averagd_ibi,c='k',s= 70)
    ax22.scatter(ipm_ref, normalizd_averagd_ibi_1st_qtl, marker='_',c='k',s= 70)
    ax22.scatter(ipm_ref, normalizd_averagd_ibi_3rd_qtl, marker='_',c='k',s= 70)

    nrow = {'gP': gP, 'ipm':ipm_ref, 'bd_median':averagd_bd, 'bd_1st':averagd_bd_1st_qtl, 'bd_3rd':averagd_bd_3rd_qtl, 'bd_norm_mean':normalizd_averagd_bd, 'bd_norm_1st': normalizd_averagd_bd_1st_qtl,
    'bd_norm_3rd':normalizd_averagd_bd_3rd_qtl, 'ibi_median':averagd_ibi, 'ibi_1st':averagd_ibi_1st_qtl , 'ibi_3rd':averagd_ibi_3rd_qtl , 'ibi_norm_mean':normalizd_averagd_ibi ,'ibi_norm_1st:':normalizd_averagd_ibi_1st_qtl ,'ibi_norm_3rd':normalizd_averagd_ibi_3rd_qtl,
    'cluster':1}
    tab1 = tab1.append(nrow, ignore_index=True)


ax11.plot(tab1.ipm, tab1.bd_median, c='k', marker = 's', linewidth=6)
ax21.plot(tab1.ipm, tab1.ibi_median, c='k', marker = 's', linewidth=6)


# ax11.set_title('pooled data')
ax11.set_ylabel('BD(s)', fontsize=25, fontweight='bold')
ax21.set_ylabel('IBI (s)', fontsize=25, fontweight='bold')

ax11.set_xlabel(r'I$^{pump}_{max}$', fontsize=25, fontweight='bold')
ax21.set_xlabel(r'I$^{pump}_{max}$', fontsize=25, fontweight='bold')

ax20.set_ylabel('BD normalized', fontsize=25, fontweight='bold')
ax22.set_ylabel('IBI normalized', fontsize=25, fontweight='bold')
ax20.set_xlabel(r'I$^{pump}_{max}$', fontsize=25, fontweight='bold')
ax22.set_xlabel(r'I$^{pump}_{max}$', fontsize=25, fontweight='bold')

ax11.tick_params(axis='both', length=10, width=4, direction='inout')
ax21.tick_params(axis='both', length=10, width=4, direction='inout')

ax20.tick_params(axis='both', length=10, width=4, direction='inout')
ax22.tick_params(axis='both', length=10, width=4, direction='inout')


ax11.plot(tab1.ipm, tab1.bd_median, c='k', linewidth=2, alpha=0.5)
ax21.plot(tab1.ipm, tab1.ibi_median, c='k', linewidth=2, alpha=0.5)

ax20.plot(tab1.ipm, tab1.bd_norm_mean, c='k', linewidth=2, alpha=0.5)
ax22.plot(tab1.ipm, tab1.ibi_norm_mean, c='k', linewidth=2, alpha=0.5)

# str_title = r'$\bar{g}_P$='+str(5.0)+' (nS)'
str_title = r'$\bar{g}_P$='+str(6.0)+' (nS)'

ax11.set_ylim(0,4.5)
ax21.set_ylim(0,4.5)

ax11.set_yticks(np.arange(1,5,1))
ax21.set_yticks(np.arange(1,5,1))

f11.suptitle(str_title, fontweight='bold', fontsize=30)
h2.suptitle(str_title, fontweight='bold', fontsize=30)
f11.savefig(dir10+'dependencies1.png')
h2.savefig(dir10+'dependencies2.png')





###################################################################################################
###################################################################################################
protokol = 'nm'
# # polyfit regression on BD data: experiment-by-experiment
# plt.close('all')
h6 = plt.figure(figsize=(10,10))
ax16 = h6.add_subplot(111)
# ax26 = h6.add_subplot(122)

# polydeg = 2     # degree of polynomial to fit

# ax16, p0 = inference_plots.polyfit_bd_ibi(tab6.ipm, tab6.bd_median, polydeg, ax16, -1, colores, tab6.cluster, 1)
# ax26, p1 = inference_plots.polyfit_bd_ibi(tab6.ipm, tab6.ibi_median, polydeg, ax26, -1, colores, tab6.cluster, 1)

# inference_plots.fancy_trends(ax16, ax26, f16, dir4, 6.0, protokol, 'robust', fileID, len(plist))

# # # generate BD vs IBI characteristics
# ax12, ax22, countr2 = inference_plots.bd_ibi_stats(list10, cluster_ID, ax12, ax22, plist, fh, minimum_spikes, colores, countr2)


# # # LINEAR regression
# r0, ax16 = inference_plots.regress_bd_ibi_g(protokol, 5.0, tab1.bd_median, tab1.ibi_median, ax16, h6,colores, -1, 10, 'dataframe_5_r', 'BD(s)','IBI(s)' , tab1.cluster ,dir3, 1)
# r0, ax16 = inference_plots.regress_bd_ibi_g(protokol, 5.0, tab1.bd_norm_mean, tab1.ibi_norm_mean, ax16, h6,colores, -1, 10, 'dataframe_5_n', 'BD(normalized)','IBI(normalized)' , tab1.cluster ,dir3, 1)

r0, ax16 = inference_plots.regress_bd_ibi_g(protokol, 6.0, tab1.bd_median, tab1.ibi_median, ax16, h6,colores, -1, 10, 'dataframe_6_r', 'BD(s)','IBI(s)' , tab1.cluster ,dir3, 1)
# r0, ax16 = inference_plots.regress_bd_ibi_g(protokol, 6.0, tab1.bd_norm_mean, tab1.ibi_norm_mean, ax16, h6,colores, -1, 10, 'dataframe_6_n', 'BD(normalized)','IBI(normalized)' , tab1.cluster ,dir3, 1)

# regress_bd_ibi_g(protocol, gP, predictor_x, predicted_y, ax, f, colores, kolor, sweep_length, sample_label, xlabl, ylabl, cluster_ID,dir3, choice):
# 

pd.DataFrame.to_excel(tab1 , dir0+'tab1'+str(gP)+'.xlsx')


