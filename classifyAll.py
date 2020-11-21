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
from scipy.signal import savgol_filter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import OrderedDict
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import linregress
from scipy.stats.stats import pearsonr
from mpl_toolkits import mplot3d
from scipy.stats import ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests 
import pandas as pd

import get_data
import inference_plots
import custom_analysis

import pandas as pd
import os
import sys

cmaps = OrderedDict()

##############################################################################
####	graphix stuff
mpl.rcParams['axes.linewidth']=10


markers = ['1', '2', 'x0', '1','4','3','p','+','h','|','.', '2', 'x0', '1','4','3','p','+','h','|']
colores = ['lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k','lightcoral','firebrick','red','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','darkorchid','fuchsia', 'k']
colores_binario = ['red','blue']
# colores = ['purple','lawngreen','gold','steelblue','darkorange','darkgoldenrod','gold','lawngreen','forestgreen','turquoise','steelblue','navy','firebrick','fuchsia']


font = {'weight' : 'bold',
        'size'   : 30}
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
    '20o14001'
	])
np.save('experiment_list',experiment_list)


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

############## IMPORTANT : minimum 3 successful experiment repeats to consider the data for subsequent analysis
os.makedirs('data_analyses/classifyAll', exist_ok=True)

dir0 = 'data_analyses/classifyAll/'
compute_data = 1

if compute_data==True:
    countr0= 1
    countr1= 1
    countr2= 1

    # {'gp':gp_vec, 'ipm': IpumpMax_vec, 'Nai':Nai_excursion, 'bd_median':bd_median, 'f_median':median_f , 'V_ex': Vm_excursion}
    tab1 = pd.DataFrame(columns = ['fileID', 'gp', 'ipm', 'Nai_ex', 'bd_median' , 'f_median' , 'f_mean', 'V_ex', 'ibi_median'])

    protokol = 'nm'
    # protokol = 'cs'

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
                if list10[-4]==protokol:

                    print(fileID,list10, countr0, countr1)
                    lab1, title1, constant = get_data.labelfun(list10)

                    plist = list10[0:-5]
                    km=10

                    # # generate classic statistics plots of trends. Upregulation of IPM
                    # gp_vec,IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_vec, bd_vec_std, ibi_vec, ibi_vec_std, hz_vec, period_vec, burst_vm_excursion, cytNa_excursion, Nai_peaks, Nai_troughs,fID, CoV = inference_plots.classic_stats(list10, ax10, ax20,plist, fh, minimum_spikes, colores, countr0)
                    # ax10, ax20, countr0 = inference_plots.classic_trend(IpumpMax_vec, bd_vec, bd_vec_std,ibi_vec, ibi_vec_std, colores, countr0, km, ax10, ax20, plist)
                    
                    # # generate robust statistics plots of trends. Upregulation of IPM
                    # gp_vec, IpumpMax_vec, cytNa_mean, cytNa_std, cytNa_max, cytNa_min, bd_median, bd_first_quartile, bd_third_quartile, ibi_median, ibi_first_quartile, ibi_third_quartile, mean_f, median_f, Vm_excursion, fID, ax11, ax21 = inference_plots.robust_stats(list10, ax11, ax21,plist, fh, minimum_spikes, colores, countr1)
                    gp_vec, IpumpMax_vec, Nai_excursion, bd_median, bd_first_quartile, bd_third_quartile, ibi_median, ibi_first_quartile, ibi_third_quartile, mean_f, median_f, Vm_excursion_vec, cluster_vec ,fID, ax11, ax21, normalization_constant_BD, normalization_constant_IBI=inference_plots.robust_stats(list10, ax11, ax21,plist, fh, minimum_spikes, colores, countr1)
                    # ax11, ax21, countr1 = inference_plots.robust_trend(cluster_ID, IpumpMax_vec, bd_median, bd_first_quartile, bd_third_quartile ,ibi_median, ibi_first_quartile, ibi_third_quartile, colores, countr1, km, ax11, ax21, plist)

                    ## The new row of data that will be added to the dataFrame with all pooled data
                    for ii in range(0,len(gp_vec)):
                        nrow = {'fileID':fileID, 'gp':gp_vec[ii], 'ipm': IpumpMax_vec[ii], 'Nai_ex':Nai_excursion[ii], 'bd_median':bd_median[ii], 'f_median':median_f[ii] , 'f_mean':mean_f[ii],'V_ex': Vm_excursion_vec[ii], 'ibi_median': ibi_median[ii]}
                        tab1 = tab1.append(nrow, ignore_index=True)

    ######################################################################################################
    ## classic stats plots
    # inference_plots.fancy_trends(ax10, ax20, f10, dir0, 10, protokol, 'classic', 'Cumulative', len(bd_vec))

    ## robust stats plots
    # inference_plots.fancy_trends(ax11, ax21, f11, dir0, 10, protokol, 'robust', 'Cumulative', len(bd_median))


    # # plt.close('all')


    # plt.ion()
    # plt.show()

    plt.close('all')

    plt.hist(tab1.Nai_ex, label='Nai', alpha=0.5)
    plt.hist(tab1.V_ex, label='Vm', alpha=0.5)
    plt.hist(tab1.bd_median, label='bd', alpha=0.5)
    plt.hist(tab1.ibi_median, label='ibi', alpha=0.5)
    plt.hist(tab1.f_median, label='f', alpha=0.5)

    p1 = tab1.f_median>90
    print(tab1[p1])


    # plt.ion()
    # plt.legend()
    # plt.show()

    plt.savefig(dir0+'histograms.png')


    pd.DataFrame.to_excel(tab1 , dir0+'data_table01.xlsx')


tab1 = pd.read_excel(dir0+'data_table01.xlsx')

# clulab = custom_analysis.clustering_KMeans(tab1.V_ex, tab1.f_median, tab1.bd_median, 2)
# IpumpMax, gp, x_all, time,V, cytNa, interspike_tolerance, burst_duration, interburst_interval, cycle_period, interspike_tolerance = get_data.experimental_variables(control, protocol, coupling, list10, list1, 0, '18914000/18914000_')

print(tab1)


