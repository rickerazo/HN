####################################################################
# Ricardo Erazo
# Neuroscience Institute
# Georgia State University
# Emory University
# All rights reserved
####################################################################

# Import numpy output from traces.py
import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import find_peaks
import os
import matplotlib
## graphs stuff
matplotlib.rcParams['axes.linewidth']=0.0
font = {'weight' : 'bold',
        'size'   : 60}
plt.rc('font', **font)
ft2 = 45
# potential
# plt.figure(1,figsize=(45,30))
# # pump
# plt.figure(2,figsize=(45,30))
# # sodium
# plt.figure(3,figsize=(45,30))

cwd = os.getcwd()

fileID= cwd[-8:len(cwd)]



######################################	TIME BINS 	########################################################################################################
# select a time range to plot (timebin), all the data can be overwhelming since it's longer than 10 minutes. Carefully discard ~10 cycles to avoid transient membrane potential response to gp
list1 = list(np.load(fileID+'_list1.npy'))
list2 = list(np.load(fileID+'_list2.npy'))
list3 = list(np.load(fileID+'_list3.npy'))

param0 = ['list1','list2']#,'list3']
np.save(fileID+'_param',param0)
timebin_plot = 37

def compute(list10):
	control = list10[-1]
	list0 = list10[0:-2]
	if control == 'I':
		IpumpMax = list10[-2]
		# list0 = np.arange(gp_min,gp_max,stepSize) #create list of gp
		list0=list0[::-1] #reverse the order of the list, for plotting purposes

		######################################################################## Plot loop
	# np.save(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),I_mem)
		for j in range(0,np.size(list0)):
			i=list0[j]
			gp = float(i)
			V = np.load(fileID+'_V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			t = np.load(fileID+'_time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			Na_in = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			I_mem = np.load(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy')
			# fh1 = fileID+'_timebin_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.npy'
			# timebin_all = np.load(fh1)
			t_ini = t[0]
			t_end = t[0]+timebin_plot
		######################################## TIMEBINS 
			time1 = t - t[0]	
			# time2 = t[np.nonzero(t<t_end)]
			# time2 = time1
			# pump = Ipump[np.nonzero(t<t_end)]
			# pump = pump/IpumpMax
			# cyt_na = Na_in[np.nonzero(t<t_end)]
			

			# Ip = I_mem- Ipump-0.1
			# # Ip = Gp*mp*(V-Ep) -> mp = Ip/(Gp*(V-ENa))
			# # compute reversial potential first
			# Na_o = 0.115 #M is a constant
			# # Na_in in nM must be converted to SI: (nano = 1e-9)
			# ENa = 0.02526*np.log(Na_o/(Na_in))
			# mp = Ip /(gp*(V-ENa))
			# Ip = I_total[np.nonzero(t<t_end)]

			x1 = -70
			
			# print(np.max(Ipump),np.min(Ipump))
			# print(np.mean(Ipump))
			if np.mean(Ipump)>0.1:
				referenceLine = 0.1
				referenceStr = ' 0.1 nA'
				x1 = -70
				# print('more')
			else:
				referenceLine = 0.0
				referenceStr = ' 0.0 nA'
				x1 = -35
				# print('less')
			x2 = x1-45
			x3 = x2 - 40
			x4 = x3 - 70
			# x5 = x4 - 20

		########################  Current calculations

			Ip = I_mem+ Ipump+0.1
			Ip = -Ip
			# print(j)
			V1 = V/1000 # V is in miliVolts,V1 is in Volts

			Na_o = 0.115*1000 #M is a constant
			ENa = 0.02526*np.log(Na_o/(Na_in))
			mp = Ip /(gp*(V1-ENa))
			np.save(fileID+'_mp_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),mp)
			np.save(fileID+'_Ip_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),Ip)


		##################################### VOLTAGE TRACES
		## plot the data
			# print(i,j)
			plt.figure(i,figsize=(50,65))
			# Ipump
			shift1 = -60
			scale_factor1 = 150
			Ipump_plot = Ipump*scale_factor1+shift1
			# Cyt Na
			shift2 = -120
			scale_factor2 = 6
			CytNa_plot = Na_in*scale_factor2+shift2
			## mp
			shift3 = -60
			scale_factor3 = 30
			mp_plot = mp*scale_factor3+shift3
			## Ip
			scale_factor4 = 200
			shift4 = 0

			Ip_plot = Ip*scale_factor4+shift4
			lw = 4
			## traces
			plt.plot(time1,V,color='black',linewidth=lw)
			plt.plot(time1,x1+Ipump_plot,color='black',linewidth=lw)
			plt.plot(time1,x2+CytNa_plot,color='black')
			plt.plot(time1,x3+mp_plot,color='black')
			plt.plot(time1,x4+Ip_plot,color='black')

			## 50 mV dotted line
			plt.plot([time1[0], time1[-1]],[-50, -50],'--',color='black')
			## 20 mV bars
			plt.plot([-3.5, -3.5],[-50,-30],color='black',linewidth=10)
			
			## 0.1 nA dotted bar
			plt.plot([time1[0],time1[-1]],[x1+referenceLine*scale_factor1+shift1,x1+ referenceLine*scale_factor1+shift1],'--',color='black')
			## 0.1 nA bar
			plt.plot([-3.5,-3.5],[x1+referenceLine*scale_factor1+shift1,x1+ (referenceLine+0.1)*scale_factor1+shift1],color='black',linewidth=8)
			## 10 nM dotted bar
			plt.plot([time1[0],time1[-1]],[x2+10*scale_factor2+shift2,x2+10*scale_factor2+shift2],'--',color='black')
			## 5 nM bar
			plt.plot([-3.5,-3.5],[x2+10*scale_factor2+shift2,x2+15*scale_factor2+shift2],linewidth=10,color='black')
			## Ip bar
			plt.plot([-3.5,-3.5],[x4-0.1*scale_factor4+shift4,x4-0.3*scale_factor4+shift4],linewidth=10,color='black')
			

			# 0 dashed
			plt.plot([time1[0],time1[-1]],[x3+0*scale_factor3+shift3,x3+0*scale_factor3+shift3],'--',color='black')
			# 1 dashed
			plt.plot([time1[0],time1[-1]],[x3+1*scale_factor3+shift3,x3+1*scale_factor3+shift3],'--',color='black')
			# -0.1 nA Ip dashed line
			plt.plot([time1[0],time1[-1]],[x4-0.1*scale_factor4+shift4,x4-0.1*scale_factor4+shift4],'--',color='black')

			##data label
			plt.text(-9,-10,r'I$_{pump}^{max}$='+str(IpumpMax)+' nA',fontsize=90,weight='bold')
			plt.text(-9,-10,r'G$_p$='+str(gp)+' nS',fontsize=85,weight='bold')
			plt.text(-7,-35,r'V',fontsize=80,weight='bold')
			plt.text(-9,(0.08+referenceLine)*scale_factor1+shift1 + x1,r'I$_{pump}$',fontsize=80,weight='bold')
			plt.text(-8,x2+13*scale_factor2+shift2,r'[Na]$_{in}$',fontsize=80,weight='bold')
			plt.text(-7,x3+0.5*scale_factor3+shift3,r'm$_P$',fontsize=80,weight='bold')
			plt.text(-7,x4-0.1*scale_factor4+shift4,r'I$_p$',fontsize=80,weight='bold')
		###################################################################### Graph config
		####################################### VOLTAGE TRACES
			plt.figure(i)
			time1 = t - t[0]
			# time = time + timebin*(i-gp_min)

			plt.tick_params(
			    axis='both',          # changes apply to the x-axis
			    which='both',      # both major and minor ticks are affected
			    bottom=False,      # ticks along the bottom edge are off
			    top=False,         # ticks along the top edge are off
			    left=False,
			    right=False,
			    labelleft=False,
			    labelbottom=False) # labels along the bottom edge are off

			# ten seconds: [0,10]
			# 8*scale_factor2+shift2 +x2
			plt.plot([25,35],[-0.15*scale_factor3+shift3 +x3,-0.15*scale_factor3+shift3 +x3],linewidth=10,color='black')
			plt.text(27,-0.35*scale_factor3+shift3 +x3,'10 seconds')
			# labels for bars
			plt.text(-4.5,-36,'20 mV',fontsize=ft2,rotation=90)
			plt.text(-4.5,(0.1+referenceLine)*scale_factor1+shift1 + x1,'0.1 nA', fontsize=ft2,rotation=90)
			plt.text(-4.5,13*scale_factor2+shift2 + x2,'5 nM', fontsize=ft2,rotation=90)
			plt.text(-4.5,x4-0.17*scale_factor4+shift4,'0.2 nA', fontsize=ft2,rotation=90)

			# labels for dashed lines
			plt.text(-2.5,-55,'-50 mV', fontsize=ft2)
			plt.text(-2,-4+referenceLine*scale_factor1+shift1 +x1,referenceStr, fontsize=ft2)
			plt.text(-2,9.25*scale_factor2+shift2 +x2,'10 nM',fontsize=ft2)
			plt.text(-2,-0.13*scale_factor4+shift4+x4,'-0.1 nA',fontsize=ft2)
			plt.text(-1,x3+0.9*scale_factor3+shift3,'1',fontsize=ft2)
			plt.text(-1,x3-0.1*scale_factor3+shift3,'0',fontsize=ft2)

			#plt.title(r'Membrane potential $HN_7$ $I_{pump}^{max}$='+str(IpumpMax)+' nA')
			# plt.ylabel(r'$V_m$(mV)')
			#plt.xlabel('Time(s)')
			plt.axis([-6,timebin_plot,-0.7*scale_factor4+x4+shift4,15])
			plt.savefig(fileID+'_traces_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'.png')
	return I_mem, Na_in, Ip

# plt.close('all')
I_mem, Na_in, Ip = compute(list1)

# Ip = I_mem+ Ipump+0.1
# Ip = -Ip
# take it wth sign minus and 

# Na_o = 0.115 #M is a constant
# ENa = 0.02526*np.log(Na_o/(Na_in))
# mp = Ip /(gp*(mV-ENa))
# plt.plot(mp,label=r'Na$_o$='+str(Na_o))
# print(ENa)

# plt.figure(12)
# plt.plot(Ip,label=r'Na$_o$='+str(Na_o))

#V1 is variable in Volts
# gp in nS
#Ip in nA

# plt.figure(11)
# Na_o = 0.115*1000 #M is a constant
# ENa = 0.02526*np.log(Na_o/(Na_in))
# mp = Ip /(gp*(V1-ENa))
# # print(ENa)
# plt.plot(mp,label=r'Na$_o$='+str(Na_o))

# plt.title(r'm$_p$')
# plt.legend()

# plt.figure(12)
# plt.plot(Ip,label=r'Na$_o$='+str(Na_o))
# plt.plot([0,len(Ip)],[0,0])
# plt.title(r'I$_{p}$')
# plt.ylabel('nA')
# plt.legend()

plt.ion()
plt.show()

# Important note
# Ipump has model notation
# Itotal has notation the way its injected