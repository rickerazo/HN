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
ft2 = 55
# potential
# plt.figure(1,figsize=(45,30))
# # pump
# plt.figure(2,figsize=(45,30))
# # sodium
# plt.figure(3,figsize=(45,30))

cwd = os.getcwd()

fileID= cwd[-8:len(cwd)]
fileID = fileID


######################################	TIME BINS 	########################################################################################################
# select a time range to plot (timebin), all the data can be overwhelming since it's longer than 10 minutes. Carefully discard ~10 cycles to avoid transient membrane potential response to gp

# param0 = ['list1','list2']#,'list3']
# np.save(fileID+'param',param0)
timebin_plot = 37

def compute(list10):
	control = list10[-1]
	list0 = list10[0:-5]
	protocol = list10[-4]
	coupling = list10[-3]
	I_bias = float(list10[-2])
	#####################################################################################
	#####################################################################################
	if control == 'G':
		gp = float(list10[-5])
		# list0 = np.arange(gp_min,gp_max,stepSize) #create list of gp
		list0=list0[::-1] #reverse the order of the list, for plotting purposes

		######################################################################## Plot loop
	# np.save(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),I_mem)
		for j in range(0,np.size(list0)):
			i=list0[j]
			IpumpMax = float(i)
			V = np.load(fileID+'_V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			t = np.load(fileID+'_time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			Na_in = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			I_mem = np.load(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
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
			x2 = x1 - 50
			x3 = x2 - 70
			x4 = x3 - 70
			# x5 = x4 - 20
			
		########################  Current calculations
			# I_bias = np.load(fileID+'_Ibias.npy')
			Ip = I_mem+ Ipump+I_bias
			Ip = -Ip
			# print(j)
			V1 = V/1000 # V is in miliVolts,V1 is in Volts

			Na_o = 0.115*1000 #M is a constant
			ENa = 0.02526*np.log(Na_o/(Na_in))
			# print(gp,len(V1),len(ENa))
			mp = Ip /(gp*(V1-ENa))
			np.save(fileID+'_mp_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling),mp)
			np.save(fileID+'_Ip_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling),Ip)


		##################################### VOLTAGE TRACES
		## plot the data
			# print(i,j)
			plt.figure(i,figsize=(40,30))
			# Ipump
			shift1 = -90
			scale_factor1 = 250
			Ipump_plot = Ipump*scale_factor1+shift1
			# Cyt Na
			shift2 = -150
			scale_factor2 = 5
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
			# plt.plot(time1,x3+mp_plot,color='black')
			plt.plot(time1,x4+Ip_plot,color='black')

			## 50 mV dotted line
			plt.plot([time1[0], time1[-1]],[-50, -50],'--',color='black',linewidth=5)
			## 20 mV bars
			plt.plot([-3.25, -3.25],[-50,-30],color='black',linewidth=20)
			
			## 0.1 nA dotted bar
			plt.plot([time1[0],time1[-1]],[x1+referenceLine*scale_factor1+shift1,x1+ referenceLine*scale_factor1+shift1],'--',color='black',linewidth=5)
			## 0.1 nA bar
			plt.plot([-3.25, -3.25],[x1+referenceLine*scale_factor1+shift1,x1+ (referenceLine+0.1)*scale_factor1+shift1],color='black',linewidth=20)
			## 10 nM dotted bar
			plt.plot([time1[0],time1[-1]],[x2+10*scale_factor2+shift2,x2+10*scale_factor2+shift2],'--',color='black',linewidth=5)
			## 5 nM bar
			plt.plot([-3.25, -3.25],[x2+10*scale_factor2+shift2,x2+15*scale_factor2+shift2],linewidth=20,color='black')
			## Ip bar
			plt.plot([-3.25, -3.25],[x4-0.1*scale_factor4+shift4,x4-0.3*scale_factor4+shift4],linewidth=20,color='black')
			

			# 0 dashed
			# plt.plot([time1[0],time1[-1]],[x3+0*scale_factor3+shift3,x3+0*scale_factor3+shift3],'--',color='black')
			# 1 dashed
			# plt.plot([time1[0],time1[-1]],[x3+1*scale_factor3+shift3,x3+1*scale_factor3+shift3],'--',color='black')
			# -0.1 nA Ip dashed line
			plt.plot([time1[0],time1[-1]],[x4-0.1*scale_factor4+shift4,x4-0.1*scale_factor4+shift4],'--',color='black',linewidth=5)

			f1 = 100
			##data label
			plt.text(-12, 50,r'I$_{pump}^{max}$='+str(IpumpMax)+' nA',fontsize=100,weight='bold')
			plt.text(-12,15,r'$\bar{g}_P$='+str(gp)+' nS',fontsize=100,weight='bold')
			plt.text(-10,-35,r'$V_m$',fontsize=f1,weight='bold')
			plt.text(-12,(0.08+referenceLine)*scale_factor1+shift1 + x1,r'I$_{pump}$',fontsize=f1,weight='bold')
			plt.text(-11,x2+13*scale_factor2+shift2,r'[Na]$_i$',fontsize=f1,weight='bold')
			# plt.text(-7,x3+0.5*scale_factor3+shift3,r'm$_P$',fontsize=80,weight='bold')
			plt.text(-10,x4-0.1*scale_factor4+shift4,r'I$_P$',fontsize=f1,weight='bold')
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
			plt.plot([10,20],[15,15],linewidth=20,color='black')
			plt.text(12,17,'10 seconds')
			# labels for bars
			plt.text(-4.5,-50,'20 mV',fontsize=ft2,rotation=90)
			x1+referenceLine*scale_factor1+shift1,x1
			plt.text(-4.5,(referenceLine)*scale_factor1+shift1 + x1,'0.1 nA', fontsize=ft2,rotation=90)
			plt.text(-4.5,x2+10*scale_factor2+shift2,'5 nM', fontsize=ft2,rotation=90)
			plt.text(-4.5,x4-0.3*scale_factor4+shift4,'0.2 nA', fontsize=ft2,rotation=90)

			if timebin_plot>time1[-1]:
				x_label = time1[-1]+0.5
			else:
				x_label = timebin_plot+0.5
			# labels for dashed lines
			plt.text(x_label,-50,'-50 mV', fontsize=ft2)
			plt.text(x_label,-4+referenceLine*scale_factor1+shift1 +x1,referenceStr, fontsize=ft2)
			plt.text(x_label,9.25*scale_factor2+shift2 +x2,'10 nM',fontsize=ft2)
			plt.text(x_label,-0.13*scale_factor4+shift4+x4,'-0.1 nA',fontsize=ft2)			# plt.text(-1,x3+0.9*scale_factor3+shift3,'1',fontsize=ft2)
			# plt.text(-1,x3-0.1*scale_factor3+shift3,'0',fontsize=ft2)

			plt.text(x_label, 50, protocol, fontsize=100 )
			#plt.title(r'Membrane potential $HN_7$ $I_{pump}^{max}$='+str(IpumpMax)+' nA')
			# plt.ylabel(r'$V_m$(mV)')
			#plt.xlabel('Time(s)')
			plt.axis([-6,timebin_plot,-0.7*scale_factor4+x4+shift4,15])
			plt.savefig(fileID+'_traces_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.png')	


	if control == 'I':
		IpumpMax = list10[-5]
		# list0 = np.arange(gp_min,gp_max,stepSize) #create list of gp
		list0=list0[::-1] #reverse the order of the list, for plotting purposes

		######################################################################## Plot loop
	# np.save(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp),I_mem)
		for j in range(0,np.size(list0)):
			i=list0[j]
			gp = float(i)
			V = np.load(fileID+'_V_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			t = np.load(fileID+'_time_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			Ipump = np.load(fileID+'_pump_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			Na_in = np.load(fileID+'_cytNa_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
			I_mem = np.load(fileID+'_Itot_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.npy')
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
			x2 = x1 - 50
			x3 = x2 - 70
			x4 = x3 - 70
			# x5 = x4 - 20

		########################  Current calculations
			# I_bias = np.load(fileID+'_Ibias.npy')
			Ip = I_mem+ Ipump+I_bias
			Ip = -Ip
			# print(j)
			V1 = V/1000 # V is in miliVolts,V1 is in Volts

			Na_o = 0.115*1000 #M is a constant
			ENa = 0.02526*np.log(Na_o/(Na_in))
			mp = Ip /(gp*(V1-ENa))
			np.save(fileID+'_mp_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling),mp)
			np.save(fileID+'_Ip_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling),Ip)


		##################################### VOLTAGE TRACES
		## plot the data
			# print(i,j)
			plt.figure(i,figsize=(40,30))
			# Ipump
			shift1 = -90
			scale_factor1 = 250
			Ipump_plot = Ipump*scale_factor1+shift1
			# Cyt Na
			shift2 = -150
			scale_factor2 = 5
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
			# plt.plot(time1,x3+mp_plot,color='black')
			plt.plot(time1,x4+Ip_plot,color='black')

			## 50 mV dotted line
			plt.plot([time1[0], time1[-1]],[-50, -50],'--',color='black',linewidth=5)
			## 20 mV bars
			plt.plot([-3.25, -3.25],[-50,-30],color='black',linewidth=20)
			
			## 0.1 nA dotted bar
			plt.plot([time1[0],time1[-1]],[x1+referenceLine*scale_factor1+shift1,x1+ referenceLine*scale_factor1+shift1],'--',color='black',linewidth=5)
			## 0.1 nA bar
			plt.plot([-3.25, -3.25],[x1+referenceLine*scale_factor1+shift1,x1+ (referenceLine+0.1)*scale_factor1+shift1],color='black',linewidth=20)
			## 10 nM dotted bar
			plt.plot([time1[0],time1[-1]],[x2+10*scale_factor2+shift2,x2+10*scale_factor2+shift2],'--',color='black',linewidth=5)
			## 5 nM bar
			plt.plot([-3.25, -3.25],[x2+10*scale_factor2+shift2,x2+15*scale_factor2+shift2],linewidth=20,color='black')
			## Ip bar
			plt.plot([-3.25, -3.25],[x4-0.1*scale_factor4+shift4,x4-0.3*scale_factor4+shift4],linewidth=20,color='black')
			

			# 0 dashed
			# plt.plot([time1[0],time1[-1]],[x3+0*scale_factor3+shift3,x3+0*scale_factor3+shift3],'--',color='black')
			# 1 dashed
			# plt.plot([time1[0],time1[-1]],[x3+1*scale_factor3+shift3,x3+1*scale_factor3+shift3],'--',color='black')
			# -0.1 nA Ip dashed line
			plt.plot([time1[0],time1[-1]],[x4-0.1*scale_factor4+shift4,x4-0.1*scale_factor4+shift4],'--',color='black',linewidth=5)

			f1 = 100
			##data label
			plt.text(-12, 50,r'I$_{pump}^{max}$='+str(IpumpMax)+' nA',fontsize=100,weight='bold')
			plt.text(-12,15,r'$\bar{g}_P$='+str(gp)+' nS',fontsize=100,weight='bold')
			plt.text(-10,-35,r'$V_m$',fontsize=f1,weight='bold')
			plt.text(-12,(0.08+referenceLine)*scale_factor1+shift1 + x1,r'I$_{pump}$',fontsize=f1,weight='bold')
			plt.text(-11,x2+13*scale_factor2+shift2,r'[Na]$_i$',fontsize=f1,weight='bold')
			# plt.text(-7,x3+0.5*scale_factor3+shift3,r'm$_P$',fontsize=80,weight='bold')
			plt.text(-10,x4-0.1*scale_factor4+shift4,r'I$_P$',fontsize=f1,weight='bold')
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
			plt.plot([10,20],[15,15],linewidth=20,color='black')
			plt.text(12,17,'10 seconds')
			# labels for bars
			plt.text(-4.5,-50,'20 mV',fontsize=ft2,rotation=90)
			plt.text(-4.5,(referenceLine)*scale_factor1+shift1 + x1,'0.1 nA', fontsize=ft2,rotation=90)
			plt.text(-4.5,x2+10*scale_factor2+shift2,'5 nM', fontsize=ft2,rotation=90)
			plt.text(-4.5,x4-0.3*scale_factor4+shift4,'0.2 nA', fontsize=ft2,rotation=90)

			if timebin_plot>time1[-1]:
				x_label = time1[-1]+0.5
			else:
				x_label = timebin_plot+0.5
			# labels for dashed lines
			plt.text(x_label,-50,'-50 mV', fontsize=ft2)
			plt.text(x_label,-4+referenceLine*scale_factor1+shift1 +x1,referenceStr, fontsize=ft2)
			plt.text(x_label,9.25*scale_factor2+shift2 +x2,'10 nM',fontsize=ft2)
			plt.text(x_label,-0.13*scale_factor4+shift4+x4,'-0.1 nA',fontsize=ft2)			# plt.text(-1,x3+0.9*scale_factor3+shift3,'1',fontsize=ft2)

			plt.text(x_label, 50, protocol, fontsize=100 )
			#plt.title(r'Membrane potential $HN_7$ $I_{pump}^{max}$='+str(IpumpMax)+' nA')
			# plt.ylabel(r'$V_m$(mV)')
			#plt.xlabel('Time(s)')
			plt.axis([-6,timebin_plot,-0.7*scale_factor4+x4+shift4,15])
			plt.savefig(fileID+'_traces_IpumpMax='+str(IpumpMax)+'_gp='+str(gp)+'_protocol='+str(protocol)+'_coupling='+str(coupling)+'.png')



params = list(np.load(fileID+'_param.npy'))
for k in range(0,len(params)):# import file
	list20 = params[k]
	han2 = fileID+'_'+str(list20)+'.npy'
	list10 =  list(np.load(han2))
	print(fileID,list10)

	compute(list10)
	plt.close('all')

# for k in range(0,len(params)):# import file
# k = 2
# list20 = params[k]
# han2 = fileID+'_'+str(list20)+'.npy'
# list10 =  list(np.load(han2))
# print(fileID,list10)
# # I_mem, Na_in, Ip = compute(list10)
# compute(list10)
# plt.close('all')

# plt.ion()
# plt.show()

# Important note
# Ipump has model notation
# Itotal has notation the way its injected