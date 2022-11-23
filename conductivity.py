import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd

#file='test/40_out.xyz'
file='../test/lys/lys_100.xyz'


#time step in picoseconds
dt=0.1

class End_of_Loop(Exception): pass #To being able of closing the loop outside the loop functions

def read_line(f):  
	#Is needed to check when the file finish and the exception need to being handled
	line=f.readline()
	if not line:
		raise End_of_Loop
	else:
		return line

def process_line(line):
	#simple rasing function
	line=[float(l) for l in line.strip('\n').split()]
	#print(line)
	return line[0],line[1:4]


def process_frame(f):
	#take file object and return position of all atoms in 1 single frame

	n_mol=read_line(f)
	n_mol=int(n_mol) 		#FIRST LINE is number of molecules in the frame
	f.readline().strip('\n')  		#SECOND LINE is the blanck line of xyz frame
	#FROM THIRD LINE
	i=0
	#initialize variables
	q_frame, x_frame, y_frame,z_frame= [],[],[],[]
	#cicle over molecules
	for i in range(n_mol):
		charge , coord = process_line(read_line(f))	
		q_frame.append(charge)
		x_frame.append(coord[0])
		y_frame.append(coord[1])
		z_frame.append(coord[2])
	return q_frame,x_frame,y_frame,z_frame



#-------------------------------------------------
#DATA
# travis_data=pd.read_csv("../test/lys/lys_c_travis.csv",sep='; ',header=0,engine='python')
# travis_data.columns = ['t', 'msd', 'derivative']
# travis_data_a=pd.read_csv("../test/lys/lys_a_travis.csv",sep='; ',header=0,engine='python')
# travis_data_a.columns = ['t', 'msd', 'derivative']
# vmd_data_c=pd.read_csv("../test/lys/vmd_c.dat",sep=' ',header=None)
# vmd_data_c.columns = ['t', 'msd']
# vmd_data_a=pd.read_csv("../test/lys/vmd_a.dat",sep=' ',header=None)
# vmd_data_a.columns = ['t', 'msd']
# vmd_data_a.t=vmd_data_a.t*1e3
# vmd_data_a.msd=vmd_data_a.msd*1e4
# vmd_data_c.t=vmd_data_c.t*1e3
# vmd_data_c.msd=vmd_data_c.msd*1e4
#DATA 
#----------------------------




#TRAJECTORY LOADING---------------------------------
try:
	with open(file,'r') as f:
		q,x,y,z=[],[],[],[]
		print("LOADING TRAJECTORY FILE")
		while True:
			q_new,x_new,y_new,z_new=process_frame(f)
			if q!=q_new and q!=[]:
				#check that the order of the atom didn't changed from frame to frame
				raise Exception("ERROR: the order of molecules has changed in the dynamics. Check your file please")
			q=q_new
			x.append(x_new)
			y.append(y_new)
			z.append(z_new)
except End_of_Loop:
	print("End of File")
	pass
x=np.array(x)
y=np.array(y)
z=np.array(z)
q=np.array(q)
#select the indexes of anion and cations
#is a tuple of array, one for axis
cation_idx=np.where(q == 1)
anion_idx=np.where(q == -1)


#PRINTING INFO ON TRAJECTORY
print("LOADING TRAJECTORY FILE : COMPLETE")
if x.shape == y.shape == z.shape and q.shape[0] == x.shape[1] == y.shape[1] == z.shape[1]:
	print("Total Number of Frame : ",x.shape[0])
	print("Total length of trajectory : ",x.shape[0]*dt, "ps")
	print("Total number of Molecules : ",x.shape[1])
else:
	raise Exception("SOMETHING WRONG WITH THE TRAJECTORY")
#---------------------------
#END OF TRAJECTORY LOADING
#---------------------------

def regression(msd,t,scaling=0.3):
	#
	#Use skitlearn tools for making the linear regression
	#
	idx=int(len(t)*scaling)
	t_pred=t[idx:].reshape((-1,1))
	msd_subset=msd[idx:]
	model_c = LinearRegression().fit(t_pred,msd_subset)
	print("LINEAR REGRESSION ON THE LAST ",100-scaling*100," %")
	print(f"slope : {model_c.coef_} ")
	print(f"Intercept: {model_c.intercept_} ")
	print(f"D : {model_c.coef_/6} pm^2/ps")
	msd_pred = model_c.predict(t_pred)
	return msd_pred, t_pred

def get_t(msd,dt):
	#create the proper time vector
	return np.arange(0,msd.shape[0]*dt,dt)


def plotting(msd_list):
	#
	# Simple plotting of raw MSD and linear regression line
	# - require a list of msd [msd_cation,msd_anion]
	#
	fig, ax = plt.subplots(1,2)
	for i,msd in enumerate(msd_list):
		t=get_t(msd,dt)
		msd_pred,t_pred = regression(msd,t)
		ax[i].plot(t,msd,linewidth=1.3,label=r'msd',color='red',zorder=2)
		ax[i].plot(t_pred,msd_pred,label='linear regression',linewidth=1,linestyle='dashed',color='blue',zorder=3)
		ax[i].legend()
		ax[i].set_xlabel(r'time / ps')
		ax[i].set_ylabel(r'MSD / pm^2')
		ax[i].set_box_aspect(1)
	fig.suptitle("Cation and Anion MSD")
	plt.show()

def self_product(dx,dy,dz):
	#Product for selfdiffusion i*i
	r2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
	#mean over N molecules
	r2=np.sum(r2,axis=1)/r2.shape[1]
	#UNIT CHANGES FROM A^2/ps -> pm^2/ps
	unit_conversion=1e4
	r2=r2*unit_conversion
	return r2


def get_self_msd(x,y,z,depth=0.3):
	#CORRELATION DEPTH
	#is the percentage of the trajectory in which the correlation between ions is took in account
	
	#loop over subsets
	print("Correlation depth : ",depth*100," %")
	subset_idx = np.arange(0,int(depth*x.shape[0]))
	max_origin_index = x.shape[0]-len(subset_idx)
	print("Number of intervals : ",max_origin_index)
	msd = np.zeros(len(subset_idx))
	i=0
	#LOOP OVER INTERVALS
	while subset_idx[0]<max_origin_index :
		print(i)
		if i%100 == 0:
			print("Intervals processed : ",i)
		x_tmp=x[subset_idx,:]
		y_tmp=y[subset_idx,:]
		z_tmp=z[subset_idx,:]

		#Deviations respect to the reference t0 of the interval
		dx=x_tmp[:,:]-x_tmp[0,:]
		dy=y_tmp[:,:]-y_tmp[0,:]
		dz=z_tmp[:,:]-z_tmp[0,:]

		r2=self_product(dx,dy,dz)
		msd=msd+r2
		i+=1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/i
	return msd





# print("-----------------------")
# print("Cation Self-diffusion   ")
# print("-----------------------")
#msd1=get_self_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.99)
# print("-----------------------")
# print("Anion Self-diffusion   ")
# print("-----------------------")
# msd2=get_self_msd(x[:,anion_idx[0]],y[:,anion_idx[0]],z[:,anion_idx[0]],depth=0.99)
# plotting([msd1,msd2])

def comb(array):
	#simple function i need to get all the combination of ij inner products without repetitions
	return [np.multiply(array[:,0],i) for i in array[:,1:].T]

import time

def inter_same_product(dx,dy,dz):
		start = time.time()
		#f=lambda a:[np.multiply(a[:,0],i) for i in a[:,1:].T]
		dx2,dy2,dz2=[],[],[]
		for a,b,c in zip(dx.T,dy.T,dz.T):
			dx2=dx2+comb(dx)
			dy2=dy2+comb(dy)
			dz2=dz2+comb(dz)
			dx=np.delete(dx,0,axis=1)
			dy=np.delete(dy,0,axis=1)
			dz=np.delete(dz,0,axis=1)
		#Then re-shape the product in the correct shape
		dx2=np.array(dx2).T
		dy2=np.array(dy2).T
		dz2=np.array(dz2).T
		r2=dx2+dy2+dz2
		r2=np.sum(r2,axis=1)/r2.shape[1]
		unit_conversion=1e4
		r2=r2*unit_conversion
		end = time.time()
		print("TIME : ",end - start)
		return r2

def get_inter_msd(x,y,z,depth=0.3):
	#CORRELATION DEPTH
	#is the percentage of the trajectory in which the correlation between ions is took in account
	
	#loop over subsets
	print("Correlation depth : ",depth*100," %")
	subset_idx = np.arange(0,int(depth*x.shape[0]))
	max_origin_index = x.shape[0]-len(subset_idx)
	print("Number of intervals : ",max_origin_index)
	msd = np.zeros(len(subset_idx))
	count=0
	#LOOP OVER INTERVALS
	while subset_idx[0]<max_origin_index :
		print(count)	
		if count%100 == 0:
			print("Intervals processed : ",count)
		x_tmp=x[subset_idx,:]
		y_tmp=y[subset_idx,:]
		z_tmp=z[subset_idx,:]

		#Deviations respect to the reference t0 of the interval
		dx=x_tmp[:,:]-x_tmp[0,:]
		dy=y_tmp[:,:]-y_tmp[0,:]
		dz=z_tmp[:,:]-z_tmp[0,:]
		#SPECIAL PHASE IF i!=j
		r2=inter_same_product(dx,dy,dz)
		#END SPECIAL PHASE

		msd=msd+r2
		count=count+1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/count
	return msd

msd = get_inter_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.99)
#plotting([msd])





