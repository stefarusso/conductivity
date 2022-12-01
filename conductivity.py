import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd
import scipy

#Trajectory file FROM USER INPUT
#file='test/40_out.xyz'
file='../test/lys/lys_100.xyz'


#time step in picoseconds FROM USER INPUT
dt=0.1

class End_of_Loop(Exception): pass #To being able of closing the loop outside the loop functions

def read_line(f):  
	#Is needed to check when the file finish and the exception need to being handled
	line=f.readline()
	#check is line is empty meaning that the file end is reach
	if not line:
		raise End_of_Loop
	else:
		return line

def process_line(line):
	#simple rasing function
	line=[float(l) for l in line.strip('\n').split()]
	#every line is composed by charge , X , Y , Z
	return line[0],line[1:4]


def process_frame(trajectory_file_object):
	#take file object and return position of all atoms in 1 single frame
	n_molecules=read_line(trajectory_file_object)
	n_molecules=int(n_molecules) 						#FIRST LINE is number of molecules in the frame
	trajectory_file_object.readline().strip('\n')  		#SECOND LINE is the blanck line of xyz frame
														#FROM THIRD LINE there are the coordinates
	i=0
	#initialize variables
	q_frame, x_frame, y_frame,z_frame= [],[],[],[]
	#cicle over molecules
	for i in range(n_molecules):
		charge , coordinates = process_line(read_line(trajectory_file_object))	
		q_frame.append(charge)
		x_frame.append(coordinates[0])
		y_frame.append(coordinates[1])
		z_frame.append(coordinates[2])
	return q_frame,x_frame,y_frame,z_frame



#-------------------------------------------------
#DATA
travis_data=pd.read_csv("../test/lys/lys_c_travis.csv",sep='; ',header=0,engine='python')
travis_data.columns = ['t', 'msd', 'derivative']
travis_data_a=pd.read_csv("../test/lys/lys_a_travis.csv",sep='; ',header=0,engine='python')
travis_data_a.columns = ['t', 'msd', 'derivative']
vmd_data_c=pd.read_csv("../test/lys/vmd_c.dat",sep=' ',header=None)
vmd_data_c.columns = ['t', 'msd']
vmd_data_a=pd.read_csv("../test/lys/vmd_a.dat",sep=' ',header=None)
vmd_data_a.columns = ['t', 'msd']
vmd_data_a.t=vmd_data_a.t*1e3
vmd_data_a.msd=vmd_data_a.msd*1e4
vmd_data_c.t=vmd_data_c.t*1e3
vmd_data_c.msd=vmd_data_c.msd*1e4
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

#CHECKING PROPER DIMENSIONS OF X Y Z 
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
	#Use skitlearn tools for making the linear regression
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


def plotting(msd_list,filename='msd.csv'):
	# Simple plotting of raw MSD and linear regression line
	# - require a list of msd [msd_cation,msd_anion]
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
		pd.DataFrame({"t":t,"msd":msd}).to_csv(filename,header=["t","msd"],index=None)
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
#msd1=get_self_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.50)
# print("-----------------------") 
# print("Anion Self-diffusion   ")
# print("-----------------------")
#msd2=get_self_msd(x[:,anion_idx[0]],y[:,anion_idx[0]],z[:,anion_idx[0]],depth=0.50)
# plotting([msd1],"cation_selfdiffusion.csv")
# plotting([msd2],"anion_selfdiffusion.csv")




import time

#i!=j Cation-Cation and Anion-Anion
def inter_same_product(dx,dy,dz):
	#it take dx dy and dz each with dimension [N_frame,N_molecules] and calculate the product mediated over uniques combinations
	
	#bin_coeff is the total number of unique products
	bin_coeff=scipy.special.comb(dx.shape[1],2,exact=True)
	dx_cs=dx[:,::-1].cumsum(axis=1)[:,::-1]-dx
	dx2=np.einsum('ij,ji->i',dx,dx_cs.T)/bin_coeff
	dy_cs=dy[:,::-1].cumsum(axis=1)[:,::-1]-dy
	dy2=np.einsum('ij,ji->i',dy,dy_cs.T)/bin_coeff
	dz_cs=dz[:,::-1].cumsum(axis=1)[:,::-1]-dz
	dz2=np.einsum('ij,ji->i',dz,dz_cs.T)/bin_coeff
	r2=dx2+dy2+dz2
	unit_conversion=1e4
	r2=r2*unit_conversion
	return r2


def product(delta1,delta2):
	#product function for inter ions inter-diffusion 
	#cation_i-anion_j with al possible product ij (N_i*N_j possibilities) 
	tot_number_combination=delta1.shape[1]*delta2.shape[1]
	#create a matrix [NxF] with N times repeated the array of summation over N molecules (axis 1) 
	delta1=np.tile(np.sum(delta1,axis=1),(delta2.shape[1],1))
	#we need only the diagonal component of the matrix multiplication
	delta_square=np.einsum('ij,ji->i',delta2,delta1)
	#mediated over total number combination N_cations*N_anions
	delta_square=delta_square/tot_number_combination
	return delta_square

#i!=j Cation-anion
def inter_product(dx,dy,dz,cation_anion_idx):
	#it require the index of cations and anions for selecting the right subdata from dx dy and dz
	#it calculate the sum of all possible product ij, where i=cations and j=anions	
	cation_idx,anion_idx = cation_anion_idx
	dx2=product(dx[:,cation_idx[0]],dx[:,anion_idx[0]])
	dy2=product(dy[:,cation_idx[0]],dx[:,anion_idx[0]])
	dz2=product(dz[:,cation_idx[0]],dx[:,anion_idx[0]])
	r2=dx2+dy2+dz2
	unit_conversion=1e4
	r2=r2*unit_conversion
	return r2

def get_inter_msd(x,y,z,depth=0.3,cation_anion_idx=None):
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
		if count%50 == 0:
			print("Intervals processed : ",count,"/",max_origin_index)
		#Deviations respect to the reference t0 of the interval
		dx=x[subset_idx,:]-x[subset_idx[0],:]
		dy=y[subset_idx,:]-y[subset_idx[0],:]
		dz=z[subset_idx,:]-z[subset_idx[0],:]

		#SPECIAL PHASE IF i!=j
		start=time.time()
		if cation_anion_idx:
			r2=inter_product(dx,dy,dz,cation_anion_idx)
		else:
			r2=inter_same_product(dx,dy,dz)
		end=time.time()
		#print("time : ",end-start)
		#END SPECIAL PHASE
		msd=msd+r2
		count=count+1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/count
	return msd


# #INTER-IONS  <-----
# msd1 = get_inter_msd(x,y,z,depth=0.7,cation_anion_idx=[cation_idx,anion_idx])

# fig, ax = plt.subplots()
# t=get_t(msd1,dt)
# msd1_pred,t_pred = regression(msd1,t)
# ax.plot(t,msd1,linewidth=1.3,label=r'msd',color='red',zorder=2)
# ax.plot(t_pred,msd1_pred,label='linear regression',linewidth=1,linestyle='dashed',color='blue',zorder=3)
# ax.legend()
# ax.set_xlabel(r'time / ps')
# ax.set_ylabel(r'MSD / pm^2')
# ax.set_box_aspect(1)
# #pd.DataFrame({"t":t,"msd":msd}).to_csv(filename,header=["t","msd"],index=None)
# fig.suptitle("Cation and Anion inter-diffusion MSD")
# plt.show()



#INTERDIFFUSION cation-cation and anion-anion
msd1 = get_inter_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.9)
msd2= get_inter_msd(x[:,anion_idx[0]],y[:,anion_idx[0]],z[:,anion_idx[0]],depth=0.9)


fig, ax = plt.subplots(1,2)
t=get_t(msd1,dt)
msd1_pred,t_pred = regression(msd1,t)
ax[0].plot(t,msd1,linewidth=1.3,label=r'msd',color='red',zorder=2)
ax[0].plot(t_pred,msd1_pred,label='linear regression',linewidth=1,linestyle='dashed',color='blue',zorder=3)
ax[0].legend()
ax[0].set_xlabel(r'time / ps')
ax[0].set_ylabel(r'MSD / pm^2')
ax[0].set_box_aspect(1)
msd2_pred,t_pred = regression(msd2,t)
ax[1].plot(t,msd2,linewidth=1.3,label=r'msd',color='red',zorder=2)
ax[1].plot(t_pred,msd2_pred,label='linear regression',linewidth=1,linestyle='dashed',color='blue',zorder=3)
ax[1].legend()
ax[1].set_xlabel(r'time / ps')
ax[1].set_ylabel(r'MSD / pm^2')
ax[1].set_box_aspect(1)
#pd.DataFrame({"t":t,"msd":msd}).to_csv(filename,header=["t","msd"],index=None)
fig.suptitle("Cation and Anion inter-diffusion MSD")
plt.show()



#plotting([msd])



