import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd
import scipy
import time

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
	q_frames, x_frames, y_frames,z_frames= [],[],[],[]
	#cicle over molecules
	for i in range(n_molecules):
		charge , coordinates = process_line(read_line(trajectory_file_object))	
		q_frames.append(charge)
		x_frames.append(coordinates[0])
		y_frames.append(coordinates[1])
		z_frames.append(coordinates[2])
	return q_frames,x_frames,y_frames,z_frames





def load_trajectory(filename):
	#TRAJECTORY LOADING
	if 'dt' not in globals():
		global dt
		dt=float(input("Input the physical time-step in ps : "))
	try:
		with open(filename,'r') as trajectory_file_object:
			q,x,y,z=[],[],[],[]
			print("LOADING TRAJECTORY FILE")
			while True:
				q_new,x_new,y_new,z_new=process_frame(trajectory_file_object)
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
	#ALL TRAJECTORY LOADED

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

	return x,y,z,q,cation_idx,anion_idx


def regression(msd,t,scaling=0.3):
	#Use skitlearn tools for making the linear regression

	#it takes the subset of time vector and MSD to perform linear regression
	idx=int(len(t)*scaling)
	t_pred=t[idx:].reshape((-1,1))
	msd_subset=msd[idx:]
	#LINEAR REGRESSION
	model_c = LinearRegression().fit(t_pred,msd_subset)
	slope=model_c.coef_
	D=model_c.coef_/6
	intercept=model_c.intercept_
	print("LINEAR REGRESSION ON THE LAST ",100-scaling*100," %")
	print(f"slope : {model_c.coef_} ")
	print(f"Intercept: {model_c.intercept_} ")
	print(f"D : {model_c.coef_/6} pm^2/ps")
	#generate msd_predition point for plotting with the same spacing and interval of t_prediction
	msd_pred = model_c.predict(t_pred)
	return [msd_pred,t_pred],[slope,intercept]

def get_t(msd):
	#create the proper time vector
	return np.arange(0,msd.shape[0]*dt,dt)


def plotting(msd_list,filename='msd.csv'):
	# Simple plotting of raw MSD and linear regression line
	# !!!! require a list of msd_list = [msd_cation,msd_anion] !!!!
	fig, ax = plt.subplots(1,2)
	for i,msd in enumerate(msd_list):
		t=get_t(msd)
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

#----------------

def self_product(dx,dy,dz):
	#Product for selfdiffusion i*i cation-cation or anion-anion
	#dx,dy and dz are matrix [F,N]
	# return the DeltaR_i*DeltaR_i product
	#
	# we can skip the root calculation in this case
	#DR = square_root(DX^2 + DY^2 + DZ^2)
	#DR2 = square_root(DX^2 + DY^2 + DZ^2)^2
	#DR2=DX^2 + DY^2 + DZ^2
	dr2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
	#mean over N molecules
	dr2=np.sum(dr2,axis=1)/dr2.shape[1]
	#UNIT CHANGES FROM A^2/ps -> pm^2/ps
	unit_conversion=1e4
	dr2=dr2*unit_conversion
	return dr2


def get_selfdiffusion_msd(x,y,z,depth=0.3):
	#CORRELATION DEPTH
	#is the percentage of the trajectory in which the correlation between ions is took in account
	
	#loop over subsets
	print("Correlation depth : ",depth*100," %")
	#subset_idx contains all the indexes for the subset of lenght given by the correlation depth
	subset_idx = np.arange(0,int(depth*x.shape[0]))
	max_origin_index = x.shape[0]-len(subset_idx)
	print("Number of intervals : ",max_origin_index)
	#initialize the msd sum
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
		dr2=self_product(dx,dy,dz)
		msd=msd+dr2
		count+=1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/count
	return msd


################################### IF WE USE DISTANCE R*R= (sqrt(X1^2 +Y^2+ Z^2))^2
#i!=j Cation-Cation and Anion-Anion
def DISTANCE_interdiffusion_same_product(dx,dy,dz):
	#i!=j Cation-Cation or Anion,Anion
	#
	#it take dx dy and dz each with dimension [N_frame,N_molecules]
	#get the distance DeltaR and make the product DR_i*DR_j=DR^2
	#the product is mediated over uniques combinations (binomial coefficient)
	
	#bin_coeff is the total number of unique products
	#it is required to mediate over all the possible unique combinations
	dr2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
	dr=np.sqrt(dr2)

	bin_coeff=scipy.special.comb(dx.shape[1],2,exact=True)

	dr_cs=dr[:,::-1].cumsum(axis=1)[:,::-1]-dr
	dr=np.einsum('ij,ji->i',dr,dr_cs.T)/bin_coeff
	unit_conversion=1e4
	dr=dr*unit_conversion
	return dr


def DISTANCE_inter_product(delta1,delta2):
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


def DISTANCE_interdiffusion_inter_product(dx,dy,dz,cation_anion_idx):
	#i!=j Cation-Anion
	#
	#it require the index of cations and anions for selecting the right subdata from dx dy and dz
	#it calculate the sum of all possible product ij, where i=cations and j=anions	
	cation_idx,anion_idx = cation_anion_idx
	#single coordinates already mediated over the total number of products (N_cations * N_anions)
	dr2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
	dr=np.sqrt(dr2)
	dr=DISTANCE_inter_product(dr[:,cation_idx[0]],dr[:,anion_idx[0]])
	unit_conversion=1e4
	dr=dr*unit_conversion
	return dr



#################################### IF WE USE PRODUCT R*R=X1*X2+Y1*Y2+Z1*Z2)

#i!=j Cation-Cation and Anion-Anion
def interdiffusion_same_product(dx,dy,dz):
	#i!=j Cation-Cation or Anion,Anion
	#
	#it take dx dy and dz each with dimension [N_frame,N_molecules]
	#get the distance DeltaR and make the product DR_i*DR_j=DR^2
	#the product is mediated over uniques combinations (binomial coefficient)
	
	#bin_coeff is the total number of unique products
	#it is required to mediate over all the possible unique combinations
	bin_coeff=scipy.special.comb(dx.shape[1],2,exact=True)

	dx_cs=dx[:,::-1].cumsum(axis=1)[:,::-1]-dx
	dx2=np.einsum('ij,ji->i',dx,dx_cs.T)/bin_coeff
	dy_cs=dy[:,::-1].cumsum(axis=1)[:,::-1]-dy
	dy2=np.einsum('ij,ji->i',dy,dy_cs.T)/bin_coeff
	dz_cs=dz[:,::-1].cumsum(axis=1)[:,::-1]-dz
	dz2=np.einsum('ij,ji->i',dz,dz_cs.T)/bin_coeff
	dr2=dx2+dy2+dz2
	unit_conversion=1e4
	dr2=dr2*unit_conversion
	return dr2


def inter_product(delta1,delta2):
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


def interdiffusion_inter_product(dx,dy,dz,cation_anion_idx):
	#i!=j Cation-Anion
	#
	#it require the index of cations and anions for selecting the right subdata from dx dy and dz
	#it calculate the sum of all possible product ij, where i=cations and j=anions	
	cation_idx,anion_idx = cation_anion_idx
	dx2=inter_product(dx[:,cation_idx[0]],dx[:,anion_idx[0]])
	dy2=inter_product(dy[:,cation_idx[0]],dx[:,anion_idx[0]])
	dz2=inter_product(dz[:,cation_idx[0]],dx[:,anion_idx[0]])
	dr2=dx2+dy2+dz2
	unit_conversion=1e4
	dr2=dr2*unit_conversion
	return dr2

#INTERDIFFUSION WHERE i!=j
#two possible cases :
#-same ion type Cation-Cation and Anion-Anion
#-inter ion type Cation-Anion
def get_interdiffusion_msd(x,y,z,depth=0.3,cation_anion_idx=None):
	#CORRELATION DEPTH
	#is the percentage of the trajectory in which the correlation between ions is took in account	
	#
	#This function calculates the MSD when i!=j 
	#In the case same-ion (cation-cation and anion-anion) inter-ion (cation-anion)

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
			dr2=interdiffusion_inter_product(dx,dy,dz,cation_anion_idx)
		else:
			dr2=interdiffusion_same_product(dx,dy,dz)
		end=time.time()
		#print("time : ",end-start)
		#END SPECIAL PHASE
		msd=msd+dr2
		count=count+1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/count
	return msd


#COLLECTIVE, ALL TERMS
#-same ion type Cation-Cation and Anion-Anion
#-inter ion type Cation-Anion
def get_collective_msd(x,y,z,depth=0.3):
	#CORRELATION DEPTH
	#is the percentage of the trajectory in which the correlation between ions is took in account	
	#
	#This function calculates  all the MSD terms

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
		dr2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
		dr=np.sqrt(dr2)
		dr=inter_product(dr,dr)
		unit_conversion=1e4
		dr=dr*unit_conversion
		end=time.time()
		#print("time : ",end-start)
		#END SPECIAL PHASE
		msd=msd+dr
		count=count+1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/count
	return msd

#INTERDIFFUSION WHERE i!=j
#two possible cases :
#-same ion type Cation-Cation and Anion-Anion
#-inter ion type Cation-Anion
def DISTANCE_get_interdiffusion_msd(x,y,z,depth=0.3,cation_anion_idx=None):
	#CORRELATION DEPTH
	#is the percentage of the trajectory in which the correlation between ions is took in account	
	#
	#This function calculates the MSD when i!=j 
	#In the case same-ion (cation-cation and anion-anion) inter-ion (cation-anion)

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
			dr2=DISTANCE_interdiffusion_inter_product(dx,dy,dz,cation_anion_idx)
		else:
			dr2=DISTANCE_interdiffusion_same_product(dx,dy,dz)
		end=time.time()
		#print("time : ",end-start)
		#END SPECIAL PHASE
		msd=msd+dr2
		count=count+1
		subset_idx=subset_idx+1
	print("All Interval processed")
	#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
	msd=msd/count
	return msd



def all(filename):
	#loading
	x,y,z,q,cation_index,anion_index = load_trajectory(filename)
	print("_________INTER CAT-CAT MSD________ ")
	# #INTERDIFFUSION same ion : cation-cation and anion-anion
	inter_1 = get_interdiffusion_msd(x[:,cation_index[0]],y[:,cation_index[0]],z[:,cation_index[0]],depth=0.7)
	t = get_t(inter_1)
	pd.DataFrame({"t":t,"msd":inter_1}).to_csv("inter1.csv",header=["t","msd"],index=None)
	print("_________INTER ANI-ANI MSD________ ")
	inter_2 = get_interdiffusion_msd(x[:,anion_index[0]],y[:,anion_index[0]],z[:,anion_index[0]],depth=0.7)
	pd.DataFrame({"t":t,"msd":inter_2}).to_csv("inter2.csv",header=["t","msd"],index=None)
	#INTER-IONS
	print("_________INTER-INTER CAT-ANI MSD________")
	inter_inter = get_interdiffusion_msd(x,y,z,depth=0.7,cation_anion_idx=[cation_index,anion_index])
	pd.DataFrame({"t":t,"msd":inter_inter}).to_csv("inter_inter.csv",header=["t","msd"],index=None)
	
	print("_________SELF-Cation MSD________ ")
	cation = get_selfdiffusion_msd(x[:,cation_index[0]],y[:,cation_index[0]],z[:,cation_index[0]],depth=0.70)
	pd.DataFrame({"t":t,"msd":cation}).to_csv("cat.csv",header=["t","msd"],index=None)
	print("_________SELF-anion MSD________ ")
	anion = get_selfdiffusion_msd(x[:,anion_index[0]],y[:,anion_index[0]],z[:,anion_index[0]],depth=0.70)
	pd.DataFrame({"t":t,"msd":anion}).to_csv("anion.csv",header=["t","msd"],index=None)
	
	print("CATION REGRESSION:")
	[pred_1,t_1],[slope_1,intercept_1] = regression(cation,t,scaling=0.3)
	print("Anion REGRESSION:")
	[pred_2,t_2],[slope_2,intercept_2] = regression(anion,t,scaling=0.3)
	
	print("Inter1 REGRESSION:")
	[pred_3,t_3],[slope_3,intercept_3] = regression(inter1,t,scaling=0.3)
	print("Inter2 REGRESSION:")
	[pred_4,t_4],[slope_4,intercept_4] = regression(inter2,t,scaling=0.3)
	print("Inter_inter REGRESSION:")
	[pred_5,t_5],[slope_5,intercept_5] = regression(inter_inter,t,scaling=0.3)
	
	tot_coll=cation+anion+inter1+inter2+inter_inter	
	print("SUM cat+anion+inter1+inter2+inter_inter REGRESSION:")
	[pred_6,t_6],[slope_6,intercept_6] = regression(tot_coll,t,scaling=0.3)
	
	fig, ax = plt.subplots()
	ax.plot(t,cation,linewidth=1.3,label=r'msd cation',color='blue',zorder=2)
	ax.plot(t_1,pred_1,label=f'linear regression, y=x {slope_1}  +({intercept_1}) pm^2/ps \n D: {slope_1/6} ',linewidth=1,linestyle='dashed',color='blue',zorder=3)
	ax.plot(t,anion,linewidth=1.3,label=r'msd anion',color='green',zorder=2)
	ax.plot(t_2,pred_2,label=f'linear regression, y=x {slope_2}  +({intercept_2}) pm^2/ps \n D: {slope_2/6} ',linewidth=1,linestyle='dashed',color='green',zorder=3)
	ax.plot(t,tot_coll,linewidth=1.3,label=r'msd sum',color='purple',zorder=2)
	ax.plot(t_6,pred_6,label=f'linear regression, y=x {slope_6}  +({intercept_6}) pm^2/ps \n D: {slope_6/6} ',linewidth=1,linestyle='dashed',color='purple',zorder=3)
	ax.plot(t,inter_1,linewidth=1.3,label=r'msd cation-cation',color='brown',zorder=2)
        ax.plot(t_3,pred_3,label=f'linear regression, y=x {slope_3}  +({intercept_3}) pm^2/ps \n D: {slope_3/6} ',linewidth=1,linestyle='dashed',color='brown',zorder=3)
	ax.plot(t,inter_2,linewidth=1.3,label=r'msd anion-anion',color='orange',zorder=2)
        ax.plot(t_4,pred_4,label=f'linear regression, y=x {slope_4}  +({intercept_4}) pm^2/ps \n D: {slope_4/6} ',linewidth=1,linestyle='dashed',color='orange',zorder=3)
	ax.plot(t,inter_inter,linewidth=1.3,label=r'msd cation-anion',color='red',zorder=2)
        ax.plot(t_5,pred_5,label=f'linear regression, y=x {slope_5}  +({intercept_5}) pm^2/ps \n D: {slope_5/6} ',linewidth=1,linestyle='dashed',color='red',zorder=3)
	ax.legend()
	ax.set_xlabel(r'time / ps')
	ax.set_ylabel(r'MSD / pm^2')
	ax.set_box_aspect(1)
	fig.savefig('plot.pdf', format='pdf')
	plt.show()





if __name__ == "__main__":
	#TESTING
	travis_data=pd.read_csv("../test/lys/lys_c_travis.csv",sep='; ',header=0,engine='python')
	travis_data.columns = ['t', 'msd', 'derivative']
	vmd_data_c=pd.read_csv("../test/lys/vmd_c.dat",sep=' ',header=None)
	vmd_data_c.columns = ['t', 'msd']
	vmd_data_c.t=vmd_data_c.t*1e3
	vmd_data_c.msd=vmd_data_c.msd*1e4
	filename='../test/lys/lys_100.xyz'
	dt=0.1
	x,y,z,q,cation_idx,anion_idx=load_trajectory(filename)
	#SELFDIFFUSION
	# CATION
	msd_self_1=get_selfdiffusion_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.70)
	t=get_t(msd_self_1,dt)
	# ANION
	msd_self_2=get_selfdiffusion_msd(x[:,anion_idx[0]],y[:,anion_idx[0]],z[:,anion_idx[0]],depth=0.70)

	# #INTERDIFFUSION same ion : cation-cation and anion-anion
	msd_inter_1 = get_interdiffusion_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.7)
	msd_inter_2 = get_interdiffusion_msd(x[:,anion_idx[0]],y[:,anion_idx[0]],z[:,anion_idx[0]],depth=0.7)

	# # #INTER-IONS  <-----
	msd_inter_inter = get_interdiffusion_msd(x,y,z,depth=0.7,cation_anion_idx=[cation_idx,anion_idx])

	# #COLLECTIVE  <-----
	msd_collective = get_collective_msd(x,y,z,depth=0.7)

	fig, ax = plt.subplots()
	[msd_collective_pred,t_coll_pred],[slope,intercept] = regression(msd_collective,t)
	ax.plot(t,msd_self_1,linewidth=1.3,label=r'msd_self_cation',color='blue',zorder=2)
	ax.plot(t,msd_self_2,linewidth=1.3,label=r'msd_self_anion',color='green',zorder=2)
	ax.plot(t,msd_inter_1,linewidth=1.3,label=r'msd_cation-cation',color='purple',zorder=2)
	ax.plot(t,msd_inter_2,linewidth=1.3,label=r'msd_anion-anion',color='grey',zorder=2)
	ax.plot(t,msd_inter_inter,linewidth=1.3,label=r'msd_cation-anion',color='orange',zorder=2)
	ax.plot(t,msd_collective,linewidth=1.3,label=r'msd_coll',color='red',zorder=2)
	ax.plot(t_coll_pred,msd_collective_pred,label=f'linear regression, y=x {slope}  +({intercept}) pm^2/ps \n D: {slope/6} ',linewidth=1,linestyle='dashed',color='black',zorder=3)
	ax.legend()
	ax.set_xlabel(r'time / ps')
	ax.set_ylabel(r'MSD / pm^2')
	ax.set_box_aspect(1)
	#pd.DataFrame({"t":t,"msd":msd}).to_csv(filename,header=["t","msd"],index=None)
	fig.suptitle("Cation and Anion inter-diffusion MSD")
	plt.show()
else:
	print("Usable Packages: ")
	print("x,y,z,q,cation_index,anion_index = conductivity.load_trajectory(\"filename\")")
	print("msd = conductivity.get_selfdiffusion_msd(x[:,cation_idx[0]],y[:,cation_idx[0]],z[:,cation_idx[0]],depth=0.70)")
	print("msd_collective = conductivity.get_collective_msd(x,y,z,depth=0.7)")
	print("t = get_t(msd)")
	print("msd_prediction,t_predtion = conductivity.regression(msd,t)")
	print("conductivity.plotting([msd1,msd2])")
	print("conductivity.collective(filename)")

