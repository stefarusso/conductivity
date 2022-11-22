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

def regression(msd,t,scaling=0.3):
	#
	#Use skitlearn tools for making the linear regression
	#
	idx=int(len(t)*scaling)
	t_pred=t[idx:].reshape((-1,1))
	msd_subset=msd[idx:]
	print("SUB_T : " ,t_pred.shape)
	print("SUB_MSD : ",msd_subset.shape)
	model_c = LinearRegression().fit(t_pred,msd_subset)
	print(f"slope CATION : {model_c.coef_} ")
	print(f"Intercept: {model_c.intercept_} ")
	print(f"D CATION : {model_c.coef_/6} pm^2/ps")
	msd_pred = model_c.predict(t_pred)
	return msd_pred, t_pred

def plotting(t,msd,t_pred,msd_pred):
	#
	# Simple plotting of raw MSD and linear regression line
	#
	fig, ax = plt.subplots()
	ax.plot(t,msd,linewidth=1.3,label=r'msd',color='red',zorder=2)
	ax.plot(t_pred,msd_pred,label='linear regression',linewidth=1,linestyle='dashed',color='blue',zorder=3)
	ax.legend()
	ax.set_title("Cation MSD")
	ax.set_xlabel(r'time / ps')
	ax.set_ylabel(r'MSD / pm^2')
	plt.show()



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
if x.shape == y.shape == z.shape and q.shape[0] == x.shape[1] == y.shape[1] == z.shape[1]:
	print("Total Number of Frame : ",x.shape[0])
	print("Total length of trajectory : ",x.shape[0]*dt, "ps")
	print("Total number of Molecules : ",x.shape[1])
else:
	raise Exception("SOMETHING WRONG WITH THE TRAJECTORY")
#---------------------------
#END OF TRAJECTORY LOADING
#---------------------------






#-------------------------------------------------------------
#CORRELATION DEPTH
#is the percentage of the trajectory in which the correlation between ions is took in account
depth=0.97

print("-----------------------")
print("TEST CORRELATION DEPTH : ",depth*100," %")
print("-----------------------")

#CATION SUBDATA
x_c=x[:,cation_idx[0]]
y_c=y[:,cation_idx[0]]
z_c=z[:,cation_idx[0]]

#loop over subsets
subset_idx = np.arange(0,int(depth*x_c.shape[0]))
max_origin_index = x_c.shape[0]-len(subset_idx)
print("Number of intervals : ",max_origin_index)

#
msd = np.zeros(len(subset_idx))
i=0
#LOOP OVER INTERVALS
while subset_idx[0]<max_origin_index :
	if i%100 == 0:
		print("Intervals processed : ",i)
	x=x_c[subset_idx,:]
	y=y_c[subset_idx,:]
	z=z_c[subset_idx,:]

	#Deviations respect to the reference t0 of the interval
	dx=x[:,:]-x[0,:]
	dy=y[:,:]-y[0,:]
	dz=z[:,:]-z[0,:]

	#Product for selfdiffusion i*i
	r2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
	#mean over N molecules
	r2=np.sum(r2,axis=1)/r2.shape[1]
	#UNIT CHANGES FROM A^2/ps -> pm^2/ps
	unit_conversion=1e4
	r2=r2*unit_conversion
	msd=msd+r2
	i+=1
	subset_idx=subset_idx+1
print("All Interval processed")
#MSD MEAN OVER TOTAL NUMBER OF INTERVALS
msd=msd/i

#T VECTOR
t_interval=np.arange(0,msd.shape[0]*dt,dt)
msd_pred,t_pred = regression(msd,t_interval)
plotting(t_interval,msd,t_pred,msd_pred)


# #PLOTTING
# fig, ax = plt.subplots(1,2)
# ax[0].plot(t_interval,msd,linewidth=1.5,label=r'msd mean',color='red',zorder=2)
# ax[0].plot(t_subset,msd_pred,label='linear regression',linewidth=1,linestyle='dashed',color='orange',zorder=3)
# ax[0].plot(travis_data.t,travis_data.msd,linewidth=1,linestyle='dotted',label=r'msd travis',color='blue',zorder=3)
# ax[0].plot(vmd_data_c.t,vmd_data_c.msd,linewidth=1.4,linestyle='dotted',label=r'msd vmd',color='purple',zorder=3)
# ax[0].legend()
# ax[0].set_title("Cation MSD")
# ax[0].set_xlabel(r'time / ps')
# ax[0].set_ylabel(r'MSD / pm^2')
# ax[0].set_box_aspect(1)
# plt.show()

