import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


#file='test/40_out.xyz'
file='../test/200.xyz'

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



#Final coord vectors
x=np.array(x)
y=np.array(y)
z=np.array(z)
q=np.array(q)
print("LOADING TRAJECTORY FILE : COMPLETE")
print("X : ",x.shape)
print("Y : ",y.shape)
print("Z : ",z.shape)
print("q : ",q.shape)


#select the indexes of anion and cations
#is a tuple of array, one for axis
cation_idx=np.where(q == 1)
anion_idx=np.where(q == -1)

#Deviations respect initial t
dx=x[:,:]-x[0,:]
dy=y[:,:]-y[0,:]
dz=z[:,:]-z[0,:]

print(dz)
#Product for selfdiffusion i*i
r2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)

#UNIT CHANGES FROM A^2/ps -> m^2/s
unit_conversion=1e-8
#r2=r2*unit_conversion
print("|R-R0|^2 : ",r2.shape)

#SUBDATA FOR ANION AND CATION
r2_cation=r2[:,cation_idx[0]]
print("CATION SUB DATA : ",r2_cation.shape)
r2_anion=r2[:,anion_idx[0]]
print("ANION SUB DATA : ",r2_cation.shape)

#T VECTOR
t=np.arange(0,r2.shape[0]*dt,dt)
print("T : ", t.shape )

#SIMPLE MEAN OF MSD 
msd_cation=np.sum(r2_cation,axis=1)/r2_cation.shape[1]
print("MSD_CATION : ",msd_cation.shape)
msd_anion=np.sum(r2_anion,axis=1)/r2_anion.shape[1]
print("MSD_ANION : ",msd_anion.shape)


#REGRESSION
#linear regression is done over 60% of the msd
#CATION
scaling=0.3
idx=int(len(t)*scaling)
t_subset=t[idx:].reshape((-1,1))
msd_subset_c=msd_cation[idx:]
print("SUB_T : " ,t_subset.shape)
print("SUB_MSD : ",msd_subset_c.shape)
model_c = LinearRegression().fit(t_subset,msd_subset_c)
print(f"D CATION : {model_c.coef_*6} m^2/s")
msd_pred_c = model_c.predict(t_subset)
#ANION
msd_subset_a=msd_anion[idx:]
print("SUB_T : " ,t_subset.shape)
print("SUB_MSD : ",msd_subset_a.shape)
model_a = LinearRegression().fit(t_subset,msd_subset_a)
print(f"D ANION : {model_a.coef_*6} m^2/s")
msd_pred_a = model_a.predict(t_subset)


#PLOTTING
fig, ax = plt.subplots(1,2)
ax[0].plot(t,msd_cation,linewidth=1.5,label=r'msd mean',color='red',zorder=2)
for mol in r2_cation.T:
	 ax[0].plot(t,mol,linewidth=0.5,alpha=0.2,color='black',zorder=1)
ax[0].plot(t_subset,msd_pred_c,label='linear regression',linewidth=1,linestyle='dashed',color='orange',zorder=3)
ax[0].legend()
ax[0].set_title("Cation MSD")
ax[0].set_xlabel(r'time / ps')
ax[0].set_ylabel(r'MSD / m^2/s')
#ax[0].set_ylim([0,max(msd_cation)*3])
ax[0].set_box_aspect(1)

ax[1].plot(t,msd_anion,linewidth=1.5,label=r'msd mean',color='red',zorder=2)
for mol in r2_anion.T:
	 ax[1].plot(t,mol,linewidth=0.5,alpha=0.2,color='black',zorder=1)
ax[1].plot(t_subset,msd_pred_a,label='linear regression',linewidth=1,linestyle='dashed',color='orange',zorder=3)
ax[1].legend()
ax[1].set_title("Anion MSD")
#ax[1].set_ylim([0,max(msd_anion)*3])
ax[1].set_xlabel(r'time / ps')
ax[1].set_ylabel(r'MSD / m^2/s')
ax[1].set_box_aspect(1)
fig.tight_layout()
plt.show()








# #correlation depth 1/3
# depth=0.33
# #CATION
# print("TEST CATION SELF DIFFUSION")
# x_c=x[:,cation_idx[0]]
# print("cation X : ",x_c.shape)
# y_c=y[:,cation_idx[0]]
# print("cation Y : ",y_c.shape)
# z_c=z[:,cation_idx[0]]
# print("cation Z : ",z_c.shape)




# #loop over subsets
# subset=np.arange(0,int(depth*x_c.shape[0]))



# x=x_c[subset,:]
# y=y_c[subset,:]
# z=z_c[subset,:]
# print("subset X : ",x.shape)
# print("subset Y : ",y.shape)
# print("subset Z : ",z.shape)


# #Deviations respect initial t
# dx=x[:,:]-x[0,:]
# dy=y[:,:]-y[0,:]
# dz=z[:,:]-z[0,:]

# #Product for selfdiffusion i*i
# r2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
# print("|R-R0|^2 : ",r2.shape)
# #mean over molecules
# r2=np.sum(r2,axis=1)/r2.shape[1]
# #UNIT CHANGES FROM A^2/ps -> m^2/s
# unit_conversion=1e-8
# r2=r2*unit_conversion
# print("|R-R0|^2 MEAN : ",r2.shape)









# #T VECTOR
# t=np.arange(0,r2.shape[0]*dt,dt)
# print("T : ", t.shape )




