import numpy as np
#file='test/40_out.xyz'
file='test/200.xyz.xyz'

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
	for i in range(2000):
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


x=np.array(x)
y=np.array(y)
z=np.array(z)
q=np.array(q)

#is a tuple of array, one for axis
cation_idx=np.where(q == 1)
print(cation_idx[0].shape)
anion_idx=np.where(q == 1)

print("LOADING TRAJECTORY FILE : COMPLETE")
print("X : ",x.shape)
print("Y : ",y.shape)
print("Z : ",z.shape)
print("q : ",q.shape)

dx=x[:,:]-x[0,:]
dy=y[:,:]-y[0,:]
dz=z[:,:]-z[0,:]

r2=np.multiply(dx,dx)+np.multiply(dy,dy)+np.multiply(dz,dz)
#UNIT CHANGES FROM A^2/ps -> m^2/s
r2=r2*1e-08
print("|R-R0|^2 : ",r2.shape)

r2_cation=r2[:,cation_idx[0]]
print("CATION SUB DATA : ",r2_cation.shape)
r2_anion=r2[:,anion_idx[0]]
print("ANION SUB DATA : ",r2_cation.shape)

t=np.arange(0,r2.shape[0]*dt,dt)

msd_cation=np.sum(r2_cation,axis=1)/r2_cation.shape[1]
print("MSD_CATION : ",msd_cation.shape)
msd_anion=np.sum(r2_anion,axis=1)/r2_anion.shape[1]
print("MSD_ANION : ",msd_anion.shape)

import matplotlib
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(t,msd_cation,linewidth=1.5       ,label=r'msd',color='red',zorder=2)
for mol in r2_cation.T:
	 ax.plot(t,mol,linewidth=0.5,alpha=0.2,label=r'msd',color='black',zorder=1)
plt.show()


