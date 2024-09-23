#!/Users/stefano/anaconda3/bin/python
class trajectory:
    def __init__(self):
        self.traj = [] #è un vettore [x,y,z,q,atom]
    def load(self):
        #trajectory loading routine
        q = 1.000
        atom = "C"
        x = 1.20
        y = 3.0
        z = -2.78
        self.traj = [x,y,z,q,atom]


class Logger():
    def __init__(self):
        import logging
        self.logger = logging.getLogger('log')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(message)s')
        if not len(self.logger.handlers):
            consoleHandler = logging.StreamHandler()
            consoleHandler.setLevel(logging.INFO)
            consoleHandler.setFormatter(formatter)
            file_handler = logging.FileHandler('logs.log')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
            self.logger.addHandler(consoleHandler)
    def print(self,message):
        self.logger.info(message)

class End_of_Loop(Exception): pass #To being able of closing the loop outside the loop functions

def parsline(file):
    line=file.readline()
    if not line:
        raise End_of_Loop
    else:
        line = line.strip('\n')
        #gromacs line format "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
        res_num = int(line[0:5])
        res_name, atom_name  = [ line[c:c+5].strip()  for c in range(5,5*2+1,5)]
        atom_name = ''.join([ i for i in atom_name if not i.isdigit()]) #remove digits from the atom name
        atom_num = int(line[5*3:5*3+5]) 
        x, y, z = [ float(line[c:c+8])  for c in range(5*4,5*4+8*3,8)]
        #end gromacs
        return x, y, z, res_num, res_name, atom_name


def read_line(f):
	#Is needed to check when the file finish and the exception need to being handled
	line=f.readline()
	#check is line is empty meaning that the file end is reach
	if not line:
		raise End_of_Loop
	else:
		return line

def load_frame(filename):
    import numpy as np
    try:
        with open(filename, 'r') as file:
            #repeat for every frame up to the EOF
            CM = [] #easier solution for appending np.arrays to an empty array without defining dimension first-hand.
            while True:
                read_line(file) #first line is a comment
                n_atoms = int(read_line(file).strip())
                sum_cm = np.zeros(4)  #x,y,z, mass
                res_num_count = 1
                cm = np.empty((0,3),float) #x,y,z, res_name
                res_names = []
                for i in range(0,n_atoms):
                    x, y, z, res_num, res_name, atom_name = parsline(file)
                    if res_num_count == res_num: #same molecule
                        sum_cm[0] += x * atomic_mass[atom_name]
                        sum_cm[1] += y * atomic_mass[atom_name]
                        sum_cm[2] += z * atomic_mass[atom_name]
                        sum_cm[3] += atomic_mass[atom_name]
                        if i == n_atoms-1 : #last line before end of frame
                            cm = np.append(cm, np.reshape(sum_cm[0:3] / sum_cm[3], (1, 3)), axis=0)
                            res_names.append(res_name)
                        else:
                            pass
                    else: #new molecule
                        res_num_count = res_num     #update the counter
                        cm = np.append(cm, np.reshape(sum_cm[0:3] / sum_cm[3], (1,3)),axis=0)                     #save center of mass for the molecule
                        sum_cm = np.array([x* atomic_mass[atom_name], y* atomic_mass[atom_name], z* atomic_mass[atom_name], atomic_mass[atom_name]])                   #start new molecule
                        res_names.append(res_name)
                read_line(file)  # last line is a cell dimensions
                CM.append(cm)
                return np.array(CM), res_names
        return np.array(CM), res_names
    except End_of_Loop:
        pass #all trajectory processed

atomic_mass = { 'H':1.008, 'C':12.011, 'O':15.999, 'N':14.0067, 'Cl':35.453, 'Al':26.9815, 'F':19.9984 }

import numpy as np
cm,res_names=load_frame("test_files/test.gro")

