## Conductivity Calculation
Simple script for calculate the diagonal terms of conductivity.

Conductivity can be calculated with:

![img](tex2img.png)

which can be calculated by the slope of mean square deviation taking into account even the off-diagonal elements i!=j:

![img2](tex2img_2.png)


#### Input data

data need to be altrady processed with Travis analyzer and has to have this formate where every line of the frame is the center of mass of one ion and the first column is the total ion charge of the molecules.
Movement of the center of mass of the cell is removed and is ceneterd on the center of mass in [0.0.0] 


Example of input data:

charge   ,X_coord      ,Y_coord      ,Z-coord
```
4
  
1     10.94772714   13.34369366   12.45011788
-1    -32.80060706  -27.17604627  -13.32500798
1     17.16800035  -14.08190288  -24.77136572
-1    -22.58886859  -28.62449835   14.05828311
```

#### Usage
python cond.py input.xyz

## How it works
#### It's just for me, it's for remembering what I have written in the code (you can ignore it)

The trajectory is extracted and stored in 4 array:
- x [F,N]
- y [F,N]
- z [F,N]
- q [1,N] (the order of molecules don't change in the trajectory)

F = number of frames

N = number of molecules in the cell

With this structure is easy to find the deviation:

![img](http://www.sciweavers.org/tex2img.php?eq=%20x%3D%0A%5Cbegin%7Bbmatrix%7D%0Ax_%7Bo%7D%5E%7B0%7D%20%26%20x_%7B1%7D%5E%7B0%7D%20%26%20..%20%26%20x_%7BN-1%7D%5E%7B0%7D%20%26%20x_%7BN%7D%5E%7B0%7D%5C%5C%20%0Ax_%7B0%7D%5E%7B1%7D%26%20%20%26%20%20%26%20%20%26%20%5C%5C%0A%5Cvdots%26%26%26%26%20%5C%5C%0Ax_%7B0%7D%5E%7BF-1%7D%26%26%26%26%20%5C%5C%0Ax_%7B0%7D%5E%7BF%7D%26%26%26%26%20%5C%5C%0A%5Cend%7Bbmatrix%7D%20%0A%5C%2C%20%0A%5CDelta%20x%20%3D%20%0A%5Cbegin%7Bbmatrix%7D%0Ax_%7Bo%7D%5E%7B0%7D%20-%20x_%7Bo%7D%5E%7B0%7D%20%26%20x_%7B1%7D%5E%7B0%7D%20-x_%7B1%7D%5E%7B0%7D%20%26%20..%20%26%20x_%7BN-1%7D%5E%7B0%7D%20-%20x_%7BN-1%7D%5E%7B0%7D%20%26%20x_%7BN%7D%5E%7B0%7D%20-%20x_%7BN%7D%5E%7B0%7D%5C%5C%20%0Ax_%7B0%7D%5E%7B1%7D%20-%20x_%7Bo%7D%5E%7B0%7D%26%20%20%26%20%20%26%20%20%26%20%5C%5C%0A%5Cvdots%26%26%26%26%20%5C%5C%0Ax_%7B0%7D%5E%7BF-1%7D%20-%20x_%7Bo%7D%5E%7B0%7D%26%26%26%26%20%5C%5C%0Ax_%7B0%7D%5E%7BF%7D%20-%20x_%7Bo%7D%5E%7B0%7D%26%26%26%26%20%5C%5C%0A%5Cend%7Bbmatrix%7D&bc=White&fc=Black&im=jpg&fs=12&ff=modern&edit=0)

from the deviation is possibile to calculate MSD which is linked to the self-diffusion coefficient:

![img](https://bit.ly/3X7713Z)




