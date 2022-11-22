## Conductivity Calculation
Simple script for calculate the diagonal terms of conductivity.

Conductivity can be calculated with:  [1]

![img](/img/tex2img.png)

which can be calculated by the slope of mean square deviation taking into account even the off-diagonal elements i!=j: [2]

![img2](/img/tex2img_2.png)


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
#### Self-diffusion (i=j)

It's easy to calculate and it's the only property from which is still available a reference value from other programs like Travis (http://www.travis-analyzer.de/) and VMD (https://github.com/giorginolab/vmd_diffusion_coefficient). MSD is taken with the simple product between the deviation of the same ions (i=j).

The trajectory is extracted and stored in 4 array:
- x [F,N]
- y [F,N]
- z [F,N]
- q [1,N] (the order of molecules don't change in the trajectory)

F = number of frames

N = number of molecules in the cell

With this structure is easy to find the deviation:

![img](https://github.com/stefarusso/conducivity/blob/master/img/tex2img_4.png?raw=true)

from the deviation is possibile to calculate MSD which is linked to the self-diffusion coefficient: [3,4]


![img](https://github.com/stefarusso/conducivity/blob/master/img/tex2img_3.png?raw=true)

This can be usefull in case like ionic solutes in neutral solvents where there is very little correlation of motion between the ions. In case of ionic liquids this can lend to over-estimation of diffusion coefficient because th of-diagonal terms D_{ij} where the MSD is calculated by the product of the deviation of different ions can be not close to zero (as the solute exemple). This balance the value too positive of D_i giving a more valueble result with the collective D_coll

#### Inter molecular correlations (i!=j)

This is the quantity where ionic liquids start to have deviations from neutral species. In water this quantity is close to zero because solutes dynamics are uncorrelated and is often ignored. In ionic liquids ionic couples can move togheter and the diffusion of single ion can be higher than expected if we don't take into account the correlated movements that happen with the other ions in the simulations.

##### Cation-Cation and Anion-Anion 
i!=j.
the procedure is the same at least up to the deviation matrix DeltaX,DeltaY,DeltaZ. This time there is the need of a loop over the ions. There are N*(N-1) products compared to the N * N in the i=j case

