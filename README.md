## Conductivity Calculation
Simple script for calculate the diagonal terms of conductivity.

Conductivity can be calculated with:

[![img](https://bit.ly/3fY21xL)](#)

which can be calculated by the slope of mean square deviation taking into account even the off-diagonal elements i!=j:

[![img](https://bit.ly/3WVaVwt)](#)


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
