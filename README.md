## Conductivity Calculation
Simple script for calculate the diagonal terms of conductivity.

Conductivity can be calculated with: 

```math
\sigma= \frac{e^2 }{V_{cell} \, K_b \, T } \lim_{t\to\infty} \frac{1}{6t} \, \frac{1}{N} \sum_{i,j}^{N} z_iz_j < [ R_i(t)-R_i(t_0) ] [ R_j(t)-R_j(t_0) ]>_t      \qquad[1]  \\

``` 

Where $V_{cell}$ is the volume of the simulation cell in angstrom, $T$ is temperature in kelvin, $t$ is the simulation time, $N$ is the total number of combinations ij and $R_i(t)$ is the position of center of mass of molecule $i$ at time $t$. 

which can be calculated by the slope of mean square deviation taking into account even the off-diagonal elements $i \neq j $:

```math
\begin{align}
\sigma &= \frac{e^2 }{V_{cell} \, K_b \, T } \, D_{coll}         &\qquad[2]  \\
D_{coll}&= \lim_{t\to\infty} \frac{1}{6t} \, \frac{1}{N} \sum_{i,j}^{N} z_iz_j < [ R_i(t)-R_i(t_0) ] [ R_j(t)-R_j(t_0) ]>_t      &\qquad[3]  \\
D_{coll}&=D_{j} + D_{i} +D_{ij}  &\qquad[4]
\end{align}
```

The average over $t$ is done by sectioning the trajectory in smaller sample which size is defined by the $correlation\ depth$. The same quantity is calcolated for all the samples and than is mediated over the number of total samples, that is different from $N$ which is the total number of combination between i and j double sum.

#### Input data

data need to be altrady processed with Travis analyzer and has to have this formate where every line of the frame is the center of mass of one ion and the first column is the total ion charge of the molecules.
Movement of the center of mass of the cell is removed and the center of mass is centered on [0.0.0] 


Example of input data:

$n_{molecules}$

$charge, \qquad X_{coord}, \qquad Y_{coord}, \qquad Z_{coord} $. 
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
#### Self-diffusion $\boldsymbol{\sigma_i}$ (i=j)

It's easy to calculate and it's the only property from which is still available a reference value from other programs like Travis (http://www.travis-analyzer.de/) and VMD (https://github.com/giorginolab/vmd_diffusion_coefficient). MSD is taken with the simple product between the deviation of the same ions (i=j).

```math
\begin{align}
D_{i}&= \lim_{t\to\infty} \frac{1}{6t} \, \frac{1}{N} \sum_{i}^{N} z_i^2 < [ R_i(t)-R_i(t_0) ]^2 >_t      &\qquad[5]  
\end{align}
```

The trajectory is extracted and stored in 4 array:
- x [F,N]
- y [F,N]
- z [F,N]
- q [1,N] (the order of molecules don't change in the trajectory)

F = number of frames

N = number of molecules in the cell

With this structure is easy to find the deviation:

```math
X=
\begin{bmatrix}
x_0^0 & x_1^0 & \dots \dots &x_{N-1}^0 & x_{N}^0 \\
x_0^1&&&&\\
\vdots&&&&\\
\vdots&&&&\\
x_0^{F-1}&&&&\\
x_0^{F}&&&&\\
\end{bmatrix}

\qquad 

\Delta X=
\begin{bmatrix}
x_0^0-x_0^0 & x_1^0-x_1^0 & \dots \dots &x_{N-1}^0-x_{N-1}^0 & x_{N}^0-x_{N}^0 \\
x_0^1-x_0^0&&&&\\
\vdots&&&&\\
\vdots&&&&\\
x_0^{F-1}x_0^0&&&&\\
x_0^{F}x_0^0&&&&\\
\end{bmatrix}
```

from the deviation is possibile to calculate MSD which is linked to the self-diffusion coefficient $D_i$:

```math
\begin{align}
MSD(t)   = \frac{1}{N} \sum_{i}^{N} z_i^2 < [ R_i(t+\Delta t)-R_i(t) ]^2 >               &\qquad[6] \\ \\
\Delta R_i(t) = (\Delta X_i(t)^2 +\Delta Y_i(t)^2 + \Delta Z_i(t)^2)^{\frac{1}{2}}      \\
\Delta R_i(t)^2 = \Delta X_i(t)^2 +\Delta Y_i(t)^2 + \Delta Z_i(t)^2                    &\qquad[7]\\ \\
MSD(t)=\frac{1}{N}*(\Delta X(t)^2+ \Delta Y(t)^2+\Delta Z(t)^2)    &\qquad[8] \\
D_i=\lim_{t\to\infty}\frac{z_i^2\;MSD(t)}{6t}.    &\qquad[9]
\end{align}
```

This can be usefull in case like ionic solutes in neutral solvents where there is very little correlation of motion between the ions. In case of ionic liquids this can lend to over-estimation of diffusion coefficient because th off-diagonal terms $D_{ij}$ where the MSD is calculated by the product of the deviation of different ions can be far from zero (as the solute exemple). This balance the positive value of self-diffusion $D_i$ giving a more reliable result with the collective $D_{coll}$.

the program can reproduce the same results of Travis and VMD as you can see here

![img](/img/selfdiffusion.png)

#### Inter molecular correlations ( $\boldsymbol{i \neq j}$ )

This is the quantity where ionic liquids start to have deviations from neutral species. In water this quantity is close to zero because solutes dynamics are uncorrelated and is often ignored. In ionic liquids ionic couples can move togheter and the diffusion of single ion can be higher than expected if we don't take into account the correlated movements that happen with the other ions in the simulations.

##### Cation-Cation and Anion-Anion $\boldsymbol{i \neq j}$
The procedure is the same at least up to the deviation matrix $\Delta X$ , $\Delta Y$ , $\Delta Z$ 

```math
\begin{align}
\frac{1}{N} \sum_{i,j}^{N} [ R_{i}(t)-R_{i}(t_0)][ R_{j}(t)-R_{j}(t_0)] &= \frac{1}{N} \sum_{i,j}^{N} \Delta R_i(t) * \Delta R_j(t) \\ \\
\Delta R_i(t)&=\sqrt{ \Delta X_{i}(t)^2 +\Delta Y_{i}(t)^2 + \Delta Z_{i}(t)^2} \\ \\
\Delta R_i(t)*\Delta R_j(t)&=\sqrt{ [\Delta X_{i}(t)^2 +\Delta Y_{i}(t)^2 + \Delta Z_{i}(t)^2]*[\Delta X_{i}(t)^2 +\Delta Y_{i}(t)^2 + \Delta Z_{i}(t)^2]} \\
&=\sqrt{ \Delta I*\Delta J} =\sqrt{ \Delta I}*\sqrt{\Delta J} \\ \\

\frac{1}{N} \sum_{i,j}^{N} [ R_{i}(t)-R_{i}(t_0)][ R_{j}(t)-R_{j}(t_0)] &= \frac{1}{N} \sum_{i,j}^{N} \sqrt{ \Delta I*\Delta J}     \\
\end{align}
```

This time there is the need of a two nested loops over the ions i-j. This make the double sum over products really slow  $O(n^2)$  because there are  ${N\choose 2}=\frac{N!}{2(N-2)!}$ unique combinations to sum over. One full iteration of two loops can take up to 6seconds in the worst case, since this sum has to be repeated over all the sub-sumple of frames (can be quite big 1000-100000) it's unrealistically too slow.
It's possible to rearrange the double sum as a matrix multiplication and take advantage of the processor efficency of numpy. Here you can see the procedure in a simple 3x1 array, in the program this has being extended in 2D:

```math
\begin{align}
A&=[a_0,a_1,a_2] \\
\sum_{i\neq j}^3 a_i * a_j & = a_0 a_1 + a_0 a_2 + a_1 a_2 \\  
&=a_0(a_1+a_2)+a_1(a_2)+a_2*(0)\\ \\
A.cumsum()-A&=[0\,,a_2\,,a_2+a_1]\\ \\
\sum_{i\neq j}^3 a_i*a_j &=A*[A.cumsum()-A]_{inverted}^T= [a_0,a_1,a_2]\begin{bmatrix}a_1+a_2\\a_2\\0\end{bmatrix}
\end{align}
```
This granted $O(n)$ gaining a 100x speedup (now it takes around 0.02seconds for iteration)

##### Cation-Anion $\boldsymbol{i \neq j}$
The same procedure is done with the only exception that here all the product ij ( $N_i*N_j$ ) need to be done compared with the cation-cation where the same ion product need to be excluded ${N\choose 2}$. So it's a bit more simple since the matrix multiplication ij is enough to probe all the deviation products.
