# About


Solves the coupled (nonlinear) Schrödinger equation either in real or imaginary time. The implemented algorithm is outlined here: https://www.sciencedirect.com/science/article/abs/pii/S0010465509001131

This code can be used to study static and dynamic properties of two-component Bose-Bose quantum droplets.
https://arxiv.org/abs/2104.09102


# Theory

Simulates the two-component $\psi_1$ and $\psi_2$ system evolving according to the eGPE equations

$$
\begin{aligned}
i \hbar \frac{\partial \psi_{1}}{\partial t}=&\left(-\frac{\hbar^{2}}{2 m} \nabla^{2}+\frac{4 \pi \hbar^{2} a_{11}}{m} \rho_{1}+\frac{4 \pi \hbar^{2} a_{12}}{m} \rho_{2}+\frac{128 \sqrt{\pi} \hbar^{2} a_{11}}{3 m}\left(a_{11} \rho_{1}+a_{22} \rho_{2}\right)^{3 / 2}\right) \psi_{1} \\
i \hbar \frac{\partial \psi_{2}}{\partial t}=&\left(-\frac{\hbar^{2}}{2 m} \nabla^{2}+\frac{4 \pi \hbar^{2} a_{22}}{m} \rho_{2}+\frac{4 \pi \hbar^{2} a_{12}}{m} \rho_{1}+\frac{128 \sqrt{\pi} \hbar^{2} a_{22}}{3 m}\left(a_{11} \rho_{1}+a_{22} \rho_{2}\right)^{3 / 2}\right) \psi_{2}
\end{aligned}
$$


# Requirements


```
pip install requirements.txt
```

# Documentation


Run the following code for help.
```
python egpe_solver.py -h 
```


# Files

[egpe_solver.py](egpe_solver.py) solves the coupled equations in imaginary or real time 

[plot_utils/2D_plot.py](plot_utils/2D_plot.py) makes 2D plots of integrated density in form of the .png files . 

[plot_utils/make_movie.sh](plot_utils/make_movie.sh) generates a movie out of a collection of .png files. 



# Example use cases

## Static properties of droplets within the MF+LHY theory