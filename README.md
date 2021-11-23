# About


Solves the coupled (nonlinear) Schrödinger equation either in real or imaginary time. The implemented algorithm is outlined here: https://www.sciencedirect.com/science/article/abs/pii/S0010465509001131

This code can be used to study static and dynamic properties of two-component Bose-Bose quantum droplets.
https://arxiv.org/abs/2104.09102

[egpe_solver.py](egpe_solver.py) solves the coupled equations in imaginary or real time 

[plot_utils/2D_plot.py](plot_utils/2D_plot.py) makes 2D plots of integrated density in form of the .png files . 

[plot_utils/make_movie.sh](plot_utils/make_movie.sh) generates a movie out of a collection of .png files. 


# Requirements

Just `numpy`. You can install it with
```
sudo apt-get install python-numpy
```

# Documentation


Run the following code for help.
```
python -h
```

# Example use cases

## Static properties of droplets within the MF+LHY theory