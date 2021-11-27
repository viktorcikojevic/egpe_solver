""" usage: train.py [-h] 

"""

import argparse
import os

def main():

    # Initiate argument parser
    parser = argparse.ArgumentParser(description=
    "Performs eGPE calculation.\n\
    Units of length given in a_11\n\
    Units of energy \hbar^2 / (m_1 a_{11}^2)\n\
    Units of time m_1 a_{11}^2 / \hbar",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a11', help='Value of a_11 (Bohr radius)', type=float,  required=True)
    parser.add_argument('-a22', help='Value of a_22 (Bohr radius)', type=float,  required=True)
    parser.add_argument('-a12', help='Value of a_12 (Bohr radius)', type=float,  required=True)
    parser.add_argument('-m1', help='Mass of component 1 (atomic mass unit)', type=float, default=None)
    parser.add_argument('-m2', help='Mass of component 2 (atomic mass unit)', type=float,  required=True)
    parser.add_argument('-nxyz',  help='Number of grid points in xyz direction', type=int, required=True)
    parser.add_argument('--grid-size',nargs='+', help='Size of box (xyz)', type=int, required=True)
    parser.add_argument('--absorbing-buffer',nargs='+', help='what percentage of x, y and z are absorbing buffers',
                        type=float, default=0., required=False)
    
    parser.add_argument('-dtau', help='Imaginary time-step (in miliseconds)',type=float, required=False)
    parser.add_argument('--entol', help='Energy tolerance to end the simulation (defaults to 1.E-08)',type=float, default=1.E-08, required=False)
    parser.add_argument('-dt', help='Real time-step (in miliseconds)',type=float, required=False)
    parser.add_argument('-sim-time', help='Simulation time (in miliseconds)',type=float, required=False)
    parser.add_argument('--printfreq', help='Frequency of output',type=int, default=10, required=False)
   
   
    args = parser.parse_args()
    print(args)
    
    from setup import SimulationSetup
    
    simulation = SimulationSetup(a11=args.a11,
                                 a22=args.a22,
                                 a12=args.a12,
                                 m1=args.m1,
                                 m2=args.m2,
                                 nxyz=args.nxyz,
                                 L_grid=args.grid_size,
                                 L_grid_nobuffer_percentage=[1-gs for gs in args.absorbing_buffer],
                                 dtau=args.dtau,
                                 dt=args.dt,
                                 entol=args.entol,
                                 printfreq=args.printfreq
                                 )
    print(simulation)
    print(args)
    


if __name__ == '__main__':
    main()