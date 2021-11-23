""" usage: train.py [-h] 

"""

import argparse
import os

def main():

    # Initiate argument parser
    parser = argparse.ArgumentParser(description="Performs eGPE calculation.\n \
                                                 Units of length given in a_11, \
                                                 Units of energy \hbar^2 / (2 m_1)    ",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-a11', help='Value of a_11', type=float, default=None
    )
    parser.add_argument(
        '-a22', help='Value of a_22', type=float, default=None
    )
    parser.add_argument(
        '-a12', help='Value of a_12', type=float, default=None
    )
    
    args = parser.parse_args()
    
    
    from sim_params import SimulationSetup
    
    simulation = SimulationSetup(a11=args.a11,
                                 a22=args.a22,
                                 a12=args.a12)
    print(simulation)
    print(args)
    


if __name__ == '__main__':
    main()