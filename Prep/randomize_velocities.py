#!/usr/bin/env python

from __future__ import print_function, division

import argparse, sys
import numpy, numpy.random
import scipy, scipy.constants

description = \
"""
A script to rewrite the velocities in a gromacs .gro file so that a simulation is guaranteed
not to have the Boltzmann distribution built in.
"""

"""
From the GROMACS manual:
Lines contain the following information (top to bottom):

        title string (free format string, optional time in ps after t=)
        number of atoms (free format integer)
        one line for each atom (fixed format, see below)
        box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), the last 6 values may be omitted (they will be set to zero). GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.

"""
GRO_ATOM_FORMAT = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
# 20 intro chars, 24 position chars, 24 velocity chars

def vprint(msg):
    global args
    if args.verbose:
        print(msg,file=sys.stderr if args.output is sys.stdout else sys.stdout)

parser = argparse.ArgumentParser(description=description)
parser.add_argument('-i', '--input', type=argparse.FileType('rt'), default=sys.stdin,
                    help='Input file. (Default: standard input.)')
parser.add_argument('-o', '--output', type=argparse.FileType('wt'), default=sys.stdout,
                    help='Output file. (Default: standard output.)')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Report extra information')
parser.add_argument('temperature', type=float,
                    help='Temperature (in K)')
parser.add_argument('mass', type=float,
                    help='Mass (in g/mol)')
            
args = parser.parse_args()

# E_K = (3/2)kBT = (1/2)mv^2 so v = sqrt(3kBT/m)
m = args.mass * scipy.constants.atomic_mass
vprint('Mass = {} g/mol, {} kg/particle'.format(args.mass,m))
vmag = numpy.sqrt(3*scipy.constants.Boltzmann*args.temperature/m)
vprint('Magnitude of velocity: {} m/s'.format(vmag))
vmag /= 1000.0 # convert from m/s to nm/ps (or km/s)

# Copy the title string
args.output.write(args.input.readline())

# Copy and interpret the number of atoms
n_atoms_line = args.input.readline()
args.output.write(n_atoms_line)
n_atoms = int(n_atoms_line)

for iatom in range(n_atoms):
    atom_line = args.input.readline()
    line_header = atom_line[:44]
    args.output.write(line_header)
    
    # Draw random numbers for each of vx, vy, or vz
    vs = numpy.random.uniform(-1.0, 1.0, size=(3,))
    
    # Normalize
    vs /= (vs**2).sum()**(0.5)
    
    # Scale to desired magnitude
    vs *= vmag
    
    # Write new velocities
    args.output.write('%8.4f%8.4f%8.4f\n' % tuple(vs))
    
    
# Copy the box vectors
args.output.write(args.input.readline())


    
