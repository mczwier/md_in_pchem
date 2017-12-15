#!/usr/bin/env python

from __future__ import division, print_function

import argparse, sys, os
import numpy as np

description="""\
A script to create histograms from a GROMACS-produced XVG file that can
then be read into Excel, which lacks a good histogram function. Some
might find this easier than using the Analysis Toolkit. By default,
a space separated file is written, but a CSV file may be specified.
"""

def read_xvg(filename):
    '''A little auxiliary function to read a GROMACS-produced XVG
    file, automatially skipping over the default xmgr formatting.'''
    
    skiplines = 0
    with open(filename, 'rt') as xvgfile:
        for line in xvgfile:
            if line[0] in ('#','@'):
                skiplines += 1
                continue
            else:
                break
        return np.loadtxt(xvgfile)

def vprint(msg):
    global args
    if args.verbose:
        print(msg,file=sys.stderr if args.output is sys.stdout else sys.stdout)

parser = argparse.ArgumentParser(description=description)
parser.add_argument('input', 
                    help='Input file. (Required)')
parser.add_argument('-o', '--output', 
                    help='Output file. (Default: based on input filename)')
parser.add_argument('-n', '--nbins', type=int, default=50,
                    help='Number of bins for the histogram. (Default: %(default)s)')
parser.add_argument('-c', '--csv', action='store_true',
                    help='Write a CSV file instead of a space-delimited file.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Report extra information')

args = parser.parse_args()

if args.output == '-':
    args.output = sys.stdout
elif not args.output:
    (basename, junk) = os.path.splitext(args.input)
    if args.csv:
        args.output = open(basename + '.csv', 'wb')
    else:
        args.output = open(basename + '.dat', 'wt')
else:
    if args.csv:
        args.output = open(args.output, 'wb')
    else:
        args.output = open(args.output, 'wt')
vprint('Output file: {}'.format(args.output.name))

veldata = read_xvg(args.input)
vprint('Velocity data has shape {} and type {}'.format(veldata.shape, veldata.dtype))
# Convert from nm/ps (km/s) to m/s
vx = veldata[:,1:]*1000
(vx_f, vx_bounds) = np.histogram(vx, args.nbins, normed=True)
vx_centers = (vx_bounds[:-1] + vx_bounds[1:])/2.0    
del veldata

if args.csv:
    import csv
    writer = csv.writer(args.output)
    writer.writerow(['Bin Center', 'Probability Density'])
    for (v,p) in zip(vx_centers, vx_f):
        writer.writerow([v,p])
else:
    args.output.write('{:24s}  {:24s}\n'.format('Bin Center', 'Probability Density'))
    for (v,p) in zip(vx_centers, vx_f):
        args.output.write('{:24s}  {:24s}\n'.format(repr(v), repr(p)))
