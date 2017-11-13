#!/bin/bash

# It is assumed that the GROMACS binaries are on the $PATH

# This disables GROMACS auto-backup for the duration of this script
export GMX_MAXBACKUP=-1

# Clean up previous runs

set -e

rm -fv \#*
rm -fv *.{gro,xtc,trr,log,edr,cpt,xvg,tpr}

cp ../Prep/equilPT_100.gro initial.gro
cp ../Prep/argon.top .

# Equilibrate to the initial conditions we need for a nice smooth
# compression and cooling.
# The instructor may wish to do this for the students, but it's not
# necessary either way, as long as it gets done.

gmx grompp -f equil.mdp -c initial.gro -p argon.top -o equil.tpr
gmx mdrun -nt 1 -nb cpu -s equil.tpr -c equil.gro -e equil.edr \
    -g equil.log -x equil.xtc -cpo equil.cpt

# Then run dynamics for compression and cooling
gmx grompp -f md.mdp -c equil.gro -t equil.cpt -p argon.top -o md.tpr
gmx mdrun -nt 1 -nb cpu -s md.tpr -c md.gro -e md.edr -g md.log -x md.xtc

# Extract all the information we need
echo  7 | gmx energy -f md.edr -s md.tpr -o temperature.xvg
echo 13 | gmx energy -f md.edr -s md.tpr -o volume.xvg
echo 14 | gmx energy -f md.edr -s md.tpr -o density.xvg
echo  4 | gmx energy -f md.edr -s md.tpr -o pot_energy.xvg
echo  5 | gmx energy -f md.edr -s md.tpr -o kin_energy.xvg
echo  6 | gmx energy -f md.edr -s md.tpr -o tot_energy.xvg

rm -fv state.cpt state_prev.cpt mdout.mdp
