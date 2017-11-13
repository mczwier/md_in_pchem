#!/bin/bash

# It is assumed that the GROMACS binaries are on the $PATH

# This disables GROMACS auto-backup for the duration of this script
export GMX_MAXBACKUP=-1

# This causes this script to abort if any command fails
set -e

# Remove any copies of previous files
rm -fv *.{trr,xtc,tpr,edr,log}
rm -fv equilPT_100.gro equilT_200.gro equilT_300.gro
rm -fv uniform_{100,200,300}.gro

# The first equilibration involves both a thermostat and barostat, and
# ensures that the simulation volume is consistent with the target
# temperature and pressure. Later equilibrations for real gas
# simulations DO NOT USE velocities from this run, so the Boltzmann
# distribution is never actually "built in" to the simulation.

gmx grompp -f equilPT_100.mdp -c argon100.gro -p argon.top -o equilPT_100.tpr
gmx mdrun -nt 1 -nb cpu -s equilPT_100.tpr -c equilPT_100.gro -e equilPT_100.edr \
          -g equilPT_100.log -x equilPT_100.xtc

# The next equilibrations adjust the temperature from 100 to 200 or
# 300 K for use in the ideal gas simulations.

gmx grompp -f equilT_200.mdp -c equilPT_100.gro -p argon.top -o equilT_200.tpr
gmx mdrun -nt 1 -nb cpu -s equilT_200.tpr -c equilT_200.gro -e equilT_200.edr \
          -g equilT_200.log -x equilT_200.xtc

gmx grompp -f equilT_300.mdp -c equilT_200.gro -p argon.top -o equilT_300.tpr
gmx mdrun -nt 1 -nb cpu -s equilT_300.tpr -c equilT_300.gro -e equilT_300.edr \
          -g equilT_300.log -x equilT_300.xtc


# Prepare initial velocities with random directions but velocities equal
# to the expected (mean) velocity at the given temperature
./randomize_velocities.py -i equilPT_100.gro -o uniform_100.gro 100.0 39.948
./randomize_velocities.py -i equilT_200.gro  -o uniform_200.gro 200.0 39.948
./randomize_velocities.py -i equilT_300.gro  -o uniform_300.gro 300.0 39.948

# Now perform a 5ns equilibration including COM velocity subtraction
# to allow collisions to restore the velocity distribution to something
# realistic
gmx grompp -f equilE.mdp -c uniform_100.gro -p argon.top -o equilE_UniVel_100.tpr -maxwarn 1
gmx mdrun -nt 1 -nb cpu -s equilE_UniVel_100.tpr -c equilE_UniVel_100.gro \
    -e equilE_UniVel_100.edr -g equilE_UniVel_100.log -x equilE_UniVel_100.xtc

gmx grompp -f equilE.mdp -c uniform_200.gro -p argon.top -o equilE_UniVel_200.tpr -maxwarn 1
gmx mdrun -nt 1 -nb cpu -s equilE_UniVel_200.tpr -c equilE_UniVel_200.gro \
    -e equilE_UniVel_200.edr -g equilE_UniVel_200.log -x equilE_UniVel_200.xtc

gmx grompp -f equilE.mdp -c uniform_300.gro -p argon.top -o equilE_UniVel_300.tpr -maxwarn 1
gmx mdrun -nt 1 -nb cpu -s equilE_UniVel_300.tpr -c equilE_UniVel_300.gro \
    -e equilE_UniVel_300.edr -g equilE_UniVel_300.log -x equilE_UniVel_300.xtc

rm -f state.cpt state_prev.cpt mdout.mdp
