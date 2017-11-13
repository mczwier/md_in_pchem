#!/bin/bash

# It is assumed that the GROMACS binaries are on the $PATH

# This disables GROMACS auto-backup for the duration of this script
export GMX_MAXBACKUP=-1

# Clean up previous runs

set -e

rm -fv \#*
rm -fv *.{gro,xtc,trr,log,edr,cpt,xvg,tpr}

# Link in previously-generated trajectory data
# and generate velocity histograms

for TEMP in 100 200 300; do
    ln -s ../Activity1_Ideal_Real_Gases/real_${TEMP}.{tpr,trr} .
    echo 0 | gmx traj -s real_$TEMP.tpr -f real_$TEMP.trr -ov vels_$TEMP.xvg -x -noy -noz    
done

rm -fv state.cpt state_prev.cpt mdout.mdp
