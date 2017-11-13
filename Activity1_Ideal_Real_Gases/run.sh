#!/bin/bash

# It is assumed that the GROMACS binaries are on the $PATH

# This disables GROMACS auto-backup for the duration of this script
export GMX_MAXBACKUP=-1

# Clean up previous runs

set -e

rm -fv \#*
rm -fv *.{gro,xtc,trr,log,edr,cpt,xvg}

# Copy in the equilibrated configurations with randomized and
# re-equilibrated velocities.

cp -v ../Prep/argon.top ../Prep/ideal.top .

for TEMP in 100 200 300; do 
    cp -v ../Prep/equilE_UniVel_$TEMP.gro equil_$TEMP.gro
    
    gmx grompp -f md_ideal.mdp -c equil_$TEMP.gro -p ideal.top -o ideal_$TEMP.tpr -maxwarn 1
    gmx mdrun -nt 1 -nb cpu -s ideal_$TEMP.tpr -c ideal_$TEMP.gro \
              -e ideal_$TEMP.edr  -g ideal_$TEMP.log -o ideal_$TEMP.trr
    gmx rdf -s ideal_$TEMP.tpr -f ideal_$TEMP.trr -ref 2 -sel 2 -o rdf_${TEMP}_ideal.xvg
    echo 6 5 4 1 7 9 | gmx energy -f ideal_$TEMP.edr -o energy_${TEMP}_ideal.xvg >& energy_${TEMP}_ideal.log

    gmx grompp -f md_real.mdp -p argon.top -c equil_$TEMP.gro -o real_$TEMP.tpr -maxwarn 1
    gmx mdrun -nt 1 -nb cpu -s real_$TEMP.tpr -c real_$TEMP.gro \
              -e real_$TEMP.edr -g real_$TEMP.log -o real_$TEMP.trr
    gmx rdf -s real_$TEMP.tpr -f real_$TEMP.trr -ref 2 -sel 2 -o rdf_${TEMP}_real.xvg
    echo 6 5 4 1 7 9 | gmx energy -f real_$TEMP.edr -o energy_${TEMP}_real.xvg >& energy_${TEMP}_real.log
done

rm -fv state.cpt state_prev.cpt mdout.mdp
