# MD for P-Chem

## Introduction

This distribution contains the following directories:

* `Documents/` — Student notes for each activity, as LaTeX sources, ODT, and DOCX files
* `Examples/` — Example output and analysis for each activity
* `Activity1_Ideal_Real_Gases/` — MD simulation to illustrate the
  differences between ideal and real gases. Also introduces radial
  distribution functions and the connection between molecular
  interaction energy and the parameters of the van der Waals equation.
* `Activity2_Meaning_Of_Beta/` — An inquiry-based activity to uncover
  the meaning of the constant (beta) in the Boltzmann distribution
  equation from the speed distribution of a real gas.
* `Activity3_Phase_Change/` — MD simulation of the compression of
  argon from gas to liquid.
  
For those less familiar with molecular simulations, the "Getting
started" section below provides instructions for how to obtain the
necessary software and use the files provided here.

## Getting started

The following instructions are for Mac OS. Instructions for other
operating systems (especially Windows) would be a welcome contribution
from the community.  The most straightforward way to get things going
is to install MacPorts (https://www.macports.org/ ). If you do not
have MacPorts installed already, follow the installation instructions
on the MacPorts website; use the Package (.pkg) installer for maximum
convenience. From there, start a terminal (the Terminal application
located under `/Applications/Utilities`).  In the terminal, use
MacPorts to install the GROMACS MD engine:

`sudo port install gromacs`

To download the files necessary to perform the activities, clone this
Git repository. Also in the terminal:

`git clone https://github.com/mczwier/md_in_pchem.git`

This creates a directory named `md_in_pchem`. This is the starting
point for all of the simulations.

Finally, download the VMD viewer from
http://www.ks.uiuc.edu/Research/vmd/ . This requires a free
registration.

## Running the simulations

Data files and scripts have been provided to run the simulations and
extract the data necessary for further analysis. Only the software
described above (GROMACS and VMD) is required. The scripts provided
that each activity is run in sequence. This is not necessarily
required from a pedagogical standpoint, but if you choose to skip an
activity, you will nonetheless need to run the simulations for each
prior activity before moving on. From the `md_in_pchem`
directory created above, do the following:

**For the first activity (ideal and real gases)**: To run the
simulations for the first activity, run the following commands from
the terminal in the `md_in_pchem` directory:

```
cd Activity1_Ideal_Real_Gases 
./run.sh
```

This will take several minutes. For example, on a 2015-era MacBook
Pro, this takes about five minutes. This runs the simulations and
extracts energies and radial distribution functions, leaving them in
the current directory (as text files ending with `.xvg`, which can be
loaded into Excel or another analysis tool). A Finder window can be
opened to interact with these files by running `open .` (note trailing
period) in the terminal. When done with this activity, run `cd ..` to
return to the `md_in_pchem` directory.




## Where to go next

### Running the Jupyter analysis notebooks

If you are interested in exploring the Jupyter notebooks for analysis,
download the Anaconda Python Distribution from
(https://www.anaconda.com/download/). This is not necessary to run the
simulations or to analyze them, but may be useful. Choose the Python 3
version (not the Python 2.7 version). 

To run one of the analysis notebooks, run `jupyter notebook` from the
Terminal in the `md_in_pchem` directory. This will start a web browser
and display a listing of the contents of the `md_in_pchem` directory,
which should include directories for each of the activities. Click on
the directory for the activity whose analysis you want to run
(e.g. `Activity1_Ideal_Real_Gases`) and click on the Jupyter notebook
(extension `.ipynb`, e.g., `Activity1_Analysis.ipynb`) to start the
analysis notebook. These notebooks run from top to bottom (like a
Mathematica notebook, for instance); pressing Shift-Enter executes a
cell and moves to the next one. To run the example analysis, just keep
pressing Shift-Enter.

### Equilibrating the simulations

Pre-equilibrated simulations (at the correct temperature and density
and every atom having an appropriate velocity) have been
provided. However, to illustrate variability from simulation to
simulation, it is necessary to re-equilibrate (generate new
velocities). This requires Python and the Numpy library (provided, for
example, by the Anaconda distribution described above). To
re-equilibrate, execute the following **prior** to running simulations
for each of the activities. From the `md_in_pchem` directory, run the
following: 

```
cd Prep
./run_equil.sh
```

