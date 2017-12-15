# Molecular Dynamics Simulations for Physical Chemistry Classes

## Introduction

This distribution contains the following directories. Simulations are
"ready to run" according to the instructions contained in this file or
the documents in the `Documents/` directory. Some software packages
(notably the GROMACS molecular dynamics engine and the VMD
visualization program) are required. Information on how to install
this software is provided in the "Getting Started" section below.

* `Activity1_Ideal_Real_Gases/` — MD simulation to illustrate the
  differences between ideal and real gases. Also introduces radial
  distribution functions and the connections between molecular size
  and interaction energy and the parameters of the van der Waals
  equation. This directory (and other activity directories) also
  contain example output and analysis. 
* `Activity2_Meaning_Of_Beta/` — An inquiry-based activity to uncover
  the meaning of the constant (beta) in the Boltzmann distribution
  equation from the speed distribution of a real gas.
* `Activity3_Phase_Change/` — MD simulation of the compression of
  argon from gas to liquid, illustrating the essential features of
  phase transitions, the structure of the liquid state, and the
  calculation of energy and enthalpy of vaporization.
* `Documents/` — Student notes for Activities 1 and 3, as LaTeX
  sources, ODT, and DOCX files. (Activity 2 is a guided inquiry, and
  so only has instructors' notes.)
* `Prep/` — Files and scripts necessary to prepare fully-independent
  simulations for different students or situations. Python and the
  numpy and scipy libraries are required for this; the version that
  ships with recent version of MacOS suffices.
  
Each simulation directory contains an example analysis notebook.
These analysis notebooks (for example `Activity1_Analysis.ipynb` in
`Activity1_Ideal_Real_Gases/`) can be viewed through GitHub even
without installing or running Python.
  
The "Getting Started" section below provides instructions for how to
obtain the necessary software and use the files provided here.

## Getting Started (MacOS)

The following instructions are for MacOS. Instructions for other
operating systems (especially Windows) would be a welcome contribution
from the community.

### Installing the necessary software
The most straightforward way to get things going is to install
[Homebrew](https://brew.sh/). Follow the (remarkably short)
installation instructions as given on the Homebrew web page.  Then, in
the terminal, run the following:
```
brew install gromacs
```

If you wish to run the Python-based analysis notebooks, also run
```
brew tap homebrew/science
brew install python3 numpy scipy ipython jupyter
brew install matplotlib --with-python3
```
or install the Anaconda Python Distribution from 
<https://www.anaconda.com/download/>

Finally, download the VMD viewer from
<http://www.ks.uiuc.edu/Research/vmd/>. This requires a free
registration.

### Downloading the files necessary to run simulations and analysis

To download the files necessary to perform the activities, clone this
Git repository. Also in the terminal:

`git clone https://github.com/mczwier/md_in_pchem.git`

This creates a directory named `md_in_pchem`. This is the starting
point for all of the simulations.


## Running the simulations

Data files and scripts have been provided to run the simulations and
extract the data necessary for further analysis. Only the software
described above (GROMACS and VMD) is required. The scripts provided
assume that each activity is run in sequence. This is not necessarily
required from a pedagogical standpoint, but if you choose to skip an
activity, you will nonetheless need to run the simulations for each
prior activity before moving on. From the `md_in_pchem` directory
created above, do the following to run the simulations for a given activity.

### The first activity (ideal and real gases)
To run the simulations for the first activity, run the following
commands from the terminal in the `md_in_pchem` directory:

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
return to the `md_in_pchem` directory. **Alternatively**, see
`Act01.pdf`, `Act01.odt`, or `Act01.docx` in the `Documents/`
directory for an equivalent step-by-step, hands-on approach.

### The second activity (the meaning of thermodynamic β)

The second activity — determining the physical meaning of β in the
Boltzmann distribution exp(–βE) — is an analysis-only activity, and
requires the trajectories (simulations) from the first activity. Thus,
if you want to use this activity, first perform the simulations for
the first activity (as described above). To prepare data for analysis,
execute the following, again from the `md_in_pchem` directory:

```
cd Activity2_Meaning_Of_Beta
./run.sh
```

This creates the files `vels_100.xvg`, `vels_200.xvg`, and
`vels_300.xvg`, which contain the x velocities for each of the 100
argon atoms in the simulation as a function of simulation time. These
files can then be loaded into your favorite analysis tool to construct
histograms. For those wishing to use Excel but not wishing to use the
Analysis ToolPak, a Python script (`make_vel_histogram.py`) has been
provided to read atomic velocities and produce histograms in either a
plain text table or CSV format. Python and Numpy are required, and the
version of Python shipped with recent versions of MacOS suffices.

### The third activity (phase changes)

The third activity — the nature of phase changes and calculation of
the heat of vaporization — does not require any preparatory
simulations to be performed. To run this simulation, do the following
from the `md_in_pchem` directory:

```
cd Activity3_Phase_Change/
./run.sh
```
**Alternatively**, see `Act03.pdf`, `Act03.odt`, or `Act03.docx` in
the `Documents` directory for step-by-step instructions for both
running and analyzing the simulation.

## Where to go next

### Running the Jupyter analysis notebooks

If you are interested in exploring the Jupyter notebooks for analysis,
install Python and the Jupyter notebook system either using Homebrew
(see "Getting Started" above) or download the Anaconda Python
Distribution from <https://www.anaconda.com/download/>. Choose the
Python 3 version (not the Python 2.7 version).

Once Python is installed, run `jupyter notebook` from the Terminal in
the `md_in_pchem` directory. This will start a web browser and display
a listing of the contents of the `md_in_pchem` directory, which should
include directories for each of the activities. Click on the directory
for the activity whose analysis you want to run
(e.g. `Activity1_Ideal_Real_Gases`) and click on the Jupyter notebook
(extension `.ipynb`, e.g., `Activity1_Analysis.ipynb`) to start the
analysis notebook. These notebooks run from top to bottom (like a
Mathematica notebook, for instance); pressing Shift-Enter executes a
cell and moves to the next one. To run the example analysis, just keep
pressing Shift-Enter.

### (Re)equilibrating the simulations

Pre-equilibrated simulations (at the correct temperature and density
and every atom having an appropriate velocity) have been
provided. However, to illustrate variability from simulation to
simulation (or equivalently to ensure that different students have
different data to work with), it is necessary to re-equilibrate
(generate new velocities). This requires Python and the Numpy library
(provided, for example, by the Anaconda distribution described above,
or by MacOS). To re-equilibrate, execute the following **prior** to
running simulations for each of the activities. From the `md_in_pchem`
directory, run the following:

```
cd Prep
./run_equil.sh
```

