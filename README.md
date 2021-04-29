### **+-----------------------------------                       -----------------------------------+**  
### **!+----------------------------              **READ ME**               ----------------------------+!**  
### **+-----------------------------------                       -----------------------------------+**  

*Author:* Stéphane Michoulier  
*Email:* <stephane.michoulier@gmail.com>

*Eden:* version 1.3

#### **DESCRIPTION**:

Eden is a 1D code primarily wrote to study the evolution of grains porosity within protoplanetary disc considering different physical process such as growth, fragmentation, drift, pressure bump, snowline (Vericel & Gonzalez 2020) etc.  
This code allow you to follow the evolution of a given number of particles from a set of initials conditions chosen by the user in a stationary non self-gravitating gas disc and was mainly developed to understand the behaviour of porous grains, but also to test algorithms for future implementation in 3D Hydrodynamics code.  

FEW WORDS ON INPUT FILE OPTION (Respect the blank space, or expect unforeseen consequences).  
2 models can be used to study the evolution of the porosity. 
One is based on the increase of mass (dm/dt) and use A. Garcia's reformulated algorithms (Garcia 2018),(Garcia & Gonzalez 2020), the other one is based on the increase of size (ds/dt) and use my algorithms.  
The dmdt model support the bounce process, but not the dsdt model.  
The degeneracy due to the filling factor in the dsdt model is taken into account.  

You are allowed to create artificially a pressure bump, by adding a gaussian profile to the gas surface density. To create the desired pressure bump, the easiest way is to set ngrains to 0, profile to 1 and plot Sigma vs R, Pg vs R and dustfrac vs R in order to visualise your bump. Then, proceed by iteration.

Be aware that using a pressure bump with radial drift and fragmentation is way slower than any other setup due to small time step.

For porous grains, physical properties of the material have to be changed if the intrinsic density is changed.  

#### **LAUNCHING THE PROGRAM**

Set your SYSTEM environment to be able to build programs in C++: export SYSTEM=g++  
To use all functionality of the program, be sure to  have the latest C++ version installed.  

If the program is not yet built, follow the instruction below:
> - Open a new terminal
> - Use the command cd to reach "eden" directory: cd 'Yourdirectory'/eden
> - Once in the directory, build the code with: make
> - Then, source file will be compiled in "source" folder and the program will be built
> - The simulation can then be launch. To do this, go inside the "build" folder where the executable is: cd build
> - First, the input file need to be created. In order to write it, write: ./eden setup
> - Then, modify the input file to make your own simulation. Don't forget to save your setup
> - Finally, run Eden: ./eden  

Once the simulation is finished, multiple output files are written into the output folder depending on the option:
> - Initials_Conditions_*Rref*.txt summarises your setup
> - output_*setup name*.out: one output file for each particle. The file name depends on the model used (dmdt/dsdt and porous/compact grains) 
> - disk_profiles.out contains different profiles (surface density, pressure, etc).
> - disrupt_param_*model*.out contains parameters of disrupted grains.

A script is available to plot directly some useful quantities. To do so, write in your terminal (in the build folder):
> - ./ploteden help to display how to plot things
All python scripts are in the python folder, you can modify them as you wish to suit your data.

#### **HOW THE PROGRAM IS BUILT ?**

The program contains multiple files in the source folder "source" and the header folder "header":

1. Source files (.CPP and .H)

- a. main.cpp

	main.cpp contains all variables declarations and time loop

- b. disc.cpp & disc.h
	
	Contain general functions describing gas disks such as gas density, pressure, etc.
    The full list is given in the .h file

- c. dust.cpp & dust.h
	
	Contain general functions describing properties of dust grains as kinetic energy, relative velocity or Stokes number.
    The full list is given in the .h file

- d. porosity.cpp & porosity.h

    Contain all the algorithms driving the evolution of the filling factor depending on the mass/size and the different processes: growth, frag, compression by collision/gravity etc 

- e. evol.cpp & evol.h

    Contain all the algorithms linked to time stepping: growth rate, drift velocity, adaptative time step etc.

- f. readwrite.cpp & readwrite.h
   
	Contain reading and writing functions

- g. constantandconversion.h

	Contains the declaration of physical constants and conversion functions between some useful  units

#### **Reference**
Arnaud Vericel and Jean-François Gonzalez. Self-induced dust traps around snow lines in pro- toplanetary discs. MNRAS, 492(1) :210–222, February 2020. doi : 10.1093/mnras/stz3444.

Anthony Garcia. Evolution of grain porosity during growth : a solution to planetary for- mation barriers? Phd thesis, Université de Lyon, September 2018. URL https://tel. archives-ouvertes.fr/tel-01977317.

Anthony J. L. Garcia and Jean-François Gonzalez. Evolution of porous dust grains in pro- toplanetary discs - I. Growing grains. MNRAS, 493(2) :1788–1800, April 2020. doi : 10.1093/mnras/staa382.

J. D. Hunter. Matplotlib : A 2d graphics environment. Computing in Science & Engineering, 9(3) :90–95, 2007. doi : 10.1109/MCSE.2007.55.
