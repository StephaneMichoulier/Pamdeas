**+-----------------------------------                     -----------------------------------+**
**!+----------------------------             **READ ME**             ----------------------------+!**
**+-----------------------------------                     -----------------------------------+**

*Author:* Stéphane Michoulier
*email:* <stephane.michoulier@univ-lyon1.fr>

*Eden:* version 1.2

**DESCRIPTION**:

Eden is a 1D code primarily wrote to study the evolution of grains porosity within protoplanetary disc considering different physical process such as growth, fragmentation, drift, pressure bump, etc.  
This code allow you to follow the evolution of a given number of particles from a set of initials conditions chosen by the user in a stationary non self-gravitating gas disk.  
This code is mainly developed to understand the behaviour of porous grains constrain by numerous physical processes and to test the algorithms for a future implementation in 3D Hydrodynamics code.  

FEW WORDS ON INPUT FILE OPTION (Respect the blank space, or expect unforeseen consequences).  
2 models can be used to study the evolution of the porosity. 
One is based on the increase of mass (dm/dt) and use A. Garcia's algorithms, the other one is based on the increase of size (ds/dt) and use my algorithms.  
The dmdt model support the bounce process, but not the dsdt model.  
The degeneracy due to the filling factor in the dsdt model is taken into account.  

You are allowed to create artificially a pressure bump, by adding a gaussian profile to the gas surface density. To create the desired pressure bump, the easiest way is to set ngrains to 0, profile to 1 and plot Sigma vs R, Pg vs R and dustfrac vs R in order to visualise your bump. Then, proceed by iteration.

Be aware that using a pressure bump with radial drift and fragmentation is way slower than any other setup due to small time step.

For porous grains, physical properties of the material have to be changed if the intrinsic density is changed.  

**LAUNCHING THE PROGRAM**

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
> - output_header.txt gives you the corresponding columns in output files.
> - disk_profiles.out contains different profiles (surface density, pressure, etc).
> - disk_profiles_header.txt gives you the corresponding columns in the diskprofiles.out files.
> - disrupt_param_<model>.out contains parameters of disrupted grains.
> - disrupt_param_header.txt gives you the corresponding columns in the disrupt_param.out  files.

**HOW THE PROGRAM IS BUILT ?**

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

