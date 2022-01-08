#Installation istruction
In order to make the code portable a docker virtual machine has been used.

1. So the first thing to do is checking if docker is installed. 
2. Clone the repo.
2. Run the **run.sh** file.

Running **run.sh** file, these things are done:

1. The Dockerfile is executed and so all the project is compiled, run and saved in a docker image. 
2. The solutions folders of all the test-cases are exported from the docker container to the local computer in order to give the chance to open the **.png** image files that show the test-case graph solutions. This step is done thanks to a cycle that loops over all the case folders. 
3. The docker image can be access through a suitable container typing:
         docker container run -it pacs_project /bin/bash
 
In order to see how the Dockerfile has been run check the appendix at the end of the report. 

#Content of the git repo 
The Project folder contains all the test-cases folders.

In particular:

##case0_example
This test case implements a simple differential problem as example for solving a PDE.
 
##case1_darcy
This test case implements a convergence analysis for the Darcy system code.
 
##case2_darcy
This test case implements and solves the Darcy system with the permability parameter K which is not continous in space.
 
##case3_transport
This test case implements and solves the transport problem of a single tracer giving the chance to choose between an implicit or explicit scheme.
 
##case4_linear_decay
This test case implements and solves the transport and linear decay problem. Tranpsort is treated in an implicit way, instead the linear decay is treated explicitely.
 
##case5_2reagents
This test case implements and solves two coupled transport and reaction PDE. The coupling is due to the reaction term. In particular, it treats the $CaSiO_3$ dissolution with the consequent generation of $Ca^{2+}$ ions.
 
##case6_all
This test case implements and solves the overall problem with all the reagents implicated in the $CO_2$ storage problem. 

Then we have all the folder containing the material needed by these test-cases to be run:

##include
This folder contains all the header files that declare functions and types common for all the test cases.

##src
This folder contains all the source files that define the functions declared in the header files of the include folder.

##external
This folder contains all the material used that has not be written by the author of this project but comes from external sources.



Each test-case folder presents this structure:

##src folder
This folder contains all the header and source files specific for the test-case with the main.cpp file and the data.pot one. 

##build folder
In this folder the CMakeLists.txt file has to be run.

##A CMakeLists.txt file
Running this file in the build folder (typing in the build folder "cmake ..") a suitable Makefile is generated in order to compile the test-case. 







 


 
 
