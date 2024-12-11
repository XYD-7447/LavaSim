# LavaSim

The simulation requires six input files, located in 'Input/' directory: 
ctl file containing model super-parameters, envir file specifying lunar environmental conditions, lava file defining lava properties, neibor file establishing neighborhood relationships, volcano file indicating vent locations and effusion rates, and zdata file representing the underground topography.

For execution, the codes have been optimized to operate with the CUDA module version 11.2, tailored for deployment on the GTX 1080 graphics processing unit.

To initiate the program, please enter the following commands into the Terminal:

'make'

'./build'
