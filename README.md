# LavaSim

The simulation requires six input files, located in 'Input/' directory: 
ctl file containing model super-parameters, envir file specifying lunar environmental conditions, lava file defining lava properties, neibor file establishing neighborhood relationships, volcano file indicating vent locations and effusion rates, and zdata file representing the underground topography.

The codes can be running using CUDA module 11.2 on the GTX 1080 card

Type the following code in the Command Prompt to run the program, LavaSim:

make
./build
