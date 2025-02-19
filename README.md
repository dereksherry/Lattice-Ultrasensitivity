This contains the primary files and code related to the LU Paper. Files are named according to which figure they correspond to, with a few exceptions.

Code for running simulations of the LU model
LUSim.m defines the class describing a single instance of the LU model
LUSimulator.m is the function which actually runs the simulations
parameterGeneratorLU.m is a function which outputs the initial conditions and default parameters for the simulations


Code to make figure creation easier:
axxytofigxy.m was written by someone else and taken from the MATLAB repositories to convert from axis coordinates to figure coordinates, and makes the generation of some figure annotations easier.
drawShadedRectangle.m is function I wrote to make the graphic design for Figure 1 easier

