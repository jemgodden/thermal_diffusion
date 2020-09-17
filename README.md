# thermal_diffusion

This program has the ability to model how the temperature of a medium changes as a nuclear waste rod diffuses heat into its surroundings, making use of the cylindrical symmetry of the rod. This is predominantly an attempt to recreate a method used by Dr Olsen-Kettle in a 2011 paper.

This program can also model how temperature changes radially in an egg as it is surrounded by boiling water, making use of the modelled spherical symmetry of the egg. This uses the same method as used for modelling a nuclear waste rod, but is adapted for this scenario.

# Downloading

The easiest way to download the code is to press the download button on the home page of the GitHub repository.

# Requirements

To run the program a C compiler is required, along with the following libraries: stdio.h, stdlib.h, string.h, stdbool.h, math.h, getopt.h and various modules from the gsl library.

# Use

This program will run from the command line, after being compiled. It can be compiled using the prewrriten compile function, via the command line: ./compile

The command line arguments take the form:

./a.out --model (--time_step val1) (--resolution val2) (--mass val3) (--temperature val4) (--output_file str1)

The models are outline at the top of the file, being either egg or rod.

Time_step and resolution can be used to specify the parameters of the simulation for either model.

Mass and temperature will only apply to the egg model.

All of these have model-specific default values if none are stated.

The output_file that data is printed to can be given a certain name, but if it is not then the file will be given the default name.

Graphs showing how the temperature of the object in the model changes over time can be printed after the simulation's conclusion.

For the egg model, the time taken for it to soft-boil can also be printed.

# Log

Initial version uploaded to GitHub.
