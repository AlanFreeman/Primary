# Primary

This repository contains the code for a signal processing model of the upstream primate visual system. The model simulates the macaque midget/parvocellular pathway from cones to layer 4Cbeta of primary visual cortex.

The code runs in Matlab. Ensure all .m and .mat files are in the current folder or in a folder on the Matlab path.

Run runPrim.m. This should plot a map of neuronal arrays.

Now edit runPrim.m. The main switch statement controls the analysis. You just ran case struct.x.y: try some other cases.

All cases set metadata used by stream, which executes the analyses. See "User guide to stream" and streamTute.m for short tutorials on stream.

The file setPrim.m contains code for setting model parameters, and is only required if changing parameters. After unzipping the data file Primary.zip, run setPrim.m.
