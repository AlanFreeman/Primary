# Primary

This repository contains the code for a signal processing model of the upstream primate visual system. The model simulates the macaque midget/parvocellular pathway from cones to layer 4Cbeta of primary visual cortex.

The code runs in Matlab. Ensure all .m and .mat files are on the Matlab path.

Run runPrim. This should plot a map of neuronal arrays.

Now edit runPrim. The main switch statement controls the analysis. You just ran case array.x.y: try some other cases.

All cases set metadata used by stream, which executes the analyses. See Guide to stream and streamTute.m for short tutorials on stream.

The file setPrim contains code for setting model parameters. After unzipping the data file Colour.zip, run setCol to fit the ROG to empirical data.
