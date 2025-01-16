[![View simple-BLDC-model on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/97082-simple-bldc-model)
# simple-BLDC-model
A simple (and probably inaccurate) brushless DC motor modelling script.

This software is entirely written in MATLAB and comes with two main files (a standard MATLAB script and a MATLAB live script) and a MATLAB function. The MATLAB live script includes also the mathematical theory behind the software and it has been already exported as PDF to be used as user manual.

Both the standard and the live script call the same MATLAB function to calculate the motor performance and both show the same charts in the end. The difference is that the MATLAB live script is interactive and shows live updates when changing the input parameters, while the MATLAB standard script can be easily expanded with other features or be included in other programs.

![example-output-charts](examples/example.png)

## Additional features
Actually, the `main.m` script also plots the motor map. The additional `propMatchCustom.m` script provides matching between motor and propeller, where the latter's data are loaded from external file.

![example-output-motor-map](examples/motor-map-example.png)

![example-output-motor-map](examples/prop-match-example.png)
