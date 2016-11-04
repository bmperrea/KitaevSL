# KitaevSL

These are the Matlab codes for computing finite correlation 
functions of the Kitaev spin liquid, an exotic many body physics model,
including finite temperature correlation functions.

For details about the Geyer initseq functions see [this wiki page](wiki/Geyer-stats-in-Matlab).

### Usage and function

One can begin by running the following to compute the Raman spectra
for a finite honeycomb model at zero temperature.
~~~~
rmax=5; b=10; s=0; p=0; T=0;
[I,En] = stretch_2D_6(rmax,b,s,p,T);
~~~~
Also try rmax = 25;

Here rmax is the radius of the honeycomb flake in layers of sites. b is the grunesien 
parameter representing the magnetic response to strain s (max fractional strain). 
Here p represents the flux sector (number of fluxes, randomly placed) and T represents 
the temperature relative to J = J^x = J^y = J^z (within that flux sector; used for the
spinons only). 

The function creates the matrices and plots spectra. It calls on dos1_loop, which does
the diagonalization. dos1_loop uses histwv_loop, a homemade vectorized code to
compute weighted histograms. The ouput I is a cell of histograms corresponding to 
DOS and Raman spectra. They are also plotted by the function (see the code for
which parts of I are which). En is simply the ground state energy of the flux sector.

The finite temperature plots are computed in a script called flux_controller_final,
which makes use of all of the other codes. This one probes flux sectors keeping a running 
best guess of the desired spectra (Raman spectra and DOS in this case). It continues to do
so until the expected error is lowenough. The code that creates the matrices for this is 
primarily stretch_flux6. However, in the beginning a rough estimate of the partition
function is made using initialize_pv2 -> find_pdist -> stretch_2D_6. 

Since this is a script, the parameters are set in the first few lines. You can change
them by changing the script itself, or commenting out the lines and setting them
outside the script first. The output is set in a .mat file, whose name is dynamically
created to contain the important parameters of the current calculation. There are
also a few plots of the current expected error that are continuously updated 
to allow the user to follow the progress.
 
flux_controller_final diagonalizes
matrices in parallel. If you don't have the parallel toolbox you can just replace 'parfor'
by 'for' and it will all run the same. 

Good luck.

bmperrea
