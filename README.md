# KitaevSL

These are the Matlab codes for computing finite correlation 
functions of the Kitaev spin liquid, an exotic many body physics model,
including finite temperature correlation functions.

For details about the Geyer initseq functions see [this wiki page](https://github.com/bmperrea/KitaevSL/wiki/Geyer-stats-in-Matlab).

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
the diagonalization to compute Raman correlation functions. dos1_loop uses histwv_loop, 
a homemade vectorized code to compute weighted histograms. The ouput I is a cell of 
histograms corresponding to DOS and Raman spectra. They are also plotted by the 
function (see the code forwhich parts of I are which). En is simply the ground state 
energy of the flux sector.

The finite temperature plots are computed in a script called flux_controller_MCMC,
which makes use of many other codes. This one runs a typical Markoc Chain Monte 
Carlo algorithm where the different configurations are different flux sectors. 
The primary endpoint to run the calculation over multiple temperatures is the 
script do_some_MCMC

Since this is a script, the parameters are set in the first few lines. You can change
them by changing the script itself, or commenting out the lines and setting them
outside the script first. The output is set in a .mat file, whose name is dynamically
created to contain the important parameters of the current calculation. There are
also a few plots of the current expected error that are continuously updated 
to allow the user to follow the progress.
 
do_some_MCMC does work in on up to 12 parallel cores (the limit as of 2016).
If you don't have the parallel toolbox you can just replace 'parfor'
by 'for' and it will all run the same. 

The codes called by flux_controller_MCMC are in the following three folders:
- work_functions
- raman_stuff
- geyer_stats

The supercomputer has an example of a pbs script that one can run on the 
University of Minnesota's msi. This is done by first ssh-ing into the 
login node at msi `ssh login` (login.msi.umn.edu) -- see msi.umn.edu for support.
If you can do that, then just do `ssh mesabi` and then go to the directory
where you put all of the codes (I recommend all in one folder on msi since 
setting the path is difficult on the supercomputer). Then submit the job like:
`qsub -q small setup.pbs`

If you screwed up you can `qsig -s SIGINT something` where 'something' 
is the name of the job returned after you started (also in your email).
Otherwise you can just watch the .o and .e (output and error) files
until it's done and then you collect the data.

Good luck.

bmperrea
