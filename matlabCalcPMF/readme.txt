===========================
	*README*
===========================

This is a package of Matlab Code for calculating the Potential of Mean Force 
from forword and reverse trajectories.

===========================
	*setup*	
===========================

Set the spring constant (ks), starting point (bA) and ending point (bB) in the 
matlab scripts (fr_pmf.m, je_pmf.m, mlm_pmf.m).

===========================
	*usage*
===========================

All the matlab code files should be copied to the folder which contains the
trajectory files. 

(Trajectory file name starting with 'F-' indicates that it is a forward 
trajectory and trajectory file name starting with 'R-' indicates that 
it is a reverse trajectory.) 

Then, in the matlab command window, simply 1) run the command je_pmf to calculate
the PMF using Jarzynski Equality (cumulant approximation) with forward 
trajectoris and reverse trajectories, respectly. 2) Run mlm_pmf to calculate the
PMF using ML method. 3) Run fr_pmf to calculate the PMF using FR method. In addition,
this command will give the work distributions for both forward and reverse trajectories.

A subfolder named "data" will be created and all the results will be stored in this folder.




Jiong Zhang
Aug 23, 2011
