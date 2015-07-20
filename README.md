Model Based Adaptive Optimal Design: Simulation of Adult to Children Bridging Study
======

These scripts was used for the project "Simulated model based adaptive optimal design of adult to children
bridging study using FDA stopping criteria which was presented at PAGE 2015
http://www.page-meeting.eu/default.asp?abstract=3614

By running the code in the page_results_analysis.R script, the plots can be replicated using the generated results
in the Page_results folder.

The code which runs MBAOD simulations for different levels of misspecification is available in the R-script file
Bridging_1_variable_weight_xoptim_log.R



## Installation

You need to have R, and the MBAOD package with all its requirements installed.  Download the latest version of R from www.r-project.org.

See the github page for the MBAOD package for more information regarding installation of this package:
https://github.com/andrewhooker/MBAOD

To allow for the utilization of the stopping critera, the file mbaod_simulate.R in the MBAOD-master/R -directory has to be replaced with the mbaod_simulate.R 
file provided with this download.

## Stopping criteria
The stopping criteria is implemented in the R-script file stop_critX.R.
In this version, the stopping criteria requires lineralizable scaling models as described by Wang et al.
In the example presented at PAGE, the scaling model used was a size and maturation model without hill coefficient and uncertianty on the TM50 parameter. 
This can however be altered within the stop_critX script to other models, as long as they are linearizable.

The next version will include a simulation based evaluation of the uncertianty around the estimated Clearance and volume parameters to allow for any scaling model