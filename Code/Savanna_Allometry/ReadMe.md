Savanna_Allometry
========================================================

General Information
-------------------------

Code for reproducing results from "Allometric convergence of savanna trees and implications for the use of plant scaling models in variable ecosystems" by Andrew T. Tredennick, Lisa Patrick Bentley, and Niall P. Hanan 2013 published in PLoS ONE, doi: ######. Code was written by Andrew T. Tredennick (JAGS code for hierarchical Bayesian model highly informed by code available as supplementary material for "Evaluating scaling models in biology using hierarchical Bayesian approaches" by Charles A. Price, Kiona Ogle, Ethan P. White, and Joshua S. Weitz published in Ecology Letters, doi: 10.1111/j.1461-0248.2009.01316.x)

The code and data in this repository allow a user to run the hierarchical Bayesian model described in the paper as well as reproduce all the figures.

Requirements: JAGS (**J**ust **A**nother **G**IBBS **S**ampler, http://mcmc-jags.sourceforge.net/), R 2.x (http://www.r-project.org/) and the following R packges:
* rjags
* coda
* rdyrad
* ggplot2
* scales
* gridExtra

The R code (``SavannaAllometry_CallHBmodel.R``) includes lines that will automatically install these packages. These lines can be commented out (``#``).

Running the HB model will take considerable time at 1,000,000 iterations of three MCMC chains. As such, it is highly recommended you change the variable ``sample.n`` in the file ``SavannaAllometry_CallHBmodel.R`` until you are sure the model is running properly. However, you should note that with fewer iterations the HB model will likely not converge.

Data use: Data is stored in a Dryad repository (link) for replication purposes. If you wish to use the data for further research please access the data from Dryad and cite appropriately. You can also optionally contact Andrew Tredennick at atredenn@gmail.com or Andrew.Tredennick@colostate.edu if you would like to collaborate.

Included Files
-------------------------
* ``SavannaAllometry_CallHBmodel.R`` -- gets data from Dryad repository, cleans the data, calls the HB model in the JAGS code (``JAGS_HBmodel_SavannaAllometry.r``)
* ``JAGS_HBmodel_SavannaAllometry.r`` -- contains code for the hiearchical Bayesian model in JAGS format

License
-------------------------


Contact Information
-------------------------
Andrew Tredennick's emails: atredenn@gmail.com or Andrew.Tredennick@colostate.edu

Andrew's website: http://warnercnr.colostate.edu/~atredenn/Welcome.html
