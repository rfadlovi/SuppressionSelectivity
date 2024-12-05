# SuppressionSelectivity
This is the repository for the code and data associated with the paper Selectivity of invasive species suppression efforts influences control efficacy in the Journal of Applied Ecology.

Associated Dryad repository DOI: 10.5061/dryad.jdfn2z3mw (Fadlovich et al. 2024).

SuppressionSelectivity.R contains all of the main code for running the model and creating the figures.

selectivity.bug is the JAGS model code. It is called in SuppressionSelectivity.R and does not need to be opened unless you are interested in modifying.  

The inputs folder contains data files that are read into the model. These are the files in Dryad, plus one file that has been published as part of an annual report cited in the text. 

The outputs folder contains data files that are produced while running this code. The creation of these data files has been commented out in the code, but can be uncommented to allow for easy retention of alternative model runs. 

The code is divided into 7 main sections 

* About this code - general description

* What data is needed to run this code - an overview of data requirements (age or stage based population estimates and catch or selectivity data) and the loading of the Utah Lake data files. 

* Selectivity Estimates - Calculating the observed selectivity for Utah Lake and 

* Simulations - the selectivity simulations discussed in the paper. Practitioners interested in tinkering with the selectivity and effort scenarios and not the underlying Utah Lake specific selectivity calculations can start here 

* Figures - plotting of figure 1 through 4 

* Additional Calculations - additional calculations that appear in the manuscript  

* Figures (Supplement) - plotting of supplemental figures 