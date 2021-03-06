# inlaVP

This package contains functions to fit Variance Partiotioning (VP) disease mapping models in R-INLA; these models are described in the paper by 
Franco-Villoria M et al. (2021) presenting a new approach to model space-time data. The main idea is to specify space-time interaction models in a way that the contribution of space, time and their interaction, to the total variance can be estimated. One relevant application is spatio-temporal disease mapping, but the proposed models may be adopted in all those applications where intrinsic GMRFs are meant as tools to perform smoothing in more than one dimension; 
for instance in the analysis of grid-data such as those arising from agricultural field trials or spatio-temporal data from environmental studies and ecological surveys.

The long-term plan is to produce some extra material on how to use the package, wignettes etc. 

In case you find the package useful and you want to use it and refer to it, please quote:

*Franco-Villoria M, Ventrucci M, Rue H. (2021). 'Variance partitioning in spatio-temporal disease mapping models'. arXiv:2109.13374}*

If you have any issue about the package (or any ideas on how to improve the material in there) please send an email to massimo.ventrucci@unibo.it.
