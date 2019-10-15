# OSA-Delirium
This repository holds code for the paper "Obstructive sleep apnea does not predispose to postoperative ICU delirium: a retrospective single-center analysis." The project is entirely in R. 
The project can be re-run using the `osa_analysis_flow.R` script. There are a small number of places where what was, in reality, an interactive calulation are simply hard-coded, such as caliper selection for matching and comparison of the sparse versus nonsparse BART variants.
Because the data are protected and difficult to de-identify, the osa_data folder is missing. It is entirely secondary data; if you want to use it you will have to obtain permission from Washington University's IRB.

## bartbcf
The primary analysis uses the methods of Hahn [1], which when we were creating the project did not have a reference implementation available. We modified the BART 1.7 R package to include the Bayesian Causal Forest calculations. The `bartbcf` folder contains the result. We validated the calculation on the examples of Hill [2], but don't intend to maintain this package as Hahn is actively developing his BCF package to accomplish the same goal.

[1] Hahn, P. R., Murray, J. S. & Carvalho, C. Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects. ArXiv170609523 Stat (2017).

[2] Hill, J. L. Bayesian Nonparametric Modeling for Causal Inference. J. Comput. Graph. Stat. 20, 217–240 (2011).
