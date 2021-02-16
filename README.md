## Supporting information for "Still no evidence for disruption of global patterns of nest predation in shorebirds"

by Martin Bulla, Mihai Valcu and Bart Kempenaers

### **Overview**

Data, codes and outputs of the analyses. 

Complement the original [Bulla et al. 2019 comment in Science](https://science.sciencemag.org/content/364/6445/eaaw8529) in showing the lack of evidence for disruption of global patterns of nest predation in shorebirds by demonstrating issues with Kubelka et al.'s response ([2019, Science](https://science.sciencemag.org/content/364/6445/eaaw9893)) and with their initial study ([2018, Science](https://science.sciencemag.org/content/362/6415/680)).

### **Folders and files**

[Data](Data/): all data used in the analyses are (for convenience copied) from [Kubelka et al. 2018](https://doi.org/10.5061/dryad.45g90h4) or [Kubelka et al. 2019](https://osf.io/46bt3/), except for 'Rebuttal_data_with_DPR_by_yr.csv' and 'sources.xlsx', which were copied from [Bulla et al. 2019](https://osf.io/x8fs6/)

[R](R/): all r-scripts used in the analyses
- Constants_Functions.R loads functions and packages used in the other R-scripts (needs to be loaded before running the other scripts)
- Prepare_data_from_Bulla_et_al_2019.R script from [Bulla et al. 2019](https://osf.io/x8fs6/) prepares the datasets (but only after Constants_Functions.R is loaded) for Analyses_Note_6_Figure&Table_Beintema.R
- Analyses scripts generate the all outputs for the main text, figures, tables and notes

[Output](Output/): outputs of the Figures and Table scripts

Still_no_evidence.sublime-project - sublime project file

LICENSE - terms of reuse

.gitignore - inclusion and exclusion criteria for the repository