This R pipeline is designed for application on EHR datasets to analyze the physiological changes of
postoperative adhesion patients across multiple encounters.

The pipeline first filters patients to get a relevant study cohort, then uses LCA on a select patient
subset to define latent subgroups. Posterior probability will be used to sort each relevant patient encounter
into the defined subgroups, and patient physiological changes over encounters will be analyzed using
linear mixed-effects models.

The analysis ultimately aims to inform about patient physiological shifts, providing insights on applications
such as risk stratification, treatment, and preventative biomaterial design.

The code is currently linked to basic synthetic data .R files. To run the pipeline in an R environment, load all the
files and run file '19.R'. It will automatically source all the preceding files. The required R libraries as of right
now are: clinfun, corpcor, dplyr, ggalluvial, ggplot2, ggrepel, lme4, lmerTest, LMest, lubridate, reshape2, stats,
tidyLPA, and tidyverse.

If you have any suggestions, concerns, or would like to know more about the specifics of the methodology,
I can be contacted at [nguyen.minh201023@gmail.com]. The project is currently being managed and developed individually,
but I am highly open to collaboration. If you plan to use this pipeline or adapt it to your own research, please credit
this repository or reach out so we can collaborate!

Some disclaimers and notes:
- The code is still a work in progress. Some features have not been added, completed, or updated yet.
- I currently do not have access to real EHR schemas yet, so the way the data is joined from synthetic data files to
  the analysis files may not be fully applicable on real data structures.
- Since the code is built based on synthetic data files, I added code to omit or ignore NA values at some
  parts to ignore missing values that weren't included in the data generation. To add features like imputation methods,
  these parts must first be changed.
- The analysis variables are not final and are subject to change depending on availability and relevance.
- File 05.5.R is currently outdated and broken but it is only an optional side-analysis and will be updated later.
- The code is currently designed for a two-dataset study, but can be easily modified to accommodate
  a single-dataset analysis.
- The numbers in the .R file names usually signify the order in the pipeline.
