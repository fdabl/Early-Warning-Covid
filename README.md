# Early Warning Signals of the Second COVID-19 Wave
This repository contains code to reproduce all analyses and figures in Dablander, F., Heesterbeek, H., Borsboom, D., & Drake, J. M. ([2022](https://royalsocietypublishing.org/doi/10.1098/rspb.2021.1809)). Overlapping Time Scales Obscure Early Warning Signals of the Second COVID-19 Wave. *Proceedings of the Royal Society B, 289*(1968), 20211809.

This repository is structured as follows:

  - **Code/helpers.R**: Contains helper functions for the empirical analysis and the simulation study.
  - **Code/setup-data.R**: Gets data from the WHO, estimates R(t), saves output into **Data/** (one file is too big for this repository).
  - **Code/methodology.R**: Creates Figure 1 describing the methodology.
  - **Code/empirical-analysis.R**: Analyses the data and creates Table 1.
  - **Code/figures.R**: Creates Figure 2 & 3 showing the data and summarizing the empirical results.
  - **Code/illustration.R**: Creates (raw) Figure 4 showing when indicators do and do not work.
  - **Code/simulation.R**: Code for running the simulation study, saves output under **Results/**, also creates (raw) Figure 5.
  - **Code/analysis.R**: Code for the empirical analysis.
