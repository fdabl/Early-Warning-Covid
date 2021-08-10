# Early Warning Signals of the Second COVID-19 Wave
This repository contains code to reproduce all analyses and figures in Dablander, F., Heesterbeek, H., Borsboom, D., & Drake, J. M. ([2021](https://www.medrxiv.org/content/10.1101/2021.07.27.21261226v2)). Overlapping Time Scales Obscure Early Warning Signals of the Second COVID-19 Wave.

This repository is structured as follows:

  - **Code/helpers.R**: Contains helper functions for the empirical analysis and the simulation study.
  - **Code/setup-data.R**: Gets data from the WHO, estimates R(t), saves output into **Data/** (some files are too big for this repository).
  - **Code/methodology.R**: Creates Figure 1 describing the methodology.
  - **Code/figures.R**: Creates Figure 2 & 3 showing the data and summarizing the empirical results.
  - **Code/illustration.R**: Creates Figures 4 & 5 showing when indicators do and do not work.
  - **Code/simulation.R**: Code for running the simulation study, saves output under **Results/**, also creates Figures 6 & 7.
  - **Code/analysis.R**: Code for the empirical analysis.