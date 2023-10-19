# An efficient approach to nowcasting the time-varying reproduction number
Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes

## Overview of Repository
Contained within this repository is R code pertaining to a simulation study discussed in the paper titled "An efficient approach to nowcasting the time-varying reproduction number" authored by Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes.

## Simulation codes
All necessary R scripts for conducting the simulation are located within this repository. These scripts are organized into two subdirectories: **Rnowcast** and **estimR**.

1. **Rnowcast**: This folder contains R codes for estimating the instantaneous reproduction number using the proposed methodology in the paper.
2. **estimR**: Uses reported cases and nowcasted values as inputs for the *estimR()* function within the EpiLPS package to estimate the reproduction number. 

For both approach the simulation codes are provided for a specific delay probabilities and nowcast date.
