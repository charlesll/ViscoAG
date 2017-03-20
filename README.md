# ViscoAG software

Copyright (C) Charles Le Losq, Daniel R. Neuville

GNU GPL V3 license, see the LICENSE file.

This software allows calculating the viscosity of Na2O-K2O-SiO2 silicate melts through the knowledge of their structure and thermodynamic properties, using the Adam and Gibbs (1965) framework.

The software is written in the Julia programming language, see http://julialang.org/

Details are provided in the article:

Le Losq C., Neuville D. R. (2017) Molecular structure, configurational entropy and viscosity of silicate melts: link through the Adam and Gibbs theory of viscous flow. Journal of Non Cristalline Solids 463:175-188.

Until the 8th of May, the article is available for free at: https://authors.elsevier.com/a/1UkW854fAxj9K

# Installation

You need a running version of Julia to be able to run the script. See http://julialang.org/ for your system.

The required Julia packages are listed in the REQUIREMENT file. Please consult the documentation of the different packages for any specific installation instructions, as well as the Julia documentation for general information on Julia.

The script can be runned either directly from the terminal, or using an IDE as the Atom/Juno IDE. See http://junolab.org/ for further information.

# Descriptions of the files and folders

Dataset.ods : .ods spreadsheet containing the viscosity and Qn distribution datasets used in the model;

./data/ : folder, contains the data in .csv that are used in the various julia scripts;

./figures/: folder, contains the figures generated by the scripts;

./model_outputs/: folder, contains the outputs generated by model_fitting.jl. Do not erase the .jld files as they contain the best model as well as all the bootstrap results.

./src/ : folder, contains the following Julia scripts

  - model_functions.jl

    This file contains several functions used in the different Julia scripts. See file comments for further details.

  - model_fitting.jl

    This script allows fitting the model to the data. It also generate the figure 6 (Figure6.pdf). See file comments for further details.

  - model_predictions.jl

    This script allows generating the predictions of the model associated with the figures 7 to 10 in the Le Losq and Neuville paper. Please see file comments for further details, and the Le Losq and Neuville (2017) paper for the legends and captions.

  - model_boot_glueresults.jl

    This scripts allows gluing the results from bootstrap. Indeed, running 40,000 resampling and fitting instances in one loop tends to saturate the computer memory. To avoid that, we runned 10 times 4000 resampling and fitting instances. The results are stick together with this code. See file comments for further details.

  - model_example_NS_KS_4_3.jl

    This script is an example of how predictions can be generated with using the functions from model_functions.jl. It generates examples with using data for KS3, KS4, NS3 and NS4 melts. See file comments for further details.

# Running the codes

With Julia and the required packages on your system, you can simply run them by typing "julia" followed by their name in a terminal. For instance:

  julia model_fitting.jl

will fit the model to the dataset. As said above, you can also use an IDE such as the Juno one.
