################################################################################
#
# Copyright (C) Charles Le Losq, Daniel R. Neuville
#
# This code is a supplementary material of the article
# Le Losq C., Neuville D. R. (2017) Molecular structure, configurational entropy and viscosity of silicate melts: link through the Adam and Gibbs theory of viscous flow. Journal of Non Cristalline Solids XX:XXX-XXX.
#
# The following code is under GNU GPL V3 license, see the LICENSE file.
#
# This file allows gluing the results from the bootstrap,
# recorded in 10 different files.
#
################################################################################

################################################################################
# LIBRARY LOADING
################################################################################

using JLD

################################################################################
# DATA LOADING
################################################################################

boot1 = load("../model_outputs/boot_results_1.jld","params_boot")
boot2 = load("../model_outputs/boot_results_2.jld","params_boot")
boot3 = load("../model_outputs/boot_results_3.jld","params_boot")
boot4 = load("../model_outputs/boot_results_4.jld","params_boot")
boot5 = load("../model_outputs/boot_results_5.jld","params_boot")
boot6 = load("../model_outputs/boot_results_6.jld","params_boot")
boot7 = load("../model_outputs/boot_results_7.jld","params_boot")
boot8 = load("../model_outputs/boot_results_8.jld","params_boot")
boot9 = load("../model_outputs/boot_results_9.jld","params_boot")
boot10 = load("../model_outputs/boot_results_10.jld","params_boot")

boot = [boot1;boot2;boot3;boot4;boot5;boot6;boot7;boot8;boot9;boot10]

################################################################################
# SIZE OF DATASET = BOOTSTRAP SIZE
################################################################################

println("Size of bootstrap array is $(size(boot,1))")

################################################################################
# SAVING THE GLUED BOOTSTRAP RESULTS
################################################################################

save("../model_outputs/boot_results.jld","params_boot",boot)
println("Bootstrap results glued together. Done.")
