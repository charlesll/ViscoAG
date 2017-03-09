################################################################################
#
# Copyright (C) Charles Le Losq, Daniel R. Neuville
#
# This code is a supplementary material of the article
# Le Losq C., Neuville D. R. (2017) Molecular structure, configurational entropy and viscosity of silicate melts: link through the Adam and Gibbs theory of viscous flow. Journal of Non Cristalline Solids XX:XXX-XXX.
#
# The following code is under GNU GPL V3 license, see the LICENSE file.
#
# This file shows an example of how to generate predictions from the model
# with using the model_direct function, example is made with the KS3, NS3, KS4,
# and NS4 melts.
#
################################################################################

################################################################################
# LIBRARY LOADING
################################################################################

# We need the following libraries
using PyPlot
using Spectra
using Dierckx
using StatsBase
using JLD

# we also include the following code that contains the function allowing us
# to make predictions from the model
include("model_functions.jl")

################################################################################
# DATA IMPORT
################################################################################
NS4_data = readcsv("../data/NS4.csv",skipstart=1)
KS4_data = readcsv("../data/KS4.csv",skipstart=1)

NS3_data = readcsv("../data/NS3.csv",skipstart=1)
KS3_data = readcsv("../data/KS3.csv",skipstart=1)

################################################################################
# QN SPLINE IMPORTATION
################################################################################
spl_q2_na = load("../data/data_model.jld","spl_q2_na")
spl_q3_na = load("../data/data_model.jld","spl_q3_na")
spl_q4_na = load("../data/data_model.jld","spl_q4_na")
spl_q2_k = load("../data/data_model.jld","spl_q2_k")
spl_q3_k = load("../data/data_model.jld","spl_q3_k")
spl_q4_k = load("../data/data_model.jld","spl_q4_k")

################################################################################
# MODEL BEST AND BOOTSTRAP PARAMS IMPORT
################################################################################
params_best = load("../model_outputs/best_results.jld","params")
params_boot = load("../model_outputs/boot_results.jld","params_boot")
mean_boot = mean(params_boot,1)
println("Done.")

################################################################################
# ANALYSIS OF PARAMETER DISTRIBUTION
################################################################################
# See the documentation of Spectra.jl for using this function and its options.
#mean_boot, std_boot = bootperf(params_boot, plotting = "False", parameter = 13, feature = 0, savefigures = "False")

# constructing the parameter array
# the 2.28% and 97.72% quantiles are taken and corrected from any bias between the best
# parameters values and the mean parameters values from the bootstrap.
params_withquantile = ones(size(params_best)[2],3)
params_withquantile[:,1] = params_best[:]
for i = 1:size(params_best,2)
    bias = mean_boot[i] - params_best[i] # we suppress any bias
    params_withquantile[i,2] = (quantile(params_boot[:,i],0.0228)-bias)
    params_withquantile[i,3] = (quantile(params_boot[:,i],0.9772)-bias)
end

################################################################################
# ESTIMATES FOR THE KS-NS COMPOSITIONS
################################################################################

# We are going to use the function with a given SiO2 content
# We provide a known Tg to the model, to allow a better estimation
T = collect(500:1.0:2300) # temperature for the estimates

n_mod_NS4, Tg_NS4, ScTg_NS4, ~, ~, ~, ~, Be_NS4, q2_ns4, q3_ns4, q4_ns4 = model_direct(80.0,0.0,T,params_withquantile,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k,Tg=753.0)
n_mod_KS4, Tg_KS4, ScTg_KS4, ~, ~, ~, ~, Be_KS4, q2_ks4, q3_ks4, q4_ks4 = model_direct(80.0,1.0,T,params_withquantile,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k,Tg=760.0)
n_mod_NS3, Tg_NS3, ScTg_NS3, ~, ~, ~, ~, Be_NS3, q2_ns3, q3_ns3, q4_ns3 = model_direct(75.0,0.0,T,params_withquantile,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k,Tg=737.0)
n_mod_KS3, Tg_KS3, ScTg_KS3, ~, ~, ~, ~, Be_KS3, q2_ks3, q3_ks3, q4_ks3 = model_direct(75.0,1.0,T,params_withquantile,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k,Tg=750.0)

################################################################################
# FIGURE GENERATION
################################################################################

fig = figure()
scatter(10000./NS4_data[:,5],NS4_data[:,6],color="black",label="NS4")
scatter(10000./KS4_data[:,5],KS4_data[:,6],color="red",label="KS4")
scatter(10000./NS3_data[:,5],NS3_data[:,6],color="blue",label="NS3")
scatter(10000./KS3_data[:,5],KS3_data[:,6],color="cyan",label="KS3")
plot(10000./T,n_mod_NS4,color="black")
plot(10000./T,n_mod_KS4,color="red")
plot(10000./T,n_mod_NS3,color="blue")
plot(10000./T,n_mod_KS3,color="cyan")

# axis stuffs
xlim(5,14)
ylim(0,12.5)
xlabel(L"1000/T, K^{-1}", fontname = "Sans", fontsize = 18)
ylabel("Viscosity, Pa s", fontname = "Sans", fontsize = 18)

# legend
legend(fancybox=true,loc="best")
annotate("Points are data \nLines are \nmodel predictions",xy=(0.75,0.2),xycoords="axes fraction",fontsize=18,fontname="Sans",ha="center")

savefig("../figures/KNS4_KNS3_example.pdf")
