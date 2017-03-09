################################################################################
#
# Copyright (C) Charles Le Losq, Daniel R. Neuville
#
# This code is a supplementary material of the article
# Le Losq C., Neuville D. R. (2017) Molecular structure, configurational entropy and viscosity of silicate melts: link through the Adam and Gibbs theory of viscous flow. Journal of Non Cristalline Solids XX:XXX-XXX.
#
# The following code is under GNU GPL V3 license, see the LICENSE file.
#
# This file generates the figures of the paper related to the model results.
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
# QN SPLINE IMPORTATION
################################################################################
spl_q2_na = load("../data/data_model.jld","spl_q2_na")
spl_q3_na = load("../data/data_model.jld","spl_q3_na")
spl_q4_na = load("../data/data_model.jld","spl_q4_na")
spl_q2_k = load("../data/data_model.jld","spl_q2_k")
spl_q3_k = load("../data/data_model.jld","spl_q3_k")
spl_q4_k = load("../data/data_model.jld","spl_q4_k")

################################################################################
# MODEL and BOOTSTRAP RESULTS IMPORT
################################################################################
params_best = load("../model_outputs/best_results.jld","params")
params_boot = load("../model_outputs/boot_results.jld", "params_boot")
println("Importations done.")

################################################################################
# ANALYSIS OF PARAMETER DISTRIBUTION
################################################################################

mean_boot = mean(params_boot,1) # the mean value of the parameter distribution = parameter value from bootstrap

# To look at all the parameter probability density we use the following line to plot all the distributions:
name_parameters=["Ae" "K1" "K2" "Be Q2 Naenv" "Be Q3 Naenv" "Be Q2 Kenv" "Be Q3 Kenv" "Be Q4" "Sconf(Tg) Q2 Naenv" "Sconf(Tg) Q3 Naenv" "Sconf(Tg) Q2 Kenv" "Sconf(Tg) Q3 Kenv" "Sconf(Tg) Q4"]
for i = 1:13;    figure();    PyPlot.plt[:hist](params_boot[:,i],1000);    xlabel("Parameter "*name_parameters[i]);    ylabel("Number of bootstrap");end

# Printing the best parameter values and their 95% confidence intervals
for i = 1:13
    bias = mean_boot[i] - params_best[i] # we suppress any bias
    println("Parameter "*name_parameters[i]*", best value $(round(params_best[i],2)), boot bias $(round(bias,2)), low and high bounds are $(round(quantile(params_boot[:,i],0.0228)-bias,2)) and $(round(quantile(params_boot[:,i],0.9772)-bias,2)),respectively.")
    println("")
end
println("\nWith bonds in a -/+ range form:")
for i = 1:13
    bias = mean_boot[i] - params_best[i] # we suppress any bias
    println("Parameter "*name_parameters[i]*", best value $(round(params_best[i],2)), low and high bounds are -$(round(params_best[i]-(quantile(params_boot[:,i],0.0228)-bias),2)) and +$(round((quantile(params_boot[:,i],0.9772)-bias)-params_best[i],2)),respectively.")
    println("")
end

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
#### GENERATING THE VALUES FOR THE ENDMEMBERS
################################################################################

sio2_sc = collect(60:0.2:100)
entropy_NS = collect(60:0.2:100)
entropy_KS = collect(60:0.2:100)
entropy_NS_lin = collect(60:0.2:100)
entropy_KS_lin = collect(60:0.2:100)
entropy_NS_QnMix = collect(60:0.2:100)
entropy_KS_QnMix = collect(60:0.2:100)
Be_NS = collect(60:0.2:100)
Be_KS = collect(60:0.2:100)
BeSc_NS = collect(60:0.2:100)
BeSc_KS = collect(60:0.2:100)
tg_NS = collect(60:0.2:100)
tg_KS = collect(60:0.2:100)
qn_NS = ones(size(tg_NS,1),3)
qn_KS = ones(size(tg_NS,1),3)

# For the errors:
entropy_NS_error = ones(size(entropy_NS,1),2)
Be_NS_error = ones(size(entropy_NS,1),2)
BeSc_NS_error = ones(size(entropy_NS,1),2)
entropy_KS_error = ones(size(entropy_KS,1),2)
Be_KS_error = ones(size(entropy_KS,1),2)
BeSc_KS_error = ones(size(entropy_NS,1),2)
tg_NS_error = ones(size(entropy_NS,1),2)
tg_KS_error = ones(size(entropy_KS,1),2)

for i = 1:size(entropy_NS,1)
    #### FOR NA2O-SiO2
    n_mod, Tg, ScTg, Sc_lin, Sc_qn, ~, ~, Be, q2_ns, q3_ns, q4_ns = model_direct(sio2_sc[i],0.0,[1000],params_withquantile[:,1],spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    entropy_NS[i] = ScTg
    entropy_NS_lin[i] = Sc_lin
    entropy_NS_QnMix[i] = Sc_qn
    Be_NS[i] = Be
    BeSc_NS[i] = Be/ScTg
    tg_NS[i] = Tg_calculation(params_withquantile[1,1],Be,ScTg)
    qn_NS[i,1] = q2_ns
    qn_NS[i,2] = q3_ns
    qn_NS[i,3] = q4_ns

    error_n, tg_NS_error[i,:], entropy_NS_error[i,:], Be_NS_error[i,:],BeSc_NS_error[i,:] = error_propagation(params_best,params_boot,sio2_sc[i],0.0,1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    #### FOR K2O-SiO2
    n_mod2, Tg2, ScTg2, Sc_lin, Sc_qn, ~, Cpc, Be2,q2_ks, q3_ks, q4_ks  = model_direct(sio2_sc[i],1.0,[1000],params_withquantile[:,1],spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    entropy_KS[i] = ScTg2
    entropy_KS_lin[i] = Sc_lin
    entropy_KS_QnMix[i] = Sc_qn
    Be_KS[i] = Be2
    BeSc_KS[i] = Be2/ScTg2
    tg_KS[i] = Tg_calculation(params_withquantile[1,1],Be2,ScTg2)
    qn_KS[i,1] = q2_ks
    qn_KS[i,2] = q3_ks
    qn_KS[i,3] = q4_ks

    error_n, tg_KS_error[i,:], entropy_KS_error[i,:], Be_KS_error[i,:],BeSc_KS_error[i,:] = error_propagation(params_best,params_boot,sio2_sc[i],1.0,1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

end

################################################################################
#### GENERATING THE VALUES FOR THE Na-K MIX at 80 and 66% SiO2
################################################################################

xk_model = collect(0.:0.01:1.0)
entropy_KNS = ones(size(xk_model,1),2) # two columns for values at 80 and 66 mol% SiO2
entropy_KNS_lin = ones(size(xk_model,1),2)
entropy_KNS_QnMix = ones(size(xk_model,1),2)
entropy_KNS_NaKMix = ones(size(xk_model,1),2)

Be_KNS = ones(size(xk_model,1),2)
BeSc_KNS = ones(size(xk_model,1),2)
tg_KNS = ones(size(xk_model,1),2)
qn_KNS = ones(size(xk_model,1),3,2)

# For the error propagation:
entropy_KNS_error = ones(size(xk_model,1),2,2)
Be_KNS_error = ones(size(xk_model,1),2,2)
BeSc_KNS_error = ones(size(xk_model,1),2,2)
tg_KNS_error = ones(size(xk_model,1),2,2)

for i = 1:size(xk_model,1)
    #### For tetrasilicates, last dimension = 1 in arrays
    n_mod, Tg, ScTg, Sc_lin, Sc_qn, ScAlka, ~, Be, q2_ns, q3_ns, q4_ns = model_direct(80.0,xk_model[i],[1000],params_best,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    entropy_KNS[i,1] = ScTg
    entropy_KNS_lin[i,1] = Sc_lin
    entropy_KNS_QnMix[i,1] = Sc_qn
    entropy_KNS_NaKMix[i,1] = ScAlka

    Be_KNS[i,1] = Be
    BeSc_KNS[i,1] = Be/ScTg
    tg_KNS[i,1] = Tg_calculation(params_withquantile[1,1],Be,ScTg)
    qn_KNS[i,1,1] = q2_ns
    qn_KNS[i,2,1] = q3_ns
    qn_KNS[i,3,1] = q4_ns

    error_n, tg_KNS_error[i,:,1], entropy_KNS_error[i,:,1], Be_KNS_error[i,:,1],BeSc_KNS_error[i,:,1] = error_propagation(params_best,params_boot,80.0,xk_model[i],1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    #### FOR disilicates, last dimension = 1 in arrays
    n_mod, Tg, ScTg, Sc_lin, Sc_qn, ScAlka, Cpc, Be, q2_ns, q3_ns, q4_ns = model_direct(66.0,xk_model[i],[1000],params_best,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    entropy_KNS[i,2] = ScTg
    entropy_KNS_lin[i,2] = Sc_lin
    entropy_KNS_QnMix[i,2] = Sc_qn
    entropy_KNS_NaKMix[i,2] = ScAlka

    Be_KNS[i,2] = Be
    BeSc_KNS[i,2] = Be/ScTg
    tg_KNS[i,2] =  Tg_calculation(params_withquantile[1,1],Be,ScTg)
    qn_KNS[i,1,2] = q2_ns
    qn_KNS[i,2,2] = q3_ns
    qn_KNS[i,3,2] = q4_ns

    error_n, tg_KNS_error[i,:,2], entropy_KNS_error[i,:,2], Be_KNS_error[i,:,2],BeSc_KNS_error[i,:,2] = error_propagation(params_best,params_boot,66.0,xk_model[i],1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

end

################################################################################
#### FIGURE 7
################################################################################

fig_paper = figure(figsize=(12,7))
PyPlot.rc("xtick", labelsize=14)
PyPlot.rc("ytick", labelsize=14)

fig_paper[:text](0.5, 0.02, L"Alkali content, as mol% M$_2$O in SiO$_2$", ha="center",fontsize = 16,fontname="Arial")
fig_paper[:text](0.5, 0.96, L"Alkali content, as mol% M$_4$O$_2$ in SiO$_2$", ha="center",fontsize = 16,fontname="Arial")

# transformation de l'axe X en M4O2 au lieu de M2O
function tick_function(X::Array{Float64})
    return round(Int,X./2./(X./2+(100-X)) .*100)
end

ax1 = subplot(1,2,1)

# MODEL VALUES
ax1[:plot](100.0-sio2_sc,entropy_NS,color="r",linestyle="--",linewidth=2.5)
ax1[:plot](100.0-sio2_sc, entropy_NS_error[:,1],color="r",linestyle="--",linewidth=0.5)
ax1[:plot](100.0-sio2_sc, entropy_NS_error[:,2],color="r",linestyle="--",linewidth=0.5)

# DATA
ax1[:plot](100.0.-[80. 66.6 60.], [9.67 8.02 6.72],linestyle="",marker="s",markersize=14,color="r",alpha=0.4) # NEUVILLE 2006 DATA
ax1[:plot](100.0.-[75.], [7.71],linestyle="",marker="d",markersize=14,color="r",alpha=0.4) # LE LOSQ 2014 DATA
ax1[:plot](100.0.-[80. 75. 66.66 60.0], [10.92 9.40 9.54 9.75],marker="o",markersize=14,color="r",alpha=0.4) # RICHET 1984

# MODEL VALUES
ax1[:plot](100.0-sio2_sc,entropy_KS,color="b",linestyle="-",linewidth=2.5)
ax1[:plot](100.0-sio2_sc, entropy_KS_error[:,1],color="b",linestyle="-",linewidth=0.5)
ax1[:plot](100.0-sio2_sc, entropy_KS_error[:,2],color="b",linestyle="-",linewidth=0.5)

# DATA
ax1[:plot](100.0.-[75.], [9.42],marker="o",markersize=14,color="b",alpha=0.4) # RICHET 1984 DATA

ax1[:set_xlim](0.,40.)
ax1[:set_ylim](5,14)

ax1[:set_ylabel](L"S$^{conf}(T_{g})$, J mol$^{-1}$ K$^{-1}$",fontsize = 16, fontname="Arial")
ax1[:annotate]("A)",xy=(0.1,0.9),xycoords="axes fraction",fontsize = 20, fontname="Arial")

ax1_b = ax1[:twiny]()
ax1_b[:set_xlim](ax1[:get_xlim]())
new_tick_locations=[9.53, 18.18, 26.08, 33.33, 40.]
ax1_b[:set_xticks](new_tick_locations)
ax1_b[:set_xticklabels](tick_function(new_tick_locations))

ax2 = subplot(1,2,2)

model_NS = ax2[:plot](100.0-sio2_sc,Be_NS,color="r",linewidth=2.5,linestyle="--",label=L"Model Na$_2$O-SiO$_2$")
ax2[:plot](100.0-sio2_sc, Be_NS_error[:,1],color="r",linestyle="--",linewidth=0.5)
ax2[:plot](100.0-sio2_sc, Be_NS_error[:,2],color="r",linestyle="--",linewidth=0.5)
N2006 = ax2[:plot](100.0.-[80. 66.6 60.], [101600 79630 63570],linestyle="",marker="s",markersize=14,color="r",alpha=0.4,label="N2006") # NEUVILLE 2006 DATA
LL2014 = ax2[:plot](100.0.-[75.], [76976],linestyle="",marker="d",markersize=14,color="r",alpha=0.4,label="LL2014") # LE LOSQ 2014 DATA
R1984 = ax2[:plot](100.0.-[80. 75. 66.66], [113711. 94685. 93560.],marker="o",markersize=14,color="r",alpha=0.4,label="R1984") # RICHET 1984

model_KS = ax2[:plot](100.0-sio2_sc,Be_KS,color="b",linewidth=2.5,linestyle="-",label=L"Model K$_2$O-SiO$_2$")
ax2[:plot](100.0-sio2_sc, Be_KS_error[:,1],color="b",linestyle="-",linewidth=0.5)
ax2[:plot](100.0-sio2_sc, Be_KS_error[:,2],color="b",linestyle="-",linewidth=0.5)
ax2[:plot](100.0.-[75.], [96600],marker="o",markersize=14,color="b",alpha=0.4) # RICHET 1984

ax2[:yaxis][:set_ticks_position]("right")
ax2[:yaxis][:set_label_position]("right")
ax2[:set_ylabel](L"B$_e$, J mol$^{-1}$",fontsize = 16, fontname="Arial")

ax2[:set_ylim]([40000,170000])
ax2[:set_xlim]([0.,40.])

ax2[:annotate]("B)",xy=(0.8,0.9),xycoords="axes fraction",fontsize = 20, fontname="Arial")
# for whatever reason the legend does not work anymore there... see the paper for the legends.
#ax2[:legend](handles=[model_NS, model_KS, N2006, LL2014, R1984],frameon=false,numpoints=1,bbox_to_anchor=(0.65, 0.3))


ax2_b = ax2[:twiny]()
ax2_b[:set_xlim](ax2[:get_xlim]())
#already defined upstairs => new_tick_locations=[9.53, 18.18, 26.08, 33.33, 40.]
ax2_b[:set_xticks](new_tick_locations)
ax2_b[:set_xticklabels](tick_function(new_tick_locations))

show(ax2)

savefig("../figures/Figure7.pdf")

################################################################################
#### FIGURE 8
################################################################################

fig_mix = figure(figsize=(12,7))
PyPlot.rc("xtick", labelsize=14)
PyPlot.rc("ytick", labelsize=14)

fig_mix[:text](0.5, 0.02, L"$X_K = \frac{K}{K+Na}$", ha="center",fontsize = 16,fontname="Arial")

ax1 = subplot(1,2,1)

# MODEL VALUES
ax1[:plot](xk_model, entropy_KNS[:,1],color="black",linestyle="-",linewidth=2.5)
ax1[:plot](xk_model, entropy_KNS_error[:,1,1],color="black",linestyle="-",linewidth=0.5)
ax1[:plot](xk_model, entropy_KNS_error[:,2,1],color="black",linestyle="-",linewidth=0.5)

# MODEL VALUES
ax1[:plot](xk_model, entropy_KNS[:,2],color="green",linestyle="--",linewidth=2.5)
ax1[:plot](xk_model, entropy_KNS_error[:,1,2],color="green",linestyle="--",linewidth=0.5)
ax1[:plot](xk_model, entropy_KNS_error[:,2,2],color="green",linestyle="--",linewidth=0.5)

ax1[:set_ylim](5,15)

ax1[:set_ylabel](L"S$^{conf}(T_{g})$, J mol$^{-1}$ K$^{-1}$",fontsize = 16, fontname="Arial")
ax1[:annotate]("A)",xy=(0.1,0.9),xycoords="axes fraction",fontsize = 20, fontname="Arial")

ax1[:annotate](L"80 mol% SiO$_2$",xy=(0.6,14.2),xycoords="data",fontsize = 16, fontname="Arial",color="black",ha="center")
ax1[:annotate](L"66 mol% SiO$_2$",xy=(0.6,11.5),xycoords="data",fontsize = 16, fontname="Arial",color="green",ha="center")

ax2 = subplot(1,2,2)

# MODEL VALUES
ax2[:plot](xk_model,Be_KNS[:,1],color="black",linestyle="-",linewidth=2.5)
ax2[:plot](xk_model, Be_KNS_error[:,1,1],color="black",linestyle="-",linewidth=0.5)
ax2[:plot](xk_model, Be_KNS_error[:,2,1],color="black",linestyle="-",linewidth=0.5)

ax2[:plot](xk_model,Be_KNS[:,2],color="green",linestyle="--",linewidth=2.5)
ax2[:plot](xk_model, Be_KNS_error[:,1,2],color="green",linestyle="--",linewidth=0.5)
ax2[:plot](xk_model, Be_KNS_error[:,2,2],color="green",linestyle="--",linewidth=0.5)

ax2[:yaxis][:set_ticks_position]("right")
ax2[:yaxis][:set_label_position]("right")
ax2[:set_ylabel](L"B$_e$, J mol$^{-1}$",fontsize = 16, fontname="Arial")

ax2[:set_ylim]([40000,160000])
#ax2[:set_xlim]([0.,40.])

ax2[:annotate](L"80 mol% SiO$_2$",xy=(0.6,145000),xycoords="data",fontsize = 16, fontname="Arial",color="black",ha="center")
ax2[:annotate](L"66 mol% SiO$_2$",xy=(0.6,100000),xycoords="data",fontsize = 16, fontname="Arial",color="green",ha="center")

ax2[:annotate]("B)",xy=(0.8,0.9),xycoords="axes fraction",fontsize = 20, fontname="Arial")

show(ax2)

savefig("../figures/Figure8.pdf")

################################################################################
#### FIGURE 9 WITH INSERT
################################################################################

# Tg data (see supplementary dataset)
Na_Tg_data_x = [100.;85.;82.5;80.;75.;73.;70.;66.6;65.;60.] #SiO2 in Na2O-SiO2
Na_Tg_data_y = [1480.;777.;762.;753.;737.;737.;728.;726.;714.;698.5] #associated Tg

K_Tg_data_x = [100.;92.;89.;82.;80.;75.;71.;66.;60.] #SiO2 in K2O-SiO2
K_Tg_data_y = [1480.;828.;792.;769.;760.;750.;727.;721.;698.] #associated Tg

NaK_Tg_data = readcsv("../data/tg_kns_mix.csv",skipstart=1) # for mixed composition

#### WE FIRST CALCULATE THE VALUES RETURNED BY THE MODEL TO COMPARE THEM DIRECTLY WITH DATA

# generating arrays
Tg_NS_calc = ones(size(Na_Tg_data_x,1),3)
Tg_KS_calc = ones(size(K_Tg_data_x,1),3)
Tg_KNS_calc = ones(size(NaK_Tg_data,1),3)

# calculating the Tgs as well as their associated errors with the appropriate functions.
for i = 1:size(Na_Tg_data_x,1) # FOR NA2O-SiO2 GLASSES
    ~, Tg, ~, ~, ~, ~, ~, ~, ~, ~, ~ = model_direct(Na_Tg_data_x[i],0.0,[1000],params_withquantile[:,1],spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    ~, ese_Tg, ~, ~,~ = error_propagation(params_best,params_boot,Na_Tg_data_x[i],0.0,1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    Tg_NS_calc[i,1] = Tg
    Tg_NS_calc[i,2:3] = ese_Tg
end

for i = 1:size(K_Tg_data_x,1) # FOR K2O-SiO2 GLASSES
    ~, Tg, ~, ~, ~, ~, ~, ~, ~, ~, ~ = model_direct(K_Tg_data_x[i],1.0,[1000],params_withquantile[:,1],spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    ~, ese_Tg, ~, ~,~ = error_propagation(params_best,params_boot,K_Tg_data_x[i],1.0,1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    Tg_KS_calc[i,1] = Tg
    Tg_KS_calc[i,2:3] = ese_Tg
end

for i = 1:size(NaK_Tg_data,1) # FOR MIXED GLASSES
    ~, Tg, ~, ~, ~, ~, ~, ~, ~, ~, ~ = model_direct(NaK_Tg_data[i,1],NaK_Tg_data[i,2],[1000],params_withquantile[:,1],spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
    ~, ese_Tg, ~, ~,~ = error_propagation(params_best,params_boot,NaK_Tg_data[i,1],NaK_Tg_data[i,2],1000.,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    Tg_KNS_calc[i,1] = Tg
    Tg_KNS_calc[i,2:3] = ese_Tg
end

# NOW WE DO THE FIGURE
fig_paper2 = figure(figsize=(10,8))
PyPlot.rc("xtick", labelsize=14)
PyPlot.rc("ytick", labelsize=14)
plot(100.0-sio2_sc,tg_NS,color="r",linestyle="--",linewidth=2.5,label=L"Na$_2$O-SiO$_2$")
plot(100.0 - sio2_sc, tg_NS_error[:,1],color="r",linestyle="--",linewidth=0.5)
plot(100.0 - sio2_sc, tg_NS_error[:,2],color="r",linestyle="--",linewidth=0.5)

plot(100.0-sio2_sc,tg_KS,color="b",linestyle="-",linewidth=2.5,label=L"K$_2$O-SiO$_2$")
plot(100.0 - sio2_sc, tg_KS_error[:,1],color="b",linestyle="-",linewidth=0.5)
plot(100.0 - sio2_sc, tg_KS_error[:,2],color="b",linestyle="-",linewidth=0.5)

ax = gca()
xlabel(L"Alkali content, as mol% M$_2$O in SiO$_2$", ha="center",fontsize = 16,fontname="Arial")

# We generated a second Y axis to be clear about Tg from model and Tg values
ax[:yaxis][:set_ticks_position]("left")
ax[:yaxis][:set_label_position]("left")
ylabel(L"B$_e$/S$^{conf}(T_{g})$ * 1/(12-A$_e$) = $T_{g}$, K",fontsize = 16, fontname="Arial")
legend(frameon=false,numpoints=1,bbox_to_anchor=(1, 1))

# we plot the data
ax2 = ax[:twinx]()
ax2[:plot](100.0.-Na_Tg_data_x, Na_Tg_data_y,marker="s",markersize=8,color="r",linestyle="") # NA
ax2[:plot](100.0.-K_Tg_data_x, K_Tg_data_y,marker="D",markersize=8,color="b",linestyle="") # K

# setting up some stuff for legend and limits
ax2[:set_ylabel](L"$T_g$, K",fontsize = 16, fontname="Arial")
ax[:set_ylim]([500,1600])
ax[:set_xlim]([0,40])
ax2[:set_xlim](ax[:get_xlim]())
ax2[:set_ylim](ax[:get_ylim]())

# we generate a top X axis showing the mol% M4O2
ax2_b = ax[:twiny]()
ax2_b[:set_xlim](ax[:get_xlim]())
#already defined upstairs => new_tick_locations=[9.53, 18.18, 26.08, 33.33, 40.]
ax2_b[:set_xticks](new_tick_locations)
ax2_b[:set_xticklabels](tick_function(new_tick_locations))
ax2_b[:set_xlabel](L"Alkali content, as mol% M$_4$O$_2$ in SiO$_2$", ha="center",fontsize = 16,fontname="Arial")

savefig("../figures/Figure9.pdf")

# NOW WE DO ANOTHER FIGURE, SMALL, WILL BE AN INSERT FOR THE PAST FIGURE
fig_paper3 = figure(figsize=(6,4))

# Plotting
errorbar(Na_Tg_data_y[:],Tg_NS_calc[:,1],yerr=Tg_NS_calc[:,3]-Tg_NS_calc[:,1],linestyle="",marker="s",color="r")
errorbar(K_Tg_data_y[:],Tg_KS_calc[:,1],yerr=Tg_KS_calc[:,3]-Tg_KS_calc[:,1],linestyle="",marker="D",color="b")
errorbar(NaK_Tg_data[:,3],Tg_KNS_calc[:,1],yerr=Tg_KNS_calc[:,3]-Tg_KNS_calc[:,1],linestyle="",marker="o",color = "g")

plot([500,1500],[500,1500],linestyle="--",color="black")

xlim([500,1000])
ylim([500,1000])

xlabel(L"Measured viscous T$_g$",fontsize = 16,fontname="Arial")
ylabel(L"Calculated T$_g$",fontsize = 16,fontname="Arial")

tight_layout()

savefig("../figures/Figure9_insert.pdf")

std_tg_na = rmsd(Tg_NS_calc[2:end,1],Na_Tg_data_y[2:end]) # excluding sio2 in the rmsd calculation
std_tg_k = rmsd(Tg_KS_calc[2:end,1],K_Tg_data_y[2:end]) # excluding sio2 in the rmsd calculation
std_tg_nak = rmsd(Tg_KNS_calc[:,1],NaK_Tg_data[:,3])

println("The std for Tg of Na melts is $(std_tg_na), \n that for Tg of K melts is $(std_tg_k), and that for Tg of Na-K melts is $(std_tg_nak).")

################################################################################
#### FIGURE 10
################################################################################

fig_10 = figure(figsize=(11,5))
PyPlot.rc("xtick", labelsize=14)
PyPlot.rc("ytick", labelsize=14)

fig_10[:text](0.30, 0.97, L"mol% M$_4$O$_2$ in SiO$_2$", ha="center",fontsize = 16,fontname="Arial")

# CONTRIBUTIONS KS
ax12 = subplot(1,2,1)
ax12[:plot](100.0-sio2_sc,entropy_KS_lin,color="b",linestyle="--",linewidth=2.0)
ax12[:plot](100.0-sio2_sc,entropy_KS_QnMix,color="b",linestyle="-",linewidth=2.0)

# CONTRIBUTIONS NS
ax12[:plot](100.0-sio2_sc,entropy_NS_lin,color="r",linestyle="--",linewidth=2.0)
ax12[:plot](100.0-sio2_sc,entropy_NS_QnMix,color="r",linestyle="-",linewidth=2.0)
ax12[:set_ylabel](L"J mol$^{-1}$ K$^{-1}$",fontsize = 16, fontname="Arial")
ax12[:set_ylim]([0,10])

ax12[:set_xlabel]( L"mol% M$_2$O in SiO$_2$", ha="center",fontsize = 16,fontname="Arial")

ax12[:annotate](L"A) SiO$_2$-M$_2$O binary",xy=(0.05,0.9),xycoords="axes fraction",fontsize = 20, fontname="Arial")

# An additional X axis
ax12_b = ax12[:twiny]()
ax12_b[:set_xlim](ax12[:get_xlim]())
#already defined upstairs => new_tick_locations=[9.53, 18.18, 26.08, 33.33, 40.]
ax12_b[:set_xticks](new_tick_locations)
ax12_b[:set_xticklabels](tick_function(new_tick_locations))

#### DURING Na-K MIXING
ax13 = subplot(1,2,2)
ax13[:plot](xk_model,entropy_KNS_lin[:,1],color="black",linestyle = "--",linewidth = 2.0)
ax13[:plot](xk_model,entropy_KNS_QnMix[:,1],color="black",linestyle = "-",linewidth = 2.0)
ax13[:plot](xk_model,entropy_KNS_NaKMix[:,1],color="black",linestyle = "-.",linewidth = 2.0)

ax13[:plot](xk_model,entropy_KNS_lin[:,2],color="green",linestyle = "--",linewidth = 2.0)
ax13[:plot](xk_model,entropy_KNS_QnMix[:,2],color="green",linestyle = "-",linewidth = 2.0)
ax13[:plot](xk_model,entropy_KNS_NaKMix[:,2],color="green",linestyle = "-.",linewidth = 2.0)

ylim(0,8)

ax13[:annotate]("B) Na-K mixing",xy=(0.05,0.9),xycoords="axes fraction",fontsize = 20, fontname="Arial")
ax13[:set_xlabel]( L"$X_K$ = K/(K+Na)", ha="center",fontsize = 16,fontname="Arial")

tight_layout()

savefig("../figures/Figure10.pdf")
