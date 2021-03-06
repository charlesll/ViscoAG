################################################################################
#
# Copyright (C) Charles Le Losq, Daniel R. Neuville
#
# This code is a supplementary material of the article
# Le Losq C., Neuville D. R. (2017) Molecular structure, configurational entropy and viscosity of silicate melts: link through the Adam and Gibbs theory of viscous flow. Journal of Non Cristalline Solids XX:XXX-XXX.
#
# The following code is under GNU GPL V3 license, see the LICENSE file.
#
# This file performs the fitting of the model to the dataset. It outputs the
# model parameters in ../model_outputs/best_results.jld
#
# 4000 model parameter sets are also generated by bootstrapping the data and
# fitting the model again and again. Those results are in
# ../model_outputs/boot_results_x.jld. Several bootstrap results
#  can be glued together with the model_boot_glueresults.jl script.
#
# A switch to turn off the bootstrap calculation has been added, see line 368.
#
################################################################################

################################################################################
# LIBRARY LOADING
################################################################################

# For doing the plots, we will use the Matplotlib Julia API (see http://matplotlib.org/ and https://github.com/JuliaPy/PyPlot.jl)
using PyPlot
PyPlot.rc("xtick", labelsize=14)
PyPlot.rc("ytick", labelsize=14)

# for doing the calculations, we use the JuMP optimisation framework with the Ipopt solver
# Spline fitting is done with the Julia wrapper of the Dierckx FORTRAN library
using JuMP
using Ipopt
using Dierckx
using StatsBase
using JLD
using Distributions

# We need Spectra if we read the results from the bootstrap at the end
using Spectra

# we also include the following code that contains the functions allowing us
# to make various calculations
include("model_functions.jl")

################################################################################
# DATA IMPORT
################################################################################
# file with columns mol%SiO2, mol%Na2O, mol%K2O Tg(K) T(K) and viscosity (Pa s)
data_tot = readcsv("../data/visco_silicates_Julia.csv",skipstart=1)

# Qn distribution in K2O-SiO2 and Na2O-SiO2
qn_k = readcsv("../data/qn_ks.csv",skipstart=1)
qn_na = readcsv("../data/qn_ns.csv",skipstart=1)

################################################################################
# Qn SPLINE FITTING
################################################################################

spl_q2_na = Spline1D(qn_na[:,2],qn_na[:,3],k=4,s=20)
spl_q3_na = Spline1D(qn_na[:,2],qn_na[:,4],k=4,s=20)
spl_q4_na = Spline1D(qn_na[:,2],qn_na[:,5],k=4,s=30)

spl_q2_k = Spline1D(qn_k[:,2],qn_k[:,3],k=4,s=20)
spl_q3_k = Spline1D(qn_k[:,2],qn_k[:,4],k=4,s=20)
spl_q4_k = Spline1D(qn_k[:,2],qn_k[:,5],k=4,s=30)

# saving the splines as .jld files (HDF5 containers) for future uses
save("../data/data_model.jld","spl_q2_na", spl_q2_na, "spl_q3_na", spl_q3_na,"spl_q4_na",spl_q4_na,"spl_q2_k",spl_q2_k,"spl_q3_k",spl_q3_k,"spl_q4_k",spl_q4_k)

################################################################################
# QN UNITS AND CHEMISTRY CALCULATIONS
################################################################################

# Explicit attributions for chemical composition to simplify things
sio2 = data_tot[:,1]
na2o = data_tot[:,2]
k2o = data_tot[:,3]
xk = k2o./(k2o+na2o)

# Reading the Tg, T and viscosity data
Tg = data_tot[:,4]
t = data_tot[:,5]
n = data_tot[:,6]

# Normalization to any small deviation to 100 of the starting compositions
sio2 = sio2./(sio2+na2o+k2o)*100
na2o = na2o./(sio2+na2o+k2o)*100
k2o = k2o./(sio2+na2o+k2o)*100

# creating arrays that will be used for storing the calculated Qn values
q2 = ones(length(sio2),1)
q3 = ones(length(sio2),1)
q4 = ones(length(sio2),1)

si = ones(length(sio2),1)
o = ones(length(sio2),1)
k = ones(length(sio2),1)
na = ones(length(sio2),1)
q2 = ones(length(sio2),1)
q2_na = ones(length(sio2),1)
q2_k = ones(length(sio2),1)
q3 = ones(length(sio2),1)
q3_na = ones(length(sio2),1)
q3_k = ones(length(sio2),1)
q4 = ones(length(sio2),1)
logQ2 = ones(length(sio2),1)
logQ2_na = ones(length(sio2),1)
logQ2_k = ones(length(sio2),1)
logQ3 = ones(length(sio2),1)
logQ3_na = ones(length(sio2),1)
logQ3_k = ones(length(sio2),1)
logQ4 = ones(length(sio2),1)
xna = ones(length(sio2),1)
logXk = ones(length(sio2),1)
logXna = ones(length(sio2),1)

# we start a loop to calculate the fraction of Qn units
# regarding the input chemical composition
for i = 1:size(sio2,1)
    si[i], o[i], k[i], na[i], q2[i], q2_na[i], q2_k[i], q3[i], q3_na[i], q3_k[i], q4[i], logQ2[i], logQ2_na[i], logQ2_k[i], logQ3[i], logQ3_na[i], logQ3_k[i], logQ4[i], xk[i], xna[i], logXk[i], logXna[i] = qn_calculation(sio2[i],na2o[i],k2o[i],spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)
end

################################################################################
# HEAT CAPACITY CALCULATION
################################################################################
Cpg, aCpl, bCpl, ap, b = cp_model(sio2,k2o,na2o,Tg)

################################################################################
# PROBLEM DEFINITION AND OPTIMIZATION
################################################################################

# We use the JuMP framework to write the optimization problem.
# It is solved using the Ipopt solver.
mod = Model(solver=IpoptSolver(acceptable_tol = 1e-10))
nb_data = size(t,1) # number of data points

#### Declaration of the model variables
# Kqn and Knk respectively correspond to K1 and K2 in the paper.
# I kept those names here as they are more explicit.
# Boundaries are applied but very large. Of course, we forbid negative values
# for Sc and Be.

@variable(mod, -10 <= Ae <= 3, start = -2)

@variable(mod, -10000<= Kqn <= 100000 , start = 1)
@variable(mod, -100000 <= Knk <= 1000000, start = 10000 )

@variable(mod, 0 <= BeQ2_Na <= 400000, start = 10555 )
@variable(mod, 0 <= BeQ3_Na <= 400000, start = 80000)
@variable(mod, 0 <= BeQ2_K <= 400000, start = 150000)
@variable(mod, 0 <= BeQ3_K <= 400000, start = 100)
@variable(mod, 0 <= BeQ4 <= 400000, start = 150000)

@variable(mod, 0 <= ScTgQ2_Na <= 30, start = 5)
@variable(mod, 0 <= ScTgQ3_Na <= 30, start = 5)
@variable(mod, 0 <= ScTgQ2_K <= 30, start = 5)
@variable(mod, 0 <= ScTgQ3_K <= 30, start = 5)
@variable(mod, 0 <= ScTgQ4 <= 30, start = 2.0)

# The expressions need to be explicit in JuMP and need to work with scalar values, as JuMP performs automatic differentiation
@NLexpression(mod, Be[j=1:nb_data], (q2_na[j] * BeQ2_Na + q2_k[j] * BeQ2_K + q3_na[j] * BeQ3_Na + q3_k[j] * BeQ3_K + q4[j] * BeQ4) + Kqn * (-8.314  * si[j]/o[j]*2. * (logQ2[j] + logQ3[j] + logQ4[j])) + Knk * ( -8.314 * (k[j]+na[j])/o[j]*2. * (logXk[j] + logXna[j])))
@NLexpression(mod, ScTg[j=1:nb_data], (q2_na[j] * ScTgQ2_Na + q2_k[j] * ScTgQ2_K + q3_na[j] * ScTgQ3_Na + q3_k[j] * ScTgQ3_K + q4[j] * ScTgQ4) + (-8.314  * si[j]/o[j]*2. * (logQ2[j] + logQ3[j] + logQ4[j]))  + (-8.314 * (k[j]+na[j])/o[j]*2. * (logXk[j] + logXna[j])))
@NLexpression(mod, n_model[j=1:nb_data], Ae + Be[j] / (t[j] * (ScTg[j] + (ap[j] * (log(t[j])-log(Tg[j])) + b[j] * (t[j]-Tg[j])))))

# The objective function to solve, a simple least-square model
@NLobjective(mod,Min,sum{(n_model[j] - n[j])^2, j=1:nb_data})

println("Model constructed.")

################################################################################
# FIRST LOOK AT THE STARTING MODEL
################################################################################

# parameter extractions
Ae_mod = getvalue(Ae)

Kqn_mod = getvalue(Kqn)
Knk_mod = getvalue(Knk)

BeQ2_Na_mod = getvalue(BeQ2_Na)
BeQ3_Na_mod = getvalue(BeQ3_Na)
BeQ2_K_mod = getvalue(BeQ2_K)
BeQ3_K_mod = getvalue(BeQ3_K)
BeQ4_mod = getvalue(BeQ4)

ScTgQ2_Na_mod = getvalue(ScTgQ2_Na)
ScTgQ3_Na_mod = getvalue(ScTgQ3_Na)
ScTgQ2_K_mod = getvalue(ScTgQ2_K)
ScTgQ3_K_mod = getvalue(ScTgQ3_K)
ScTgQ4_mod = getvalue(ScTgQ4)

# calculation of the parameters
Be = Be_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,BeQ2_Na_mod,BeQ2_K_mod,BeQ3_Na_mod,BeQ3_K_mod,BeQ4_mod,Kqn_mod,Knk_mod)
ScTg, ~,~,~ = Sc_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,ScTgQ2_Na_mod,ScTgQ2_K_mod,ScTgQ3_Na_mod,ScTgQ3_K_mod,ScTgQ4_mod)
model_n = n_calculation(Ae_mod,Be,ScTg,ap,b,t,Tg)

# Plotting the initial guess model and recording it as Start_model_example1.pdf
# As it can be seen, this is far from being ideal...
fig1 = figure()
scatter(10000./t,n, color="black")
scatter(10000./t,model_n,color="red",marker="x")
xlabel(L"10000/T, K$^{-1}$", fontsize = 16, fontname="Arial")
ylabel("Viscosity, log Pa s", fontsize = 16, fontname="Arial")
#show() # activate to show the plot in a X11 window
savefig("../figures/Start_model_example1.pdf")

################################################################################
# SOLVING THE MODEL
################################################################################

# Solve for the control and state
println("Solving...")
status = solve(mod)

# Display results
println("Solver status: ", status)
println("Objective value: ", getobjectivevalue(mod))

################################################################################
# EXTRACTING THE CALCULATED PARAMETERS
################################################################################

Ae_mod = getvalue(Ae)

Kqn_mod = getvalue(Kqn)
Knk_mod = getvalue(Knk)

BeQ2_Na_mod = getvalue(BeQ2_Na)
BeQ3_Na_mod = getvalue(BeQ3_Na)
BeQ2_K_mod = getvalue(BeQ2_K)
BeQ3_K_mod = getvalue(BeQ3_K)
BeQ4_mod = getvalue(BeQ4)

ScTgQ2_Na_mod = getvalue(ScTgQ2_Na)
ScTgQ3_Na_mod = getvalue(ScTgQ3_Na)
ScTgQ2_K_mod = getvalue(ScTgQ2_K)
ScTgQ3_K_mod = getvalue(ScTgQ3_K)
ScTgQ4_mod = getvalue(ScTgQ4)

# calculation of Be, Sconf(Tg) and viscosity
Be = Be_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,BeQ2_Na_mod,BeQ2_K_mod,BeQ3_Na_mod,BeQ3_K_mod,BeQ4_mod,Kqn_mod,Knk_mod)
ScTg, ~,~,~ = Sc_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,ScTgQ2_Na_mod,ScTgQ2_K_mod,ScTgQ3_Na_mod,ScTgQ3_K_mod,ScTgQ4_mod)
model_n = n_calculation(Ae_mod,Be,ScTg,ap,b,t,Tg)

# calculation of rmsd: we use the function from StatsBase
rmsd_model = round(rmsd(n,model_n),3)

################################################################################
# FIGURE 6 OF THE PAPER
################################################################################

fig2 = figure(figsize=(9,7))

# We are going to distinguish Na2O-SiO2, K2O-SiO2 and Na2O-K2O-SiO2 data in the figure
# For doing so, we need to read the chemical composition and to plot things accordingly

subplot(121) # First subplot
for i = 1:size(t,1)
    if sio2[i] == 100.0
        scatter(10000./t[i],n[i], edgecolor="black",facecolor="none",marker="h") # pure SiO2 Data
        scatter(10000./t[i],model_n[i],color="black",marker="x") # pure SiO2 Model
    elseif na2o[i] == 0.0
        scatter(10000./t[i],n[i], edgecolor="blue",facecolor="none",marker="D") # K2O-SiO2 Data
        scatter(10000./t[i],model_n[i],color="blue",marker="x") # K2O-SiO2 Model
    elseif k2o[i] == 0.0
        scatter(10000./t[i],n[i], edgecolor="red",facecolor="none",marker="o") # K2O-SiO2 Data
        scatter(10000./t[i],model_n[i],color="red",marker="x") # K2O-SiO2 Model
    else
        scatter(10000./t[i],n[i], edgecolor="green",facecolor="none",marker="s") # K2O-SiO2 Data
        scatter(10000./t[i],model_n[i],color="green",marker="x") # K2O-SiO2 Model
    end
end

# Defining the axes limits
xlim(3,16)
ylim(0, 15)

# Setting the labels
xlabel(L"10000/T, K$^{-1}$", fontsize = 16, fontname="Arial")
ylabel("Viscosity, log Pa s", fontsize = 16, fontname="Arial")

# Some Annotation
annotate("A)", xy=(0.1, 0.9), xycoords="axes fraction",fontsize = 20, fontname="Arial")

# Creating a quick legend for indicating the model symbols
# A little trick is to create fake data with the right symbols and to use handles.
model_dt = scatter([],[],color="black",marker="x",label="Model values")
legend(loc=4,frameon=false,handles=[model_dt])

subplot(122)

for i = 1:size(t,1)
    if sio2[i] == 100.0
        scatter(n[i],model_n[i], edgecolor="black",facecolor="none",marker="h") # K2O-SiO2 Data
    elseif na2o[i] == 0.0
        scatter(n[i],model_n[i],edgecolor="blue",facecolor="none",marker="D") # K2O-SiO2
    elseif k2o[i] == 0.0
        scatter(n[i],model_n[i],edgecolor="red",facecolor="none",marker="o") # Na2O-SiO2
    else
        scatter(n[i],model_n[i],edgecolor="green",facecolor="none",marker="s") # Na2O-SiO2
    end
end

plot([0,15],[0,15],linestyle="--",color="black")

# Setting the labels
xlabel("Measured viscosity, log Pa s", fontsize = 16, fontname="Arial")
ylabel("Modeled viscosity, log Pa s", fontsize = 16, fontname="Arial")

# Defining the axes limits
xlim(0,15)
ylim(0,15)

# Some Annotation
annotate(L"$\sigma$ = "*string(rmsd_model), xy=(2., 11.), xycoords="data",fontsize = 18, fontname="Arial")
annotate("B)", xy=(0.1, 0.9), xycoords="axes fraction",fontsize = 20, fontname="Arial")

# We also need to create our custom legend
# A little trick is to create fake data with the right symbols and to use handles.
sio2_dt = scatter([],[],edgecolor="black",facecolor="none",marker="h",label=L"SiO$_2$")
k_dt = scatter([], [], edgecolor="blue",facecolor="none", marker="D", label=L"K$_2$O-SiO$_2$")
na_dt = scatter([], [], color="red", marker="o", label=L"Na$_2$O-SiO$_2$")
nak_dt = scatter([], [], color="green", marker="s", label=L"Na$_2$O-K$_2$O-SiO$_2$")

legend(loc=4,frameon=false,handles=[sio2_dt,na_dt,k_dt,nak_dt])

#show()

savefig("../figures/Figure6.pdf")

################################################################################
# PRINTING AND SAVING BEST PARAMETERS
################################################################################

println("Model standard deviation is: $(rmsd_model)")

println("Ae is $(Ae_mod)")

println("Kqn (K1) is $(Kqn_mod)")
println("Knk (K2) is $(Knk_mod)")

println("BeQ2_Na is $(BeQ2_Na_mod)")
println("BeQ3_Na is $(BeQ3_Na_mod)")
println("BeQ2_K is $(BeQ2_K_mod)")
println("BeQ3_K is $(BeQ3_K_mod)")
println("BeQ4 is $(BeQ4_mod)")

println("ScTgQ2_Na is $(ScTgQ2_Na_mod)")
println("ScTgQ3_Na is $(ScTgQ3_Na_mod)")
println("ScTgQ2_K is $(ScTgQ2_K_mod)")
println("ScTgQ3_K is $(ScTgQ3_K_mod)")
println("ScTgQ4 is $(ScTgQ4_mod)")

save("../model_outputs/best_results.jld","params",[Ae_mod Kqn_mod Knk_mod BeQ2_Na_mod BeQ3_Na_mod BeQ2_K_mod BeQ3_K_mod BeQ4_mod ScTgQ2_Na_mod ScTgQ3_Na_mod ScTgQ2_K_mod ScTgQ3_K_mod ScTgQ4_mod])

################################################################################
# ERROR CALCULATIONS
################################################################################
# Errors are calculated by bootstrapping.
#
# We resample with replacement the entire viscosity dataset and then perform
# nb_boot new model optimisations.
#
# This allows to estimate the errors without any assumptions on their distribution
# See the references in the paper for further details (Efron and related references)

# The software performs 4000 bootstrap, it may take some time.
# I leave a switch there that allows to turn on or off this procedure.
# The 40000 bootstrap fits used in the paper to calculate the errors have been
# performed running 10 different sessions and gluing the results together.
# See the script model_boot_glueresults.jl
#
# Please note that the random seed was not fixed: this may give very small differences
# between different runs, negigible compared to the errors affecting the parameters.

# if you want to do the error calculation, turn boot_switch to 1
boot_switch = 0

if boot_switch == 1
  nb_boot = 4000

  params_boot = Array{Float64}(nb_boot,13)
  compteur = 100
  for i = 1:nb_boot

      #resampling data with replacement = non-parametric bootstrapping
      vect = collect(1:size(t,1))
      idx = sample(vect,size(vect,1),replace=true)

      si_b = si[idx]
      o_b = o[idx]
      k_b = k[idx]
      na_b = na[idx]
      q2_na_b = q2_na[idx]
      q2_k_b = q2_k[idx]
      q3_na_b = q3_na[idx]
      q3_k_b = q3_k[idx]
      q4_b = q4[idx]
      logXk_b = logXk[idx]
      logXna_b = logXna[idx]
      logQ2_b = logQ2[idx]
      logQ3_b = logQ3[idx]
      logQ4_b = logQ4[idx]
      t_b = t[idx]
      Tg_b = Tg[idx]
      ap_b = ap[idx]
      b_b = b[idx]
      n_b = n[idx]

      mod2 = Model(solver=IpoptSolver(print_level=0)) # to be safe we re-define the model

      # JuMP variables, expressions and objective functions
      @variable(mod2, -10 <= Ae_b <= 3, start = -1.81)

      @variable(mod2, -10000<= Kqn_b <= 100000 , start = 300)
      @variable(mod2, -100000 <= Knk_b <= 1000000, start = 6000 )

      @variable(mod2, 0 <= BeQ2_Na_b <= 400000, start = 87000 )
      @variable(mod2, 0 <= BeQ3_Na_b <= 400000, start = 61000)
      @variable(mod2, 0 <= BeQ2_K_b <= 400000, start = 103000)
      @variable(mod2, 0 <= BeQ3_K_b <= 400000, start = 104000)
      @variable(mod2, 0 <= BeQ4_b <= 400000, start = 155000)

      @variable(mod2, 0 <= ScTgQ2_Na_b <= 30, start = 5)
      @variable(mod2, 0 <= ScTgQ3_Na_b <= 30, start = 1)
      @variable(mod2, 0 <= ScTgQ2_K_b <= 30, start = 5)
      @variable(mod2, 0 <= ScTgQ3_K_b <= 30, start = 5)
      @variable(mod2, 0 <= ScTgQ4_b <= 30, start = 8)

      @NLexpression(mod2, Be_b[j=1:nb_data], (q2_na_b[j] * BeQ2_Na_b + q2_k_b[j] * BeQ2_K_b + q3_na_b[j] * BeQ3_Na_b + q3_k_b[j] * BeQ3_K_b + q4_b[j] * BeQ4_b) + Kqn_b * (-8.314  * si_b[j]/o_b[j]*2. * (logQ2_b[j] + logQ3_b[j] + logQ4_b[j])) + Knk_b * ( -8.314 * (k_b[j]+na_b[j])/o_b[j]*2 * (logXk_b[j] + logXna_b[j]) ))
      @NLexpression(mod2, ScTg_b[j=1:nb_data], (q2_na_b[j] * ScTgQ2_Na_b + q2_k_b[j] * ScTgQ2_K_b + q3_na_b[j] * ScTgQ3_Na_b + q3_k_b[j] * ScTgQ3_K_b + q4_b[j] * ScTgQ4_b) + (-8.314  * si_b[j]/o_b[j]*2. * (logQ2_b[j] + logQ3_b[j] + logQ4_b[j]))  + (-8.314 * (k_b[j]+na_b[j])/o_b[j]*2 * (logXk_b[j] + logXna_b[j])))
      @NLexpression(mod2, n_model_b[j=1:nb_data], Ae_b + Be_b[j] / (t_b[j] * (ScTg_b[j] + (ap_b[j] * (log(t_b[j])-log(Tg_b[j])) + b_b[j] * (t_b[j]-Tg_b[j])))))

      @NLobjective(mod2,Min,sum{(n_model_b[j] - n_b[j])^2, j=1:nb_data}) # The objective function to solve

      status = solve(mod2) # Solve

      if compteur == i # Display status each 100 iterations
          println("Iteration $(i), Solver status: ", status)
          compteur = compteur + 100
      end

      params_boot[i,1] = getvalue(Ae_b)
      params_boot[i,2] = getvalue(Kqn_b)
      params_boot[i,3] = getvalue(Knk_b)
      params_boot[i,4] = getvalue(BeQ2_Na_b)
      params_boot[i,5] = getvalue(BeQ3_Na_b)
      params_boot[i,6] = getvalue(BeQ2_K_b)
      params_boot[i,7] = getvalue(BeQ3_K_b)
      params_boot[i,8] = getvalue(BeQ4_b)
      params_boot[i,9] = getvalue(ScTgQ2_Na_b)
      params_boot[i,10] = getvalue(ScTgQ3_Na_b)
      params_boot[i,11] = getvalue(ScTgQ2_K_b)
      params_boot[i,12] = getvalue(ScTgQ3_K_b)
      params_boot[i,13] = getvalue(ScTgQ4_b)
  end

  save("../model_outputs/boot_results_10.jld","params_boot",params_boot)

  println("Done.")

end
