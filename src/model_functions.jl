################################################################################
#
# Copyright (C) Charles Le Losq, Daniel R. Neuville
#
# This code is a supplementary material of the article
# Le Losq C., Neuville D. R. (2017) Molecular structure, configurational entropy and viscosity of silicate melts: link through the Adam and Gibbs theory of viscous flow. Journal of Non Cristalline Solids XX:XXX-XXX.
#
# The following code is under GNU GPL V3 license, see the LICENSE file.
#
# This file defines several functions that are used in the model and for
# making the figures.
#
################################################################################

"""
Function cp_model(sio2,k2o,na2o,Tg)

calculates the glass and liquid heat capacity from the models of Richet, Chem. Geol. 62, 1987, 111-124 & Richet and Bottinga, Geochem. Cosmochem. Acta 49, 1985, 471.

you can either provide all the inputs as Floats or as Arrays, to have the corresponding outputs.

Inputs
  sio2: Float64 OR Array{Float64}, the sio2 content in mol%;
  k2o: Float64 OR Array{Float64}, the k2o content in mol%;
  na2o: Float64 OR Array{Float64}, the na2o content in mol%;
  Tg: Float64 OR Array{Float64}, the glass transition temperature content in kelvin;

Outputs
  Cpg: Float64 OR Array{Float64}, the glass heat capacity at Tg, in J mol^-1;
  aCpl: Float64 OR Array{Float64}, the constant term of the liquid heat capacity;
  bCpl: Float64 OR Array{Float64}, the temperature-dependent term of the liquid heat capacity;
  ap: Float64 OR Array{Float64}, the constant term of the configurational heat capacity;
  b: Float64 OR Array{Float64}, the temperature-dependent term of the configurational heat capacity.
"""
function cp_model(sio2,k2o,na2o,Tg)
    Cpg = sio2./100.0.*(127.2 - 0.010777.*Tg + 431270.0./Tg.^2 -1463.9./Tg.^0.5) + k2o./100.0.*(84.323 +0.000731.*Tg -829800.0./Tg.^2) + na2o./100.0.*(70.884 +0.02611.*Tg -358200.0./Tg.^2)

    aCpl = 81.37.*sio2./100 + 100.6.*na2o./100 + 50.13.*k2o./100 + sio2./100.0.*(k2o./100.0.*k2o./100).*151.7
    bCpl = 0.01578.*k2o./100

    ap = aCpl - Cpg
    b = bCpl

    return Cpg, aCpl, bCpl, ap, b
end

"""
Function Be_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,BeQ2_Na_mod,BeQ2_K_mod,BeQ3_Na_mod,BeQ3_K_mod,BeQ4_mod,Kqn_mod,Knk_mod)

calculates the Be term of the Adam and Gibbs equation.

Inputs
  q2_na: Float64 OR Array{Float64}, the fraction of Q2(Naenv) in the glass;
  q2_k: Float64 OR Array{Float64}, the fraction of Q2(Kenv) in the glass;
  q3_na: Float64 OR Array{Float64}, the fraction of Q3(Naenv) in the glass;
  q3_k: Float64 OR Array{Float64}, the fraction of Q3(Kenv) in the glass;
  q4: Float64 OR Array{Float64}, the fraction of Q4 in the glass;
  si: Float64 OR Array{Float64}, the fraction of Si in the glass;, at%
  o: Float64 OR Array{Float64}, the fraction of O in the glass;, at%
  na: Float64 OR Array{Float64}, the fraction of Na in the glass;, at%
  k: Float64 OR Array{Float64}, the fraction of K in the glass;, at%
  logQ2: Float64 OR Array{Float64}, the logarithm of the fraction of Q2;
  logQ3: Float64 OR Array{Float64}, the logarithm of the fraction of Q3;
  logQ4: Float64 OR Array{Float64}, the logarithm of the fraction of Q4;
  logXk: Float64 OR Array{Float64}, the logarithm of the K/(K+Na) ratio;
  logXna: Float64 OR Array{Float64}, the logarithm of the Na/(K+Na) ratio;
  BeQ2_Na_mod: Float64 OR Array{Float64}, the partial molar Be term associated to Q2(Naenv) units;
  BeQ2_K_mod: Float64 OR Array{Float64}, the partial molar Be term associated to Q2(Kenv) units;
  BeQ3_Na_mod: Float64 OR Array{Float64}, the partial molar Be term associated to Q3(Naenv) units;
  BeQ3_K_mod: Float64 OR Array{Float64}, the partial molar Be term associated to Q3(Kenv) units;
  BeQ4_mod: Float64 OR Array{Float64}, the partial molar Be term associated to Q4 units;
  Kqn_mod: Float64 OR Array{Float64}, the K1 scaling term, see Le Losq and Neuville, JNCS 2017;
  Knk_mod: Float64 OR Array{Float64}, the K2 scaling term, see Le Losq and Neuville, JNCS 2017.

Outputs
  Be: the constant of the Adam and Gibbs equation proportional to the energy barriers opposed to viscous flow, J mol^-1.
"""
function Be_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,BeQ2_Na_mod,BeQ2_K_mod,BeQ3_Na_mod,BeQ3_K_mod,BeQ4_mod,Kqn_mod,Knk_mod)
    return (q2_na .* BeQ2_Na_mod + q2_k .* BeQ2_K_mod + q3_na .* BeQ3_Na_mod + q3_k .* BeQ3_K_mod+ q4 .* BeQ4_mod) + Kqn_mod .* (-8.314  .* si./o .*2. .* (logQ2+logQ3+logQ4)) + Knk_mod .* ( -8.314 .* (k+na)./o .*2 .*  (logXk+logXna))
end

"""
Function Sc_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,ScTgQ2_Na_mod,ScTgQ2_K_mod,ScTgQ3_Na_mod,ScTgQ3_K_mod,ScTgQ4_mod)

calculates the Be term of the Adam and Gibbs equation.

Inputs
  q2_na: Float64 OR Array{Float64}, the fraction of Q2(Naenv) in the glass;
  q2_k: Float64 OR Array{Float64}, the fraction of Q2(Kenv) in the glass;
  q3_na: Float64 OR Array{Float64}, the fraction of Q3(Naenv) in the glass;
  q3_k: Float64 OR Array{Float64}, the fraction of Q3(Kenv) in the glass;
  q4: Float64 OR Array{Float64}, the fraction of Q4 in the glass;
  si: Float64 OR Array{Float64}, the fraction of Si in the glass, at%;
  o: Float64 OR Array{Float64}, the fraction of O in the glass, at%;
  na: Float64 OR Array{Float64}, the fraction of Na in the glass, at%;
  k: Float64 OR Array{Float64}, the fraction of K in the glass, at%;
  logQ2: Float64 OR Array{Float64}, the logarithm of the fraction of Q2;
  logQ3: Float64 OR Array{Float64}, the logarithm of the fraction of Q3;
  logQ4: Float64 OR Array{Float64}, the logarithm of the fraction of Q4;
  logXk: Float64 OR Array{Float64}, the logarithm of the K/(K+Na) ratio;
  logXna: Float64 OR Array{Float64}, the logarithm of the Na/(K+Na) ratio;
  ScTgQ2_Na_mod: Float64 OR Array{Float64}, the partial molar Scong(Tg) term associated to Q2(Naenv) units;
  ScTgQ2_K_mod: Float64 OR Array{Float64}, the partial molar Scong(Tg) term associated to Q2(Kenv) units;
  ScTgQ3_Na_mod: Float64 OR Array{Float64}, the partial molar Scong(Tg) term associated to Q3(Naenv) units;
  ScTgQ3_K_mod: Float64 OR Array{Float64}, the partial molar Scong(Tg) term associated to Q3(Kenv) units;
  ScTgQ4_mod: Float64 OR Array{Float64}, the partial molar Scong(Tg) term associated to Q4 units;

Outputs
  ScTg: the configurational entropy at the glass transition temperature Tg, J mol^-1 K^-1;
  Sc_lin: the topological contribution to ScTg;
  Sc_qn: the chemical contribution due to mixing of Si between Q2, Q3 qnd Q4 units;
  Sc_alka: the chemical contribution due to mixing of Na and K.
"""
function Sc_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,ScTgQ2_Na_mod,ScTgQ2_K_mod,ScTgQ3_Na_mod,ScTgQ3_K_mod,ScTgQ4_mod)
    Sc_lin = (q2_na .* ScTgQ2_Na_mod + q2_k .* ScTgQ2_K_mod + q3_na .* ScTgQ3_Na_mod + q3_k .* ScTgQ3_K_mod + q4 .* ScTgQ4_mod)
    Sc_qn = (-8.314  .* si./o.*2. .* (logQ2+logQ3+logQ4))
    Sc_alka = (-8.314 .* ((k+na)./o) .*2 .*  (logXk+logXna))

    ScTg =  Sc_lin+ Sc_qn + Sc_alka

    return ScTg, Sc_lin, Sc_qn, Sc_alka
end

"""
Function n_calculation(Ae,Be,ScTg,ap,b,T,Tg)

Inputs
  Ae: Float64 OR Array{Float64}, the pre-exponential term of the AG equation, Pa s;
  Be: Float64 OR Array{Float64}, the Be term of the AG equation, J mol^-1;
  ScTg: Float64 OR Array{Float64}, the configurational entropy at Tg, J mol^-1 K^-1;
  ap: Float64 OR Array{Float64}, the constant term of the configurational heat capacity;
  b: Float64 OR Array{Float64}, the temperature-dependent term of the configurational heat capacity;
  T: Float64 OR Array{Float64}, the melt temperature, K;
  Tg: Float64 OR Array{Float64}, the glass transition temperature, K;

Outputs
  viscosity in Pa s, Float64 OR Array{Float64}.
"""
function n_calculation(Ae,Be,ScTg,ap,b,T,Tg)
    return Ae .+ Be ./ (T .* (ScTg .+ (ap .* (log(T).-log(Tg)) .+ b .* (T.-Tg))))
end

"""
Function Tg_calculation(Ae,Be,ScTg)

calculates the glass transition temperature from the Ae, Be and ScTg inputs.

Inputs
  Ae: Float64 OR Array{Float64}, the pre-exponential term of the AG equation, Pa s;
  Be: Float64 OR Array{Float64}, the Be term of the AG equation, J mol^-1;
  ScTg: Float64 OR Array{Float64}, the configurational entropy at Tg, J mol^-1 K^-1;

Outputs
  The glass transition temperature Tg in K, Float64 OR Array{Float64}.
"""
function Tg_calculation(Ae,Be,ScTg)
    return Be./((12-Ae).*ScTg)
end

"""
Function qn_calculation(sio2::Float64,na2o::Float64,k2o::Float64,spl_q2_na::Spline1D,spl_q3_na::Spline1D,spl_q4_na::Spline1D,spl_q2_k::Spline1D,spl_q3_k::Spline1D,spl_q4_k::Spline1D)

calculates the structure of the melt in term of Qn units, as well as several other chemical parameters used in the model.

warning: this functions works only with single Float64 for sio2, na2o and k2o, and will return Float64 values.

Inputs
  sio2: Float64, the mol% of sio2 in the melt;
  na2o: Float64, the mol% of na2o in the melt;
  k2o: Float64,  the mol% of k2o in the melt;
  spl_q2_na: Spline1D, contains the cubic spline describing the Q2 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q3_na: Spline1D, contains the cubic spline describing the Q3 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q4_na: Spline1D, contains the cubic spline describing the Q4 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q2_k: Spline1D, contains the cubic spline describing the Q2 vs % K2O relationship in K2O-SiO2 glasses;
  spl_q3_k: Spline1D, contains the cubic spline describing the Q3 vs % K2O relationship in K2O-SiO2 glasses;
  spl_q4_k: Spline1D, contains the cubic spline describing the Q4 vs % K2O relationship in K2O-SiO2 glasses;

Outputs
  si: Float64, the fraction of Si in the glass, at%;
  o: Float64, the fraction of O in the glass, at%;
  k: Float64, the fraction of K in the glass, at%;
  na: Float64, the fraction of Na in the glass, at%;
  q2: Float64, the fraction of Q2 in the glass;
  q2_na: Float64, the fraction of Q2(Naenv) in the glass;
  q2_k: Float64, the fraction of Q2(Kenv) in the glass;
  q3: Float64, the fraction of Q3 in the glass;
  q3_na: Float64, the fraction of Q3(Naenv) in the glass;
  q3_k: Float64, the fraction of Q3(Kenv) in the glass;
  q4: Float64, the fraction of Q4 in the glass;
  logQ2: Float64, the logarithm of the fraction of Q2;
  logQ2_na: Float64, the logarithm of the fraction of Q2(Naenv);
  logQ2_k: Float64, the logarithm of the fraction of Q2(Kenv);
  logQ3: Float64, the logarithm of the fraction of Q3;
  logQ3_na: Float64, the logarithm of the fraction of Q3(Naenv);
  logQ3_k: Float64, the logarithm of the fraction of Q3(Kenv);
  logQ4: Float64, the logarithm of the fraction of Q4;
  xk: Float64, the K/(K+Na) ratio;
  xna: Float64, the Na/(K+Na) ratio;
  logXk: Float64, the logarithm of the K/(K+Na) ratio;
  logXna: Float64, the logarithm of the Na/(K+Na) ratio.
"""
function qn_calculation(sio2::Float64,na2o::Float64,k2o::Float64,spl_q2_na::Spline1D,spl_q3_na::Spline1D,spl_q4_na::Spline1D,spl_q2_k::Spline1D,spl_q3_k::Spline1D,spl_q4_k::Spline1D)

    # To be safe, re-normalisation to any deviation to 100
    sio2 = sio2./(sio2+na2o+k2o)*100
    na2o = na2o./(sio2+na2o+k2o)*100
    k2o = k2o./(sio2+na2o+k2o)*100
    alka_t = na2o + k2o # total alkali content

    # calcul of unitary formula and chemical atom molar fractions
    si::Float64 = sio2
    na::Float64 = 2.*na2o
    k::Float64 = 2.*k2o
    o::Float64 = 2.*sio2 + na2o + k2o
    totions::Float64 = si+na+k+o

    fraction_qn::Float64 = (si+o)./(si+na+k+o)
    xk::Float64 = k./(na+k)
    # we need to be careful as the above will yield NaN values when K2O and Na2O will be both equal to 0
    # so we add a safeguard there
    if isnan(xk) == true
      xk =0.0
    end
    xna::Float64 = 1-xk

    # we calculate the fraction of Q2, Q3 and Q4 units for the sodic endmember
    q3_na_peak::Float64 = evaluate(spl_q3_na,alka_t)
    q2_na_peak::Float64 = evaluate(spl_q2_na,alka_t)
    q4_na_peak::Float64 = 100.0 - q3_na_peak - q2_na_peak # q4 by difference

    # now the same for the potassic endmember
    q3_k_peak::Float64 = evaluate(spl_q3_k,alka_t)
    q2_k_peak::Float64 = evaluate(spl_q2_k,alka_t)
    q4_k_peak::Float64 = 100.0 - q3_k_peak - q2_k_peak # q4 by difference

    # we need to define some rules as the Q2 splines are behaving not perfectly for low alkali contents (they do not tend to 0 as they can explore negative values).
    # For K melts, we assume regarding the existing data that there is no more Q2 below 20 mol% K2O and a smooth transition toward the spline value between 20 and 25% K2O.
    if alka_t < 20.0
        q2_k_peak = 0.0
    elseif alka_t < 25.
        q2_k_peak = (evaluate(spl_q2_k,25.0)-0.0)/5.0 .* (alka_t-20.)
    else # If [SiO2] < 75 mol%, we use the Q2 spline to calculate the fraction of Q2 units
        q2_k_peak = evaluate(spl_q2_k,alka_t)
    end

    # For Na melts, we assume regarding the existing data that there is no more Q2 below 10 mol% Na2O and a smooth transition toward the spline value between 10 and 15% Na2O.
    if alka_t < 10.0
        q2_na_peak = 0.0
    elseif alka_t < 15.
        q2_na_peak = (evaluate(spl_q2_na,15.0)-0.0)/5.0 .* (alka_t-10.)
    else # If [SiO2] <= 85 mol%, we simply use the Q2 spline to calculate the fraction of Q2 units
        q2_na_peak = evaluate(spl_q2_na,alka_t)
    end

    # Now we calculate the real fraction of Q2, Q3 and Q4 units in the system, with assuming a linear mixing between the Na and K endmembers for
    # calculating the fractions of Q2 and Q3 units. This is a simplification that has been made following the Raman analysis (see paper discussion)
    q2::Float64 = (xk * q2_k_peak + (1-xk) * q2_na_peak)/100.0
    q3::Float64 = (xk * q3_k_peak + (1-xk) * q3_na_peak)/100.0
    q4::Float64 = (xk * q4_k_peak + (1-xk) * q4_na_peak)/100.0

    # Final normalisation of the fractions for safety
    q2 = q2/(q2+q3+q4)
    q3 = q3/(q2+q3+q4)
    q4 = q4/(q2+q3+q4)

    # In the case of pure SiO2 the above does not work because of NaN values, so we define the qn values here
    if sio2 == 100
        q2 = 0.0
        q3 = 0.0
        q4 = 1.0 # q2 and q3 stay at zero, as well as na2o and k2o
    end

    # we calculate the fractions of Q2Naenv and Q2Kenv
    q2_na::Float64 = xna * q2
    q2_k::Float64  = xk * q2

    q3_na::Float64 = xna * q3
    q3_k::Float64  = xk * q3

    # we calculate the xi log(xi) values, with being careful with 0 values...
    q2_na > 0.0 ? logQ2_na::Float64 = q2_na*log(q2_na) : logQ2_na = 0.0
    q3_na > 0.0 ? logQ3_na::Float64 = q3_na*log(q3_na) : logQ3_na = 0.0

    q2_k > 0.0 ? logQ2_k::Float64 = q2_k*log(q2_k) : logQ2_k = 0.0
    q3_k > 0.0 ? logQ3_k::Float64 = q3_k*log(q3_k) : logQ3_k = 0.0

    q2 > 0.0 ? logQ2::Float64 = q2 * log(q2) : logQ2 = 0.0
    q3 > 0.0 ? logQ3::Float64 = q3 * log(q3) : logQ3 = 0.0
    q4 > 0.0 ? logQ4::Float64 = q4 * log(q4) : logQ4 = 0.0

    xk > 0.0 ? logXk::Float64 = xk * log(xk) : logXk = 0.0
    xna > 0.0 ? logXna::Float64 = xna * log(xna) : logXna = 0.0

    # returning the results
    return si, o, k, na, q2, q2_na, q2_k, q3, q3_na, q3_k, q4, logQ2, logQ2_na, logQ2_k, logQ3, logQ3_na, logQ3_k, logQ4, xk, xna, logXk, logXna
end

"""
function model_direct(sio2::Float64,xk::Float64,T::AbstractArray,params_withstuff::AbstractArray,spl_q2_na::Spline1D,spl_q3_na::Spline1D,spl_q4_na::Spline1D,spl_q2_k::Spline1D,spl_q3_k::Spline1D,spl_q4_k::Spline1D;Tg::Float64=0.0)

forward calculation of the viscosity with using the model.

warning: this functions works only with single Float64 for sio2 and xk. The calcul is done for a given composition.

Inputs
  sio2: Float64, the SiO2 content of the melt in mol%;
  xk: Float64, the K/(K+Na) ratio of the melt, calculated with K2O and Na2O in mol%;
  T: AbstractArray, temperatures at which the values are wanted;
  params: 1D AbstractArray, contains the parameters with their 0.05 and 0.95 quantile values in first, second and third columns, respectively;
  spl_q2_na: Spline1D, contains the cubic spline describing the Q2 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q3_na: Spline1D, contains the cubic spline describing the Q3 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q4_na: Spline1D, contains the cubic spline describing the Q4 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q2_k: Spline1D, contains the cubic spline describing the Q2 vs % K2O relationship in K2O-SiO2 glasses;
  spl_q3_k: Spline1D, contains the cubic spline describing the Q3 vs % K2O relationship in K2O-SiO2 glasses;
  spl_q4_k: Spline1D, contains the cubic spline describing the Q4 vs % K2O relationship in K2O-SiO2 glasses;

Optional inputs
  Tg: Float64=0.0, the glass transition temperature in K can be provided;

Outputs
  model_n: Float64, the melt viscosity, Pa s;
  Tg: Float64, the melt Tg , K; if provided as optional input, it corresponds to the provided value, otherwise it corresponds to the value calculated using the model;
  ScTg_opt: Float64, the configurational entropy at the glass transition temperature, J mol^-1 K^-1;
  Sc_lin: Float64, the topological contribution to ScTg;
  Sc_qn: Float64, the chemical contribution to ScTg due to mixing of Si between Q2, Q3 qnd Q4 units;
  Sc_alka: Float64, the chemical contribution to ScTg due to mixing of Na and K;
  Cpc_Tg: Float64, the configurational heat capacity at the glass transition temperature, J mol-1 K-1;
  Be_opt: Float64, the constant of the Adam and Gibbs equation proportional to the energy barriers opposed to viscous flow, J mol^-1;
  q2: Float64, the fraction of Q2 species in the melt;
  q3: Float64, the fraction of Q3 species in the melt;
  q4: Float64, the fraction of Q4 species in the melt;
"""
function model_direct(sio2::Float64,xk::Float64,T::AbstractArray,params::AbstractArray,spl_q2_na::Spline1D,spl_q3_na::Spline1D,spl_q4_na::Spline1D,spl_q2_k::Spline1D,spl_q3_k::Spline1D,spl_q4_k::Spline1D;Tg::Float64=0.0)

    # calcul of nao and k2o quantities
    nak2o::Float64 = 100.-sio2 # the total sum of the alkalis
    k2o::Float64 = xk*nak2o
    na2o::Float64 = 100.-sio2-k2o

    # using qn_calculation for calculating the structure
    si, o, k, na, q2, q2_na, q2_k, q3, q3_na, q3_k, q4, logQ2, logQ2_na, logQ2_k, logQ3, logQ3_na, logQ3_k, logQ4, xk, xna, logXk, logXna = qn_calculation(sio2,na2o,k2o,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    # parameter extractions
    Ae_mod = params[1]

    Kqn_mod = params[2]
    Knk_mod = params[3]

    BeQ2_Na_mod = params[4]
    BeQ3_Na_mod = params[5]
    BeQ2_K_mod = params[6]
    BeQ3_K_mod = params[7]
    BeQ4_mod = params[8]

    ScTgQ2_Na_mod = params[9]
    ScTgQ3_Na_mod = params[10]
    ScTgQ2_K_mod = params[11]
    ScTgQ3_K_mod = params[12]
    ScTgQ4_mod = params[13]

    Be_opt = Be_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,BeQ2_Na_mod,BeQ2_K_mod,BeQ3_Na_mod,BeQ3_K_mod,BeQ4_mod,Kqn_mod,Knk_mod)
    ScTg_opt, Sc_lin, Sc_qn, Sc_alka = Sc_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,ScTgQ2_Na_mod,ScTgQ2_K_mod,ScTgQ3_Na_mod,ScTgQ3_K_mod,ScTgQ4_mod)

    # In case Tg is not provided, we calculate it
    if Tg == 0
        Tg = Tg_calculation(Ae_mod,Be_opt,ScTg_opt)
    end

    # Calculation of Cp
    Cpg, aCpl, bCpl, ap, b = cp_model(sio2,k2o,na2o,Tg) # glass and melt Cp
    Cpc_Tg::Float64 = ap + b.*Tg # configuration heat capacity at Tg

    model_n = n_calculation(Ae_mod,Be_opt,ScTg_opt,ap,b,T,Tg)

    return model_n, Tg, ScTg_opt, Sc_lin, Sc_qn, Sc_alka, Cpc_Tg, Be_opt, q2, q3, q4
end

"""
Function error_propagation(params_best::AbstractArray,params_boot::AbstractArray,sio2::Float64,xk::Float64,T::Float64,spl_q2_na::Spline1D,spl_q3_na::Spline1D,spl_q4_na::Spline1D,spl_q2_k::Spline1D,spl_q3_k::Spline1D,spl_q4_k::Spline1D;Tg_input::Float64=0.0)

forward calculation with error prediction of the melt viscosity, entropy, and Be terms.

warning: this functions works only with single Float64 for sio2, xk and T. The calcul is done for a given composition and temperature.

Inputs
  params_best: AbstractArray, contains the best parameters;
  params_boot: AbstractArray, contains the parameters in columns from the n bootstrap fits (n = number of lines);
  sio2: Float64, the SiO2 content of the melt in mol%;
  xk: Float64, the K/(K+Na) ratio of the melt, calculated with K2O and Na2O in mol%;
  T: AbstractArray, temperatures at which the values are wanted. Warning: T is a Float64 number in this function, contrary to the above model_direct function that expects an Array of temperature;
  spl_q2_na: Spline1D, contains the cubic spline describing the Q2 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q3_na: Spline1D, contains the cubic spline describing the Q3 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q4_na: Spline1D, contains the cubic spline describing the Q4 vs % Na2O relationship in Na2O-SiO2 glasses;
  spl_q2_k: Spline1D, contains the cubic spline describing the Q2 vs % K2O relationship in K2O-SiO2 glasses;
  spl_q3_k: Spline1D, contains the cubic spline describing the Q3 vs % K2O relationship in K2O-SiO2 glasses;
  spl_q4_k: Spline1D, contains the cubic spline describing the Q4 vs % K2O relationship in K2O-SiO2 glasses;

Optional inputs
  Tg_input: Float64=0.0, the glass transition temperature in Kelvin can be provided;

Outputs

"""
function error_propagation(params_best::AbstractArray,params_boot::AbstractArray,sio2::Float64,xk::Float64,T::Float64,spl_q2_na::Spline1D,spl_q3_na::Spline1D,spl_q4_na::Spline1D,spl_q2_k::Spline1D,spl_q3_k::Spline1D,spl_q4_k::Spline1D;Tg_input::Float64=0.0)

    # calcul of nao and k2o quantities
    nak2o::Float64 = 100.-sio2
    k2o::Float64 = xk*nak2o
    na2o::Float64 = 100.-sio2-k2o

    # calling qn_calculation for calculating the structure
    si, o, k, na, q2, q2_na, q2_k, q3, q3_na, q3_k, q4, logQ2, logQ2_na, logQ2_k, logQ3, logQ3_na, logQ3_k, logQ4, xk, xna, logXk, logXna = qn_calculation(sio2,na2o,k2o,spl_q2_na,spl_q3_na,spl_q4_na,spl_q2_k,spl_q3_k,spl_q4_k)

    # We take the parameters from the bootstrap array
    # Values are corrected from slight biases, calculated using the best value of the model
    # The generated vectors are also shuffled to remove covariance arising from the fitting protocol and promote randomness
    biases = mean(params_boot,1) - params_best

    Ae_mod = shuffle(params_boot[:,1] - biases[1])

    Kqn_mod = shuffle(params_boot[:,2] - biases[2])
    Knk_mod = shuffle(params_boot[:,3] - biases[3])

    BeQ2_Na_mod = shuffle(params_boot[:,4] - biases[4])
    BeQ3_Na_mod = shuffle(params_boot[:,5] - biases[5])
    BeQ2_K_mod = shuffle(params_boot[:,6] - biases[6])
    BeQ3_K_mod = shuffle(params_boot[:,7] - biases[7])
    BeQ4_mod = shuffle(params_boot[:,8] - biases[8])

    ScTgQ2_Na_mod = shuffle(params_boot[:,9] - biases[9])
    ScTgQ3_Na_mod = shuffle(params_boot[:,10] - biases[10])
    ScTgQ2_K_mod = shuffle(params_boot[:,11] - biases[11])
    ScTgQ3_K_mod = shuffle(params_boot[:,12] - biases[12])
    ScTgQ4_mod = shuffle(params_boot[:,13] - biases[13])

    Be_opt_b = Be_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,BeQ2_Na_mod,BeQ2_K_mod,BeQ3_Na_mod,BeQ3_K_mod,BeQ4_mod,Kqn_mod,Knk_mod)
    ScTg_opt_b, Sc_lin_b, Sc_qn_b, Sc_alka_b = Sc_calculation(q2_na,q2_k,q3_na,q3_k,q4,si,o,na,k,logQ2,logQ3,logQ4,logXk,logXna,ScTgQ2_Na_mod,ScTgQ2_K_mod,ScTgQ3_Na_mod,ScTgQ3_K_mod,ScTgQ4_mod)

    BeSc_b = Be_opt_b./ScTg_opt_b

    Tg_b = Tg_calculation(Ae_mod,Be_opt_b,ScTg_opt_b)

    # In case Tg is provided, we calculate the viscosity with it, else we use the calculated value
    if Tg_input == 0.0
        Cpg, aCpl, bCpl, ap, b = cp_model(sio2,k2o,na2o,Tg_b)

        model_n_b = n_calculation(Ae_mod,Be_opt_b,ScTg_opt_b,ap,b,T,Tg_b)
    else
        Cpg, aCpl, bCpl, ap, b = cp_model(sio2,k2o,na2o,Tg_input)

        model_n_b = n_calculation(Ae_mod,Be_opt_b,ScTg_opt_b,ap,b,T,Tg_input)
    end

    # 95% confidence intervals
    error_n = [quantile(model_n_b[:],0.0228) quantile(model_n_b[:],0.9772)]
    error_Tg = [quantile(Tg_b[:],0.0228) quantile(Tg_b[:],0.9772)]
    error_ScTg_opt = [quantile(ScTg_opt_b[:],0.0228) quantile(ScTg_opt_b[:],0.9772)]
    error_Be_opt = [quantile(Be_opt_b[:],0.0228) quantile(Be_opt_b[:],0.9772)]
    error_Be_Sc_opt = [quantile(BeSc_b[:],0.0228) quantile(BeSc_b[:],0.9772)]

    return error_n, error_Tg, error_ScTg_opt, error_Be_opt, error_Be_Sc_opt

end
