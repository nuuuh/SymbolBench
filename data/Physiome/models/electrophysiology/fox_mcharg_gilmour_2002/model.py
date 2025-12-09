# Size of variable arrays:
sizeAlgebraic = 47
sizeStates = 13
sizeConstants = 55
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "R in component membrane (joule_per_mole_kelvin)"
    legend_constants[1] = "T in component membrane (kelvin)"
    legend_constants[2] = "F in component membrane (coulomb_per_millimole)"
    legend_algebraic[20] = "i_Na in component fast_sodium_current (microA_per_microF)"
    legend_algebraic[39] = "i_Ca in component L_type_Ca_current (microA_per_microF)"
    legend_algebraic[40] = "i_CaK in component L_type_Ca_current (microA_per_microF)"
    legend_algebraic[26] = "i_Kr in component rapid_activating_delayed_rectifiyer_K_current (microA_per_microF)"
    legend_algebraic[27] = "i_Ks in component slow_activating_delayed_rectifiyer_K_current (microA_per_microF)"
    legend_algebraic[28] = "i_to in component transient_outward_potassium_current (microA_per_microF)"
    legend_algebraic[24] = "i_K1 in component time_independent_potassium_current (microA_per_microF)"
    legend_algebraic[30] = "i_Kp in component plateau_potassium_current (microA_per_microF)"
    legend_algebraic[33] = "i_NaCa in component Na_Ca_exchanger (microA_per_microF)"
    legend_algebraic[32] = "i_NaK in component sodium_potassium_pump (microA_per_microF)"
    legend_algebraic[34] = "i_p_Ca in component sarcolemmal_calcium_pump (microA_per_microF)"
    legend_algebraic[36] = "i_Ca_b in component calcium_background_current (microA_per_microF)"
    legend_algebraic[37] = "i_Na_b in component sodium_background_current (microA_per_microF)"
    legend_algebraic[10] = "i_Stim in component membrane (microA_per_microF)"
    legend_constants[3] = "stim_start in component membrane (millisecond)"
    legend_constants[4] = "stim_end in component membrane (millisecond)"
    legend_constants[5] = "stim_period in component membrane (millisecond)"
    legend_constants[6] = "stim_duration in component membrane (millisecond)"
    legend_constants[7] = "stim_amplitude in component membrane (microA_per_microF)"
    legend_constants[50] = "E_Na in component fast_sodium_current (millivolt)"
    legend_constants[8] = "g_Na in component fast_sodium_current (milliS_per_microF)"
    legend_constants[9] = "Na_o in component standard_ionic_concentrations (millimolar)"
    legend_constants[10] = "Na_i in component standard_ionic_concentrations (millimolar)"
    legend_states[1] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[2] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_states[3] = "j in component fast_sodium_current_j_gate (dimensionless)"
    legend_algebraic[11] = "alpha_m in component fast_sodium_current_m_gate (per_millisecond)"
    legend_algebraic[21] = "beta_m in component fast_sodium_current_m_gate (per_millisecond)"
    legend_algebraic[0] = "E0_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[1] = "alpha_h in component fast_sodium_current_h_gate (per_millisecond)"
    legend_algebraic[12] = "beta_h in component fast_sodium_current_h_gate (per_millisecond)"
    legend_constants[11] = "shift_h in component fast_sodium_current_h_gate (millivolt)"
    legend_algebraic[2] = "alpha_j in component fast_sodium_current_j_gate (per_millisecond)"
    legend_algebraic[13] = "beta_j in component fast_sodium_current_j_gate (per_millisecond)"
    legend_constants[12] = "shift_j in component fast_sodium_current_j_gate (millivolt)"
    legend_constants[13] = "g_K1 in component time_independent_potassium_current (milliS_per_microF)"
    legend_constants[14] = "K_mK1 in component time_independent_potassium_current (millimolar)"
    legend_constants[51] = "E_K in component rapid_activating_delayed_rectifiyer_K_current (millivolt)"
    legend_constants[15] = "K_o in component standard_ionic_concentrations (millimolar)"
    legend_algebraic[23] = "K1_infinity in component time_independent_potassium_current_K1_gate (dimensionless)"
    legend_constants[16] = "g_Kr in component rapid_activating_delayed_rectifiyer_K_current (milliS_per_microF)"
    legend_algebraic[25] = "R_V in component rapid_activating_delayed_rectifiyer_K_current (dimensionless)"
    legend_constants[17] = "K_i in component standard_ionic_concentrations (millimolar)"
    legend_states[4] = "X_kr in component rapid_activating_delayed_rectifiyer_K_current_X_kr_gate (dimensionless)"
    legend_algebraic[3] = "X_kr_inf in component rapid_activating_delayed_rectifiyer_K_current_X_kr_gate (dimensionless)"
    legend_algebraic[14] = "tau_X_kr in component rapid_activating_delayed_rectifiyer_K_current_X_kr_gate (millisecond)"
    legend_constants[18] = "g_Ks in component slow_activating_delayed_rectifiyer_K_current (milliS_per_microF)"
    legend_constants[52] = "E_Ks in component slow_activating_delayed_rectifiyer_K_current (millivolt)"
    legend_states[5] = "X_ks in component slow_activating_delayed_rectifiyer_K_current_X_ks_gate (dimensionless)"
    legend_algebraic[15] = "tau_X_ks in component slow_activating_delayed_rectifiyer_K_current_X_ks_gate (millisecond)"
    legend_algebraic[4] = "X_ks_infinity in component slow_activating_delayed_rectifiyer_K_current_X_ks_gate (dimensionless)"
    legend_constants[19] = "g_to in component transient_outward_potassium_current (milliS_per_microF)"
    legend_states[6] = "X_to in component transient_outward_potassium_current_X_to_gate (dimensionless)"
    legend_states[7] = "Y_to in component transient_outward_potassium_current_Y_to_gate (dimensionless)"
    legend_algebraic[5] = "alpha_X_to in component transient_outward_potassium_current_X_to_gate (per_millisecond)"
    legend_algebraic[16] = "beta_X_to in component transient_outward_potassium_current_X_to_gate (per_millisecond)"
    legend_algebraic[6] = "alpha_Y_to in component transient_outward_potassium_current_Y_to_gate (per_millisecond)"
    legend_algebraic[17] = "beta_Y_to in component transient_outward_potassium_current_Y_to_gate (per_millisecond)"
    legend_constants[20] = "g_Kp in component plateau_potassium_current (milliS_per_microF)"
    legend_algebraic[29] = "Kp_V in component plateau_potassium_current_Kp_gate (dimensionless)"
    legend_constants[21] = "i_NaK_max in component sodium_potassium_pump (microA_per_microF)"
    legend_algebraic[31] = "f_NaK in component sodium_potassium_pump (dimensionless)"
    legend_constants[22] = "K_mNai in component sodium_potassium_pump (millimolar)"
    legend_constants[23] = "K_mKo in component sodium_potassium_pump (millimolar)"
    legend_constants[53] = "sigma in component sodium_potassium_pump (dimensionless)"
    legend_constants[24] = "K_mCa in component Na_Ca_exchanger (micromolar)"
    legend_constants[25] = "K_mNa in component Na_Ca_exchanger (millimolar)"
    legend_constants[26] = "K_NaCa in component Na_Ca_exchanger (microA_per_microF)"
    legend_constants[27] = "K_sat in component Na_Ca_exchanger (dimensionless)"
    legend_constants[28] = "eta in component Na_Ca_exchanger (dimensionless)"
    legend_states[8] = "Ca_i in component calcium_dynamics (micromolar)"
    legend_constants[29] = "Ca_o in component standard_ionic_concentrations (micromolar)"
    legend_constants[30] = "K_mpCa in component sarcolemmal_calcium_pump (micromolar)"
    legend_constants[31] = "i_pCa_max in component sarcolemmal_calcium_pump (microA_per_microF)"
    legend_constants[32] = "g_Cab in component calcium_background_current (milliS_per_microF)"
    legend_algebraic[35] = "E_Ca in component calcium_background_current (millivolt)"
    legend_constants[33] = "g_Nab in component sodium_background_current (milliS_per_microF)"
    legend_constants[34] = "P_Ca in component L_type_Ca_current (cm_per_millisecond)"
    legend_constants[35] = "P_CaK in component L_type_Ca_current (cm_per_millisecond)"
    legend_constants[36] = "i_Ca_half in component L_type_Ca_current (microA_per_microF)"
    legend_algebraic[38] = "i_Ca_max in component L_type_Ca_current (microA_per_microF)"
    legend_constants[37] = "C_sc in component L_type_Ca_current (microF_per_cm2)"
    legend_states[9] = "f in component L_type_Ca_current_f_gate (dimensionless)"
    legend_states[10] = "d in component L_type_Ca_current_d_gate (dimensionless)"
    legend_states[11] = "f_Ca in component L_type_Ca_current_f_Ca_gate (dimensionless)"
    legend_algebraic[7] = "f_infinity in component L_type_Ca_current_f_gate (dimensionless)"
    legend_algebraic[18] = "tau_f in component L_type_Ca_current_f_gate (millisecond)"
    legend_algebraic[8] = "d_infinity in component L_type_Ca_current_d_gate (dimensionless)"
    legend_algebraic[22] = "tau_d in component L_type_Ca_current_d_gate (millisecond)"
    legend_algebraic[19] = "E0_m in component L_type_Ca_current_d_gate (millivolt)"
    legend_constants[54] = "tau_f_Ca in component L_type_Ca_current_f_Ca_gate (millisecond)"
    legend_algebraic[9] = "f_Ca_infinity in component L_type_Ca_current_f_Ca_gate (dimensionless)"
    legend_constants[38] = "K_mfCa in component L_type_Ca_current_f_Ca_gate (micromolar)"
    legend_algebraic[45] = "beta_i in component calcium_dynamics (dimensionless)"
    legend_constants[39] = "K_mCMDN in component calcium_dynamics (micromolar)"
    legend_constants[40] = "CMDN_tot in component calcium_dynamics (micromolar)"
    legend_constants[41] = "V_myo in component calcium_dynamics (microlitre)"
    legend_constants[42] = "A_Cap in component calcium_dynamics (cm2)"
    legend_algebraic[43] = "J_rel in component calcium_dynamics (micromolar_per_millisecond)"
    legend_algebraic[44] = "J_leak in component calcium_dynamics (micromolar_per_millisecond)"
    legend_algebraic[41] = "J_up in component calcium_dynamics (micromolar_per_millisecond)"
    legend_states[12] = "Ca_SR in component calcium_dynamics (micromolar)"
    legend_constants[43] = "P_rel in component calcium_dynamics (per_millisecond)"
    legend_constants[44] = "P_leak in component calcium_dynamics (per_millisecond)"
    legend_constants[45] = "K_mCSQN in component calcium_dynamics (micromolar)"
    legend_constants[46] = "CSQN_tot in component calcium_dynamics (micromolar)"
    legend_constants[47] = "V_SR in component calcium_dynamics (microlitre)"
    legend_constants[48] = "V_up in component calcium_dynamics (micromolar_per_millisecond)"
    legend_constants[49] = "K_mup in component calcium_dynamics (micromolar)"
    legend_algebraic[42] = "gamma in component calcium_dynamics (dimensionless)"
    legend_algebraic[46] = "beta_SR in component calcium_dynamics (dimensionless)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[2] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[3] = "d/dt j in component fast_sodium_current_j_gate (dimensionless)"
    legend_rates[4] = "d/dt X_kr in component rapid_activating_delayed_rectifiyer_K_current_X_kr_gate (dimensionless)"
    legend_rates[5] = "d/dt X_ks in component slow_activating_delayed_rectifiyer_K_current_X_ks_gate (dimensionless)"
    legend_rates[6] = "d/dt X_to in component transient_outward_potassium_current_X_to_gate (dimensionless)"
    legend_rates[7] = "d/dt Y_to in component transient_outward_potassium_current_Y_to_gate (dimensionless)"
    legend_rates[9] = "d/dt f in component L_type_Ca_current_f_gate (dimensionless)"
    legend_rates[10] = "d/dt d in component L_type_Ca_current_d_gate (dimensionless)"
    legend_rates[11] = "d/dt f_Ca in component L_type_Ca_current_f_Ca_gate (dimensionless)"
    legend_rates[12] = "d/dt Ca_SR in component calcium_dynamics (micromolar)"
    legend_rates[8] = "d/dt Ca_i in component calcium_dynamics (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -94.7
    constants[0] = 8.314
    constants[1] = 310
    constants[2] = 96.5
    constants[3] = 50
    constants[4] = 9000
    constants[5] = 1000
    constants[6] = 1
    constants[7] = -80
    constants[8] = 12.8
    constants[9] = 138
    constants[10] = 10
    states[1] = 0.00024676
    states[2] = 0.99869
    states[3] = 0.99887
    constants[11] = 0
    constants[12] = 0
    constants[13] = 2.8
    constants[14] = 13
    constants[15] = 4
    constants[16] = 0.0136
    constants[17] = 149.4
    states[4] = 0.229
    constants[18] = 0.0245
    states[5] = 0.0001
    constants[19] = 0.23815
    states[6] = 0.00003742
    states[7] = 1
    constants[20] = 0.002216
    constants[21] = 0.693
    constants[22] = 10
    constants[23] = 1.5
    constants[24] = 1380
    constants[25] = 87.5
    constants[26] = 1500
    constants[27] = 0.2
    constants[28] = 0.35
    states[8] = 0.0472
    constants[29] = 2000
    constants[30] = 0.05
    constants[31] = 0.05
    constants[32] = 0.0003842
    constants[33] = 0.0031
    constants[34] = 0.0000226
    constants[35] = 0.000000579
    constants[36] = -0.265
    constants[37] = 1
    states[9] = 0.983
    states[10] = 0.0001
    states[11] = 0.942
    constants[38] = 0.18
    constants[39] = 2
    constants[40] = 10
    constants[41] = 0.00002584
    constants[42] = 0.0001534
    states[12] = 320
    constants[43] = 6
    constants[44] = 0.000001
    constants[45] = 600
    constants[46] = 10000
    constants[47] = 0.000002
    constants[48] = 0.1
    constants[49] = 0.32
    constants[50] = ((constants[0]*constants[1])/constants[2])*log(constants[9]/constants[10])
    constants[51] = ((constants[0]*constants[1])/constants[2])*log(constants[15]/constants[17])
    constants[52] = ((constants[0]*constants[1])/constants[2])*log((constants[15]+0.0183300*constants[9])/(constants[17]+0.0183300*constants[10]))
    constants[53] = (1.00000/7.00000)*(exp(constants[9]/67.3000)-1.00000)
    constants[54] = 30.0000
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[9] = 1.00000/(1.00000+power(states[8]/constants[38], 3.00000))
    rates[11] = (algebraic[9]-states[11])/constants[54]
    algebraic[1] = 0.135000*exp(((states[0]+80.0000)-constants[11])/-6.80000)
    algebraic[12] = 7.50000/(1.00000+exp(-0.100000*((states[0]+11.0000)-constants[11])))
    rates[2] = algebraic[1]*(1.00000-states[2])-algebraic[12]*states[2]
    algebraic[2] = (0.175000*exp(((states[0]+100.000)-constants[12])/-23.0000))/(1.00000+exp(0.150000*((states[0]+79.0000)-constants[12])))
    algebraic[13] = 0.300000/(1.00000+exp(-0.100000*((states[0]+32.0000)-constants[12])))
    rates[3] = algebraic[2]*(1.00000-states[3])-algebraic[13]*states[3]
    algebraic[3] = 1.00000/(1.00000+exp(-2.18200-0.181900*states[0]))
    algebraic[14] = 43.0000+1.00000/(exp(-5.49500+0.169100*states[0])+exp(-7.67700-0.0128000*states[0]))
    rates[4] = (algebraic[3]-states[4])/algebraic[14]
    algebraic[15] = 1.00000/((7.19000e-05*(states[0]-10.0000))/(1.00000-exp(-0.148000*(states[0]-10.0000)))+(0.000131000*(states[0]-10.0000))/(exp(0.0687000*(states[0]-10.0000))-1.00000))
    algebraic[4] = 1.00000/(1.00000+exp((states[0]-16.0000)/-13.6000))
    rates[5] = (algebraic[4]-states[5])/algebraic[15]
    algebraic[5] = 0.0451600*exp(0.0357700*states[0])
    algebraic[16] = 0.0989000*exp(-0.0623700*states[0])
    rates[6] = algebraic[5]*(1.00000-states[6])-algebraic[16]*states[6]
    algebraic[6] = (0.00541500*exp((states[0]+33.5000)/-5.00000))/(1.00000+0.0513350*exp((states[0]+33.5000)/-5.00000))
    algebraic[17] = (0.00541500*exp((states[0]+33.5000)/5.00000))/(1.00000+0.0513350*exp((states[0]+33.5000)/5.00000))
    rates[7] = algebraic[6]*(1.00000-states[7])-algebraic[17]*states[7]
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+12.5000)/5.00000))
    algebraic[18] = 30.0000+200.000/(1.00000+exp((states[0]+20.0000)/9.50000))
    rates[9] = (algebraic[7]-states[9])/algebraic[18]
    algebraic[0] = states[0]+47.1300
    algebraic[11] = (0.320000*algebraic[0])/(1.00000-exp(-0.100000*algebraic[0]))
    algebraic[21] = 0.0800000*exp(-states[0]/11.0000)
    rates[1] = algebraic[11]*(1.00000-states[1])-algebraic[21]*states[1]
    algebraic[8] = 1.00000/(1.00000+exp((states[0]+10.0000)/-6.24000))
    algebraic[19] = states[0]+40.0000
    algebraic[22] = 1.00000/((0.250000*exp(-0.0100000*states[0]))/(1.00000+exp(-0.0700000*states[0]))+(0.0700000*exp(-0.0500000*algebraic[19]))/(1.00000+exp(0.0500000*algebraic[19])))
    rates[10] = (algebraic[8]-states[10])/algebraic[22]
    algebraic[20] = constants[8]*(power(states[1], 3.00000))*states[2]*states[3]*(states[0]-constants[50])
    algebraic[38] = ((((constants[34]/constants[37])*4.00000*states[0]*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(states[8]*exp((2.00000*states[0]*constants[2])/(constants[0]*constants[1]))-0.341000*constants[29]))/(exp((2.00000*states[0]*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[39] = algebraic[38]*states[9]*states[10]*states[11]
    algebraic[40] = ((((((constants[35]/constants[37])*states[9]*states[10]*states[11])/(1.00000+algebraic[38]/constants[36]))*1000.00*states[0]*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(constants[17]*exp((states[0]*constants[2])/(constants[0]*constants[1]))-constants[15]))/(exp((states[0]*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[25] = 1.00000/(1.00000+2.50000*exp(0.100000*(states[0]+28.0000)))
    algebraic[26] = constants[16]*algebraic[25]*states[4]*(power(constants[15]/4.00000, 1.0/2))*(states[0]-constants[51])
    algebraic[27] = constants[18]*(power(states[5], 2.00000))*(states[0]-constants[52])
    algebraic[28] = constants[19]*states[6]*states[7]*(states[0]-constants[51])
    algebraic[23] = 1.00000/(2.00000+exp(((1.62000*constants[2])/(constants[0]*constants[1]))*(states[0]-constants[51])))
    algebraic[24] = ((constants[13]*algebraic[23]*constants[15])/(constants[15]+constants[14]))*(states[0]-constants[51])
    algebraic[29] = 1.00000/(1.00000+exp((7.48800-states[0])/5.98000))
    algebraic[30] = constants[20]*algebraic[29]*(states[0]-constants[51])
    algebraic[33] = (constants[26]/((power(constants[25], 3.00000)+power(constants[9], 3.00000))*(constants[24]+constants[29])*(1.00000+constants[27]*exp(((constants[28]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1])))))*(exp((constants[28]*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[10], 3.00000))*constants[29]-exp(((constants[28]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[9], 3.00000))*states[8])
    algebraic[31] = 1.00000/(1.00000+0.124500*exp((-0.100000*states[0]*constants[2])/(constants[0]*constants[1]))+0.0365000*constants[53]*exp((-states[0]*constants[2])/(constants[0]*constants[1])))
    algebraic[32] = (((constants[21]*algebraic[31])/(1.00000+power(constants[22]/constants[10], 1.50000)))*constants[15])/(constants[15]+constants[23])
    algebraic[34] = (constants[31]*states[8])/(constants[30]+states[8])
    algebraic[35] = ((constants[0]*constants[1])/(2.00000*constants[2]))*log(constants[29]/states[8])
    algebraic[36] = constants[32]*(states[0]-algebraic[35])
    algebraic[37] = constants[33]*(states[0]-constants[50])
    algebraic[10] = custom_piecewise([greater_equal(voi , constants[3]) & less_equal(voi , constants[4]) & less_equal((voi-constants[3])-floor((voi-constants[3])/constants[5])*constants[5] , constants[6]), constants[7] , True, 0.00000])
    rates[0] = -(algebraic[20]+algebraic[39]+algebraic[40]+algebraic[26]+algebraic[27]+algebraic[28]+algebraic[24]+algebraic[30]+algebraic[33]+algebraic[32]+algebraic[34]+algebraic[37]+algebraic[36]+algebraic[10])
    algebraic[42] = 1.00000/(1.00000+power(2000.00/states[12], 3.00000))
    algebraic[43] = (constants[43]*states[9]*states[10]*states[11]*(algebraic[42]*states[12]-states[8]))/(1.00000+1.65000*exp(states[0]/20.0000))
    algebraic[44] = constants[44]*(states[12]-states[8])
    algebraic[41] = constants[48]/(1.00000+power(constants[49]/states[8], 2.00000))
    algebraic[46] = 1.00000/(1.00000+(constants[46]*constants[45])/(power(constants[45]+states[12], 2.00000)))
    rates[12] = (algebraic[46]*((algebraic[41]-algebraic[44])-algebraic[43])*constants[41])/constants[47]
    algebraic[45] = 1.00000/(1.00000+(constants[40]*constants[39])/(power(constants[39]+states[8], 2.00000)))
    rates[8] = algebraic[45]*(((algebraic[43]+algebraic[44])-algebraic[41])-((constants[42]*constants[37])/(2.00000*constants[2]*constants[41]))*((algebraic[39]+algebraic[36]+algebraic[34])-2.00000*algebraic[33]))
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[9] = 1.00000/(1.00000+power(states[8]/constants[38], 3.00000))
    algebraic[1] = 0.135000*exp(((states[0]+80.0000)-constants[11])/-6.80000)
    algebraic[12] = 7.50000/(1.00000+exp(-0.100000*((states[0]+11.0000)-constants[11])))
    algebraic[2] = (0.175000*exp(((states[0]+100.000)-constants[12])/-23.0000))/(1.00000+exp(0.150000*((states[0]+79.0000)-constants[12])))
    algebraic[13] = 0.300000/(1.00000+exp(-0.100000*((states[0]+32.0000)-constants[12])))
    algebraic[3] = 1.00000/(1.00000+exp(-2.18200-0.181900*states[0]))
    algebraic[14] = 43.0000+1.00000/(exp(-5.49500+0.169100*states[0])+exp(-7.67700-0.0128000*states[0]))
    algebraic[15] = 1.00000/((7.19000e-05*(states[0]-10.0000))/(1.00000-exp(-0.148000*(states[0]-10.0000)))+(0.000131000*(states[0]-10.0000))/(exp(0.0687000*(states[0]-10.0000))-1.00000))
    algebraic[4] = 1.00000/(1.00000+exp((states[0]-16.0000)/-13.6000))
    algebraic[5] = 0.0451600*exp(0.0357700*states[0])
    algebraic[16] = 0.0989000*exp(-0.0623700*states[0])
    algebraic[6] = (0.00541500*exp((states[0]+33.5000)/-5.00000))/(1.00000+0.0513350*exp((states[0]+33.5000)/-5.00000))
    algebraic[17] = (0.00541500*exp((states[0]+33.5000)/5.00000))/(1.00000+0.0513350*exp((states[0]+33.5000)/5.00000))
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+12.5000)/5.00000))
    algebraic[18] = 30.0000+200.000/(1.00000+exp((states[0]+20.0000)/9.50000))
    algebraic[0] = states[0]+47.1300
    algebraic[11] = (0.320000*algebraic[0])/(1.00000-exp(-0.100000*algebraic[0]))
    algebraic[21] = 0.0800000*exp(-states[0]/11.0000)
    algebraic[8] = 1.00000/(1.00000+exp((states[0]+10.0000)/-6.24000))
    algebraic[19] = states[0]+40.0000
    algebraic[22] = 1.00000/((0.250000*exp(-0.0100000*states[0]))/(1.00000+exp(-0.0700000*states[0]))+(0.0700000*exp(-0.0500000*algebraic[19]))/(1.00000+exp(0.0500000*algebraic[19])))
    algebraic[20] = constants[8]*(power(states[1], 3.00000))*states[2]*states[3]*(states[0]-constants[50])
    algebraic[38] = ((((constants[34]/constants[37])*4.00000*states[0]*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(states[8]*exp((2.00000*states[0]*constants[2])/(constants[0]*constants[1]))-0.341000*constants[29]))/(exp((2.00000*states[0]*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[39] = algebraic[38]*states[9]*states[10]*states[11]
    algebraic[40] = ((((((constants[35]/constants[37])*states[9]*states[10]*states[11])/(1.00000+algebraic[38]/constants[36]))*1000.00*states[0]*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(constants[17]*exp((states[0]*constants[2])/(constants[0]*constants[1]))-constants[15]))/(exp((states[0]*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[25] = 1.00000/(1.00000+2.50000*exp(0.100000*(states[0]+28.0000)))
    algebraic[26] = constants[16]*algebraic[25]*states[4]*(power(constants[15]/4.00000, 1.0/2))*(states[0]-constants[51])
    algebraic[27] = constants[18]*(power(states[5], 2.00000))*(states[0]-constants[52])
    algebraic[28] = constants[19]*states[6]*states[7]*(states[0]-constants[51])
    algebraic[23] = 1.00000/(2.00000+exp(((1.62000*constants[2])/(constants[0]*constants[1]))*(states[0]-constants[51])))
    algebraic[24] = ((constants[13]*algebraic[23]*constants[15])/(constants[15]+constants[14]))*(states[0]-constants[51])
    algebraic[29] = 1.00000/(1.00000+exp((7.48800-states[0])/5.98000))
    algebraic[30] = constants[20]*algebraic[29]*(states[0]-constants[51])
    algebraic[33] = (constants[26]/((power(constants[25], 3.00000)+power(constants[9], 3.00000))*(constants[24]+constants[29])*(1.00000+constants[27]*exp(((constants[28]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1])))))*(exp((constants[28]*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[10], 3.00000))*constants[29]-exp(((constants[28]-1.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[9], 3.00000))*states[8])
    algebraic[31] = 1.00000/(1.00000+0.124500*exp((-0.100000*states[0]*constants[2])/(constants[0]*constants[1]))+0.0365000*constants[53]*exp((-states[0]*constants[2])/(constants[0]*constants[1])))
    algebraic[32] = (((constants[21]*algebraic[31])/(1.00000+power(constants[22]/constants[10], 1.50000)))*constants[15])/(constants[15]+constants[23])
    algebraic[34] = (constants[31]*states[8])/(constants[30]+states[8])
    algebraic[35] = ((constants[0]*constants[1])/(2.00000*constants[2]))*log(constants[29]/states[8])
    algebraic[36] = constants[32]*(states[0]-algebraic[35])
    algebraic[37] = constants[33]*(states[0]-constants[50])
    algebraic[10] = custom_piecewise([greater_equal(voi , constants[3]) & less_equal(voi , constants[4]) & less_equal((voi-constants[3])-floor((voi-constants[3])/constants[5])*constants[5] , constants[6]), constants[7] , True, 0.00000])
    algebraic[42] = 1.00000/(1.00000+power(2000.00/states[12], 3.00000))
    algebraic[43] = (constants[43]*states[9]*states[10]*states[11]*(algebraic[42]*states[12]-states[8]))/(1.00000+1.65000*exp(states[0]/20.0000))
    algebraic[44] = constants[44]*(states[12]-states[8])
    algebraic[41] = constants[48]/(1.00000+power(constants[49]/states[8], 2.00000))
    algebraic[46] = 1.00000/(1.00000+(constants[46]*constants[45])/(power(constants[45]+states[12], 2.00000)))
    algebraic[45] = 1.00000/(1.00000+(constants[40]*constants[39])/(power(constants[39]+states[8], 2.00000)))
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0, 10, 500)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.xlabel(legend_voi)
    pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)