# Size of variable arrays:
sizeAlgebraic = 62
sizeStates = 25
sizeConstants = 92
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "R in component membrane (joule_per_kilomole_kelvin)"
    legend_constants[1] = "T in component membrane (kelvin)"
    legend_constants[2] = "F in component membrane (coulomb_per_mole)"
    legend_constants[3] = "Cm in component membrane (microF)"
    legend_algebraic[35] = "i_K1 in component time_independent_potassium_current (nanoA)"
    legend_algebraic[52] = "i_to in component transient_outward_current (nanoA)"
    legend_algebraic[37] = "i_Kr in component rapid_delayed_rectifier_potassium_current (nanoA)"
    legend_algebraic[39] = "i_Ks in component slow_delayed_rectifier_potassium_current (nanoA)"
    legend_algebraic[40] = "i_K_ATP in component ATP_dependent_potassium_current (nanoA)"
    legend_algebraic[53] = "i_K_ACh in component ACh_dependent_potassium_current (nanoA)"
    legend_algebraic[45] = "i_Ca_L_K_cyt in component L_type_Ca_channel (nanoA)"
    legend_algebraic[48] = "i_Ca_L_K_ds in component L_type_Ca_channel (nanoA)"
    legend_algebraic[54] = "i_NaK in component sodium_potassium_pump (nanoA)"
    legend_algebraic[41] = "i_Na in component fast_sodium_current (nanoA)"
    legend_algebraic[43] = "i_b_Na in component sodium_background_current (nanoA)"
    legend_algebraic[42] = "i_p_Na in component persistent_sodium_current (nanoA)"
    legend_algebraic[46] = "i_Ca_L_Na_cyt in component L_type_Ca_channel (nanoA)"
    legend_algebraic[49] = "i_Ca_L_Na_ds in component L_type_Ca_channel (nanoA)"
    legend_algebraic[55] = "i_NaCa_cyt in component sodium_calcium_exchanger (nanoA)"
    legend_algebraic[56] = "i_NaCa_ds in component sodium_calcium_exchanger (nanoA)"
    legend_algebraic[44] = "i_Ca_L_Ca_cyt in component L_type_Ca_channel (nanoA)"
    legend_algebraic[47] = "i_Ca_L_Ca_ds in component L_type_Ca_channel (nanoA)"
    legend_algebraic[51] = "i_b_Ca in component calcium_background_current (nanoA)"
    legend_algebraic[8] = "i_Stim in component membrane (nanoA)"
    legend_constants[4] = "stim_start in component membrane (second)"
    legend_constants[5] = "stim_end in component membrane (second)"
    legend_constants[6] = "stim_period in component membrane (second)"
    legend_constants[7] = "stim_duration in component membrane (second)"
    legend_constants[8] = "stim_amplitude in component membrane (nanoA)"
    legend_algebraic[20] = "E_Na in component reversal_potentials (millivolt)"
    legend_algebraic[26] = "E_K in component reversal_potentials (millivolt)"
    legend_algebraic[29] = "E_Ks in component reversal_potentials (millivolt)"
    legend_algebraic[31] = "E_Ca in component reversal_potentials (millivolt)"
    legend_algebraic[33] = "E_mh in component reversal_potentials (millivolt)"
    legend_constants[9] = "P_kna in component reversal_potentials (dimensionless)"
    legend_constants[10] = "K_o in component extracellular_potassium_concentration (millimolar)"
    legend_constants[11] = "Na_o in component extracellular_sodium_concentration (millimolar)"
    legend_states[1] = "K_i in component intracellular_potassium_concentration (millimolar)"
    legend_states[2] = "Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_constants[12] = "Ca_o in component extracellular_calcium_concentration (millimolar)"
    legend_states[3] = "Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_constants[13] = "K_mk1 in component time_independent_potassium_current (millimolar)"
    legend_constants[14] = "g_K1 in component time_independent_potassium_current (microS)"
    legend_constants[15] = "shiftK1 in component time_independent_potassium_current (millivolt)"
    legend_constants[16] = "steepK1 in component time_independent_potassium_current (dimensionless)"
    legend_constants[17] = "g_Kr1 in component rapid_delayed_rectifier_potassium_current (microS)"
    legend_constants[18] = "g_Kr2 in component rapid_delayed_rectifier_potassium_current (microS)"
    legend_states[4] = "xr1 in component rapid_delayed_rectifier_potassium_current_xr1_gate (dimensionless)"
    legend_states[5] = "xr2 in component rapid_delayed_rectifier_potassium_current_xr2_gate (dimensionless)"
    legend_algebraic[0] = "alpha_xr1 in component rapid_delayed_rectifier_potassium_current_xr1_gate (per_second)"
    legend_algebraic[12] = "beta_xr1 in component rapid_delayed_rectifier_potassium_current_xr1_gate (per_second)"
    legend_algebraic[1] = "alpha_xr2 in component rapid_delayed_rectifier_potassium_current_xr2_gate (per_second)"
    legend_algebraic[13] = "beta_xr2 in component rapid_delayed_rectifier_potassium_current_xr2_gate (per_second)"
    legend_constants[19] = "g_Ks in component slow_delayed_rectifier_potassium_current (microS)"
    legend_states[6] = "xs in component slow_delayed_rectifier_potassium_current_xs_gate (dimensionless)"
    legend_algebraic[2] = "alpha_xs in component slow_delayed_rectifier_potassium_current_xs_gate (per_second)"
    legend_algebraic[14] = "beta_xs in component slow_delayed_rectifier_potassium_current_xs_gate (per_second)"
    legend_constants[20] = "g_K_ATP in component ATP_dependent_potassium_current (microS)"
    legend_constants[21] = "K_ATP in component ATP_dependent_potassium_current (millimolar)"
    legend_constants[22] = "ATP in component ATP_dependent_potassium_current (millimolar)"
    legend_algebraic[28] = "i_KNa in component sodium_activated_potassium_current (nanoA)"
    legend_constants[23] = "g_K_Na in component sodium_activated_potassium_current (microS)"
    legend_constants[24] = "K_kna in component sodium_activated_potassium_current (millimolar)"
    legend_constants[25] = "g_Na in component fast_sodium_current (microS)"
    legend_states[7] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[8] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_algebraic[15] = "alpha_m in component fast_sodium_current_m_gate (per_second)"
    legend_algebraic[23] = "beta_m in component fast_sodium_current_m_gate (per_second)"
    legend_constants[26] = "delta_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[3] = "E0_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[4] = "alpha_h in component fast_sodium_current_h_gate (per_second)"
    legend_algebraic[16] = "beta_h in component fast_sodium_current_h_gate (per_second)"
    legend_constants[27] = "shift_h in component fast_sodium_current_h_gate (millivolt)"
    legend_constants[28] = "g_pna in component persistent_sodium_current (microS)"
    legend_constants[29] = "g_bna in component sodium_background_current (microS)"
    legend_algebraic[50] = "i_Ca_L in component L_type_Ca_channel (nanoA)"
    legend_constants[30] = "P_Ca_L in component L_type_Ca_channel (nanoA_per_millimolar)"
    legend_constants[31] = "P_CaK in component L_type_Ca_channel (dimensionless)"
    legend_constants[32] = "P_CaNa in component L_type_Ca_channel (dimensionless)"
    legend_states[9] = "Ca_ds in component intracellular_calcium_concentration (millimolar)"
    legend_states[10] = "d in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_states[11] = "f in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_states[12] = "f2 in component L_type_Ca_channel_f2_gate (dimensionless)"
    legend_states[13] = "f2ds in component L_type_Ca_channel_f2ds_gate (dimensionless)"
    legend_constants[33] = "Km_f2 in component L_type_Ca_channel (millimolar)"
    legend_constants[34] = "Km_f2ds in component L_type_Ca_channel (millimolar)"
    legend_constants[35] = "R_decay in component L_type_Ca_channel (per_second)"
    legend_constants[36] = "FrICa in component L_type_Ca_channel (dimensionless)"
    legend_algebraic[17] = "alpha_d in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[24] = "beta_d in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[5] = "E0_d in component L_type_Ca_channel_d_gate (millivolt)"
    legend_constants[37] = "speed_d in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_algebraic[18] = "alpha_f in component L_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[25] = "beta_f in component L_type_Ca_channel_f_gate (per_second)"
    legend_constants[38] = "speed_f in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_constants[39] = "delta_f in component L_type_Ca_channel_f_gate (millivolt)"
    legend_algebraic[6] = "E0_f in component L_type_Ca_channel_f_gate (millivolt)"
    legend_constants[40] = "g_bca in component calcium_background_current (microS)"
    legend_constants[41] = "g_to in component transient_outward_current (microS)"
    legend_constants[42] = "g_tos in component transient_outward_current (dimensionless)"
    legend_states[14] = "s in component transient_outward_current_s_gate (dimensionless)"
    legend_states[15] = "r in component transient_outward_current_r_gate (dimensionless)"
    legend_algebraic[7] = "alpha_s in component transient_outward_current_s_gate (per_second)"
    legend_algebraic[19] = "beta_s in component transient_outward_current_s_gate (per_second)"
    legend_constants[43] = "g_KACh in component ACh_dependent_potassium_current (microS)"
    legend_constants[44] = "ACh in component ACh_dependent_potassium_current (millimolar)"
    legend_constants[45] = "K_D in component ACh_dependent_potassium_current (millimolar)"
    legend_states[16] = "x_ACh in component ACh_dependent_potassium_current_xACh_gate (dimensionless)"
    legend_constants[46] = "alpha_ACh in component ACh_dependent_potassium_current_xACh_gate (per_second)"
    legend_constants[47] = "beta_ACh in component ACh_dependent_potassium_current_xACh_gate (per_second)"
    legend_constants[48] = "i_NaK_max in component sodium_potassium_pump (nanoA)"
    legend_constants[49] = "K_mK in component sodium_potassium_pump (millimolar)"
    legend_constants[50] = "K_mNa in component sodium_potassium_pump (millimolar)"
    legend_algebraic[57] = "i_NaCa in component sodium_calcium_exchanger (nanoA)"
    legend_constants[51] = "k_NaCa in component sodium_calcium_exchanger (nanoA)"
    legend_constants[52] = "n_NaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[53] = "d_NaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[54] = "gamma in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[55] = "FRiNaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_algebraic[59] = "i_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[85] = "K_1 in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_algebraic[58] = "K_2 in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[56] = "K_cyca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[57] = "K_xcs in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_constants[58] = "K_srca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[59] = "alpha_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[60] = "beta_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_states[17] = "Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[60] = "i_trans in component calcium_translocation (millimolar_per_second)"
    legend_states[18] = "Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[61] = "i_rel in component calcium_release (millimolar_per_second)"
    legend_algebraic[9] = "VoltDep in component calcium_release (dimensionless)"
    legend_algebraic[30] = "RegBindSite in component calcium_release (dimensionless)"
    legend_algebraic[21] = "CaiReg in component calcium_release (dimensionless)"
    legend_algebraic[27] = "CadsReg in component calcium_release (dimensionless)"
    legend_algebraic[32] = "ActRate in component calcium_release (per_second)"
    legend_algebraic[34] = "InactRate in component calcium_release (per_second)"
    legend_constants[61] = "K_leak_rate in component calcium_release (per_second)"
    legend_constants[62] = "K_m_rel in component calcium_release (per_second)"
    legend_constants[63] = "K_m_Ca_cyt in component calcium_release (millimolar)"
    legend_constants[64] = "K_m_Ca_ds in component calcium_release (millimolar)"
    legend_algebraic[38] = "PrecFrac in component calcium_release (dimensionless)"
    legend_states[19] = "ActFrac in component calcium_release (dimensionless)"
    legend_states[20] = "ProdFrac in component calcium_release (dimensionless)"
    legend_algebraic[36] = "SpeedRel in component calcium_release (dimensionless)"
    legend_constants[90] = "V_i in component intracellular_calcium_concentration (micrometre3)"
    legend_states[21] = "Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_states[22] = "Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    legend_constants[65] = "Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_constants[66] = "Trop in component intracellular_calcium_concentration (millimolar)"
    legend_constants[67] = "alpha_Calmod in component intracellular_calcium_concentration (per_millimolar_second)"
    legend_constants[68] = "beta_Calmod in component intracellular_calcium_concentration (per_second)"
    legend_constants[69] = "alpha_Trop in component intracellular_calcium_concentration (per_millimolar_second)"
    legend_constants[70] = "beta_Trop in component intracellular_calcium_concentration (per_second)"
    legend_constants[71] = "radius in component intracellular_calcium_concentration (micrometre)"
    legend_constants[72] = "length in component intracellular_calcium_concentration (micrometre)"
    legend_constants[84] = "V_Cell in component intracellular_calcium_concentration (micrometre3)"
    legend_constants[88] = "V_i_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[73] = "V_ds_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[74] = "V_rel_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[75] = "V_e_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[76] = "V_up_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[77] = "Kdecay in component intracellular_calcium_concentration (per_second)"
    legend_states[23] = "light_chain in component contraction (dimensionless)"
    legend_states[24] = "cross_bridge in component contraction (dimensionless)"
    legend_constants[78] = "KCont1 in component contraction (per_second)"
    legend_algebraic[10] = "XCont2 in component contraction (dimensionless)"
    legend_algebraic[22] = "XCont1 in component contraction (dimensionless)"
    legend_constants[79] = "KCont2 in component contraction (per_second)"
    legend_constants[80] = "KCont3 in component contraction (per_second)"
    legend_constants[81] = "KCont4 in component contraction (per_second)"
    legend_constants[82] = "sarcomere_length in component contraction (micrometre)"
    legend_constants[83] = "cross_bridge_density in component contraction (per_micrometre)"
    legend_constants[86] = "tension_rest in component contraction (dimensionless)"
    legend_constants[87] = "tension_active in component contraction (dimensionless)"
    legend_constants[89] = "overlap in component contraction (micrometre)"
    legend_constants[91] = "cross_bridge_availability in component contraction (dimensionless)"
    legend_algebraic[11] = "isometric_tension in component contraction (dimensionless)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[4] = "d/dt xr1 in component rapid_delayed_rectifier_potassium_current_xr1_gate (dimensionless)"
    legend_rates[5] = "d/dt xr2 in component rapid_delayed_rectifier_potassium_current_xr2_gate (dimensionless)"
    legend_rates[6] = "d/dt xs in component slow_delayed_rectifier_potassium_current_xs_gate (dimensionless)"
    legend_rates[7] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[8] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[10] = "d/dt d in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_rates[11] = "d/dt f in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_rates[12] = "d/dt f2 in component L_type_Ca_channel_f2_gate (dimensionless)"
    legend_rates[13] = "d/dt f2ds in component L_type_Ca_channel_f2ds_gate (dimensionless)"
    legend_rates[14] = "d/dt s in component transient_outward_current_s_gate (dimensionless)"
    legend_rates[15] = "d/dt r in component transient_outward_current_r_gate (dimensionless)"
    legend_rates[16] = "d/dt x_ACh in component ACh_dependent_potassium_current_xACh_gate (dimensionless)"
    legend_rates[19] = "d/dt ActFrac in component calcium_release (dimensionless)"
    legend_rates[20] = "d/dt ProdFrac in component calcium_release (dimensionless)"
    legend_rates[2] = "d/dt Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_rates[1] = "d/dt K_i in component intracellular_potassium_concentration (millimolar)"
    legend_rates[3] = "d/dt Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_rates[21] = "d/dt Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_rates[22] = "d/dt Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    legend_rates[9] = "d/dt Ca_ds in component intracellular_calcium_concentration (millimolar)"
    legend_rates[17] = "d/dt Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_rates[18] = "d/dt Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_rates[23] = "d/dt light_chain in component contraction (dimensionless)"
    legend_rates[24] = "d/dt cross_bridge in component contraction (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -92.849333
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 9.5e-5
    constants[4] = 0.1
    constants[5] = 1000
    constants[6] = 0.5
    constants[7] = 0.003
    constants[8] = -3
    constants[9] = 0.03
    constants[10] = 4
    constants[11] = 140
    states[1] = 136.5644281
    states[2] = 7.3321223
    constants[12] = 2
    states[3] = 1.4e-5
    constants[13] = 10
    constants[14] = 0.5
    constants[15] = 10
    constants[16] = 1.25
    constants[17] = 0.0021
    constants[18] = 0.0013
    states[4] = 1.03e-5
    states[5] = 2e-7
    constants[19] = 0.0026
    states[6] = 0.001302
    constants[20] = 0
    constants[21] = 0.1
    constants[22] = 0
    constants[23] = 0
    constants[24] = 20
    constants[25] = 2.5
    states[7] = 0.0016203
    states[8] = 0.9944036
    constants[26] = 1e-5
    constants[27] = 0
    constants[28] = 0.004
    constants[29] = 0.0006
    constants[30] = 0.1
    constants[31] = 0.002
    constants[32] = 0.01
    states[9] = 1.88e-5
    states[10] = 0
    states[11] = 1
    states[12] = 0.9349197
    states[13] = 0.9651958
    constants[33] = 100000
    constants[34] = 0.001
    constants[35] = 20
    constants[36] = 1
    constants[37] = 3
    constants[38] = 0.3
    constants[39] = 0.0001
    constants[40] = 0.00025
    constants[41] = 0.005
    constants[42] = 0
    states[14] = 0.9948645
    states[15] = 0
    constants[43] = 0
    constants[44] = 5
    constants[45] = 0.00013
    states[16] = 0
    constants[46] = 0.5
    constants[47] = 0.5
    constants[48] = 0.7
    constants[49] = 1
    constants[50] = 40
    constants[51] = 0.0005
    constants[52] = 3
    constants[53] = 0
    constants[54] = 0.5
    constants[55] = 0.001
    constants[56] = 0.0003
    constants[57] = 0.4
    constants[58] = 0.5
    constants[59] = 0.4
    constants[60] = 0.03
    states[17] = 0.4531889
    states[18] = 0.4481927
    constants[61] = 0.05
    constants[62] = 250
    constants[63] = 0.0005
    constants[64] = 0.01
    states[19] = 0.0042614
    states[20] = 0.4068154
    states[21] = 0.0005555
    states[22] = 0.0003542
    constants[65] = 0.02
    constants[66] = 0.05
    constants[67] = 100000
    constants[68] = 50
    constants[69] = 100000
    constants[70] = 200
    constants[71] = 12
    constants[72] = 74
    constants[73] = 0.1
    constants[74] = 0.1
    constants[75] = 0.4
    constants[76] = 0.01
    constants[77] = 10
    states[23] = 3.32e-5
    states[24] = 8.09e-5
    constants[78] = 12000
    constants[79] = 100
    constants[80] = 60
    constants[81] = 25
    constants[82] = 2
    constants[83] = 0.05
    constants[84] = (3.14159*(power(constants[71]/1000.00, 2.00000))*constants[72])/1000.00
    constants[85] = (constants[56]*constants[57])/constants[58]
    constants[86] = 0.000200000*exp(2.00000*constants[82])
    constants[87] = custom_piecewise([greater(constants[82] , 1.00000), 1.00000-exp(-3.00000*(constants[82]-1.00000)) , True, 0.00000])
    constants[88] = ((1.00000-constants[75])-constants[76])-constants[74]
    constants[89] = custom_piecewise([greater(constants[82] , 2.00000), 1.00000-0.625000*(constants[82]-2.00000) , True, 1.00000])
    constants[90] = constants[84]*constants[88]
    constants[91] = constants[87]*constants[89]*constants[83]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[12] = 1.00000-1.00000*(states[3]/(constants[33]+states[3])+states[12])
    rates[13] = constants[35]*(1.00000-(states[9]/(constants[34]+states[9])+states[13]))
    rates[15] = 333.000*(1.00000/(1.00000+exp(-(states[0]+4.00000)/5.00000))-states[15])
    rates[16] = constants[46]*(1.00000-states[16])-constants[47]*states[16]
    rates[24] = constants[80]*states[23]*(1.00000-states[24])-constants[81]*states[24]
    algebraic[0] = 50.0000/(1.00000+exp(-(states[0]-5.00000)/9.00000))
    algebraic[12] = 0.0500000*exp(-(states[0]-20.0000)/15.0000)
    rates[4] = algebraic[0]*(1.00000-states[4])-algebraic[12]*states[4]
    algebraic[1] = 50.0000/(1.00000+exp(-(states[0]-5.00000)/9.00000))
    algebraic[13] = 0.400000*exp(-(power((states[0]+30.0000)/30.0000, 3.00000)))
    rates[5] = algebraic[1]*(1.00000-states[5])-algebraic[13]*states[5]
    algebraic[2] = 14.0000/(1.00000+exp(-(states[0]-40.0000)/9.00000))
    algebraic[14] = 1.00000*exp(-states[0]/45.0000)
    rates[6] = algebraic[2]*(1.00000-states[6])-algebraic[14]*states[6]
    algebraic[4] = 20.0000*exp(-0.125000*((states[0]+75.0000)-constants[27]))
    algebraic[16] = 2000.00/(1.00000+320.000*exp(-0.100000*((states[0]+75.0000)-constants[27])))
    rates[8] = algebraic[4]*(1.00000-states[8])-algebraic[16]*states[8]
    algebraic[7] = 0.0330000*exp(-states[0]/17.0000)
    algebraic[19] = 33.0000/(1.00000+exp(-0.125000*(states[0]+10.0000)))
    rates[14] = algebraic[7]*(1.00000-states[14])-algebraic[19]*states[14]
    algebraic[10] = states[21]/constants[65]
    algebraic[22] = states[22]/constants[66]
    rates[23] = constants[78]*(power(algebraic[22], 2.00000))*algebraic[10]*(1.00000-states[23])-constants[79]*states[23]
    algebraic[3] = states[0]+41.0000
    algebraic[15] = custom_piecewise([less(fabs(algebraic[3]) , constants[26]), 2000.00 , True, (200.000*algebraic[3])/(1.00000-exp(-0.100000*algebraic[3]))])
    algebraic[23] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    rates[7] = algebraic[15]*(1.00000-states[7])-algebraic[23]*states[7]
    algebraic[5] = (states[0]+24.0000)-5.00000
    algebraic[17] = custom_piecewise([less(fabs(algebraic[5]) , 0.000100000), 120.000 , True, (30.0000*algebraic[5])/(1.00000-exp(-algebraic[5]/4.00000))])
    algebraic[24] = custom_piecewise([less(fabs(algebraic[5]) , 0.000100000), 120.000 , True, (12.0000*algebraic[5])/(exp(algebraic[5]/10.0000)-1.00000)])
    rates[10] = constants[37]*(algebraic[17]*(1.00000-states[10])-algebraic[24]*states[10])
    algebraic[6] = states[0]+34.0000
    algebraic[18] = custom_piecewise([less(fabs(algebraic[6]) , constants[39]), 25.0000 , True, (6.25000*algebraic[6])/(exp(algebraic[6]/4.00000)-1.00000)])
    algebraic[25] = 12.0000/(1.00000+exp((-1.00000*(states[0]+34.0000))/4.00000))
    rates[11] = constants[38]*(algebraic[18]*(1.00000-states[11])-algebraic[25]*states[11])
    algebraic[21] = states[3]/(states[3]+constants[63])
    algebraic[27] = states[9]/(states[9]+constants[64])
    algebraic[30] = algebraic[21]+(1.00000-algebraic[21])*algebraic[27]
    algebraic[34] = 60.0000+500.000*(power(algebraic[30], 2.00000))
    algebraic[36] = custom_piecewise([less(states[0] , -50.0000), 5.00000 , True, 1.00000])
    rates[20] = states[19]*algebraic[36]*algebraic[34]-algebraic[36]*1.00000*states[20]
    algebraic[9] = exp(0.0800000*(states[0]-40.0000))
    algebraic[32] = 0.00000*algebraic[9]+500.000*(power(algebraic[30], 2.00000))
    algebraic[38] = (1.00000-states[19])-states[20]
    rates[19] = algebraic[38]*algebraic[36]*algebraic[32]-states[19]*algebraic[36]*algebraic[34]
    algebraic[47] = (((constants[36]*4.00000*constants[30]*states[10]*states[11]*states[13]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1]))))*(states[3]*exp((100.000*constants[2])/(constants[0]*constants[1]))-constants[12]*exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    rates[9] = (-1.00000*algebraic[47])/(2.00000*1.00000*constants[73]*constants[90]*constants[2])-states[9]*constants[77]
    algebraic[26] = ((constants[0]*constants[1])/constants[2])*log(constants[10]/states[1])
    algebraic[35] = (((constants[14]*constants[10])/(constants[10]+constants[13]))*(states[0]-algebraic[26]))/(1.00000+exp((((states[0]-algebraic[26])-constants[15])*constants[2]*constants[16])/(constants[0]*constants[1])))
    algebraic[52] = constants[41]*(constants[42]+states[14]*(1.00000-constants[42]))*states[15]*(states[0]-algebraic[26])
    algebraic[37] = (((constants[17]*states[4]+constants[18]*states[5])*1.00000)/(1.00000+exp((states[0]+9.00000)/22.4000)))*(states[0]-algebraic[26])
    algebraic[29] = ((constants[0]*constants[1])/constants[2])*log((constants[10]+constants[9]*constants[11])/(states[1]+constants[9]*states[2]))
    algebraic[39] = constants[19]*(power(states[6], 2.00000))*(states[0]-algebraic[29])
    algebraic[45] = ((((1.00000-constants[36])*constants[31]*constants[30]*states[10]*states[11]*states[12]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[1]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[10]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[48] = (((constants[36]*constants[31]*constants[30]*states[10]*states[11]*states[13]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[1]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[10]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[54] = (((constants[48]*constants[10])/(constants[49]+constants[10]))*states[2])/(constants[50]+states[2])
    rates[1] = (-1.00000/(1.00000*constants[90]*constants[2]))*((algebraic[35]+algebraic[37]+algebraic[39]+algebraic[45]+algebraic[48]+algebraic[52])-2.00000*algebraic[54])
    algebraic[33] = ((constants[0]*constants[1])/constants[2])*log((constants[11]+0.120000*constants[10])/(states[2]+0.120000*states[1]))
    algebraic[41] = constants[25]*(power(states[7], 3.00000))*states[8]*(states[0]-algebraic[33])
    algebraic[20] = ((constants[0]*constants[1])/constants[2])*log(constants[11]/states[2])
    algebraic[43] = constants[29]*(states[0]-algebraic[20])
    algebraic[42] = ((constants[28]*1.00000)/(1.00000+exp(-(states[0]+52.0000)/8.00000)))*(states[0]-algebraic[20])
    algebraic[46] = ((((1.00000-constants[36])*constants[32]*constants[30]*states[10]*states[11]*states[12]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[2]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[11]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[49] = (((constants[36]*constants[32]*constants[30]*states[10]*states[11]*states[13]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[2]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[11]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[55] = ((1.00000-constants[55])*constants[51]*(exp((constants[54]*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], constants[52]))*constants[12]-exp(((constants[54]-1.00000)*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[11], constants[52]))*states[3]))/((1.00000+constants[53]*(states[3]*(power(constants[11], constants[52]))+constants[12]*(power(states[2], constants[52]))))*(1.00000+states[3]/0.00690000))
    rates[2] = (-1.00000/(1.00000*constants[90]*constants[2]))*(algebraic[41]+algebraic[42]+algebraic[43]+3.00000*algebraic[54]+3.00000*algebraic[55]+algebraic[46]+algebraic[49])
    algebraic[40] = (constants[20]*(states[0]+80.0000))/(1.00000+power(constants[22]/constants[21], 2.00000))
    algebraic[53] = (constants[43]*(constants[10]/(constants[10]+constants[13]))*states[16]*(1.00000/(1.00000+power(constants[45]/constants[44], 2.00000)))*(states[0]-algebraic[26]))/(1.00000+exp((2.00000*constants[2]*(states[0]-(algebraic[26]+10.0000)))/(constants[0]*constants[1])))
    algebraic[56] = (constants[55]*constants[51]*(exp((constants[54]*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], constants[52]))*constants[12]-exp(((constants[54]-1.00000)*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[11], constants[52]))*states[9]))/((1.00000+constants[53]*(states[9]*(power(constants[11], constants[52]))+constants[12]*(power(states[2], constants[52]))))*(1.00000+states[9]/0.00690000))
    algebraic[44] = ((((1.00000-constants[36])*4.00000*constants[30]*states[10]*states[11]*states[12]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1]))))*(states[3]*exp((100.000*constants[2])/(constants[0]*constants[1]))-constants[12]*exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[31] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[12]/states[3])
    algebraic[51] = constants[40]*(states[0]-algebraic[31])
    algebraic[8] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    rates[0] = (-1.00000/constants[3])*(algebraic[8]+algebraic[35]+algebraic[52]+algebraic[37]+algebraic[39]+algebraic[40]+algebraic[53]+algebraic[54]+algebraic[41]+algebraic[43]+algebraic[42]+algebraic[46]+algebraic[49]+algebraic[55]+algebraic[56]+algebraic[44]+algebraic[47]+algebraic[45]+algebraic[48]+algebraic[51])
    algebraic[58] = states[3]+states[17]*constants[85]+constants[56]*constants[57]+constants[56]
    algebraic[59] = (states[3]/algebraic[58])*constants[59]-((states[17]*constants[85])/algebraic[58])*constants[60]
    algebraic[60] = 50.0000*(states[17]-states[18])
    rates[17] = (constants[88]/constants[76])*algebraic[59]-algebraic[60]
    rates[21] = constants[67]*states[3]*(constants[65]-states[21])-constants[68]*states[21]
    algebraic[61] = ((power(states[19]/(states[19]+0.250000), 2.00000))*constants[62]+constants[61])*states[18]
    rates[18] = (constants[76]/constants[74])*algebraic[60]-algebraic[61]
    rates[22] = constants[69]*states[3]*(constants[66]-states[22])-constants[70]*states[22]
    rates[3] = ((((-1.00000/(2.00000*1.00000*constants[90]*constants[2]))*(((algebraic[44]+algebraic[51])-2.00000*algebraic[55])-2.00000*algebraic[56])+states[9]*constants[73]*constants[77]+(algebraic[61]*constants[74])/constants[88])-rates[21])-rates[22])-algebraic[59]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = 50.0000/(1.00000+exp(-(states[0]-5.00000)/9.00000))
    algebraic[12] = 0.0500000*exp(-(states[0]-20.0000)/15.0000)
    algebraic[1] = 50.0000/(1.00000+exp(-(states[0]-5.00000)/9.00000))
    algebraic[13] = 0.400000*exp(-(power((states[0]+30.0000)/30.0000, 3.00000)))
    algebraic[2] = 14.0000/(1.00000+exp(-(states[0]-40.0000)/9.00000))
    algebraic[14] = 1.00000*exp(-states[0]/45.0000)
    algebraic[4] = 20.0000*exp(-0.125000*((states[0]+75.0000)-constants[27]))
    algebraic[16] = 2000.00/(1.00000+320.000*exp(-0.100000*((states[0]+75.0000)-constants[27])))
    algebraic[7] = 0.0330000*exp(-states[0]/17.0000)
    algebraic[19] = 33.0000/(1.00000+exp(-0.125000*(states[0]+10.0000)))
    algebraic[10] = states[21]/constants[65]
    algebraic[22] = states[22]/constants[66]
    algebraic[3] = states[0]+41.0000
    algebraic[15] = custom_piecewise([less(fabs(algebraic[3]) , constants[26]), 2000.00 , True, (200.000*algebraic[3])/(1.00000-exp(-0.100000*algebraic[3]))])
    algebraic[23] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    algebraic[5] = (states[0]+24.0000)-5.00000
    algebraic[17] = custom_piecewise([less(fabs(algebraic[5]) , 0.000100000), 120.000 , True, (30.0000*algebraic[5])/(1.00000-exp(-algebraic[5]/4.00000))])
    algebraic[24] = custom_piecewise([less(fabs(algebraic[5]) , 0.000100000), 120.000 , True, (12.0000*algebraic[5])/(exp(algebraic[5]/10.0000)-1.00000)])
    algebraic[6] = states[0]+34.0000
    algebraic[18] = custom_piecewise([less(fabs(algebraic[6]) , constants[39]), 25.0000 , True, (6.25000*algebraic[6])/(exp(algebraic[6]/4.00000)-1.00000)])
    algebraic[25] = 12.0000/(1.00000+exp((-1.00000*(states[0]+34.0000))/4.00000))
    algebraic[21] = states[3]/(states[3]+constants[63])
    algebraic[27] = states[9]/(states[9]+constants[64])
    algebraic[30] = algebraic[21]+(1.00000-algebraic[21])*algebraic[27]
    algebraic[34] = 60.0000+500.000*(power(algebraic[30], 2.00000))
    algebraic[36] = custom_piecewise([less(states[0] , -50.0000), 5.00000 , True, 1.00000])
    algebraic[9] = exp(0.0800000*(states[0]-40.0000))
    algebraic[32] = 0.00000*algebraic[9]+500.000*(power(algebraic[30], 2.00000))
    algebraic[38] = (1.00000-states[19])-states[20]
    algebraic[47] = (((constants[36]*4.00000*constants[30]*states[10]*states[11]*states[13]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1]))))*(states[3]*exp((100.000*constants[2])/(constants[0]*constants[1]))-constants[12]*exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[26] = ((constants[0]*constants[1])/constants[2])*log(constants[10]/states[1])
    algebraic[35] = (((constants[14]*constants[10])/(constants[10]+constants[13]))*(states[0]-algebraic[26]))/(1.00000+exp((((states[0]-algebraic[26])-constants[15])*constants[2]*constants[16])/(constants[0]*constants[1])))
    algebraic[52] = constants[41]*(constants[42]+states[14]*(1.00000-constants[42]))*states[15]*(states[0]-algebraic[26])
    algebraic[37] = (((constants[17]*states[4]+constants[18]*states[5])*1.00000)/(1.00000+exp((states[0]+9.00000)/22.4000)))*(states[0]-algebraic[26])
    algebraic[29] = ((constants[0]*constants[1])/constants[2])*log((constants[10]+constants[9]*constants[11])/(states[1]+constants[9]*states[2]))
    algebraic[39] = constants[19]*(power(states[6], 2.00000))*(states[0]-algebraic[29])
    algebraic[45] = ((((1.00000-constants[36])*constants[31]*constants[30]*states[10]*states[11]*states[12]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[1]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[10]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[48] = (((constants[36]*constants[31]*constants[30]*states[10]*states[11]*states[13]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[1]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[10]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[54] = (((constants[48]*constants[10])/(constants[49]+constants[10]))*states[2])/(constants[50]+states[2])
    algebraic[33] = ((constants[0]*constants[1])/constants[2])*log((constants[11]+0.120000*constants[10])/(states[2]+0.120000*states[1]))
    algebraic[41] = constants[25]*(power(states[7], 3.00000))*states[8]*(states[0]-algebraic[33])
    algebraic[20] = ((constants[0]*constants[1])/constants[2])*log(constants[11]/states[2])
    algebraic[43] = constants[29]*(states[0]-algebraic[20])
    algebraic[42] = ((constants[28]*1.00000)/(1.00000+exp(-(states[0]+52.0000)/8.00000)))*(states[0]-algebraic[20])
    algebraic[46] = ((((1.00000-constants[36])*constants[32]*constants[30]*states[10]*states[11]*states[12]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[2]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[11]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[49] = (((constants[36]*constants[32]*constants[30]*states[10]*states[11]*states[13]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[2]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[11]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[55] = ((1.00000-constants[55])*constants[51]*(exp((constants[54]*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], constants[52]))*constants[12]-exp(((constants[54]-1.00000)*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[11], constants[52]))*states[3]))/((1.00000+constants[53]*(states[3]*(power(constants[11], constants[52]))+constants[12]*(power(states[2], constants[52]))))*(1.00000+states[3]/0.00690000))
    algebraic[40] = (constants[20]*(states[0]+80.0000))/(1.00000+power(constants[22]/constants[21], 2.00000))
    algebraic[53] = (constants[43]*(constants[10]/(constants[10]+constants[13]))*states[16]*(1.00000/(1.00000+power(constants[45]/constants[44], 2.00000)))*(states[0]-algebraic[26]))/(1.00000+exp((2.00000*constants[2]*(states[0]-(algebraic[26]+10.0000)))/(constants[0]*constants[1])))
    algebraic[56] = (constants[55]*constants[51]*(exp((constants[54]*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], constants[52]))*constants[12]-exp(((constants[54]-1.00000)*(constants[52]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[11], constants[52]))*states[9]))/((1.00000+constants[53]*(states[9]*(power(constants[11], constants[52]))+constants[12]*(power(states[2], constants[52]))))*(1.00000+states[9]/0.00690000))
    algebraic[44] = ((((1.00000-constants[36])*4.00000*constants[30]*states[10]*states[11]*states[12]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1]))))*(states[3]*exp((100.000*constants[2])/(constants[0]*constants[1]))-constants[12]*exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[31] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[12]/states[3])
    algebraic[51] = constants[40]*(states[0]-algebraic[31])
    algebraic[8] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    algebraic[58] = states[3]+states[17]*constants[85]+constants[56]*constants[57]+constants[56]
    algebraic[59] = (states[3]/algebraic[58])*constants[59]-((states[17]*constants[85])/algebraic[58])*constants[60]
    algebraic[60] = 50.0000*(states[17]-states[18])
    algebraic[61] = ((power(states[19]/(states[19]+0.250000), 2.00000))*constants[62]+constants[61])*states[18]
    algebraic[11] = states[24]*constants[91]+constants[86]
    algebraic[28] = ((constants[23]*states[2])/(states[2]+constants[24]))*(states[0]-algebraic[26])
    algebraic[50] = algebraic[44]+algebraic[45]+algebraic[46]+algebraic[47]+algebraic[48]+algebraic[49]
    algebraic[57] = algebraic[55]+algebraic[56]
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