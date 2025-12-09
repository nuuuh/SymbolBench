# Size of variable arrays:
sizeAlgebraic = 45
sizeStates = 17
sizeConstants = 59
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
    legend_algebraic[27] = "i_K1 in component time_independent_potassium_current (nanoA)"
    legend_algebraic[38] = "i_to in component transient_outward_current (nanoA)"
    legend_algebraic[29] = "i_K in component time_dependent_potassium_current (nanoA)"
    legend_algebraic[34] = "i_Ca_L_K in component L_type_Ca_channel (nanoA)"
    legend_algebraic[30] = "i_b_K in component potassium_background_current (nanoA)"
    legend_algebraic[39] = "i_NaK in component sodium_potassium_pump (nanoA)"
    legend_algebraic[31] = "i_Na in component fast_sodium_current (nanoA)"
    legend_algebraic[32] = "i_b_Na in component sodium_background_current (nanoA)"
    legend_algebraic[35] = "i_Ca_L_Na in component L_type_Ca_channel (nanoA)"
    legend_algebraic[40] = "i_NaCa in component sodium_calcium_exchanger (nanoA)"
    legend_algebraic[33] = "i_Ca_L_Ca in component L_type_Ca_channel (nanoA)"
    legend_algebraic[37] = "i_b_Ca in component calcium_background_current (nanoA)"
    legend_algebraic[6] = "i_Stim in component membrane (nanoA)"
    legend_constants[4] = "stim_start in component membrane (second)"
    legend_constants[5] = "stim_end in component membrane (second)"
    legend_constants[6] = "stim_period in component membrane (second)"
    legend_constants[7] = "stim_duration in component membrane (second)"
    legend_constants[8] = "stim_amplitude in component membrane (nanoA)"
    legend_algebraic[15] = "E_Na in component reversal_potentials (millivolt)"
    legend_algebraic[21] = "E_K in component reversal_potentials (millivolt)"
    legend_algebraic[24] = "E_Ca in component reversal_potentials (millivolt)"
    legend_algebraic[26] = "E_mh in component reversal_potentials (millivolt)"
    legend_constants[9] = "K_o in component extracellular_potassium_concentration (millimolar)"
    legend_constants[10] = "Na_o in component extracellular_sodium_concentration (millimolar)"
    legend_states[1] = "K_i in component intracellular_potassium_concentration (millimolar)"
    legend_states[2] = "Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_constants[11] = "Ca_o in component extracellular_calcium_concentration (millimolar)"
    legend_states[3] = "Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_constants[12] = "K_mk1 in component time_independent_potassium_current (millimolar)"
    legend_constants[13] = "g_K1 in component time_independent_potassium_current (microS)"
    legend_algebraic[28] = "I_K in component time_dependent_potassium_current (nanoA)"
    legend_constants[14] = "i_K_max in component time_dependent_potassium_current (nanoA)"
    legend_states[4] = "x in component time_dependent_potassium_current_x_gate (dimensionless)"
    legend_algebraic[0] = "E0xa in component time_dependent_potassium_current_x_gate (millivolt)"
    legend_algebraic[17] = "E0xb in component time_dependent_potassium_current_x_gate (millivolt)"
    legend_algebraic[9] = "alpha_x in component time_dependent_potassium_current_x_gate (per_second)"
    legend_algebraic[23] = "beta_x in component time_dependent_potassium_current_x_gate (per_second)"
    legend_constants[15] = "g_bk in component potassium_background_current (microS)"
    legend_constants[16] = "g_Na in component fast_sodium_current (microS)"
    legend_states[5] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[6] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_algebraic[10] = "alpha_m in component fast_sodium_current_m_gate (per_second)"
    legend_algebraic[18] = "beta_m in component fast_sodium_current_m_gate (per_second)"
    legend_constants[17] = "delta_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[1] = "E0_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[2] = "alpha_h in component fast_sodium_current_h_gate (per_second)"
    legend_algebraic[11] = "beta_h in component fast_sodium_current_h_gate (per_second)"
    legend_constants[18] = "shift_h in component fast_sodium_current_h_gate (millivolt)"
    legend_constants[19] = "g_bna in component sodium_background_current (microS)"
    legend_algebraic[36] = "i_Ca_L in component L_type_Ca_channel (nanoA)"
    legend_constants[20] = "P_Ca_L in component L_type_Ca_channel (nanoA_per_millimolar)"
    legend_constants[21] = "P_CaK in component L_type_Ca_channel (dimensionless)"
    legend_constants[22] = "P_CaNa in component L_type_Ca_channel (dimensionless)"
    legend_states[7] = "d in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_states[8] = "f in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[12] = "alpha_d in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[19] = "beta_d in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[3] = "E0_d in component L_type_Ca_channel_d_gate (millivolt)"
    legend_constants[23] = "speed_d in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_algebraic[13] = "alpha_f in component L_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[20] = "beta_f in component L_type_Ca_channel_f_gate (per_second)"
    legend_constants[24] = "speed_f in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_constants[25] = "delta_f in component L_type_Ca_channel_f_gate (millivolt)"
    legend_algebraic[4] = "E0_f in component L_type_Ca_channel_f_gate (millivolt)"
    legend_constants[26] = "g_bca in component calcium_background_current (microS)"
    legend_constants[27] = "g_to in component transient_outward_current (microS)"
    legend_constants[28] = "g_tos in component transient_outward_current (dimensionless)"
    legend_states[9] = "s in component transient_outward_current_s_gate (dimensionless)"
    legend_states[10] = "r in component transient_outward_current_r_gate (dimensionless)"
    legend_algebraic[5] = "alpha_s in component transient_outward_current_s_gate (per_second)"
    legend_algebraic[14] = "beta_s in component transient_outward_current_s_gate (per_second)"
    legend_constants[29] = "i_NaK_max in component sodium_potassium_pump (nanoA)"
    legend_constants[30] = "K_mK in component sodium_potassium_pump (millimolar)"
    legend_constants[31] = "K_mNa in component sodium_potassium_pump (millimolar)"
    legend_constants[32] = "k_NaCa in component sodium_calcium_exchanger (nanoA)"
    legend_constants[33] = "n_NaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[34] = "d_NaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[35] = "gamma in component sodium_calcium_exchanger (dimensionless)"
    legend_algebraic[42] = "i_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[56] = "K_1 in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_algebraic[41] = "K_2 in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[36] = "K_cyca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[37] = "K_xcs in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_constants[38] = "K_srca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[39] = "alpha_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[40] = "beta_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_states[11] = "Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[43] = "i_trans in component calcium_translocation (millimolar_per_second)"
    legend_states[12] = "Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[44] = "i_rel in component calcium_release (millimolar_per_second)"
    legend_algebraic[7] = "VoltDep in component calcium_release (dimensionless)"
    legend_algebraic[8] = "RegBindSite in component calcium_release (dimensionless)"
    legend_algebraic[16] = "ActRate in component calcium_release (per_second)"
    legend_algebraic[22] = "InactRate in component calcium_release (per_second)"
    legend_constants[41] = "K_leak_rate in component calcium_release (per_second)"
    legend_constants[42] = "K_m_Ca in component calcium_release (millimolar)"
    legend_constants[43] = "K_m_rel in component calcium_release (per_second)"
    legend_algebraic[25] = "PrecFrac in component calcium_release (dimensionless)"
    legend_states[13] = "ActFrac in component calcium_release (dimensionless)"
    legend_states[14] = "ProdFrac in component calcium_release (dimensionless)"
    legend_constants[58] = "V_i in component intracellular_calcium_concentration (micrometre3)"
    legend_states[15] = "Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_states[16] = "Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    legend_constants[44] = "Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_constants[45] = "Trop in component intracellular_calcium_concentration (millimolar)"
    legend_constants[46] = "alpha_Calmod in component intracellular_calcium_concentration (per_millimolar_second)"
    legend_constants[47] = "beta_Calmod in component intracellular_calcium_concentration (per_second)"
    legend_constants[48] = "alpha_Trop in component intracellular_calcium_concentration (per_millimolar_second)"
    legend_constants[49] = "beta_Trop in component intracellular_calcium_concentration (per_second)"
    legend_constants[50] = "radius in component intracellular_calcium_concentration (micrometre)"
    legend_constants[51] = "length in component intracellular_calcium_concentration (micrometre)"
    legend_constants[55] = "V_Cell in component intracellular_calcium_concentration (micrometre3)"
    legend_constants[57] = "V_i_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[52] = "V_rel_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[53] = "V_e_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[54] = "V_up_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[4] = "d/dt x in component time_dependent_potassium_current_x_gate (dimensionless)"
    legend_rates[5] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[6] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[7] = "d/dt d in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_rates[8] = "d/dt f in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_rates[9] = "d/dt s in component transient_outward_current_s_gate (dimensionless)"
    legend_rates[10] = "d/dt r in component transient_outward_current_r_gate (dimensionless)"
    legend_rates[13] = "d/dt ActFrac in component calcium_release (dimensionless)"
    legend_rates[14] = "d/dt ProdFrac in component calcium_release (dimensionless)"
    legend_rates[2] = "d/dt Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_rates[1] = "d/dt K_i in component intracellular_potassium_concentration (millimolar)"
    legend_rates[3] = "d/dt Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_rates[15] = "d/dt Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_rates[16] = "d/dt Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    legend_rates[11] = "d/dt Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_rates[12] = "d/dt Ca_rel in component intracellular_calcium_concentration (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -93.7400119196694
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 9.5e-5
    constants[4] = 0.1
    constants[5] = 100
    constants[6] = 1
    constants[7] = 0.002
    constants[8] = -6
    constants[9] = 4
    constants[10] = 140
    states[1] = 136.604284305878
    states[2] = 7.50547214142684
    constants[11] = 2
    states[3] = 1.34858164771406e-5
    constants[12] = 10
    constants[13] = 1
    constants[14] = 1
    states[4] = 0.00938586574433011
    constants[15] = 0.0006
    constants[16] = 2.5
    states[5] = 0.00143405969732302
    states[6] = 0.995414125415674
    constants[17] = 1e-5
    constants[18] = 0
    constants[19] = 0.0006
    constants[20] = 0.25
    constants[21] = 0.002
    constants[22] = 0.01
    states[7] = 1.91821833548952e-8
    states[8] = 0.999999956287155
    constants[23] = 3
    constants[24] = 0.5
    constants[25] = 0.0001
    constants[26] = 0.00025
    constants[27] = 0.005
    constants[28] = 0
    states[9] = 0.997644968939185
    states[10] = 1.60424507876553e-8
    constants[29] = 0.7
    constants[30] = 1
    constants[31] = 40
    constants[32] = 0.0005
    constants[33] = 3
    constants[34] = 0
    constants[35] = 0.5
    constants[36] = 0.0003
    constants[37] = 0.4
    constants[38] = 0.5
    constants[39] = 0.4
    constants[40] = 0.03
    states[11] = 0.59333810408885
    states[12] = 0.591323137897127
    constants[41] = 0
    constants[42] = 0.0005
    constants[43] = 250
    states[13] = 0.00267040300939318
    states[14] = 0.522949441962453
    states[15] = 0.000524570960945961
    states[16] = 0.00033477224086766
    constants[44] = 0.02
    constants[45] = 0.05
    constants[46] = 100000
    constants[47] = 50
    constants[48] = 100000
    constants[49] = 200
    constants[50] = 0.012
    constants[51] = 0.074
    constants[52] = 0.1
    constants[53] = 0.4
    constants[54] = 0.01
    constants[55] = 3.14159*(power(constants[50], 2.00000))*constants[51]
    constants[56] = (constants[36]*constants[37])/constants[38]
    constants[57] = ((1.00000-constants[53])-constants[54])-constants[52]
    constants[58] = constants[55]*constants[57]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[10] = 333.000*(1.00000/(1.00000+exp(-(states[0]+4.00000)/5.00000))-states[10])
    algebraic[2] = 20.0000*exp(-0.125000*((states[0]+75.0000)-constants[18]))
    algebraic[11] = 2000.00/(1.00000+320.000*exp(-0.100000*((states[0]+75.0000)-constants[18])))
    rates[6] = algebraic[2]*(1.00000-states[6])-algebraic[11]*states[6]
    algebraic[5] = 0.0330000*exp(-states[0]/17.0000)
    algebraic[14] = 33.0000/(1.00000+exp(-0.125000*(states[0]+10.0000)))
    rates[9] = algebraic[5]*(1.00000-states[9])-algebraic[14]*states[9]
    algebraic[1] = states[0]+41.0000
    algebraic[10] = custom_piecewise([less(fabs(algebraic[1]) , constants[17]), 2000.00 , True, (200.000*algebraic[1])/(1.00000-exp(-0.100000*algebraic[1]))])
    algebraic[18] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    rates[5] = algebraic[10]*(1.00000-states[5])-algebraic[18]*states[5]
    algebraic[3] = (states[0]+24.0000)-5.00000
    algebraic[12] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (30.0000*algebraic[3])/(1.00000-exp(-algebraic[3]/4.00000))])
    algebraic[19] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (12.0000*algebraic[3])/(exp(algebraic[3]/10.0000)-1.00000)])
    rates[7] = constants[23]*(algebraic[12]*(1.00000-states[7])-algebraic[19]*states[7])
    algebraic[4] = states[0]+34.0000
    algebraic[13] = custom_piecewise([less(fabs(algebraic[4]) , constants[25]), 25.0000 , True, (6.25000*algebraic[4])/(exp(algebraic[4]/4.00000)-1.00000)])
    algebraic[20] = 50.0000/(1.00000+exp((-1.00000*(states[0]+34.0000))/4.00000))
    rates[8] = constants[24]*(algebraic[13]*(1.00000-states[8])-algebraic[20]*states[8])
    algebraic[8] = power(states[3]/(states[3]+constants[42]), 2.00000)
    algebraic[22] = 60.0000+500.000*algebraic[8]
    rates[14] = states[13]*algebraic[22]-1.00000*states[14]
    algebraic[0] = states[0]+50.0000
    algebraic[9] = (0.500000*exp(0.0826000*algebraic[0]))/(1.00000+exp(0.0570000*algebraic[0]))
    algebraic[17] = states[0]+20.0000
    algebraic[23] = (1.30000*exp(-0.0600000*algebraic[17]))/(1.00000+exp(-0.0400000*algebraic[17]))
    rates[4] = algebraic[9]*(1.00000-states[4])-algebraic[23]*states[4]
    algebraic[16] = 500.000*algebraic[8]
    algebraic[25] = (1.00000-states[13])-states[14]
    rates[13] = algebraic[25]*algebraic[16]-states[13]*algebraic[22]
    algebraic[21] = ((constants[0]*constants[1])/constants[2])*log(constants[9]/states[1])
    algebraic[27] = (((constants[13]*constants[9])/(constants[9]+constants[12]))*(states[0]-algebraic[21]))/(1.00000+exp((((states[0]-algebraic[21])-10.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[38] = constants[27]*(constants[28]+states[9]*(1.00000-constants[28]))*states[10]*(states[0]-algebraic[21])
    algebraic[28] = (constants[14]*(states[1]-constants[9]*exp((-states[0]*constants[2])/(constants[0]*constants[1]))))/140.000
    algebraic[29] = states[4]*algebraic[28]
    algebraic[34] = (((constants[21]*constants[20]*states[7]*states[8]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[1]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[9]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[30] = constants[15]*(states[0]-algebraic[21])
    algebraic[39] = (((constants[29]*constants[9])/(constants[30]+constants[9]))*states[2])/(constants[31]+states[2])
    rates[1] = (-1.00000/(1.00000*constants[58]*constants[2]))*((algebraic[27]+algebraic[29]+algebraic[34]+algebraic[38]+algebraic[30])-2.00000*algebraic[39])
    algebraic[26] = ((constants[0]*constants[1])/constants[2])*log((constants[10]+0.120000*constants[9])/(states[2]+0.120000*states[1]))
    algebraic[31] = constants[16]*(power(states[5], 3.00000))*states[6]*(states[0]-algebraic[26])
    algebraic[15] = ((constants[0]*constants[1])/constants[2])*log(constants[10]/states[2])
    algebraic[32] = constants[19]*(states[0]-algebraic[15])
    algebraic[35] = (((constants[22]*constants[20]*states[7]*states[8]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[2]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[10]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[40] = (constants[32]*(exp((constants[35]*(constants[33]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], constants[33]))*constants[11]-exp(((constants[35]-1.00000)*(constants[33]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[10], constants[33]))*states[3]))/((1.00000+constants[34]*(states[3]*(power(constants[10], constants[33]))+constants[11]*(power(states[2], constants[33]))))*(1.00000+states[3]/0.00690000))
    algebraic[33] = (((4.00000*constants[20]*states[7]*states[8]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1]))))*(states[3]*exp((100.000*constants[2])/(constants[0]*constants[1]))-constants[11]*exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[24] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[11]/states[3])
    algebraic[37] = constants[26]*(states[0]-algebraic[24])
    algebraic[6] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    rates[0] = (-1.00000/constants[3])*(algebraic[6]+algebraic[27]+algebraic[38]+algebraic[29]+algebraic[30]+algebraic[39]+algebraic[31]+algebraic[32]+algebraic[35]+algebraic[40]+algebraic[33]+algebraic[34]+algebraic[37])
    rates[2] = (-1.00000/(1.00000*constants[58]*constants[2]))*(algebraic[31]+algebraic[32]+3.00000*algebraic[39]+3.00000*algebraic[40]+algebraic[35])
    algebraic[41] = states[3]+states[11]*constants[56]+constants[36]*constants[37]+constants[36]
    algebraic[42] = (states[3]/algebraic[41])*constants[39]-((states[11]*constants[56])/algebraic[41])*constants[40]
    algebraic[43] = 50.0000*(states[11]-states[12])
    rates[11] = (constants[57]/constants[54])*algebraic[42]-algebraic[43]
    rates[15] = constants[46]*states[3]*(constants[44]-states[15])-constants[47]*states[15]
    algebraic[44] = ((power(states[13]/(states[13]+0.250000), 2.00000))*constants[43]+constants[41])*states[12]
    rates[12] = (constants[54]/constants[52])*algebraic[43]-algebraic[44]
    rates[16] = constants[48]*states[3]*(constants[45]-states[16])-constants[49]*states[16]
    rates[3] = ((((-1.00000/(2.00000*1.00000*constants[58]*constants[2]))*((algebraic[33]+algebraic[37])-2.00000*algebraic[40])+(algebraic[44]*constants[52])/constants[57])-rates[15])-rates[16])-algebraic[42]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[2] = 20.0000*exp(-0.125000*((states[0]+75.0000)-constants[18]))
    algebraic[11] = 2000.00/(1.00000+320.000*exp(-0.100000*((states[0]+75.0000)-constants[18])))
    algebraic[5] = 0.0330000*exp(-states[0]/17.0000)
    algebraic[14] = 33.0000/(1.00000+exp(-0.125000*(states[0]+10.0000)))
    algebraic[1] = states[0]+41.0000
    algebraic[10] = custom_piecewise([less(fabs(algebraic[1]) , constants[17]), 2000.00 , True, (200.000*algebraic[1])/(1.00000-exp(-0.100000*algebraic[1]))])
    algebraic[18] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    algebraic[3] = (states[0]+24.0000)-5.00000
    algebraic[12] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (30.0000*algebraic[3])/(1.00000-exp(-algebraic[3]/4.00000))])
    algebraic[19] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (12.0000*algebraic[3])/(exp(algebraic[3]/10.0000)-1.00000)])
    algebraic[4] = states[0]+34.0000
    algebraic[13] = custom_piecewise([less(fabs(algebraic[4]) , constants[25]), 25.0000 , True, (6.25000*algebraic[4])/(exp(algebraic[4]/4.00000)-1.00000)])
    algebraic[20] = 50.0000/(1.00000+exp((-1.00000*(states[0]+34.0000))/4.00000))
    algebraic[8] = power(states[3]/(states[3]+constants[42]), 2.00000)
    algebraic[22] = 60.0000+500.000*algebraic[8]
    algebraic[0] = states[0]+50.0000
    algebraic[9] = (0.500000*exp(0.0826000*algebraic[0]))/(1.00000+exp(0.0570000*algebraic[0]))
    algebraic[17] = states[0]+20.0000
    algebraic[23] = (1.30000*exp(-0.0600000*algebraic[17]))/(1.00000+exp(-0.0400000*algebraic[17]))
    algebraic[16] = 500.000*algebraic[8]
    algebraic[25] = (1.00000-states[13])-states[14]
    algebraic[21] = ((constants[0]*constants[1])/constants[2])*log(constants[9]/states[1])
    algebraic[27] = (((constants[13]*constants[9])/(constants[9]+constants[12]))*(states[0]-algebraic[21]))/(1.00000+exp((((states[0]-algebraic[21])-10.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[38] = constants[27]*(constants[28]+states[9]*(1.00000-constants[28]))*states[10]*(states[0]-algebraic[21])
    algebraic[28] = (constants[14]*(states[1]-constants[9]*exp((-states[0]*constants[2])/(constants[0]*constants[1]))))/140.000
    algebraic[29] = states[4]*algebraic[28]
    algebraic[34] = (((constants[21]*constants[20]*states[7]*states[8]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[1]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[9]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[30] = constants[15]*(states[0]-algebraic[21])
    algebraic[39] = (((constants[29]*constants[9])/(constants[30]+constants[9]))*states[2])/(constants[31]+states[2])
    algebraic[26] = ((constants[0]*constants[1])/constants[2])*log((constants[10]+0.120000*constants[9])/(states[2]+0.120000*states[1]))
    algebraic[31] = constants[16]*(power(states[5], 3.00000))*states[6]*(states[0]-algebraic[26])
    algebraic[15] = ((constants[0]*constants[1])/constants[2])*log(constants[10]/states[2])
    algebraic[32] = constants[19]*(states[0]-algebraic[15])
    algebraic[35] = (((constants[22]*constants[20]*states[7]*states[8]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))))*(states[2]*exp((50.0000*constants[2])/(constants[0]*constants[1]))-constants[10]*exp((-(states[0]-50.0000)*constants[2])/(constants[0]*constants[1])))
    algebraic[40] = (constants[32]*(exp((constants[35]*(constants[33]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(states[2], constants[33]))*constants[11]-exp(((constants[35]-1.00000)*(constants[33]-2.00000)*states[0]*constants[2])/(constants[0]*constants[1]))*(power(constants[10], constants[33]))*states[3]))/((1.00000+constants[34]*(states[3]*(power(constants[10], constants[33]))+constants[11]*(power(states[2], constants[33]))))*(1.00000+states[3]/0.00690000))
    algebraic[33] = (((4.00000*constants[20]*states[7]*states[8]*(states[0]-50.0000)*constants[2])/(constants[0]*constants[1]))/(1.00000-exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1]))))*(states[3]*exp((100.000*constants[2])/(constants[0]*constants[1]))-constants[11]*exp((-(states[0]-50.0000)*constants[2]*2.00000)/(constants[0]*constants[1])))
    algebraic[24] = ((0.500000*constants[0]*constants[1])/constants[2])*log(constants[11]/states[3])
    algebraic[37] = constants[26]*(states[0]-algebraic[24])
    algebraic[6] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    algebraic[41] = states[3]+states[11]*constants[56]+constants[36]*constants[37]+constants[36]
    algebraic[42] = (states[3]/algebraic[41])*constants[39]-((states[11]*constants[56])/algebraic[41])*constants[40]
    algebraic[43] = 50.0000*(states[11]-states[12])
    algebraic[44] = ((power(states[13]/(states[13]+0.250000), 2.00000))*constants[43]+constants[41])*states[12]
    algebraic[7] = exp(0.0800000*(states[0]-40.0000))
    algebraic[36] = algebraic[33]+algebraic[34]+algebraic[35]
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