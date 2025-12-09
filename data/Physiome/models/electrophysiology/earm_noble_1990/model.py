# Size of variable arrays:
sizeAlgebraic = 41
sizeStates = 16
sizeConstants = 54
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
    legend_constants[49] = "RTONF in component membrane (millivolt)"
    legend_constants[3] = "C_m in component membrane (microF)"
    legend_algebraic[26] = "i_b_K in component potassium_background_current (nanoA)"
    legend_algebraic[27] = "i_K1 in component time_independent_potassium_current (nanoA)"
    legend_algebraic[18] = "i_to in component transient_outward_current (nanoA)"
    legend_algebraic[22] = "i_b_Na in component sodium_background_current (nanoA)"
    legend_algebraic[24] = "i_b_Ca in component calcium_background_current (nanoA)"
    legend_algebraic[20] = "i_NaK in component sodium_potassium_pump (nanoA)"
    legend_algebraic[25] = "i_NaCa in component Na_Ca_exchanger (nanoA)"
    legend_algebraic[14] = "i_Na in component fast_sodium_current (nanoA)"
    legend_algebraic[38] = "i_Ca_L in component L_type_calcium_current (nanoA)"
    legend_algebraic[4] = "i_Stim in component membrane (nanoA)"
    legend_constants[4] = "stim_start in component membrane (second)"
    legend_constants[5] = "stim_end in component membrane (second)"
    legend_constants[6] = "stim_period in component membrane (second)"
    legend_constants[7] = "stim_duration in component membrane (second)"
    legend_constants[8] = "stim_amplitude in component membrane (nanoA)"
    legend_constants[9] = "g_Na in component fast_sodium_current (microS)"
    legend_algebraic[10] = "E_mh in component fast_sodium_current (millivolt)"
    legend_constants[10] = "Na_o in component extracellular_sodium_concentration (millimolar)"
    legend_states[1] = "Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_constants[11] = "K_c in component extracellular_potassium_concentration (millimolar)"
    legend_states[2] = "K_i in component intracellular_potassium_concentration (millimolar)"
    legend_states[3] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[4] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_algebraic[6] = "alpha_m in component fast_sodium_current_m_gate (per_second)"
    legend_algebraic[12] = "beta_m in component fast_sodium_current_m_gate (per_second)"
    legend_constants[12] = "delta_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[0] = "E0_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[1] = "alpha_h in component fast_sodium_current_h_gate (per_second)"
    legend_algebraic[7] = "beta_h in component fast_sodium_current_h_gate (per_second)"
    legend_constants[13] = "g_to in component transient_outward_current (microS)"
    legend_algebraic[16] = "E_K in component transient_outward_current (millivolt)"
    legend_constants[14] = "g_to_s in component transient_outward_current (dimensionless)"
    legend_states[5] = "r in component transient_outward_current_r_gate (dimensionless)"
    legend_states[6] = "s in component transient_outward_current_s_gate (dimensionless)"
    legend_algebraic[2] = "alpha_s in component transient_outward_current_s_gate (per_second)"
    legend_algebraic[8] = "beta_s in component transient_outward_current_s_gate (per_second)"
    legend_constants[15] = "i_NaK_max in component sodium_potassium_pump (nanoA)"
    legend_constants[16] = "K_mK in component sodium_potassium_pump (millimolar)"
    legend_constants[17] = "K_mNa in component sodium_potassium_pump (millimolar)"
    legend_algebraic[21] = "E_Na in component sodium_background_current (millivolt)"
    legend_constants[18] = "g_b_Na in component sodium_background_current (microS)"
    legend_algebraic[23] = "E_Ca in component calcium_background_current (millivolt)"
    legend_constants[19] = "g_b_Ca in component calcium_background_current (microS)"
    legend_constants[20] = "Ca_o in component extracellular_calcium_concentration (millimolar)"
    legend_states[7] = "Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_constants[21] = "k_NaCa in component Na_Ca_exchanger (nanoA)"
    legend_constants[22] = "n_NaCa in component Na_Ca_exchanger (dimensionless)"
    legend_constants[23] = "d_NaCa in component Na_Ca_exchanger (dimensionless)"
    legend_constants[24] = "gamma in component Na_Ca_exchanger (dimensionless)"
    legend_constants[25] = "g_b_K in component potassium_background_current (microS)"
    legend_constants[26] = "g_K1 in component time_independent_potassium_current (microS)"
    legend_constants[27] = "K_m_K1 in component time_independent_potassium_current (millimolar)"
    legend_algebraic[33] = "i_Ca_L_Ca in component L_type_calcium_current (nanoA)"
    legend_algebraic[34] = "i_Ca_L_K in component L_type_calcium_current (nanoA)"
    legend_algebraic[36] = "i_Ca_L_Na in component L_type_calcium_current (nanoA)"
    legend_constants[28] = "P_Ca_L in component L_type_calcium_current (nanoA_per_millimolar)"
    legend_states[8] = "d in component L_type_calcium_current_d_gate (dimensionless)"
    legend_states[9] = "f_Ca in component L_type_calcium_current_f_Ca_gate (dimensionless)"
    legend_algebraic[32] = "CaChon in component L_type_calcium_current_f_Ca_gate (dimensionless)"
    legend_algebraic[9] = "alpha_d in component L_type_calcium_current_d_gate (per_second)"
    legend_algebraic[13] = "beta_d in component L_type_calcium_current_d_gate (per_second)"
    legend_algebraic[3] = "E0_d in component L_type_calcium_current_d_gate (millivolt)"
    legend_algebraic[29] = "alpha_f_Ca in component L_type_calcium_current_f_Ca_gate (per_second)"
    legend_algebraic[30] = "beta_f_Ca in component L_type_calcium_current_f_Ca_gate (per_second)"
    legend_algebraic[31] = "CaChoff in component L_type_calcium_current_f_Ca_gate (dimensionless)"
    legend_algebraic[28] = "E0_f in component L_type_calcium_current_f_Ca_gate (millivolt)"
    legend_algebraic[37] = "i_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[51] = "K_1 in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_algebraic[35] = "K_2 in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[29] = "K_cyca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[30] = "K_xcs in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_constants[31] = "K_srca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[32] = "alpha_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[33] = "beta_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_states[10] = "Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[39] = "i_rel in component calcium_release (millimolar_per_second)"
    legend_algebraic[11] = "VoltDep in component calcium_release (dimensionless)"
    legend_algebraic[15] = "RegBindSite in component calcium_release (dimensionless)"
    legend_algebraic[17] = "ActRate in component calcium_release (per_second)"
    legend_algebraic[19] = "InactRate in component calcium_release (per_second)"
    legend_constants[34] = "K_leak_rate in component calcium_release (per_second)"
    legend_constants[35] = "K_m_rel in component calcium_release (per_second)"
    legend_algebraic[5] = "PrecFrac in component calcium_release (dimensionless)"
    legend_states[11] = "ActFrac in component calcium_release (dimensionless)"
    legend_states[12] = "ProdFrac in component calcium_release (dimensionless)"
    legend_constants[36] = "ProdFracRate in component calcium_release (per_second)"
    legend_states[13] = "Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[40] = "i_trans in component calcium_translocation (millimolar_per_second)"
    legend_constants[37] = "alpha_tr in component calcium_translocation (per_second)"
    legend_constants[53] = "V_i in component intracellular_calcium_concentration (micrometre3)"
    legend_states[14] = "Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_states[15] = "Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    legend_constants[38] = "Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_constants[39] = "Trop in component intracellular_calcium_concentration (millimolar)"
    legend_constants[40] = "alpha_Calmod in component intracellular_calcium_concentration (per_millimolar_second)"
    legend_constants[41] = "beta_Calmod in component intracellular_calcium_concentration (per_second)"
    legend_constants[42] = "alpha_Trop in component intracellular_calcium_concentration (per_millimolar_second)"
    legend_constants[43] = "beta_Trop in component intracellular_calcium_concentration (per_second)"
    legend_constants[44] = "radius in component intracellular_calcium_concentration (micrometre)"
    legend_constants[45] = "length in component intracellular_calcium_concentration (micrometre)"
    legend_constants[50] = "V_Cell in component intracellular_calcium_concentration (micrometre3)"
    legend_constants[52] = "V_i_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[46] = "V_rel_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[47] = "V_e_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[48] = "V_up_ratio in component intracellular_calcium_concentration (dimensionless)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[3] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[4] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[5] = "d/dt r in component transient_outward_current_r_gate (dimensionless)"
    legend_rates[6] = "d/dt s in component transient_outward_current_s_gate (dimensionless)"
    legend_rates[8] = "d/dt d in component L_type_calcium_current_d_gate (dimensionless)"
    legend_rates[9] = "d/dt f_Ca in component L_type_calcium_current_f_Ca_gate (dimensionless)"
    legend_rates[11] = "d/dt ActFrac in component calcium_release (dimensionless)"
    legend_rates[12] = "d/dt ProdFrac in component calcium_release (dimensionless)"
    legend_rates[1] = "d/dt Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_rates[2] = "d/dt K_i in component intracellular_potassium_concentration (millimolar)"
    legend_rates[7] = "d/dt Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_rates[14] = "d/dt Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_rates[15] = "d/dt Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    legend_rates[10] = "d/dt Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_rates[13] = "d/dt Ca_rel in component intracellular_calcium_concentration (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -91.6
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 4e-5
    constants[4] = 0.1
    constants[5] = 100000
    constants[6] = 1
    constants[7] = 0.002
    constants[8] = -1.3
    constants[9] = 0.5
    constants[10] = 140
    states[1] = 6.48
    constants[11] = 4
    states[2] = 140
    states[3] = 0.076
    states[4] = 0.015
    constants[12] = 1e-5
    constants[13] = 0.01
    constants[14] = 0
    states[5] = 0
    states[6] = 1
    constants[15] = 0.14
    constants[16] = 1
    constants[17] = 40
    constants[18] = 0.00012
    constants[19] = 5e-5
    constants[20] = 2
    states[7] = 1e-5
    constants[21] = 0.0001
    constants[22] = 3
    constants[23] = 0.0001
    constants[24] = 0.5
    constants[25] = 0.0017
    constants[26] = 0.017
    constants[27] = 10
    constants[28] = 0.05
    states[8] = 0.0011
    states[9] = 0.785
    constants[29] = 0.0003
    constants[30] = 0.4
    constants[31] = 0.5
    constants[32] = 3
    constants[33] = 0.23
    states[10] = 0.3
    constants[34] = 0
    constants[35] = 250
    states[11] = 0
    states[12] = 0
    constants[36] = 1
    states[13] = 0.3
    constants[37] = 50
    states[14] = 0.0005
    states[15] = 0.0015
    constants[38] = 0.02
    constants[39] = 0.15
    constants[40] = 100000
    constants[41] = 50
    constants[42] = 100000
    constants[43] = 200
    constants[44] = 0.01
    constants[45] = 0.08
    constants[46] = 0.1
    constants[47] = 0.4
    constants[48] = 0.01
    constants[49] = (constants[0]*constants[1])/constants[2]
    constants[50] = 3.14159*(power(constants[44], 2.00000))*constants[45]
    constants[51] = (constants[29]*constants[30])/constants[31]
    constants[52] = ((1.00000-constants[47])-constants[48])-constants[46]
    constants[53] = constants[50]*constants[52]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[5] = 333.000*(1.00000/(1.00000+exp(-(states[0]+4.00000)/5.00000))-states[5])
    algebraic[1] = 20.0000*exp(-0.125000*(states[0]+75.0000))
    algebraic[7] = 2000.00/(1.00000+320.000*exp(-0.100000*(states[0]+75.0000)))
    rates[4] = algebraic[1]*(1.00000-states[4])-algebraic[7]*states[4]
    algebraic[2] = 0.0330000*exp(-states[0]/17.0000)
    algebraic[8] = 33.0000/(1.00000+exp(-0.125000*(states[0]+10.0000)))
    rates[6] = algebraic[2]*(1.00000-states[6])-algebraic[8]*states[6]
    algebraic[0] = states[0]+41.0000
    algebraic[6] = custom_piecewise([less(fabs(algebraic[0]) , constants[12]), 2000.00 , True, (200.000*algebraic[0])/(1.00000-exp(-0.100000*algebraic[0]))])
    algebraic[12] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    rates[3] = algebraic[6]*(1.00000-states[3])-algebraic[12]*states[3]
    algebraic[3] = states[0]+19.0000
    algebraic[9] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (30.0000*algebraic[3])/(1.00000-exp(-algebraic[3]/4.00000))])
    algebraic[13] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (12.0000*algebraic[3])/(exp(algebraic[3]/10.0000)-1.00000)])
    rates[8] = algebraic[9]*(1.00000-states[8])-algebraic[13]*states[8]
    algebraic[11] = exp(0.0800000*(states[0]-40.0000))
    algebraic[15] = power(states[7]/(states[7]+0.000500000), 2.00000)
    algebraic[17] = 600.000*algebraic[11]+500.000*algebraic[15]
    algebraic[19] = 60.0000+500.000*algebraic[15]
    algebraic[5] = (1.00000-states[11])-states[12]
    rates[11] = algebraic[5]*algebraic[17]-states[11]*algebraic[19]
    rates[12] = states[11]*algebraic[19]-constants[36]*states[12]
    algebraic[28] = states[0]+34.0000
    algebraic[29] = custom_piecewise([less(fabs(algebraic[28]) , 0.000100000), 25.0000 , True, (6.25000*algebraic[28])/(exp(algebraic[28]/4.00000)-1.00000)])
    algebraic[30] = 12.0000/(1.00000+exp(-algebraic[28]/4.00000))
    algebraic[31] = states[7]/(0.00100000+states[7])
    rates[9] = (120.000*(1.00000-states[9])*algebraic[31]+(1.00000-states[9])*(1.00000-algebraic[31]))*algebraic[30]-algebraic[29]*states[9]
    algebraic[16] = constants[49]*log(constants[11]/states[2])
    algebraic[26] = constants[25]*(states[0]-algebraic[16])
    algebraic[27] = (((constants[26]*constants[11])/(constants[11]+constants[27]))*(states[0]-algebraic[16]))/(1.00000+exp((((states[0]-algebraic[16])-10.0000)*2.00000)/constants[49]))
    algebraic[18] = constants[13]*(constants[14]+states[6]*(1.00000-constants[14]))*states[5]*(states[0]-algebraic[16])
    algebraic[20] = (((constants[15]*constants[11])/(constants[16]+constants[11]))*states[1])/(constants[17]+states[1])
    algebraic[32] = (1.00000-states[9])*(1.00000-algebraic[31])
    algebraic[34] = (((0.00200000*constants[28]*states[8]*algebraic[32]*(states[0]-50.0000))/constants[49])/(1.00000-exp(-(states[0]-50.0000)/constants[49])))*(states[2]*exp(50.0000/constants[49])-constants[11]*exp(-(states[0]-50.0000)/constants[49]))
    rates[2] = (-1.00000/(1.00000*constants[53]*constants[2]))*((algebraic[27]+algebraic[34]+algebraic[18]+algebraic[26])-2.00000*algebraic[20])
    algebraic[21] = constants[49]*log(constants[10]/states[1])
    algebraic[22] = constants[18]*(states[0]-algebraic[21])
    algebraic[25] = (constants[21]*(exp((constants[24]*(constants[22]-2.00000)*states[0])/constants[49])*(power(states[1], constants[22]))*constants[20]-exp(((constants[24]-1.00000)*(constants[22]-2.00000)*states[0])/constants[49])*(power(constants[10], constants[22]))*states[7]))/((1.00000+constants[23]*(states[7]*(power(constants[10], constants[22]))+constants[20]*(power(states[1], constants[22]))))*(1.00000+states[7]/0.00690000))
    algebraic[10] = constants[49]*log((constants[10]+0.120000*constants[11])/(states[1]+0.120000*states[2]))
    algebraic[14] = constants[9]*(power(states[3], 3.00000))*states[4]*(states[0]-algebraic[10])
    algebraic[36] = (((0.0100000*constants[28]*states[8]*algebraic[32]*(states[0]-50.0000))/constants[49])/(1.00000-exp(-(states[0]-50.0000)/constants[49])))*(states[1]*exp(50.0000/constants[49])-constants[10]*exp(-(states[0]-50.0000)/constants[49]))
    rates[1] = (-1.00000/(1.00000*constants[53]*constants[2]))*(algebraic[14]+algebraic[22]+3.00000*algebraic[20]+3.00000*algebraic[25]+algebraic[36])
    algebraic[23] = 0.500000*constants[49]*log(constants[20]/states[7])
    algebraic[24] = constants[19]*(states[0]-algebraic[23])
    algebraic[33] = (((4.00000*constants[28]*states[8]*algebraic[32]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-(states[0]-50.0000)*2.00000)/constants[49])))*(states[7]*exp(100.000/constants[49])-constants[20]*exp((-(states[0]-50.0000)*2.00000)/constants[49]))
    algebraic[38] = algebraic[33]+algebraic[34]+algebraic[36]
    algebraic[4] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    rates[0] = (-1.00000/constants[3])*(algebraic[4]+algebraic[26]+algebraic[27]+algebraic[18]+algebraic[22]+algebraic[24]+algebraic[20]+algebraic[25]+algebraic[14]+algebraic[38])
    rates[14] = constants[40]*states[7]*(constants[38]-states[14])-constants[41]*states[14]
    rates[15] = constants[42]*states[7]*(constants[39]-states[15])-constants[43]*states[15]
    algebraic[35] = states[7]+states[10]*constants[51]+constants[29]*constants[30]+constants[29]
    algebraic[37] = (states[7]/algebraic[35])*constants[32]-((states[10]*constants[51])/algebraic[35])*constants[33]
    algebraic[40] = (states[10]-states[13])*constants[37]
    rates[10] = (constants[52]/constants[48])*algebraic[37]-algebraic[40]
    algebraic[39] = ((power(states[11]/(states[11]+0.250000), 2.00000))*constants[35]+constants[34])*states[13]
    rates[13] = (constants[48]/constants[46])*algebraic[40]-algebraic[39]
    rates[7] = ((((-1.00000/(2.00000*1.00000*constants[53]*constants[2]))*((algebraic[33]+algebraic[24])-2.00000*algebraic[25])+(algebraic[39]*constants[46])/constants[52])-rates[14])-rates[15])-algebraic[37]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = 20.0000*exp(-0.125000*(states[0]+75.0000))
    algebraic[7] = 2000.00/(1.00000+320.000*exp(-0.100000*(states[0]+75.0000)))
    algebraic[2] = 0.0330000*exp(-states[0]/17.0000)
    algebraic[8] = 33.0000/(1.00000+exp(-0.125000*(states[0]+10.0000)))
    algebraic[0] = states[0]+41.0000
    algebraic[6] = custom_piecewise([less(fabs(algebraic[0]) , constants[12]), 2000.00 , True, (200.000*algebraic[0])/(1.00000-exp(-0.100000*algebraic[0]))])
    algebraic[12] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    algebraic[3] = states[0]+19.0000
    algebraic[9] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (30.0000*algebraic[3])/(1.00000-exp(-algebraic[3]/4.00000))])
    algebraic[13] = custom_piecewise([less(fabs(algebraic[3]) , 0.000100000), 120.000 , True, (12.0000*algebraic[3])/(exp(algebraic[3]/10.0000)-1.00000)])
    algebraic[11] = exp(0.0800000*(states[0]-40.0000))
    algebraic[15] = power(states[7]/(states[7]+0.000500000), 2.00000)
    algebraic[17] = 600.000*algebraic[11]+500.000*algebraic[15]
    algebraic[19] = 60.0000+500.000*algebraic[15]
    algebraic[5] = (1.00000-states[11])-states[12]
    algebraic[28] = states[0]+34.0000
    algebraic[29] = custom_piecewise([less(fabs(algebraic[28]) , 0.000100000), 25.0000 , True, (6.25000*algebraic[28])/(exp(algebraic[28]/4.00000)-1.00000)])
    algebraic[30] = 12.0000/(1.00000+exp(-algebraic[28]/4.00000))
    algebraic[31] = states[7]/(0.00100000+states[7])
    algebraic[16] = constants[49]*log(constants[11]/states[2])
    algebraic[26] = constants[25]*(states[0]-algebraic[16])
    algebraic[27] = (((constants[26]*constants[11])/(constants[11]+constants[27]))*(states[0]-algebraic[16]))/(1.00000+exp((((states[0]-algebraic[16])-10.0000)*2.00000)/constants[49]))
    algebraic[18] = constants[13]*(constants[14]+states[6]*(1.00000-constants[14]))*states[5]*(states[0]-algebraic[16])
    algebraic[20] = (((constants[15]*constants[11])/(constants[16]+constants[11]))*states[1])/(constants[17]+states[1])
    algebraic[32] = (1.00000-states[9])*(1.00000-algebraic[31])
    algebraic[34] = (((0.00200000*constants[28]*states[8]*algebraic[32]*(states[0]-50.0000))/constants[49])/(1.00000-exp(-(states[0]-50.0000)/constants[49])))*(states[2]*exp(50.0000/constants[49])-constants[11]*exp(-(states[0]-50.0000)/constants[49]))
    algebraic[21] = constants[49]*log(constants[10]/states[1])
    algebraic[22] = constants[18]*(states[0]-algebraic[21])
    algebraic[25] = (constants[21]*(exp((constants[24]*(constants[22]-2.00000)*states[0])/constants[49])*(power(states[1], constants[22]))*constants[20]-exp(((constants[24]-1.00000)*(constants[22]-2.00000)*states[0])/constants[49])*(power(constants[10], constants[22]))*states[7]))/((1.00000+constants[23]*(states[7]*(power(constants[10], constants[22]))+constants[20]*(power(states[1], constants[22]))))*(1.00000+states[7]/0.00690000))
    algebraic[10] = constants[49]*log((constants[10]+0.120000*constants[11])/(states[1]+0.120000*states[2]))
    algebraic[14] = constants[9]*(power(states[3], 3.00000))*states[4]*(states[0]-algebraic[10])
    algebraic[36] = (((0.0100000*constants[28]*states[8]*algebraic[32]*(states[0]-50.0000))/constants[49])/(1.00000-exp(-(states[0]-50.0000)/constants[49])))*(states[1]*exp(50.0000/constants[49])-constants[10]*exp(-(states[0]-50.0000)/constants[49]))
    algebraic[23] = 0.500000*constants[49]*log(constants[20]/states[7])
    algebraic[24] = constants[19]*(states[0]-algebraic[23])
    algebraic[33] = (((4.00000*constants[28]*states[8]*algebraic[32]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-(states[0]-50.0000)*2.00000)/constants[49])))*(states[7]*exp(100.000/constants[49])-constants[20]*exp((-(states[0]-50.0000)*2.00000)/constants[49]))
    algebraic[38] = algebraic[33]+algebraic[34]+algebraic[36]
    algebraic[4] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    algebraic[35] = states[7]+states[10]*constants[51]+constants[29]*constants[30]+constants[29]
    algebraic[37] = (states[7]/algebraic[35])*constants[32]-((states[10]*constants[51])/algebraic[35])*constants[33]
    algebraic[40] = (states[10]-states[13])*constants[37]
    algebraic[39] = ((power(states[11]/(states[11]+0.250000), 2.00000))*constants[35]+constants[34])*states[13]
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