# Size of variable arrays:
sizeAlgebraic = 40
sizeStates = 15
sizeConstants = 55
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
    legend_algebraic[24] = "i_K1 in component time_independent_potassium_current (nanoA)"
    legend_algebraic[18] = "i_b_Na in component sodium_background_current (nanoA)"
    legend_algebraic[20] = "i_b_Ca in component calcium_background_current (nanoA)"
    legend_algebraic[23] = "i_b_K in component potassium_background_current (nanoA)"
    legend_algebraic[14] = "i_NaK in component sodium_potassium_pump (nanoA)"
    legend_algebraic[21] = "i_NaCa in component Na_Ca_exchanger (nanoA)"
    legend_algebraic[12] = "i_Na in component fast_sodium_current (nanoA)"
    legend_algebraic[35] = "i_si in component second_inward_calcium_current (nanoA)"
    legend_algebraic[3] = "i_Stim in component membrane (nanoA)"
    legend_constants[4] = "stim_start in component membrane (second)"
    legend_constants[5] = "stim_end in component membrane (second)"
    legend_constants[6] = "stim_period in component membrane (second)"
    legend_constants[7] = "stim_duration in component membrane (second)"
    legend_constants[8] = "stim_amplitude in component membrane (nanoA)"
    legend_constants[9] = "g_Na in component fast_sodium_current (microS)"
    legend_algebraic[8] = "E_mh in component fast_sodium_current (millivolt)"
    legend_constants[10] = "Na_o in component extracellular_sodium_concentration (millimolar)"
    legend_states[1] = "Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_constants[11] = "K_c in component extracellular_potassium_concentration (millimolar)"
    legend_states[2] = "K_i in component intracellular_potassium_concentration (millimolar)"
    legend_states[3] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[4] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_algebraic[5] = "alpha_m in component fast_sodium_current_m_gate (per_second)"
    legend_algebraic[10] = "beta_m in component fast_sodium_current_m_gate (per_second)"
    legend_constants[12] = "delta_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[0] = "E0_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[1] = "alpha_h in component fast_sodium_current_h_gate (per_second)"
    legend_algebraic[6] = "beta_h in component fast_sodium_current_h_gate (per_second)"
    legend_constants[13] = "i_NaK_max in component sodium_potassium_pump (nanoA)"
    legend_constants[14] = "K_mK in component sodium_potassium_pump (millimolar)"
    legend_constants[15] = "K_mNa in component sodium_potassium_pump (millimolar)"
    legend_algebraic[16] = "E_Na in component sodium_background_current (millivolt)"
    legend_constants[16] = "g_b_Na in component sodium_background_current (microS)"
    legend_algebraic[19] = "E_Ca in component calcium_background_current (millivolt)"
    legend_constants[17] = "g_b_Ca in component calcium_background_current (microS)"
    legend_states[5] = "Ca_o in component extracellular_calcium_concentration (millimolar)"
    legend_states[6] = "Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_constants[18] = "k_NaCa in component Na_Ca_exchanger (nanoA)"
    legend_constants[19] = "n_NaCa in component Na_Ca_exchanger (dimensionless)"
    legend_constants[20] = "d_NaCa in component Na_Ca_exchanger (dimensionless)"
    legend_constants[21] = "gamma in component Na_Ca_exchanger (dimensionless)"
    legend_algebraic[22] = "E_K in component time_independent_potassium_current (millivolt)"
    legend_constants[22] = "g_b_K in component potassium_background_current (microS)"
    legend_constants[23] = "g_K1 in component time_independent_potassium_current (microS)"
    legend_constants[24] = "K_m_K1 in component time_independent_potassium_current (millimolar)"
    legend_algebraic[30] = "i_siCa in component second_inward_calcium_current (nanoA)"
    legend_algebraic[31] = "i_siK in component second_inward_calcium_current (nanoA)"
    legend_algebraic[33] = "i_siNa in component second_inward_calcium_current (nanoA)"
    legend_constants[25] = "P_si in component second_inward_calcium_current (nanoA_per_millimolar)"
    legend_states[7] = "d in component second_inward_calcium_current_d_gate (dimensionless)"
    legend_states[8] = "f_Ca in component second_inward_calcium_current_f_Ca_gate (dimensionless)"
    legend_algebraic[29] = "CaChon in component second_inward_calcium_current_f_Ca_gate (dimensionless)"
    legend_algebraic[7] = "alpha_d in component second_inward_calcium_current_d_gate (per_second)"
    legend_algebraic[11] = "beta_d in component second_inward_calcium_current_d_gate (per_second)"
    legend_constants[26] = "delta_d in component second_inward_calcium_current_d_gate (millivolt)"
    legend_algebraic[2] = "E0_d in component second_inward_calcium_current_d_gate (millivolt)"
    legend_algebraic[26] = "alpha_f_Ca in component second_inward_calcium_current_f_Ca_gate (per_second)"
    legend_algebraic[27] = "beta_f_Ca in component second_inward_calcium_current_f_Ca_gate (per_second)"
    legend_algebraic[28] = "CaChoff in component second_inward_calcium_current_f_Ca_gate (dimensionless)"
    legend_constants[27] = "delta_f in component second_inward_calcium_current_f_Ca_gate (millivolt)"
    legend_algebraic[25] = "E0_f in component second_inward_calcium_current_f_Ca_gate (millivolt)"
    legend_algebraic[34] = "i_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[51] = "K_1 in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_algebraic[32] = "K_2 in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[28] = "K_cyca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[29] = "K_xcs in component sarcoplasmic_reticulum_calcium_pump (dimensionless)"
    legend_constants[30] = "K_srca in component sarcoplasmic_reticulum_calcium_pump (millimolar)"
    legend_constants[31] = "alpha_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_constants[32] = "beta_up in component sarcoplasmic_reticulum_calcium_pump (millimolar_per_second)"
    legend_states[9] = "Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[36] = "i_rel in component calcium_release (millimolar_per_second)"
    legend_algebraic[9] = "VoltDep in component calcium_release (dimensionless)"
    legend_algebraic[13] = "RegBindSite in component calcium_release (dimensionless)"
    legend_algebraic[15] = "ActRate in component calcium_release (per_second)"
    legend_algebraic[17] = "InactRate in component calcium_release (per_second)"
    legend_constants[33] = "K_leak_rate in component calcium_release (per_second)"
    legend_constants[34] = "K_m_rel in component calcium_release (per_second)"
    legend_algebraic[4] = "PrecFrac in component calcium_release (dimensionless)"
    legend_states[10] = "ActFrac in component calcium_release (dimensionless)"
    legend_states[11] = "ProdFrac in component calcium_release (dimensionless)"
    legend_states[12] = "Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_algebraic[37] = "i_trans in component calcium_translocation (millimolar_per_second)"
    legend_constants[35] = "alpha_tr in component calcium_translocation (per_second)"
    legend_constants[54] = "V_i in component intracellular_calcium_concentration (micrometre3)"
    legend_constants[36] = "Cab in component extracellular_calcium_concentration (millimolar)"
    legend_constants[37] = "K_diff in component extracellular_calcium_concentration (per_second)"
    legend_constants[53] = "Ve in component intracellular_calcium_concentration (micrometre3)"
    legend_states[13] = "Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_states[14] = "Ca_Trop in component intracellular_calcium_concentration (millimolar)"
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
    legend_algebraic[38] = "dCaCalmoddt in component intracellular_calcium_concentration (millimolar_per_second)"
    legend_algebraic[39] = "dCaTropdt in component intracellular_calcium_concentration (millimolar_per_second)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[3] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[4] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[7] = "d/dt d in component second_inward_calcium_current_d_gate (dimensionless)"
    legend_rates[8] = "d/dt f_Ca in component second_inward_calcium_current_f_Ca_gate (dimensionless)"
    legend_rates[10] = "d/dt ActFrac in component calcium_release (dimensionless)"
    legend_rates[11] = "d/dt ProdFrac in component calcium_release (dimensionless)"
    legend_rates[1] = "d/dt Na_i in component intracellular_sodium_concentration (millimolar)"
    legend_rates[5] = "d/dt Ca_o in component extracellular_calcium_concentration (millimolar)"
    legend_rates[2] = "d/dt K_i in component intracellular_potassium_concentration (millimolar)"
    legend_rates[6] = "d/dt Ca_i in component intracellular_calcium_concentration (millimolar)"
    legend_rates[9] = "d/dt Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_rates[12] = "d/dt Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_rates[13] = "d/dt Ca_Calmod in component intracellular_calcium_concentration (millimolar)"
    legend_rates[14] = "d/dt Ca_Trop in component intracellular_calcium_concentration (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -88
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 0.006
    constants[4] = 0.1
    constants[5] = 10000
    constants[6] = 1
    constants[7] = 0.002
    constants[8] = -200
    constants[9] = 50
    constants[10] = 140
    states[1] = 6.5
    constants[11] = 4
    states[2] = 140
    states[3] = 0.076
    states[4] = 0.015
    constants[12] = 1e-5
    constants[13] = 14
    constants[14] = 1
    constants[15] = 40
    constants[16] = 0.012
    constants[17] = 0.005
    states[5] = 2
    states[6] = 1e-5
    constants[18] = 0.01
    constants[19] = 3
    constants[20] = 0.0001
    constants[21] = 0.5
    constants[22] = 0.17
    constants[23] = 1.7
    constants[24] = 10
    constants[25] = 5
    states[7] = 0.0011
    states[8] = 0.785
    constants[26] = 0.0001
    constants[27] = 0.0001
    constants[28] = 0.0003
    constants[29] = 0.4
    constants[30] = 0.5
    constants[31] = 3
    constants[32] = 0.23
    states[9] = 0.3
    constants[33] = 0
    constants[34] = 250
    states[10] = 0
    states[11] = 0
    states[12] = 0.3
    constants[35] = 50
    constants[36] = 2
    constants[37] = 0.0005
    states[13] = 0.0005
    states[14] = 0.0015
    constants[38] = 0.02
    constants[39] = 0.15
    constants[40] = 100000
    constants[41] = 50
    constants[42] = 100000
    constants[43] = 200
    constants[44] = 0.08
    constants[45] = 0.08
    constants[46] = 0.1
    constants[47] = 0.4
    constants[48] = 0.01
    constants[49] = (constants[0]*constants[1])/constants[2]
    constants[50] = 3.14159*(power(constants[44], 2.00000))*constants[45]
    constants[51] = (constants[28]*constants[29])/constants[30]
    constants[52] = ((1.00000-constants[47])-constants[48])-constants[46]
    constants[53] = constants[50]*constants[47]
    constants[54] = constants[50]*constants[52]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[13] = constants[40]*states[6]*(constants[38]-states[13])-constants[41]*states[13]
    rates[14] = constants[42]*states[6]*(constants[39]-states[14])-constants[43]*states[14]
    algebraic[1] = 20.0000*exp(-0.125000*(states[0]+75.0000))
    algebraic[6] = 2000.00/(1.00000+320.000*exp(-0.100000*(states[0]+75.0000)))
    rates[4] = algebraic[1]*(1.00000-states[4])-algebraic[6]*states[4]
    algebraic[0] = states[0]+41.0000
    algebraic[5] = custom_piecewise([less(fabs(algebraic[0]) , constants[12]), 2000.00 , True, (200.000*algebraic[0])/(1.00000-exp(-0.100000*algebraic[0]))])
    algebraic[10] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    rates[3] = algebraic[5]*(1.00000-states[3])-algebraic[10]*states[3]
    algebraic[2] = (states[0]+24.0000)-5.00000
    algebraic[7] = custom_piecewise([less(fabs(algebraic[2]) , constants[26]), 120.000 , True, (30.0000*algebraic[2])/(1.00000-exp((-1.00000*algebraic[2])/4.00000))])
    algebraic[11] = custom_piecewise([less(fabs(algebraic[2]) , constants[26]), 120.000 , True, (12.0000*algebraic[2])/(exp(algebraic[2]/10.0000)-1.00000)])
    rates[7] = algebraic[7]*(1.00000-states[7])-algebraic[11]*states[7]
    algebraic[9] = exp(0.0800000*(states[0]-40.0000))
    algebraic[13] = power(states[6]/(states[6]+0.000500000), 2.00000)
    algebraic[15] = 600.000*algebraic[9]+500.000*algebraic[13]
    algebraic[17] = 60.0000+500.000*algebraic[13]
    algebraic[4] = (1.00000-states[10])-states[11]
    rates[10] = algebraic[4]*algebraic[15]-states[10]*algebraic[17]
    rates[11] = states[10]*algebraic[17]-0.600000*states[11]
    algebraic[25] = states[0]+34.0000
    algebraic[26] = custom_piecewise([less(fabs(algebraic[25]) , constants[27]), 25.0000 , True, (6.25000*algebraic[25])/(exp(algebraic[25]/4.00000)-1.00000)])
    algebraic[27] = 12.0000/(1.00000+exp(-algebraic[25]/4.00000))
    algebraic[28] = states[6]/(0.00100000+states[6])
    rates[8] = (120.000*(1.00000-states[8])*algebraic[28]+(1.00000-states[8])*(1.00000-algebraic[28]))*algebraic[27]-algebraic[26]*states[8]
    algebraic[19] = 0.500000*constants[49]*log(states[5]/states[6])
    algebraic[20] = constants[17]*(states[0]-algebraic[19])
    algebraic[21] = (constants[18]*(exp((constants[21]*(constants[19]-2.00000)*states[0])/constants[49])*(power(states[1], constants[19]))*states[5]-exp(((constants[21]-1.00000)*(constants[19]-2.00000)*states[0])/constants[49])*(power(constants[10], constants[19]))*states[6]))/((1.00000+constants[20]*(states[6]*(power(constants[10], constants[19]))+states[5]*(power(states[1], constants[19]))))*(1.00000+states[6]/0.00690000))
    algebraic[29] = (1.00000-states[8])*(1.00000-algebraic[28])
    algebraic[30] = (((4.00000*constants[25]*states[7]*algebraic[29]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-1.00000*(states[0]-50.0000)*2.00000)/constants[49])))*(states[6]*exp(100.000/constants[49])-states[5]*exp((-2.00000*(states[0]-50.0000))/constants[49]))
    rates[5] = (constants[36]-states[5])*constants[37]-(1.00000*(algebraic[30]+algebraic[21]+algebraic[20]))/(2.00000*1.00000*constants[53]*constants[2])
    algebraic[22] = constants[49]*log(constants[11]/states[2])
    algebraic[24] = (((constants[23]*constants[11])/(constants[11]+constants[24]))*(states[0]-algebraic[22]))/(1.00000+exp((((states[0]-algebraic[22])-10.0000)*2.00000)/constants[49]))
    algebraic[23] = constants[22]*(states[0]-algebraic[22])
    algebraic[14] = (((constants[13]*constants[11])/(constants[14]+constants[11]))*states[1])/(constants[15]+states[1])
    algebraic[31] = (((0.00200000*constants[25]*states[7]*algebraic[29]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[49])))*(states[2]*exp(50.0000/constants[49])-constants[11]*exp((-1.00000*(states[0]-50.0000))/constants[49]))
    rates[2] = (-1.00000/(1.00000*constants[54]*constants[2]))*((algebraic[24]+algebraic[31]+algebraic[23])-2.00000*algebraic[14])
    algebraic[16] = constants[49]*log(constants[10]/states[1])
    algebraic[18] = constants[16]*(states[0]-algebraic[16])
    algebraic[8] = constants[49]*log((constants[10]+0.120000*constants[11])/(states[1]+0.120000*states[2]))
    algebraic[12] = constants[9]*(power(states[3], 3.00000))*states[4]*(states[0]-algebraic[8])
    algebraic[33] = (((0.0100000*constants[25]*states[7]*algebraic[29]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[49])))*(states[1]*exp(50.0000/constants[49])-constants[10]*exp((-1.00000*(states[0]-50.0000))/constants[49]))
    rates[1] = (-1.00000/(1.00000*constants[54]*constants[2]))*(algebraic[12]+algebraic[18]+algebraic[14]*3.00000+algebraic[21]*3.00000+algebraic[33])
    algebraic[35] = algebraic[30]+algebraic[31]+algebraic[33]
    algebraic[3] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    rates[0] = -(algebraic[3]+algebraic[24]+algebraic[18]+algebraic[20]+algebraic[23]+algebraic[14]+algebraic[21]+algebraic[12]+algebraic[35])/constants[3]
    algebraic[32] = states[6]+states[9]*constants[51]+constants[28]*constants[29]+constants[28]
    algebraic[34] = (states[6]/algebraic[32])*constants[31]-((states[9]*constants[51])/algebraic[32])*constants[32]
    algebraic[37] = (states[9]-states[12])*constants[35]
    rates[9] = (constants[52]/constants[48])*algebraic[34]-algebraic[37]
    algebraic[36] = ((power(states[10]/(states[10]+0.250000), 2.00000))*constants[34]+constants[33])*states[12]
    rates[12] = (constants[48]/constants[46])*algebraic[37]-algebraic[36]
    algebraic[38] = constants[40]*states[6]*(constants[38]-states[13])-constants[41]*states[13]
    algebraic[39] = constants[42]*states[6]*(constants[39]-states[14])-constants[43]*states[14]
    rates[6] = ((((-1.00000/(2.00000*1.00000*constants[54]*constants[2]))*((algebraic[30]+algebraic[20])-2.00000*algebraic[21])+(algebraic[36]*constants[46])/constants[52])-algebraic[38])-algebraic[39])-algebraic[34]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = 20.0000*exp(-0.125000*(states[0]+75.0000))
    algebraic[6] = 2000.00/(1.00000+320.000*exp(-0.100000*(states[0]+75.0000)))
    algebraic[0] = states[0]+41.0000
    algebraic[5] = custom_piecewise([less(fabs(algebraic[0]) , constants[12]), 2000.00 , True, (200.000*algebraic[0])/(1.00000-exp(-0.100000*algebraic[0]))])
    algebraic[10] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    algebraic[2] = (states[0]+24.0000)-5.00000
    algebraic[7] = custom_piecewise([less(fabs(algebraic[2]) , constants[26]), 120.000 , True, (30.0000*algebraic[2])/(1.00000-exp((-1.00000*algebraic[2])/4.00000))])
    algebraic[11] = custom_piecewise([less(fabs(algebraic[2]) , constants[26]), 120.000 , True, (12.0000*algebraic[2])/(exp(algebraic[2]/10.0000)-1.00000)])
    algebraic[9] = exp(0.0800000*(states[0]-40.0000))
    algebraic[13] = power(states[6]/(states[6]+0.000500000), 2.00000)
    algebraic[15] = 600.000*algebraic[9]+500.000*algebraic[13]
    algebraic[17] = 60.0000+500.000*algebraic[13]
    algebraic[4] = (1.00000-states[10])-states[11]
    algebraic[25] = states[0]+34.0000
    algebraic[26] = custom_piecewise([less(fabs(algebraic[25]) , constants[27]), 25.0000 , True, (6.25000*algebraic[25])/(exp(algebraic[25]/4.00000)-1.00000)])
    algebraic[27] = 12.0000/(1.00000+exp(-algebraic[25]/4.00000))
    algebraic[28] = states[6]/(0.00100000+states[6])
    algebraic[19] = 0.500000*constants[49]*log(states[5]/states[6])
    algebraic[20] = constants[17]*(states[0]-algebraic[19])
    algebraic[21] = (constants[18]*(exp((constants[21]*(constants[19]-2.00000)*states[0])/constants[49])*(power(states[1], constants[19]))*states[5]-exp(((constants[21]-1.00000)*(constants[19]-2.00000)*states[0])/constants[49])*(power(constants[10], constants[19]))*states[6]))/((1.00000+constants[20]*(states[6]*(power(constants[10], constants[19]))+states[5]*(power(states[1], constants[19]))))*(1.00000+states[6]/0.00690000))
    algebraic[29] = (1.00000-states[8])*(1.00000-algebraic[28])
    algebraic[30] = (((4.00000*constants[25]*states[7]*algebraic[29]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-1.00000*(states[0]-50.0000)*2.00000)/constants[49])))*(states[6]*exp(100.000/constants[49])-states[5]*exp((-2.00000*(states[0]-50.0000))/constants[49]))
    algebraic[22] = constants[49]*log(constants[11]/states[2])
    algebraic[24] = (((constants[23]*constants[11])/(constants[11]+constants[24]))*(states[0]-algebraic[22]))/(1.00000+exp((((states[0]-algebraic[22])-10.0000)*2.00000)/constants[49]))
    algebraic[23] = constants[22]*(states[0]-algebraic[22])
    algebraic[14] = (((constants[13]*constants[11])/(constants[14]+constants[11]))*states[1])/(constants[15]+states[1])
    algebraic[31] = (((0.00200000*constants[25]*states[7]*algebraic[29]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[49])))*(states[2]*exp(50.0000/constants[49])-constants[11]*exp((-1.00000*(states[0]-50.0000))/constants[49]))
    algebraic[16] = constants[49]*log(constants[10]/states[1])
    algebraic[18] = constants[16]*(states[0]-algebraic[16])
    algebraic[8] = constants[49]*log((constants[10]+0.120000*constants[11])/(states[1]+0.120000*states[2]))
    algebraic[12] = constants[9]*(power(states[3], 3.00000))*states[4]*(states[0]-algebraic[8])
    algebraic[33] = (((0.0100000*constants[25]*states[7]*algebraic[29]*(states[0]-50.0000))/constants[49])/(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[49])))*(states[1]*exp(50.0000/constants[49])-constants[10]*exp((-1.00000*(states[0]-50.0000))/constants[49]))
    algebraic[35] = algebraic[30]+algebraic[31]+algebraic[33]
    algebraic[3] = custom_piecewise([greater_equal(voi , constants[4]) & less_equal(voi , constants[5]) & less_equal((voi-constants[4])-floor((voi-constants[4])/constants[6])*constants[6] , constants[7]), constants[8] , True, 0.00000])
    algebraic[32] = states[6]+states[9]*constants[51]+constants[28]*constants[29]+constants[28]
    algebraic[34] = (states[6]/algebraic[32])*constants[31]-((states[9]*constants[51])/algebraic[32])*constants[32]
    algebraic[37] = (states[9]-states[12])*constants[35]
    algebraic[36] = ((power(states[10]/(states[10]+0.250000), 2.00000))*constants[34]+constants[33])*states[12]
    algebraic[38] = constants[40]*states[6]*(constants[38]-states[13])-constants[41]*states[13]
    algebraic[39] = constants[42]*states[6]*(constants[39]-states[14])-constants[43]*states[14]
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