# Size of variable arrays:
sizeAlgebraic = 44
sizeStates = 15
sizeConstants = 48
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
    legend_constants[3] = "C in component membrane (microF)"
    legend_constants[43] = "RTONF in component membrane (millivolt)"
    legend_algebraic[25] = "i_f in component hyperpolarising_activated_current (nanoA)"
    legend_algebraic[27] = "i_K in component time_dependent_potassium_current (nanoA)"
    legend_algebraic[28] = "i_K1 in component time_independent_potassium_current (nanoA)"
    legend_algebraic[29] = "i_Na_b in component sodium_background_current (nanoA)"
    legend_algebraic[31] = "i_Ca_b in component calcium_background_current (nanoA)"
    legend_algebraic[32] = "i_p in component sodium_potassium_pump (nanoA)"
    legend_algebraic[33] = "i_NaCa in component Na_Ca_exchanger (nanoA)"
    legend_algebraic[35] = "i_Na in component fast_sodium_current (nanoA)"
    legend_algebraic[42] = "i_si in component second_inward_current (nanoA)"
    legend_algebraic[22] = "i_fNa in component hyperpolarising_activated_current (nanoA)"
    legend_algebraic[0] = "E_Na in component hyperpolarising_activated_current (millivolt)"
    legend_algebraic[9] = "E_K in component hyperpolarising_activated_current (millivolt)"
    legend_algebraic[24] = "i_fK in component hyperpolarising_activated_current (nanoA)"
    legend_constants[4] = "g_f_Na in component hyperpolarising_activated_current (microS)"
    legend_constants[5] = "g_f_K in component hyperpolarising_activated_current (microS)"
    legend_constants[6] = "Km_f in component hyperpolarising_activated_current (millimolar)"
    legend_states[1] = "Kc in component extracellular_potassium_concentration (millimolar)"
    legend_states[2] = "Ki in component intracellular_potassium_concentration (millimolar)"
    legend_states[3] = "Nai in component intracellular_sodium_concentration (millimolar)"
    legend_constants[7] = "Nao in component extracellular_sodium_concentration (millimolar)"
    legend_states[4] = "y in component hyperpolarising_activated_current_y_gate (dimensionless)"
    legend_algebraic[10] = "alpha_y in component hyperpolarising_activated_current_y_gate (per_second)"
    legend_algebraic[17] = "beta_y in component hyperpolarising_activated_current_y_gate (per_second)"
    legend_constants[8] = "delta_y in component hyperpolarising_activated_current_y_gate (millivolt)"
    legend_algebraic[1] = "E0_y in component hyperpolarising_activated_current_y_gate (millivolt)"
    legend_constants[9] = "speed_y in component hyperpolarising_activated_current_y_gate (dimensionless)"
    legend_algebraic[26] = "I_K in component time_dependent_potassium_current (nanoA)"
    legend_constants[10] = "i_K_max in component time_dependent_potassium_current (nanoA)"
    legend_states[5] = "x in component time_dependent_potassium_current_x_gate (dimensionless)"
    legend_algebraic[11] = "alpha_x in component time_dependent_potassium_current_x_gate (per_second)"
    legend_algebraic[18] = "beta_x in component time_dependent_potassium_current_x_gate (per_second)"
    legend_constants[11] = "delta_x in component time_dependent_potassium_current_x_gate (millivolt)"
    legend_algebraic[2] = "E0_x in component time_dependent_potassium_current_x_gate (millivolt)"
    legend_constants[12] = "g_K1 in component time_independent_potassium_current (microS)"
    legend_constants[13] = "Km_K1 in component time_independent_potassium_current (millimolar)"
    legend_constants[14] = "g_Nab in component sodium_background_current (microS)"
    legend_algebraic[30] = "E_Ca in component calcium_background_current (millivolt)"
    legend_constants[15] = "g_Cab in component calcium_background_current (microS)"
    legend_states[6] = "Cai in component intracellular_calcium_concentration (millimolar)"
    legend_constants[16] = "Cao in component extracellular_calcium_concentration (millimolar)"
    legend_constants[17] = "I_p in component sodium_potassium_pump (nanoA)"
    legend_constants[18] = "K_mK in component sodium_potassium_pump (millimolar)"
    legend_constants[19] = "K_mNa in component sodium_potassium_pump (millimolar)"
    legend_constants[20] = "n_NaCa in component Na_Ca_exchanger (dimensionless)"
    legend_constants[21] = "K_NaCa in component Na_Ca_exchanger (nanoA)"
    legend_constants[22] = "d_NaCa in component Na_Ca_exchanger (dimensionless)"
    legend_constants[23] = "gamma in component Na_Ca_exchanger (dimensionless)"
    legend_constants[24] = "g_Na in component fast_sodium_current (microS)"
    legend_algebraic[34] = "E_mh in component fast_sodium_current (millivolt)"
    legend_states[7] = "m in component fast_sodium_current_m_gate (dimensionless)"
    legend_states[8] = "h in component fast_sodium_current_h_gate (dimensionless)"
    legend_algebraic[12] = "alpha_m in component fast_sodium_current_m_gate (per_second)"
    legend_algebraic[19] = "beta_m in component fast_sodium_current_m_gate (per_second)"
    legend_constants[25] = "delta_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[3] = "E0_m in component fast_sodium_current_m_gate (millivolt)"
    legend_algebraic[4] = "alpha_h in component fast_sodium_current_h_gate (per_second)"
    legend_algebraic[13] = "beta_h in component fast_sodium_current_h_gate (per_second)"
    legend_algebraic[36] = "i_siCa in component second_inward_current (nanoA)"
    legend_algebraic[37] = "i_siK in component second_inward_current (nanoA)"
    legend_algebraic[39] = "i_siNa in component second_inward_current (nanoA)"
    legend_constants[26] = "P_si in component second_inward_current (nanoA_per_millimolar)"
    legend_states[9] = "d in component second_inward_current_d_gate (dimensionless)"
    legend_states[10] = "f in component second_inward_current_f_gate (dimensionless)"
    legend_states[11] = "f2 in component second_inward_current_f2_gate (dimensionless)"
    legend_algebraic[14] = "alpha_d in component second_inward_current_d_gate (per_second)"
    legend_algebraic[20] = "beta_d in component second_inward_current_d_gate (per_second)"
    legend_constants[27] = "delta_d in component second_inward_current_d_gate (millivolt)"
    legend_algebraic[5] = "E0_d in component second_inward_current_d_gate (millivolt)"
    legend_algebraic[15] = "alpha_f in component second_inward_current_f_gate (per_second)"
    legend_algebraic[21] = "beta_f in component second_inward_current_f_gate (per_second)"
    legend_constants[28] = "delta_f in component second_inward_current_f_gate (millivolt)"
    legend_algebraic[6] = "E0_f in component second_inward_current_f_gate (millivolt)"
    legend_constants[29] = "alpha_f2 in component second_inward_current_f2_gate (per_second)"
    legend_algebraic[7] = "beta_f2 in component second_inward_current_f2_gate (per_second)"
    legend_constants[30] = "K_mf2 in component second_inward_current_f2_gate (millimolar)"
    legend_constants[31] = "radius in component intracellular_sodium_concentration (micrometre)"
    legend_constants[32] = "length in component intracellular_sodium_concentration (micrometre)"
    legend_constants[33] = "V_e_ratio in component intracellular_sodium_concentration (dimensionless)"
    legend_constants[44] = "V_Cell in component intracellular_sodium_concentration (micrometre3)"
    legend_constants[45] = "Vi in component intracellular_sodium_concentration (micrometre3)"
    legend_constants[46] = "V_up in component intracellular_calcium_concentration (micrometre3)"
    legend_constants[47] = "V_rel in component intracellular_calcium_concentration (micrometre3)"
    legend_algebraic[38] = "i_up in component intracellular_calcium_concentration (nanoA)"
    legend_algebraic[40] = "i_tr in component intracellular_calcium_concentration (nanoA)"
    legend_algebraic[43] = "i_rel in component intracellular_calcium_concentration (nanoA)"
    legend_states[12] = "Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_states[13] = "Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_constants[34] = "Ca_up_max in component intracellular_calcium_concentration (millimolar)"
    legend_constants[35] = "K_mCa in component intracellular_calcium_concentration (millimolar)"
    legend_states[14] = "p in component intracellular_calcium_concentration (dimensionless)"
    legend_algebraic[16] = "alpha_p in component intracellular_calcium_concentration (per_second)"
    legend_algebraic[23] = "beta_p in component intracellular_calcium_concentration (per_second)"
    legend_algebraic[8] = "E0_p in component intracellular_calcium_concentration (millivolt)"
    legend_constants[36] = "tau_up in component intracellular_calcium_concentration (second)"
    legend_constants[37] = "tau_rep in component intracellular_calcium_concentration (second)"
    legend_constants[38] = "tau_rel in component intracellular_calcium_concentration (second)"
    legend_constants[39] = "rCa in component intracellular_calcium_concentration (dimensionless)"
    legend_constants[40] = "V_e in component extracellular_potassium_concentration (micrometre3)"
    legend_constants[41] = "Kb in component extracellular_potassium_concentration (millimolar)"
    legend_algebraic[41] = "i_mK in component extracellular_potassium_concentration (nanoA)"
    legend_constants[42] = "pf in component extracellular_potassium_concentration (per_second)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[4] = "d/dt y in component hyperpolarising_activated_current_y_gate (dimensionless)"
    legend_rates[5] = "d/dt x in component time_dependent_potassium_current_x_gate (dimensionless)"
    legend_rates[7] = "d/dt m in component fast_sodium_current_m_gate (dimensionless)"
    legend_rates[8] = "d/dt h in component fast_sodium_current_h_gate (dimensionless)"
    legend_rates[9] = "d/dt d in component second_inward_current_d_gate (dimensionless)"
    legend_rates[10] = "d/dt f in component second_inward_current_f_gate (dimensionless)"
    legend_rates[11] = "d/dt f2 in component second_inward_current_f2_gate (dimensionless)"
    legend_rates[3] = "d/dt Nai in component intracellular_sodium_concentration (millimolar)"
    legend_rates[14] = "d/dt p in component intracellular_calcium_concentration (dimensionless)"
    legend_rates[12] = "d/dt Ca_up in component intracellular_calcium_concentration (millimolar)"
    legend_rates[13] = "d/dt Ca_rel in component intracellular_calcium_concentration (millimolar)"
    legend_rates[6] = "d/dt Cai in component intracellular_calcium_concentration (millimolar)"
    legend_rates[1] = "d/dt Kc in component extracellular_potassium_concentration (millimolar)"
    legend_rates[2] = "d/dt Ki in component intracellular_potassium_concentration (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -60
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 0.006
    constants[4] = 6
    constants[5] = 6
    constants[6] = 45
    states[1] = 3
    states[2] = 140
    states[3] = 7.5
    constants[7] = 140
    states[4] = 0.007
    constants[8] = 1e-5
    constants[9] = 2
    constants[10] = 20
    states[5] = 0.54
    constants[11] = 0.0001
    constants[12] = 0.75
    constants[13] = 10
    constants[14] = 0.07
    constants[15] = 0.01
    states[6] = 5.8e-5
    constants[16] = 2
    constants[17] = 50
    constants[18] = 1
    constants[19] = 40
    constants[20] = 3
    constants[21] = 0.002
    constants[22] = 0.0001
    constants[23] = 0.5
    constants[24] = 1.25
    states[7] = 0.076
    states[8] = 0.015
    constants[25] = 1e-5
    constants[26] = 7.5
    states[9] = 0.0011
    states[10] = 0.785
    states[11] = 0.785
    constants[27] = 0.0001
    constants[28] = 0.0001
    constants[29] = 10
    constants[30] = 0.0005
    constants[31] = 0.08
    constants[32] = 0.08
    constants[33] = 0.1
    states[12] = 1.98
    states[13] = 0.55
    constants[34] = 5
    constants[35] = 0.002
    states[14] = 0.785
    constants[36] = 0.005
    constants[37] = 0.2
    constants[38] = 0.01
    constants[39] = 2
    constants[40] = 0.00016077
    constants[41] = 3
    constants[42] = 1
    constants[43] = (constants[0]*constants[1])/constants[2]
    constants[44] = 3.14159*(power(constants[31], 2.00000))*constants[32]
    constants[45] = constants[44]*(1.00000-constants[33])
    constants[46] = constants[45]*0.0500000
    constants[47] = constants[45]*0.0200000
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[7] = (states[6]*constants[29])/constants[30]
    rates[11] = constants[29]-states[11]*(constants[29]+algebraic[7])
    algebraic[4] = 20.0000*exp(-0.125000*(states[0]+75.0000))
    algebraic[13] = 2000.00/(320.000*exp(-0.100000*(states[0]+75.0000))+1.00000)
    rates[8] = algebraic[4]*(1.00000-states[8])-algebraic[13]*states[8]
    algebraic[1] = states[0]+52.0000
    algebraic[10] = 0.0500000*exp(-0.0670000*algebraic[1])
    algebraic[17] = custom_piecewise([less(fabs(algebraic[1]) , constants[8]), 2.50000 , True, algebraic[1]/(1.00000-1.00000*exp(-0.200000*algebraic[1]))])
    rates[4] = constants[9]*(algebraic[10]*(1.00000-states[4])-algebraic[17]*states[4])
    algebraic[2] = states[0]+22.0000
    algebraic[11] = custom_piecewise([less(fabs(algebraic[2]) , constants[11]), 2.50000 , True, (0.500000*algebraic[2])/(1.00000-exp(-algebraic[2]/5.00000))])
    algebraic[18] = custom_piecewise([less(fabs(algebraic[2]) , constants[11]), 2.50000 , True, (0.178000*algebraic[2])/(exp(algebraic[2]/15.0000)-1.00000)])
    rates[5] = algebraic[11]*(1.00000-states[5])-algebraic[18]*states[5]
    algebraic[3] = states[0]+41.0000
    algebraic[12] = custom_piecewise([less(fabs(algebraic[3]) , constants[25]), 2000.00 , True, (200.000*algebraic[3])/(1.00000-exp(-0.100000*algebraic[3]))])
    algebraic[19] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    rates[7] = algebraic[12]*(1.00000-states[7])-algebraic[19]*states[7]
    algebraic[5] = (states[0]+24.0000)-5.00000
    algebraic[14] = custom_piecewise([less(fabs(algebraic[5]) , constants[27]), 120.000 , True, (30.0000*algebraic[5])/(1.00000-exp((-1.00000*algebraic[5])/4.00000))])
    algebraic[20] = custom_piecewise([less(fabs(algebraic[5]) , constants[27]), 120.000 , True, (12.0000*algebraic[5])/(exp(algebraic[5]/10.0000)-1.00000)])
    rates[9] = algebraic[14]*(1.00000-states[9])-algebraic[20]*states[9]
    algebraic[6] = states[0]+34.0000
    algebraic[15] = custom_piecewise([less(fabs(algebraic[6]) , constants[28]), 25.0000 , True, (6.25000*algebraic[6])/(exp(algebraic[6]/4.00000)-1.00000)])
    algebraic[21] = 50.0000/(1.00000+exp((-1.00000*(states[0]+34.0000))/4.00000))
    rates[10] = algebraic[15]*(1.00000-states[10])-algebraic[21]*states[10]
    algebraic[8] = (states[0]+34.0000)--30.0000
    algebraic[16] = (0.625000*algebraic[8])/(exp(algebraic[8]/4.00000)-1.00000)
    algebraic[23] = 5.00000/(1.00000+exp((-1.00000*algebraic[8])/4.00000))
    rates[14] = algebraic[16]*(1.00000-states[14])-algebraic[23]*states[14]
    algebraic[0] = constants[43]*log(constants[7]/states[3])
    algebraic[29] = constants[14]*(states[0]-algebraic[0])
    algebraic[32] = (((constants[17]*states[1])/(constants[18]+states[1]))*states[3])/(constants[19]+states[3])
    algebraic[33] = (constants[21]*(exp((constants[23]*(constants[20]-2.00000)*states[0])/constants[43])*(power(states[3], constants[20]))*constants[16]-exp(((constants[23]-1.00000)*(constants[20]-2.00000)*states[0])/constants[43])*(power(constants[7], constants[20]))*states[6]))/((1.00000+constants[22]*(states[6]*(power(constants[7], constants[20]))+constants[16]*(power(states[3], constants[20]))))*(1.00000+states[6]/0.00690000))
    algebraic[34] = constants[43]*log((constants[7]+0.120000*states[1])/(states[3]+0.120000*states[2]))
    algebraic[35] = constants[24]*(power(states[7], 3.00000))*states[8]*(states[0]-algebraic[34])
    algebraic[22] = ((states[4]*states[1])/(states[1]+constants[6]))*constants[4]*(states[0]-algebraic[0])
    algebraic[39] = ((0.0100000*constants[26]*(states[0]-50.0000))/(constants[43]*(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[43]))))*(states[3]*exp(50.0000/constants[43])-constants[7]*exp((-1.00000*(states[0]-50.0000))/constants[43]))*states[9]*states[10]*states[11]
    rates[3] = (-1.00000*(algebraic[35]+algebraic[29]+algebraic[22]+algebraic[39]+algebraic[32]*3.00000+(algebraic[33]*constants[20])/(constants[20]-2.00000)))/(1.00000*constants[45]*constants[2])
    algebraic[38] = ((2.00000*1.00000*constants[45]*constants[2])/(1.00000*constants[36]*constants[34]))*states[6]*(constants[34]-states[12])
    algebraic[40] = ((2.00000*1.00000*constants[47]*constants[2])/(1.00000*constants[37]))*states[14]*(states[12]-states[13])
    rates[12] = (1.00000*(algebraic[38]-algebraic[40]))/(2.00000*1.00000*constants[46]*constants[2])
    algebraic[26] = (constants[10]*(states[2]-states[1]*exp(-states[0]/constants[43])))/140.000
    algebraic[27] = states[5]*algebraic[26]
    algebraic[9] = constants[43]*log(states[1]/states[2])
    algebraic[28] = (((constants[12]*states[1])/(states[1]+constants[13]))*(states[0]-algebraic[9]))/(1.00000+exp((((states[0]+10.0000)-algebraic[9])*2.00000)/constants[43]))
    algebraic[24] = ((states[4]*states[1])/(states[1]+constants[6]))*constants[5]*(states[0]-algebraic[9])
    algebraic[37] = ((0.0100000*constants[26]*(states[0]-50.0000))/(constants[43]*(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[43]))))*(states[2]*exp(50.0000/constants[43])-states[1]*exp((-1.00000*(states[0]-50.0000))/constants[43]))*states[9]*states[10]*states[11]
    algebraic[41] = (algebraic[28]+algebraic[27]+algebraic[24]+algebraic[37])-2.00000*algebraic[32]
    rates[1] = -constants[42]*(states[1]-constants[41])+(1.00000*algebraic[41])/(1.00000*constants[40]*constants[2])
    rates[2] = (-1.00000*algebraic[41])/(1.00000*constants[45]*constants[2])
    algebraic[25] = algebraic[22]+algebraic[24]
    algebraic[30] = 0.500000*constants[43]*log(constants[16]/states[6])
    algebraic[31] = constants[15]*(states[0]-algebraic[30])
    algebraic[36] = ((4.00000*constants[26]*(states[0]-50.0000))/(constants[43]*(1.00000-exp((-1.00000*(states[0]-50.0000)*2.00000)/constants[43]))))*(states[6]*exp(100.000/constants[43])-constants[16]*exp((-2.00000*(states[0]-50.0000))/constants[43]))*states[9]*states[10]*states[11]
    algebraic[42] = algebraic[36]+algebraic[37]+algebraic[39]
    rates[0] = -(algebraic[25]+algebraic[27]+algebraic[28]+algebraic[29]+algebraic[31]+algebraic[32]+algebraic[33]+algebraic[35]+algebraic[42])/constants[3]
    algebraic[43] = (((2.00000*1.00000*constants[47]*constants[2])/(1.00000*constants[38]))*states[13]*(power(states[6], constants[39])))/(power(states[6], constants[39])+power(constants[35], constants[39]))
    rates[13] = (1.00000*(algebraic[40]-algebraic[43]))/(2.00000*1.00000*constants[47]*constants[2])
    rates[6] = (-1.00000*((((algebraic[36]+algebraic[31])-(2.00000*algebraic[33])/(constants[20]-2.00000))-algebraic[43])+algebraic[38]))/(2.00000*1.00000*constants[45]*constants[2])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[7] = (states[6]*constants[29])/constants[30]
    algebraic[4] = 20.0000*exp(-0.125000*(states[0]+75.0000))
    algebraic[13] = 2000.00/(320.000*exp(-0.100000*(states[0]+75.0000))+1.00000)
    algebraic[1] = states[0]+52.0000
    algebraic[10] = 0.0500000*exp(-0.0670000*algebraic[1])
    algebraic[17] = custom_piecewise([less(fabs(algebraic[1]) , constants[8]), 2.50000 , True, algebraic[1]/(1.00000-1.00000*exp(-0.200000*algebraic[1]))])
    algebraic[2] = states[0]+22.0000
    algebraic[11] = custom_piecewise([less(fabs(algebraic[2]) , constants[11]), 2.50000 , True, (0.500000*algebraic[2])/(1.00000-exp(-algebraic[2]/5.00000))])
    algebraic[18] = custom_piecewise([less(fabs(algebraic[2]) , constants[11]), 2.50000 , True, (0.178000*algebraic[2])/(exp(algebraic[2]/15.0000)-1.00000)])
    algebraic[3] = states[0]+41.0000
    algebraic[12] = custom_piecewise([less(fabs(algebraic[3]) , constants[25]), 2000.00 , True, (200.000*algebraic[3])/(1.00000-exp(-0.100000*algebraic[3]))])
    algebraic[19] = 8000.00*exp(-0.0560000*(states[0]+66.0000))
    algebraic[5] = (states[0]+24.0000)-5.00000
    algebraic[14] = custom_piecewise([less(fabs(algebraic[5]) , constants[27]), 120.000 , True, (30.0000*algebraic[5])/(1.00000-exp((-1.00000*algebraic[5])/4.00000))])
    algebraic[20] = custom_piecewise([less(fabs(algebraic[5]) , constants[27]), 120.000 , True, (12.0000*algebraic[5])/(exp(algebraic[5]/10.0000)-1.00000)])
    algebraic[6] = states[0]+34.0000
    algebraic[15] = custom_piecewise([less(fabs(algebraic[6]) , constants[28]), 25.0000 , True, (6.25000*algebraic[6])/(exp(algebraic[6]/4.00000)-1.00000)])
    algebraic[21] = 50.0000/(1.00000+exp((-1.00000*(states[0]+34.0000))/4.00000))
    algebraic[8] = (states[0]+34.0000)--30.0000
    algebraic[16] = (0.625000*algebraic[8])/(exp(algebraic[8]/4.00000)-1.00000)
    algebraic[23] = 5.00000/(1.00000+exp((-1.00000*algebraic[8])/4.00000))
    algebraic[0] = constants[43]*log(constants[7]/states[3])
    algebraic[29] = constants[14]*(states[0]-algebraic[0])
    algebraic[32] = (((constants[17]*states[1])/(constants[18]+states[1]))*states[3])/(constants[19]+states[3])
    algebraic[33] = (constants[21]*(exp((constants[23]*(constants[20]-2.00000)*states[0])/constants[43])*(power(states[3], constants[20]))*constants[16]-exp(((constants[23]-1.00000)*(constants[20]-2.00000)*states[0])/constants[43])*(power(constants[7], constants[20]))*states[6]))/((1.00000+constants[22]*(states[6]*(power(constants[7], constants[20]))+constants[16]*(power(states[3], constants[20]))))*(1.00000+states[6]/0.00690000))
    algebraic[34] = constants[43]*log((constants[7]+0.120000*states[1])/(states[3]+0.120000*states[2]))
    algebraic[35] = constants[24]*(power(states[7], 3.00000))*states[8]*(states[0]-algebraic[34])
    algebraic[22] = ((states[4]*states[1])/(states[1]+constants[6]))*constants[4]*(states[0]-algebraic[0])
    algebraic[39] = ((0.0100000*constants[26]*(states[0]-50.0000))/(constants[43]*(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[43]))))*(states[3]*exp(50.0000/constants[43])-constants[7]*exp((-1.00000*(states[0]-50.0000))/constants[43]))*states[9]*states[10]*states[11]
    algebraic[38] = ((2.00000*1.00000*constants[45]*constants[2])/(1.00000*constants[36]*constants[34]))*states[6]*(constants[34]-states[12])
    algebraic[40] = ((2.00000*1.00000*constants[47]*constants[2])/(1.00000*constants[37]))*states[14]*(states[12]-states[13])
    algebraic[26] = (constants[10]*(states[2]-states[1]*exp(-states[0]/constants[43])))/140.000
    algebraic[27] = states[5]*algebraic[26]
    algebraic[9] = constants[43]*log(states[1]/states[2])
    algebraic[28] = (((constants[12]*states[1])/(states[1]+constants[13]))*(states[0]-algebraic[9]))/(1.00000+exp((((states[0]+10.0000)-algebraic[9])*2.00000)/constants[43]))
    algebraic[24] = ((states[4]*states[1])/(states[1]+constants[6]))*constants[5]*(states[0]-algebraic[9])
    algebraic[37] = ((0.0100000*constants[26]*(states[0]-50.0000))/(constants[43]*(1.00000-exp((-1.00000*(states[0]-50.0000))/constants[43]))))*(states[2]*exp(50.0000/constants[43])-states[1]*exp((-1.00000*(states[0]-50.0000))/constants[43]))*states[9]*states[10]*states[11]
    algebraic[41] = (algebraic[28]+algebraic[27]+algebraic[24]+algebraic[37])-2.00000*algebraic[32]
    algebraic[25] = algebraic[22]+algebraic[24]
    algebraic[30] = 0.500000*constants[43]*log(constants[16]/states[6])
    algebraic[31] = constants[15]*(states[0]-algebraic[30])
    algebraic[36] = ((4.00000*constants[26]*(states[0]-50.0000))/(constants[43]*(1.00000-exp((-1.00000*(states[0]-50.0000)*2.00000)/constants[43]))))*(states[6]*exp(100.000/constants[43])-constants[16]*exp((-2.00000*(states[0]-50.0000))/constants[43]))*states[9]*states[10]*states[11]
    algebraic[42] = algebraic[36]+algebraic[37]+algebraic[39]
    algebraic[43] = (((2.00000*1.00000*constants[47]*constants[2])/(1.00000*constants[38]))*states[13]*(power(states[6], constants[39])))/(power(states[6], constants[39])+power(constants[35], constants[39]))
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