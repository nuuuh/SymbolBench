# Size of variable arrays:
sizeAlgebraic = 64
sizeStates = 27
sizeConstants = 40
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
    legend_algebraic[56] = "i_Na in component sodium_current (nanoA)"
    legend_algebraic[36] = "i_Ca_T in component T_type_Ca_channel (nanoA)"
    legend_algebraic[34] = "i_Ca_L in component L_type_Ca_channel (nanoA)"
    legend_algebraic[59] = "i_K in component delayed_rectifying_potassium_current (nanoA)"
    legend_algebraic[39] = "i_f in component hyperpolarisation_activated_current (nanoA)"
    legend_algebraic[63] = "i_B in component linear_background_current (nanoA)"
    legend_algebraic[40] = "i_NaK in component sodium_potassium_pump (nanoA)"
    legend_algebraic[42] = "i_NaCa in component sodium_calcium_pump (nanoA)"
    legend_algebraic[41] = "i_Ca_P in component calcium_pump_current (nanoA)"
    legend_constants[4] = "P_Na in component sodium_current (mul_per_second)"
    legend_algebraic[54] = "E_Na in component reversal_potentials (millivolt)"
    legend_states[1] = "Na_c in component cleft_space_equations (millimolar)"
    legend_states[2] = "m in component sodium_current_m_gate (dimensionless)"
    legend_states[3] = "h1 in component sodium_current_h_gate (dimensionless)"
    legend_states[4] = "h2 in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[21] = "m_infinity in component sodium_current_m_gate (dimensionless)"
    legend_algebraic[27] = "tau_m in component sodium_current_m_gate (second)"
    legend_algebraic[0] = "alpha_m in component sodium_current_m_gate (per_second)"
    legend_algebraic[11] = "beta_m in component sodium_current_m_gate (per_second)"
    legend_algebraic[22] = "h1_infinity in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[33] = "h2_infinity in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[28] = "tau_h1 in component sodium_current_h_gate (second)"
    legend_algebraic[35] = "tau_h2 in component sodium_current_h_gate (second)"
    legend_algebraic[1] = "alpha_h1 in component sodium_current_h_gate (per_second)"
    legend_algebraic[12] = "beta_h1 in component sodium_current_h_gate (per_second)"
    legend_constants[5] = "g_Ca_L in component L_type_Ca_channel (microS)"
    legend_constants[6] = "E_Ca_L in component L_type_Ca_channel (millivolt)"
    legend_states[5] = "d_L in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_algebraic[32] = "d_L_infinity in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_states[6] = "f_L in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[8] = "alpha_d_L in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[19] = "beta_d_L in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[26] = "tau_d_L in component L_type_Ca_channel_d_gate (second)"
    legend_algebraic[2] = "alpha_f_L in component L_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[13] = "beta_f_L in component L_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[29] = "f_L_infinity in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[23] = "tau_f_L in component L_type_Ca_channel_f_gate (second)"
    legend_constants[7] = "g_Ca_T in component T_type_Ca_channel (microS)"
    legend_constants[8] = "E_Ca_T in component T_type_Ca_channel (millivolt)"
    legend_states[7] = "d_T in component T_type_Ca_channel_d_gate (dimensionless)"
    legend_states[8] = "f_T in component T_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[3] = "alpha_d_T in component T_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[14] = "beta_d_T in component T_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[30] = "d_T_infinity in component T_type_Ca_channel_d_gate (dimensionless)"
    legend_algebraic[24] = "tau_d_T in component T_type_Ca_channel_d_gate (second)"
    legend_algebraic[4] = "alpha_f_T in component T_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[15] = "beta_f_T in component T_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[31] = "f_T_infinity in component T_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[25] = "tau_f_T in component T_type_Ca_channel_f_gate (second)"
    legend_constants[34] = "g_K in component delayed_rectifying_potassium_current (microS)"
    legend_algebraic[58] = "E_K in component reversal_potentials (millivolt)"
    legend_constants[9] = "K_b in component cleft_space_equations (millimolar)"
    legend_states[9] = "P_a in component delayed_rectifying_potassium_current_P_a_gate (dimensionless)"
    legend_states[10] = "P_i in component delayed_rectifying_potassium_current_P_i_gate (dimensionless)"
    legend_algebraic[16] = "tau_P_a in component delayed_rectifying_potassium_current_P_a_gate (second)"
    legend_algebraic[5] = "P_a_infinity in component delayed_rectifying_potassium_current_P_a_gate (dimensionless)"
    legend_algebraic[6] = "alpha_P_i in component delayed_rectifying_potassium_current_P_i_gate (per_second)"
    legend_algebraic[17] = "beta_P_i in component delayed_rectifying_potassium_current_P_i_gate (per_second)"
    legend_algebraic[57] = "i_B_Na in component linear_background_current (nanoA)"
    legend_algebraic[62] = "i_B_Ca in component linear_background_current (nanoA)"
    legend_algebraic[60] = "i_B_K in component linear_background_current (nanoA)"
    legend_constants[10] = "g_B_Na in component linear_background_current (microS)"
    legend_constants[11] = "g_B_Ca in component linear_background_current (microS)"
    legend_constants[12] = "g_B_K in component linear_background_current (microS)"
    legend_algebraic[61] = "E_Ca in component reversal_potentials (millivolt)"
    legend_algebraic[37] = "i_f_Na in component hyperpolarisation_activated_current (nanoA)"
    legend_algebraic[38] = "i_f_K in component hyperpolarisation_activated_current (nanoA)"
    legend_constants[13] = "g_f_Na in component hyperpolarisation_activated_current (microS)"
    legend_constants[14] = "g_f_K in component hyperpolarisation_activated_current (microS)"
    legend_states[11] = "y in component hyperpolarisation_activated_current_y_gate (dimensionless)"
    legend_algebraic[7] = "y_infinity in component hyperpolarisation_activated_current_y_gate (dimensionless)"
    legend_algebraic[18] = "tau_y in component hyperpolarisation_activated_current_y_gate (second)"
    legend_constants[15] = "K_m_Na in component sodium_potassium_pump (millimolar)"
    legend_constants[16] = "K_m_K in component sodium_potassium_pump (millimolar)"
    legend_constants[17] = "i_NaK_max in component sodium_potassium_pump (nanoA)"
    legend_states[12] = "Na_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_states[13] = "K_c in component cleft_space_equations (millimolar)"
    legend_constants[18] = "i_Ca_P_max in component calcium_pump_current (nanoA)"
    legend_states[14] = "Ca_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_constants[19] = "K_NaCa in component sodium_calcium_pump (nanoA)"
    legend_constants[20] = "d_NaCa in component sodium_calcium_pump (dimensionless)"
    legend_constants[21] = "gamma in component sodium_calcium_pump (dimensionless)"
    legend_states[15] = "Ca_c in component cleft_space_equations (millimolar)"
    legend_states[16] = "K_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_states[17] = "Ca_Calmod in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_states[18] = "Ca_Trop in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_states[19] = "Ca_Mg_Trop in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_states[20] = "Mg_Mg_Trop in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_algebraic[43] = "phi_C in component intracellular_concentrations_and_buffer_equations (per_second)"
    legend_algebraic[44] = "phi_TC in component intracellular_concentrations_and_buffer_equations (per_second)"
    legend_algebraic[45] = "phi_TMgC in component intracellular_concentrations_and_buffer_equations (per_second)"
    legend_algebraic[9] = "phi_TMgM in component intracellular_concentrations_and_buffer_equations (per_second)"
    legend_algebraic[49] = "phi_B in component intracellular_concentrations_and_buffer_equations (millimolar_per_second)"
    legend_constants[22] = "Mg_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_algebraic[46] = "F_C in component intracellular_concentrations_and_buffer_equations (millimolar_per_second)"
    legend_algebraic[47] = "F_TC in component intracellular_concentrations_and_buffer_equations (millimolar_per_second)"
    legend_algebraic[48] = "F_TMgC in component intracellular_concentrations_and_buffer_equations (millimolar_per_second)"
    legend_constants[23] = "Vol in component cleft_space_equations (microLitre)"
    legend_constants[35] = "V_i in component intracellular_concentrations_and_buffer_equations (microLitre)"
    legend_algebraic[52] = "i_up in component SR_Ca_uptake_and_release (nanoA)"
    legend_algebraic[53] = "i_rel in component SR_Ca_uptake_and_release (nanoA)"
    legend_constants[24] = "Na_b in component cleft_space_equations (millimolar)"
    legend_constants[25] = "Ca_b in component cleft_space_equations (millimolar)"
    legend_constants[36] = "V_c in component cleft_space_equations (microLitre)"
    legend_constants[26] = "tau_p in component cleft_space_equations (second)"
    legend_algebraic[55] = "i_tr in component SR_Ca_uptake_and_release (nanoA)"
    legend_states[21] = "Ca_up in component SR_Ca_uptake_and_release (millimolar)"
    legend_states[22] = "Ca_rel in component SR_Ca_uptake_and_release (millimolar)"
    legend_constants[27] = "alpha_up in component SR_Ca_uptake_and_release (nanoA)"
    legend_constants[28] = "beta_up in component SR_Ca_uptake_and_release (nanoA)"
    legend_constants[29] = "alpha_rel in component SR_Ca_uptake_and_release (nanoA_per_millimolar)"
    legend_constants[37] = "K1 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_algebraic[51] = "K2 in component SR_Ca_uptake_and_release (millimolar)"
    legend_constants[30] = "k_cyca in component SR_Ca_uptake_and_release (millimolar)"
    legend_constants[31] = "k_xcs in component SR_Ca_uptake_and_release (dimensionless)"
    legend_constants[32] = "k_SRCa in component SR_Ca_uptake_and_release (millimolar)"
    legend_constants[33] = "k_rel in component SR_Ca_uptake_and_release (millimolar)"
    legend_algebraic[10] = "r_act in component SR_Ca_uptake_and_release (per_second)"
    legend_algebraic[20] = "r_inact in component SR_Ca_uptake_and_release (per_second)"
    legend_states[23] = "Ca_Calse in component SR_Ca_uptake_and_release (dimensionless)"
    legend_algebraic[50] = "phi_Calse in component SR_Ca_uptake_and_release (per_second)"
    legend_states[24] = "F1 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_states[25] = "F2 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_states[26] = "F3 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_constants[38] = "V_up in component SR_Ca_uptake_and_release (microLitre)"
    legend_constants[39] = "V_rel in component SR_Ca_uptake_and_release (microLitre)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[2] = "d/dt m in component sodium_current_m_gate (dimensionless)"
    legend_rates[3] = "d/dt h1 in component sodium_current_h_gate (dimensionless)"
    legend_rates[4] = "d/dt h2 in component sodium_current_h_gate (dimensionless)"
    legend_rates[5] = "d/dt d_L in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_rates[6] = "d/dt f_L in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_rates[7] = "d/dt d_T in component T_type_Ca_channel_d_gate (dimensionless)"
    legend_rates[8] = "d/dt f_T in component T_type_Ca_channel_f_gate (dimensionless)"
    legend_rates[9] = "d/dt P_a in component delayed_rectifying_potassium_current_P_a_gate (dimensionless)"
    legend_rates[10] = "d/dt P_i in component delayed_rectifying_potassium_current_P_i_gate (dimensionless)"
    legend_rates[11] = "d/dt y in component hyperpolarisation_activated_current_y_gate (dimensionless)"
    legend_rates[17] = "d/dt Ca_Calmod in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_rates[18] = "d/dt Ca_Trop in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_rates[19] = "d/dt Ca_Mg_Trop in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_rates[20] = "d/dt Mg_Mg_Trop in component intracellular_concentrations_and_buffer_equations (dimensionless)"
    legend_rates[12] = "d/dt Na_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_rates[16] = "d/dt K_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_rates[14] = "d/dt Ca_i in component intracellular_concentrations_and_buffer_equations (millimolar)"
    legend_rates[1] = "d/dt Na_c in component cleft_space_equations (millimolar)"
    legend_rates[13] = "d/dt K_c in component cleft_space_equations (millimolar)"
    legend_rates[15] = "d/dt Ca_c in component cleft_space_equations (millimolar)"
    legend_rates[23] = "d/dt Ca_Calse in component SR_Ca_uptake_and_release (dimensionless)"
    legend_rates[24] = "d/dt F1 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_rates[25] = "d/dt F2 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_rates[26] = "d/dt F3 in component SR_Ca_uptake_and_release (dimensionless)"
    legend_rates[21] = "d/dt Ca_up in component SR_Ca_uptake_and_release (millimolar)"
    legend_rates[22] = "d/dt Ca_rel in component SR_Ca_uptake_and_release (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -49.54105
    constants[0] = 8314.472
    constants[1] = 310
    constants[2] = 96485.3415
    constants[3] = 5.5e-5
    constants[4] = 0.00344
    states[1] = 139.9988
    states[2] = 0.250113
    states[3] = 0.001386897
    states[4] = 0.002065463
    constants[5] = 0.02115
    constants[6] = 46.4
    states[5] = 0.002572773
    states[6] = 0.98651
    constants[7] = 0.02521
    constants[8] = 45
    states[7] = 0.02012114
    states[8] = 0.1945111
    constants[9] = 5.4
    states[9] = 0.02302278
    states[10] = 0.3777728
    constants[10] = 0.00016
    constants[11] = 0.0000364
    constants[12] = 0.0000694
    constants[13] = 0.0067478
    constants[14] = 0.0128821
    states[11] = 0.09227776
    constants[15] = 5.46
    constants[16] = 0.621
    constants[17] = 0.2192
    states[12] = 9.701621
    states[13] = 5.389014
    constants[18] = 0.02869
    states[14] = 3.787018e-4
    constants[19] = 0.00001248
    constants[20] = 0.0001
    constants[21] = 0.5
    states[15] = 2.00474
    states[16] = 1.407347e2
    states[17] = 0.1411678
    states[18] = 0.07331396
    states[19] = 0.7618549
    states[20] = 0.2097049
    constants[22] = 2.5
    constants[23] = 3.497e-6
    constants[24] = 140
    constants[25] = 2
    constants[26] = 0.01
    states[21] = 16.95311
    states[22] = 16.85024
    constants[27] = 0.08
    constants[28] = 0.072
    constants[29] = 0.5
    constants[30] = 0.00005
    constants[31] = 0.9
    constants[32] = 22
    constants[33] = 0.004
    states[23] = 0.9528726
    states[24] = 0.1133251
    states[25] = 0.0007594214
    states[26] = 0.8859153
    constants[34] = 0.00693000*(power(constants[9]/1.00000, 0.590000))
    constants[35] = 0.465000*constants[23]
    constants[36] = 0.136000*constants[23]
    constants[37] = (constants[30]*constants[31])/constants[32]
    constants[38] = 0.0116600*constants[35]
    constants[39] = 0.00129600*constants[35]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[9] = 1290.00*constants[22]*(1.00000-(states[19]+states[20]))-429.000*states[20]
    rates[20] = algebraic[9]
    algebraic[10] = 240.000*exp((states[0]-40.0000)*0.0800000)+240.000*(power(states[14]/(states[14]+constants[33]), 4.00000))
    rates[24] = 0.960000*states[26]-algebraic[10]*states[24]
    algebraic[16] = 1.00000/(17.0000*exp(0.0398000*states[0])+2.11000*exp(-0.0510000*states[0]))
    algebraic[5] = 1.00000/(1.00000+exp((states[0]+5.10000)/-7.40000))
    rates[9] = (algebraic[5]-states[9])/algebraic[16]
    algebraic[6] = 100.000*exp(-0.0183000*states[0])
    algebraic[17] = 656.000*exp(0.00942000*states[0])
    rates[10] = algebraic[6]*(1.00000-states[10])-algebraic[17]*states[10]
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+72.2000)/9.00000))
    algebraic[18] = 1.00000/(1.64830*exp((states[0]+54.0600)/-24.3300)+14.0106/(0.700000+exp((states[0]+60.0000)/-5.50000)))
    rates[11] = (algebraic[7]-states[11])/algebraic[18]
    algebraic[20] = 40.0000+240.000*(power(states[14]/(states[14]+constants[33]), 4.00000))
    rates[25] = algebraic[10]*states[24]-algebraic[20]*states[25]
    rates[26] = algebraic[20]*states[25]-0.960000*states[26]
    algebraic[0] = (-824.000*(states[0]+51.9000))/(exp((states[0]+51.9000)/-8.90000)-1.00000)
    algebraic[11] = 32960.0*exp((states[0]+51.9000)/-8.90000)
    algebraic[21] = algebraic[0]/(algebraic[0]+algebraic[11])
    algebraic[27] = 1.00000/(algebraic[0]+algebraic[11])+1.50000e-05
    rates[2] = (algebraic[21]-states[2])/algebraic[27]
    algebraic[1] = 165.000*exp((states[0]+101.300)/-12.6000)
    algebraic[12] = 12360.0/(320.000*exp((states[0]+101.300)/-12.6000)+1.00000)
    algebraic[22] = algebraic[1]/(algebraic[1]+algebraic[12])
    algebraic[28] = 1.00000/(algebraic[1]+algebraic[12])
    rates[3] = (algebraic[22]-states[3])/algebraic[28]
    algebraic[32] = 1.00000/(1.00000+exp((states[0]+14.1000)/-6.00000))
    algebraic[8] = (-28.3900*(states[0]+35.0000))/(exp((states[0]+35.0000)/-2.50000)-1.00000)+(-84.9000*states[0])/(exp(-0.208000*states[0])-1.00000)
    algebraic[19] = (11.4300*(states[0]-5.00000))/(exp(0.400000*(states[0]-5.00000))-1.00000)
    algebraic[26] = 1.00000/(algebraic[8]+algebraic[19])
    rates[5] = (algebraic[32]-states[5])/algebraic[26]
    algebraic[29] = 1.00000/(1.00000+exp((states[0]+30.0000)/5.00000))
    algebraic[2] = (3.75000*(states[0]+28.0000))/(exp((states[0]+28.0000)/4.00000)-1.00000)
    algebraic[13] = 30.0000/(1.00000+exp((states[0]+28.0000)/-4.00000))
    algebraic[23] = 1.00000/(algebraic[2]+algebraic[13])
    rates[6] = (algebraic[29]-states[6])/algebraic[23]
    algebraic[30] = 1.00000/(1.00000+exp((states[0]+26.3000)/-6.00000))
    algebraic[3] = 1068.00*exp((states[0]+26.3000)/30.0000)
    algebraic[14] = 1068.00*exp((states[0]+26.3000)/-30.0000)
    algebraic[24] = 1.00000/(algebraic[3]+algebraic[14])
    rates[7] = (algebraic[30]-states[7])/algebraic[24]
    algebraic[31] = 1.00000/(1.00000+exp((states[0]+61.7000)/5.60000))
    algebraic[4] = 15.3000*exp((states[0]+61.7000)/-83.3000)
    algebraic[15] = 15.0000*exp((states[0]+61.7000)/15.3800)
    algebraic[25] = 1.00000/(algebraic[4]+algebraic[15])
    rates[8] = (algebraic[31]-states[8])/algebraic[25]
    algebraic[33] = algebraic[22]
    algebraic[35] = 20.0000*algebraic[28]
    rates[4] = (algebraic[33]-states[4])/algebraic[35]
    algebraic[43] = 129000.*states[14]*(1.00000-states[17])-307.000*states[17]
    rates[17] = algebraic[43]
    algebraic[44] = 50500.0*states[14]*(1.00000-states[18])-252.000*states[18]
    rates[18] = algebraic[44]
    algebraic[45] = 129000.*states[14]*(1.00000-(states[19]+states[20]))-4.25000*states[19]
    rates[19] = algebraic[45]
    algebraic[50] = 770.000*states[22]*(1.00000-states[23])-641.000*states[23]
    rates[23] = algebraic[50]
    algebraic[51] = states[14]+states[21]*constants[37]+constants[30]*constants[31]+constants[30]
    algebraic[52] = (constants[27]*states[14]-constants[28]*states[21]*constants[37])/algebraic[51]
    algebraic[55] = ((states[21]-states[22])*2.00000*constants[2]*constants[38])/0.0641800
    rates[21] = (algebraic[52]-algebraic[55])/(2.00000*constants[38]*constants[2])
    algebraic[53] = constants[29]*(power(states[25]/(states[25]+0.250000), 2.00000))*states[22]
    rates[22] = (algebraic[55]-algebraic[53])/(2.00000*constants[39]*constants[2])-11.4800*algebraic[50]
    algebraic[54] = ((constants[0]*constants[1])/constants[2])*log(states[1]/states[12])
    algebraic[56] = (((constants[4]*(power(states[2], 3.00000))*states[3]*states[4]*states[1]*states[0]*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(exp(((states[0]-algebraic[54])*constants[2])/(constants[0]*constants[1]))-1.00000))/(exp((states[0]*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[40] = (constants[17]*(power(states[12]/(constants[15]+states[12]), 3.00000))*(power(states[13]/(constants[16]+states[13]), 2.00000))*1.60000)/(1.50000+exp((states[0]+60.0000)/-40.0000))
    algebraic[42] = (constants[19]*((power(states[12], 3.00000))*states[15]*exp(0.0374300*states[0]*constants[21])-(power(states[1], 3.00000))*states[14]*exp(0.0374300*states[0]*(constants[21]-1.00000))))/(1.00000+constants[20]*(states[14]*(power(states[1], 3.00000))+states[15]*(power(states[12], 3.00000))))
    algebraic[57] = constants[10]*(states[0]-algebraic[54])
    algebraic[37] = constants[13]*(power(states[11], 2.00000))*(states[0]-75.0000)
    rates[12] = -(3.00000*algebraic[40]+3.00000*algebraic[42]+algebraic[57]+algebraic[37]+algebraic[56])/(constants[2]*constants[35])
    rates[1] = (constants[24]-states[1])/constants[26]+(algebraic[56]+3.00000*algebraic[42]+3.00000*algebraic[40]+algebraic[57]+algebraic[37])/(constants[2]*constants[36])
    algebraic[58] = ((constants[0]*constants[1])/constants[2])*log(states[13]/states[16])
    algebraic[59] = constants[34]*states[9]*states[10]*(states[0]-algebraic[58])
    algebraic[60] = constants[12]*(states[0]-algebraic[58])
    algebraic[38] = constants[14]*(power(states[11], 2.00000))*(states[0]+85.0000)
    rates[16] = (2.00000*algebraic[40]-(algebraic[59]+algebraic[38]+algebraic[60]))/(constants[2]*constants[35])
    rates[13] = (constants[9]-states[13])/constants[26]+(-2.00000*algebraic[40]+algebraic[59]+algebraic[60]+algebraic[38])/(constants[2]*constants[36])
    algebraic[36] = constants[7]*states[7]*states[8]*(states[0]-constants[8])
    algebraic[34] = constants[5]*(states[6]*states[5]+0.0950000*algebraic[32])*(states[0]-constants[6])
    algebraic[41] = (constants[18]*states[14])/(states[14]+0.000400000)
    algebraic[61] = ((0.500000*constants[0]*constants[1])/constants[2])*log(states[15]/states[14])
    algebraic[62] = constants[11]*(states[0]-algebraic[61])
    algebraic[46] = 0.0900000*algebraic[43]
    algebraic[47] = 0.0310000*algebraic[44]
    algebraic[48] = 0.0620000*algebraic[45]
    algebraic[49] = algebraic[46]+algebraic[47]+algebraic[48]
    rates[14] = ((2.00000*algebraic[42]+algebraic[53])-(algebraic[34]+algebraic[36]+algebraic[41]+algebraic[62]+algebraic[52]))/(2.00000*constants[35]*constants[2])-algebraic[49]
    rates[15] = (constants[25]-states[15])/constants[26]+(-2.00000*algebraic[42]+algebraic[34]+algebraic[36]+algebraic[41]+algebraic[62])/(2.00000*constants[2]*constants[36])
    algebraic[39] = algebraic[37]+algebraic[38]
    algebraic[63] = algebraic[57]+algebraic[62]+algebraic[60]
    rates[0] = -(algebraic[56]+algebraic[36]+algebraic[34]+algebraic[59]+algebraic[39]+algebraic[63]+algebraic[40]+algebraic[42]+algebraic[41])/constants[3]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[9] = 1290.00*constants[22]*(1.00000-(states[19]+states[20]))-429.000*states[20]
    algebraic[10] = 240.000*exp((states[0]-40.0000)*0.0800000)+240.000*(power(states[14]/(states[14]+constants[33]), 4.00000))
    algebraic[16] = 1.00000/(17.0000*exp(0.0398000*states[0])+2.11000*exp(-0.0510000*states[0]))
    algebraic[5] = 1.00000/(1.00000+exp((states[0]+5.10000)/-7.40000))
    algebraic[6] = 100.000*exp(-0.0183000*states[0])
    algebraic[17] = 656.000*exp(0.00942000*states[0])
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+72.2000)/9.00000))
    algebraic[18] = 1.00000/(1.64830*exp((states[0]+54.0600)/-24.3300)+14.0106/(0.700000+exp((states[0]+60.0000)/-5.50000)))
    algebraic[20] = 40.0000+240.000*(power(states[14]/(states[14]+constants[33]), 4.00000))
    algebraic[0] = (-824.000*(states[0]+51.9000))/(exp((states[0]+51.9000)/-8.90000)-1.00000)
    algebraic[11] = 32960.0*exp((states[0]+51.9000)/-8.90000)
    algebraic[21] = algebraic[0]/(algebraic[0]+algebraic[11])
    algebraic[27] = 1.00000/(algebraic[0]+algebraic[11])+1.50000e-05
    algebraic[1] = 165.000*exp((states[0]+101.300)/-12.6000)
    algebraic[12] = 12360.0/(320.000*exp((states[0]+101.300)/-12.6000)+1.00000)
    algebraic[22] = algebraic[1]/(algebraic[1]+algebraic[12])
    algebraic[28] = 1.00000/(algebraic[1]+algebraic[12])
    algebraic[32] = 1.00000/(1.00000+exp((states[0]+14.1000)/-6.00000))
    algebraic[8] = (-28.3900*(states[0]+35.0000))/(exp((states[0]+35.0000)/-2.50000)-1.00000)+(-84.9000*states[0])/(exp(-0.208000*states[0])-1.00000)
    algebraic[19] = (11.4300*(states[0]-5.00000))/(exp(0.400000*(states[0]-5.00000))-1.00000)
    algebraic[26] = 1.00000/(algebraic[8]+algebraic[19])
    algebraic[29] = 1.00000/(1.00000+exp((states[0]+30.0000)/5.00000))
    algebraic[2] = (3.75000*(states[0]+28.0000))/(exp((states[0]+28.0000)/4.00000)-1.00000)
    algebraic[13] = 30.0000/(1.00000+exp((states[0]+28.0000)/-4.00000))
    algebraic[23] = 1.00000/(algebraic[2]+algebraic[13])
    algebraic[30] = 1.00000/(1.00000+exp((states[0]+26.3000)/-6.00000))
    algebraic[3] = 1068.00*exp((states[0]+26.3000)/30.0000)
    algebraic[14] = 1068.00*exp((states[0]+26.3000)/-30.0000)
    algebraic[24] = 1.00000/(algebraic[3]+algebraic[14])
    algebraic[31] = 1.00000/(1.00000+exp((states[0]+61.7000)/5.60000))
    algebraic[4] = 15.3000*exp((states[0]+61.7000)/-83.3000)
    algebraic[15] = 15.0000*exp((states[0]+61.7000)/15.3800)
    algebraic[25] = 1.00000/(algebraic[4]+algebraic[15])
    algebraic[33] = algebraic[22]
    algebraic[35] = 20.0000*algebraic[28]
    algebraic[43] = 129000.*states[14]*(1.00000-states[17])-307.000*states[17]
    algebraic[44] = 50500.0*states[14]*(1.00000-states[18])-252.000*states[18]
    algebraic[45] = 129000.*states[14]*(1.00000-(states[19]+states[20]))-4.25000*states[19]
    algebraic[50] = 770.000*states[22]*(1.00000-states[23])-641.000*states[23]
    algebraic[51] = states[14]+states[21]*constants[37]+constants[30]*constants[31]+constants[30]
    algebraic[52] = (constants[27]*states[14]-constants[28]*states[21]*constants[37])/algebraic[51]
    algebraic[55] = ((states[21]-states[22])*2.00000*constants[2]*constants[38])/0.0641800
    algebraic[53] = constants[29]*(power(states[25]/(states[25]+0.250000), 2.00000))*states[22]
    algebraic[54] = ((constants[0]*constants[1])/constants[2])*log(states[1]/states[12])
    algebraic[56] = (((constants[4]*(power(states[2], 3.00000))*states[3]*states[4]*states[1]*states[0]*(power(constants[2], 2.00000)))/(constants[0]*constants[1]))*(exp(((states[0]-algebraic[54])*constants[2])/(constants[0]*constants[1]))-1.00000))/(exp((states[0]*constants[2])/(constants[0]*constants[1]))-1.00000)
    algebraic[40] = (constants[17]*(power(states[12]/(constants[15]+states[12]), 3.00000))*(power(states[13]/(constants[16]+states[13]), 2.00000))*1.60000)/(1.50000+exp((states[0]+60.0000)/-40.0000))
    algebraic[42] = (constants[19]*((power(states[12], 3.00000))*states[15]*exp(0.0374300*states[0]*constants[21])-(power(states[1], 3.00000))*states[14]*exp(0.0374300*states[0]*(constants[21]-1.00000))))/(1.00000+constants[20]*(states[14]*(power(states[1], 3.00000))+states[15]*(power(states[12], 3.00000))))
    algebraic[57] = constants[10]*(states[0]-algebraic[54])
    algebraic[37] = constants[13]*(power(states[11], 2.00000))*(states[0]-75.0000)
    algebraic[58] = ((constants[0]*constants[1])/constants[2])*log(states[13]/states[16])
    algebraic[59] = constants[34]*states[9]*states[10]*(states[0]-algebraic[58])
    algebraic[60] = constants[12]*(states[0]-algebraic[58])
    algebraic[38] = constants[14]*(power(states[11], 2.00000))*(states[0]+85.0000)
    algebraic[36] = constants[7]*states[7]*states[8]*(states[0]-constants[8])
    algebraic[34] = constants[5]*(states[6]*states[5]+0.0950000*algebraic[32])*(states[0]-constants[6])
    algebraic[41] = (constants[18]*states[14])/(states[14]+0.000400000)
    algebraic[61] = ((0.500000*constants[0]*constants[1])/constants[2])*log(states[15]/states[14])
    algebraic[62] = constants[11]*(states[0]-algebraic[61])
    algebraic[46] = 0.0900000*algebraic[43]
    algebraic[47] = 0.0310000*algebraic[44]
    algebraic[48] = 0.0620000*algebraic[45]
    algebraic[49] = algebraic[46]+algebraic[47]+algebraic[48]
    algebraic[39] = algebraic[37]+algebraic[38]
    algebraic[63] = algebraic[57]+algebraic[62]+algebraic[60]
    return algebraic

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