# Size of variable arrays:
sizeAlgebraic = 52
sizeStates = 15
sizeConstants = 133
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "Version in component membrane (dimensionless)"
    legend_constants[1] = "dCell in component membrane (dimensionless)"
    legend_constants[2] = "FCellConstant in component membrane (dimensionless)"
    legend_constants[110] = "FCell in component membrane (dimensionless)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[3] = "R in component membrane (millijoule_per_mole_kelvin)"
    legend_constants[4] = "T in component membrane (kelvin)"
    legend_constants[5] = "F in component membrane (coulomb_per_mole)"
    legend_constants[116] = "Cm in component membrane (microF)"
    legend_constants[6] = "CmCentre in component membrane (microF)"
    legend_constants[7] = "CmPeriphery in component membrane (microF)"
    legend_algebraic[27] = "i_Na in component sodium_current (nanoA)"
    legend_algebraic[34] = "i_Ca_L in component L_type_Ca_channel (nanoA)"
    legend_algebraic[39] = "i_Ca_T in component T_type_Ca_channel (nanoA)"
    legend_algebraic[40] = "i_to in component four_AP_sensitive_currents (nanoA)"
    legend_algebraic[41] = "i_sus in component four_AP_sensitive_currents (nanoA)"
    legend_algebraic[43] = "i_K_r in component rapid_delayed_rectifying_potassium_current (nanoA)"
    legend_algebraic[44] = "i_K_s in component slow_delayed_rectifying_potassium_current (nanoA)"
    legend_algebraic[45] = "i_f_Na in component hyperpolarisation_activated_current (nanoA)"
    legend_algebraic[46] = "i_f_K in component hyperpolarisation_activated_current (nanoA)"
    legend_algebraic[47] = "i_b_Na in component sodium_background_current (nanoA)"
    legend_algebraic[49] = "i_b_Ca in component calcium_background_current (nanoA)"
    legend_algebraic[48] = "i_b_K in component potassium_background_current (nanoA)"
    legend_algebraic[50] = "i_NaCa in component sodium_calcium_exchanger (nanoA)"
    legend_algebraic[51] = "i_p in component sodium_potassium_pump (nanoA)"
    legend_constants[132] = "i_Ca_p in component persistent_calcium_current (nanoA)"
    legend_constants[117] = "g_Na in component sodium_current (microlitre_per_second)"
    legend_constants[8] = "g_Na_Centre_Published in component sodium_current (microlitre_per_second)"
    legend_constants[9] = "g_Na_Centre_0DCapable in component sodium_current (microlitre_per_second)"
    legend_constants[10] = "g_Na_Centre_1DCapable in component sodium_current (microlitre_per_second)"
    legend_constants[11] = "g_Na_Periphery_Published in component sodium_current (microlitre_per_second)"
    legend_constants[12] = "g_Na_Periphery_0DCapable in component sodium_current (microlitre_per_second)"
    legend_constants[13] = "g_Na_Periphery_1DCapable in component sodium_current (microlitre_per_second)"
    legend_constants[111] = "E_Na in component reversal_and_equilibrium_potentials (millivolt)"
    legend_constants[14] = "Na_o in component ionic_concentrations (millimolar)"
    legend_states[1] = "m in component sodium_current_m_gate (dimensionless)"
    legend_algebraic[13] = "h in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[1] = "m_infinity in component sodium_current_m_gate (dimensionless)"
    legend_algebraic[14] = "tau_m in component sodium_current_m_gate (second)"
    legend_algebraic[0] = "F_Na in component sodium_current_h_gate (dimensionless)"
    legend_states[2] = "h1 in component sodium_current_h_gate (dimensionless)"
    legend_states[3] = "h2 in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[2] = "h1_infinity in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[15] = "h2_infinity in component sodium_current_h_gate (dimensionless)"
    legend_algebraic[16] = "tau_h1 in component sodium_current_h_gate (second)"
    legend_algebraic[28] = "tau_h2 in component sodium_current_h_gate (second)"
    legend_constants[15] = "g_Ca_L_Centre_Published in component L_type_Ca_channel (microS)"
    legend_constants[16] = "g_Ca_L_Centre_0DCapable in component L_type_Ca_channel (microS)"
    legend_constants[17] = "g_Ca_L_Centre_1DCapable in component L_type_Ca_channel (microS)"
    legend_constants[18] = "g_Ca_L_Periphery_Published in component L_type_Ca_channel (microS)"
    legend_constants[19] = "g_Ca_L_Periphery_0DCapable in component L_type_Ca_channel (microS)"
    legend_constants[20] = "g_Ca_L_Periphery_1DCapable in component L_type_Ca_channel (microS)"
    legend_constants[118] = "g_Ca_L in component L_type_Ca_channel (microS)"
    legend_constants[21] = "E_Ca_L in component L_type_Ca_channel (millivolt)"
    legend_states[4] = "d_L in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_states[5] = "f_L in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[3] = "alpha_d_L in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[17] = "beta_d_L in component L_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[35] = "d_L_infinity in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_algebraic[29] = "tau_d_L in component L_type_Ca_channel_d_gate (second)"
    legend_algebraic[4] = "alpha_f_L in component L_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[18] = "beta_f_L in component L_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[36] = "f_L_infinity in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[30] = "tau_f_L in component L_type_Ca_channel_f_gate (second)"
    legend_constants[22] = "g_Ca_T_Centre_Published in component T_type_Ca_channel (microS)"
    legend_constants[23] = "g_Ca_T_Centre_0DCapable in component T_type_Ca_channel (microS)"
    legend_constants[24] = "g_Ca_T_Centre_1DCapable in component T_type_Ca_channel (microS)"
    legend_constants[25] = "g_Ca_T_Periphery_Published in component T_type_Ca_channel (microS)"
    legend_constants[26] = "g_Ca_T_Periphery_0DCapable in component T_type_Ca_channel (microS)"
    legend_constants[27] = "g_Ca_T_Periphery_1DCapable in component T_type_Ca_channel (microS)"
    legend_constants[119] = "g_Ca_T in component T_type_Ca_channel (microS)"
    legend_constants[28] = "E_Ca_T in component T_type_Ca_channel (millivolt)"
    legend_states[6] = "d_T in component T_type_Ca_channel_d_gate (dimensionless)"
    legend_states[7] = "f_T in component T_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[5] = "alpha_d_T in component T_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[19] = "beta_d_T in component T_type_Ca_channel_d_gate (per_second)"
    legend_algebraic[37] = "d_T_infinity in component T_type_Ca_channel_d_gate (dimensionless)"
    legend_algebraic[31] = "tau_d_T in component T_type_Ca_channel_d_gate (second)"
    legend_algebraic[6] = "alpha_f_T in component T_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[20] = "beta_f_T in component T_type_Ca_channel_f_gate (per_second)"
    legend_algebraic[38] = "f_T_infinity in component T_type_Ca_channel_f_gate (dimensionless)"
    legend_algebraic[32] = "tau_f_T in component T_type_Ca_channel_f_gate (second)"
    legend_constants[29] = "g_to_Centre_Published in component four_AP_sensitive_currents (microS)"
    legend_constants[30] = "g_to_Centre_0DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[31] = "g_to_Centre_1DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[32] = "g_to_Periphery_Published in component four_AP_sensitive_currents (microS)"
    legend_constants[33] = "g_to_Periphery_0DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[34] = "g_to_Periphery_1DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[120] = "g_to in component four_AP_sensitive_currents (microS)"
    legend_constants[35] = "g_sus_Centre_Published in component four_AP_sensitive_currents (microS)"
    legend_constants[36] = "g_sus_Centre_0DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[37] = "g_sus_Centre_1DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[38] = "g_sus_Periphery_Published in component four_AP_sensitive_currents (microS)"
    legend_constants[39] = "g_sus_Periphery_0DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[40] = "g_sus_Periphery_1DCapable in component four_AP_sensitive_currents (microS)"
    legend_constants[121] = "g_sus in component four_AP_sensitive_currents (microS)"
    legend_constants[112] = "E_K in component reversal_and_equilibrium_potentials (millivolt)"
    legend_states[8] = "q in component four_AP_sensitive_currents_q_gate (dimensionless)"
    legend_states[9] = "r in component four_AP_sensitive_currents_r_gate (dimensionless)"
    legend_algebraic[7] = "q_infinity in component four_AP_sensitive_currents_q_gate (dimensionless)"
    legend_algebraic[21] = "tau_q in component four_AP_sensitive_currents_q_gate (second)"
    legend_algebraic[8] = "r_infinity in component four_AP_sensitive_currents_r_gate (dimensionless)"
    legend_algebraic[22] = "tau_r in component four_AP_sensitive_currents_r_gate (second)"
    legend_constants[41] = "g_K_r_Centre_Published in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_constants[42] = "g_K_r_Centre_0DCapable in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_constants[43] = "g_K_r_Centre_1DCapable in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_constants[44] = "g_K_r_Periphery_Published in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_constants[45] = "g_K_r_Periphery_0DCapable in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_constants[46] = "g_K_r_Periphery_1DCapable in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_constants[122] = "g_K_r in component rapid_delayed_rectifying_potassium_current (microS)"
    legend_algebraic[42] = "P_a in component rapid_delayed_rectifying_potassium_current (dimensionless)"
    legend_states[10] = "P_af in component rapid_delayed_rectifying_potassium_current_P_af_gate (dimensionless)"
    legend_states[11] = "P_as in component rapid_delayed_rectifying_potassium_current_P_as_gate (dimensionless)"
    legend_states[12] = "P_i in component rapid_delayed_rectifying_potassium_current_P_i_gate (dimensionless)"
    legend_algebraic[9] = "P_af_infinity in component rapid_delayed_rectifying_potassium_current_P_af_gate (dimensionless)"
    legend_algebraic[23] = "tau_P_af in component rapid_delayed_rectifying_potassium_current_P_af_gate (second)"
    legend_algebraic[24] = "P_as_infinity in component rapid_delayed_rectifying_potassium_current_P_as_gate (dimensionless)"
    legend_algebraic[33] = "tau_P_as in component rapid_delayed_rectifying_potassium_current_P_as_gate (second)"
    legend_algebraic[10] = "P_i_infinity in component rapid_delayed_rectifying_potassium_current_P_i_gate (dimensionless)"
    legend_constants[113] = "tau_P_i in component rapid_delayed_rectifying_potassium_current_P_i_gate (second)"
    legend_constants[47] = "g_K_s_Centre_Published in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[48] = "g_K_s_Centre_0DCapable in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[49] = "g_K_s_Centre_1DCapable in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[50] = "g_K_s_Periphery_Published in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[51] = "g_K_s_Periphery_0DCapable in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[52] = "g_K_s_Periphery_1DCapable in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[123] = "g_K_s in component slow_delayed_rectifying_potassium_current (microS)"
    legend_constants[114] = "E_K_s in component reversal_and_equilibrium_potentials (millivolt)"
    legend_states[13] = "xs in component slow_delayed_rectifying_potassium_current_xs_gate (dimensionless)"
    legend_algebraic[11] = "alpha_xs in component slow_delayed_rectifying_potassium_current_xs_gate (per_second)"
    legend_algebraic[25] = "beta_xs in component slow_delayed_rectifying_potassium_current_xs_gate (per_second)"
    legend_constants[53] = "g_f_Na_Centre_Published in component hyperpolarisation_activated_current (microS)"
    legend_constants[54] = "g_f_Na_Centre_0DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[55] = "g_f_Na_Centre_1DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[56] = "g_f_Na_Periphery_Published in component hyperpolarisation_activated_current (microS)"
    legend_constants[57] = "g_f_Na_Periphery_0DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[58] = "g_f_Na_Periphery_1DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[124] = "g_f_Na in component hyperpolarisation_activated_current (microS)"
    legend_constants[59] = "g_f_K_Centre_Published in component hyperpolarisation_activated_current (microS)"
    legend_constants[60] = "g_f_K_Centre_0DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[61] = "g_f_K_Centre_1DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[62] = "g_f_K_Periphery_Published in component hyperpolarisation_activated_current (microS)"
    legend_constants[63] = "g_f_K_Periphery_0DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[64] = "g_f_K_Periphery_1DCapable in component hyperpolarisation_activated_current (microS)"
    legend_constants[125] = "g_f_K in component hyperpolarisation_activated_current (microS)"
    legend_states[14] = "y in component hyperpolarisation_activated_current_y_gate (dimensionless)"
    legend_algebraic[12] = "alpha_y in component hyperpolarisation_activated_current_y_gate (per_second)"
    legend_algebraic[26] = "beta_y in component hyperpolarisation_activated_current_y_gate (per_second)"
    legend_constants[65] = "g_b_Na_Centre_Published in component sodium_background_current (microS)"
    legend_constants[66] = "g_b_Na_Centre_0DCapable in component sodium_background_current (microS)"
    legend_constants[67] = "g_b_Na_Centre_1DCapable in component sodium_background_current (microS)"
    legend_constants[68] = "g_b_Na_Periphery_Published in component sodium_background_current (microS)"
    legend_constants[69] = "g_b_Na_Periphery_0DCapable in component sodium_background_current (microS)"
    legend_constants[70] = "g_b_Na_Periphery_1DCapable in component sodium_background_current (microS)"
    legend_constants[126] = "g_b_Na in component sodium_background_current (microS)"
    legend_constants[71] = "g_b_K_Centre_Published in component potassium_background_current (microS)"
    legend_constants[72] = "g_b_K_Centre_0DCapable in component potassium_background_current (microS)"
    legend_constants[73] = "g_b_K_Centre_1DCapable in component potassium_background_current (microS)"
    legend_constants[74] = "g_b_K_Periphery_Published in component potassium_background_current (microS)"
    legend_constants[75] = "g_b_K_Periphery_0DCapable in component potassium_background_current (microS)"
    legend_constants[76] = "g_b_K_Periphery_1DCapable in component potassium_background_current (microS)"
    legend_constants[127] = "g_b_K in component potassium_background_current (microS)"
    legend_constants[77] = "g_b_Ca_Centre_Published in component calcium_background_current (microS)"
    legend_constants[78] = "g_b_Ca_Centre_0DCapable in component calcium_background_current (microS)"
    legend_constants[79] = "g_b_Ca_Centre_1DCapable in component calcium_background_current (microS)"
    legend_constants[80] = "g_b_Ca_Periphery_Published in component calcium_background_current (microS)"
    legend_constants[81] = "g_b_Ca_Periphery_0DCapable in component calcium_background_current (microS)"
    legend_constants[82] = "g_b_Ca_Periphery_1DCapable in component calcium_background_current (microS)"
    legend_constants[128] = "g_b_Ca in component calcium_background_current (microS)"
    legend_constants[115] = "E_Ca in component reversal_and_equilibrium_potentials (millivolt)"
    legend_constants[83] = "k_NaCa_Centre_Published in component sodium_calcium_exchanger (nanoA)"
    legend_constants[84] = "k_NaCa_Centre_0DCapable in component sodium_calcium_exchanger (nanoA)"
    legend_constants[85] = "k_NaCa_Centre_1DCapable in component sodium_calcium_exchanger (nanoA)"
    legend_constants[86] = "k_NaCa_Periphery_Published in component sodium_calcium_exchanger (nanoA)"
    legend_constants[87] = "k_NaCa_Periphery_0DCapable in component sodium_calcium_exchanger (nanoA)"
    legend_constants[88] = "k_NaCa_Periphery_1DCapable in component sodium_calcium_exchanger (nanoA)"
    legend_constants[129] = "k_NaCa in component sodium_calcium_exchanger (nanoA)"
    legend_constants[89] = "d_NaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[90] = "gamma_NaCa in component sodium_calcium_exchanger (dimensionless)"
    legend_constants[91] = "Na_i in component ionic_concentrations (millimolar)"
    legend_constants[92] = "Ca_i in component ionic_concentrations (millimolar)"
    legend_constants[93] = "Ca_o in component ionic_concentrations (millimolar)"
    legend_constants[94] = "K_m_Na in component sodium_potassium_pump (millimolar)"
    legend_constants[95] = "K_m_K in component sodium_potassium_pump (millimolar)"
    legend_constants[96] = "i_p_max_Centre_Published in component sodium_potassium_pump (nanoA)"
    legend_constants[97] = "i_p_max_Centre_0DCapable in component sodium_potassium_pump (nanoA)"
    legend_constants[98] = "i_p_max_Centre_1DCapable in component sodium_potassium_pump (nanoA)"
    legend_constants[99] = "i_p_max_Periphery_Published in component sodium_potassium_pump (nanoA)"
    legend_constants[100] = "i_p_max_Periphery_0DCapable in component sodium_potassium_pump (nanoA)"
    legend_constants[101] = "i_p_max_Periphery_1DCapable in component sodium_potassium_pump (nanoA)"
    legend_constants[130] = "i_p_max in component sodium_potassium_pump (nanoA)"
    legend_constants[102] = "K_o in component ionic_concentrations (millimolar)"
    legend_constants[103] = "i_Ca_p_max_Centre_Published in component persistent_calcium_current (nanoA)"
    legend_constants[104] = "i_Ca_p_max_Centre_0DCapable in component persistent_calcium_current (nanoA)"
    legend_constants[105] = "i_Ca_p_max_Centre_1DCapable in component persistent_calcium_current (nanoA)"
    legend_constants[106] = "i_Ca_p_max_Periphery_Published in component persistent_calcium_current (nanoA)"
    legend_constants[107] = "i_Ca_p_max_Periphery_0DCapable in component persistent_calcium_current (nanoA)"
    legend_constants[108] = "i_Ca_p_max_Periphery_1DCapable in component persistent_calcium_current (nanoA)"
    legend_constants[131] = "i_Ca_p_max in component persistent_calcium_current (nanoA)"
    legend_constants[109] = "K_i in component ionic_concentrations (millimolar)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt m in component sodium_current_m_gate (dimensionless)"
    legend_rates[2] = "d/dt h1 in component sodium_current_h_gate (dimensionless)"
    legend_rates[3] = "d/dt h2 in component sodium_current_h_gate (dimensionless)"
    legend_rates[4] = "d/dt d_L in component L_type_Ca_channel_d_gate (dimensionless)"
    legend_rates[5] = "d/dt f_L in component L_type_Ca_channel_f_gate (dimensionless)"
    legend_rates[6] = "d/dt d_T in component T_type_Ca_channel_d_gate (dimensionless)"
    legend_rates[7] = "d/dt f_T in component T_type_Ca_channel_f_gate (dimensionless)"
    legend_rates[8] = "d/dt q in component four_AP_sensitive_currents_q_gate (dimensionless)"
    legend_rates[9] = "d/dt r in component four_AP_sensitive_currents_r_gate (dimensionless)"
    legend_rates[10] = "d/dt P_af in component rapid_delayed_rectifying_potassium_current_P_af_gate (dimensionless)"
    legend_rates[11] = "d/dt P_as in component rapid_delayed_rectifying_potassium_current_P_as_gate (dimensionless)"
    legend_rates[12] = "d/dt P_i in component rapid_delayed_rectifying_potassium_current_P_i_gate (dimensionless)"
    legend_rates[13] = "d/dt xs in component slow_delayed_rectifying_potassium_current_xs_gate (dimensionless)"
    legend_rates[14] = "d/dt y in component hyperpolarisation_activated_current_y_gate (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1
    constants[1] = 0
    constants[2] = 1.0309347
    states[0] = -39.013558536
    constants[3] = 8314
    constants[4] = 310
    constants[5] = 96845
    constants[6] = 2e-5
    constants[7] = 6.5e-5
    constants[8] = 0
    constants[9] = 0
    constants[10] = 0
    constants[11] = 1.2e-6
    constants[12] = 1.204e-6
    constants[13] = 3.7e-7
    constants[14] = 140
    states[1] = 0.092361701692
    states[2] = 0.015905380261
    states[3] = 0.01445216109
    constants[15] = 0.0058
    constants[16] = 0.0057938
    constants[17] = 0.0082
    constants[18] = 0.0659
    constants[19] = 0.06588648
    constants[20] = 0.0659
    constants[21] = 46.4
    states[4] = 0.04804900895
    states[5] = 0.48779845203
    constants[22] = 0.0043
    constants[23] = 0.00427806
    constants[24] = 0.0021
    constants[25] = 0.0139
    constants[26] = 0.0138823
    constants[27] = 0.00694
    constants[28] = 45
    states[6] = 0.42074047435
    states[7] = 0.038968420558
    constants[29] = 0.00491
    constants[30] = 0.004905
    constants[31] = 0.004905
    constants[32] = 0.03649
    constants[33] = 0.036495
    constants[34] = 0.0365
    constants[35] = 6.65e-5
    constants[36] = 6.645504e-5
    constants[37] = 0.000266
    constants[38] = 0.0114
    constants[39] = 0.01138376
    constants[40] = 0.0114
    states[8] = 0.29760539675
    states[9] = 0.064402950262
    constants[41] = 0.000797
    constants[42] = 0.00079704
    constants[43] = 0.000738
    constants[44] = 0.016
    constants[45] = 0.016
    constants[46] = 0.0208
    states[10] = 0.13034201158
    states[11] = 0.46960956028
    states[12] = 0.87993375273
    constants[47] = 0.000518
    constants[48] = 0.0003445
    constants[49] = 0.000345
    constants[50] = 0.0104
    constants[51] = 0.0104
    constants[52] = 0.0104
    states[13] = 0.082293827208
    constants[53] = 0.000548
    constants[54] = 0.0005465
    constants[55] = 0.000437
    constants[56] = 0.0069
    constants[57] = 0.006875
    constants[58] = 0.0055
    constants[59] = 0.000548
    constants[60] = 0.0005465
    constants[61] = 0.000437
    constants[62] = 0.0069
    constants[63] = 0.006875
    constants[64] = 0.0055
    states[14] = 0.03889291759
    constants[65] = 5.8e-5
    constants[66] = 5.81818e-5
    constants[67] = 5.8e-5
    constants[68] = 0.000189
    constants[69] = 0.0001888
    constants[70] = 0.000189
    constants[71] = 2.52e-5
    constants[72] = 2.523636e-5
    constants[73] = 2.52e-5
    constants[74] = 8.19e-5
    constants[75] = 8.1892e-5
    constants[76] = 8.19e-5
    constants[77] = 1.32e-5
    constants[78] = 1.3236e-5
    constants[79] = 1.323e-5
    constants[80] = 4.3e-5
    constants[81] = 4.2952e-5
    constants[82] = 4.29e-5
    constants[83] = 2.7e-6
    constants[84] = 2.7229e-6
    constants[85] = 2.8e-6
    constants[86] = 8.8e-6
    constants[87] = 8.83584e-6
    constants[88] = 8.8e-6
    constants[89] = 0.0001
    constants[90] = 0.5
    constants[91] = 8
    constants[92] = 0.0001
    constants[93] = 2
    constants[94] = 5.64
    constants[95] = 0.621
    constants[96] = 0.0478
    constants[97] = 0.04782545
    constants[98] = 0.0478
    constants[99] = 0.16
    constants[100] = 0.1551936
    constants[101] = 0.16
    constants[102] = 5.4
    constants[103] = 0
    constants[104] = 0
    constants[105] = 0.0042
    constants[106] = 0
    constants[107] = 0
    constants[108] = 0.03339
    constants[109] = 140
    constants[110] = custom_piecewise([equal(constants[0] , 0.00000), (1.07000*(3.00000*constants[1]-0.100000))/(3.00000*(1.00000+0.774500*exp(-(3.00000*constants[1]-2.05000)/0.295000))) , equal(constants[0] , 1.00000), (constants[2]*constants[1])/(1.00000+0.774500*exp(-(3.00000*constants[1]-2.05000)/0.295000)) , True, (1.07000*29.0000*constants[1])/(30.0000*(1.00000+0.774500*exp(-(29.0000*constants[1]-24.5000)/1.95000)))])
    constants[111] = ((constants[3]*constants[4])/constants[5])*log(constants[14]/constants[91])
    constants[112] = ((constants[3]*constants[4])/constants[5])*log(constants[102]/constants[109])
    constants[113] = custom_piecewise([equal(constants[0] , 0.00000), 0.00200000 , equal(constants[0] , 1.00000), 0.00200000 , True, 0.00600000])
    constants[114] = custom_piecewise([equal(constants[0] , 0.00000), ((constants[3]*constants[4])/constants[5])*log((constants[102]+0.120000*constants[14])/(constants[109]+0.120000*constants[91])) , True, ((constants[3]*constants[4])/constants[5])*log((constants[102]+0.0300000*constants[14])/(constants[109]+0.0300000*constants[91]))])
    constants[115] = ((constants[3]*constants[4])/(2.00000*constants[5]))*log(constants[93]/constants[92])
    constants[116] = constants[6]+constants[110]*(constants[7]-constants[6])
    constants[117] = custom_piecewise([equal(constants[0] , 0.00000), constants[8]+constants[110]*(constants[11]-constants[8]) , equal(constants[0] , 1.00000), constants[9]+constants[110]*(constants[12]-constants[9]) , True, constants[10]+constants[110]*(constants[13]-constants[10])])
    constants[118] = custom_piecewise([equal(constants[0] , 0.00000), constants[15]+constants[110]*(constants[18]-constants[15]) , equal(constants[0] , 1.00000), constants[16]+constants[110]*(constants[19]-constants[16]) , True, constants[17]+constants[110]*(constants[20]-constants[17])])
    constants[119] = custom_piecewise([equal(constants[0] , 0.00000), constants[22]+constants[110]*(constants[25]-constants[22]) , equal(constants[0] , 1.00000), constants[23]+constants[110]*(constants[26]-constants[23]) , True, constants[24]+constants[110]*(constants[27]-constants[24])])
    constants[120] = custom_piecewise([equal(constants[0] , 0.00000), constants[29]+constants[110]*(constants[32]-constants[29]) , equal(constants[0] , 1.00000), constants[30]+constants[110]*(constants[33]-constants[30]) , True, constants[31]+constants[110]*(constants[34]-constants[31])])
    constants[121] = custom_piecewise([equal(constants[0] , 0.00000), constants[35]+constants[110]*(constants[38]-constants[35]) , equal(constants[0] , 1.00000), constants[36]+constants[110]*(constants[39]-constants[36]) , True, constants[37]+constants[110]*(constants[40]-constants[37])])
    constants[122] = custom_piecewise([equal(constants[0] , 0.00000), constants[41]+constants[110]*(constants[44]-constants[41]) , equal(constants[0] , 1.00000), constants[42]+constants[110]*(constants[45]-constants[42]) , True, constants[43]+constants[110]*(constants[46]-constants[43])])
    constants[123] = custom_piecewise([equal(constants[0] , 0.00000), constants[47]+constants[110]*(constants[50]-constants[47]) , equal(constants[0] , 1.00000), constants[48]+constants[110]*(constants[51]-constants[48]) , True, constants[49]+constants[110]*(constants[52]-constants[49])])
    constants[124] = custom_piecewise([equal(constants[0] , 0.00000), constants[53]+constants[110]*(constants[56]-constants[53]) , equal(constants[0] , 1.00000), constants[54]+constants[110]*(constants[57]-constants[54]) , True, constants[55]+constants[110]*(constants[58]-constants[55])])
    constants[125] = custom_piecewise([equal(constants[0] , 0.00000), constants[59]+constants[110]*(constants[62]-constants[59]) , equal(constants[0] , 1.00000), constants[60]+constants[110]*(constants[63]-constants[60]) , True, constants[61]+constants[110]*(constants[64]-constants[61])])
    constants[126] = custom_piecewise([equal(constants[0] , 0.00000), constants[65]+constants[110]*(constants[68]-constants[65]) , equal(constants[0] , 1.00000), constants[66]+constants[110]*(constants[69]-constants[66]) , True, constants[67]+constants[110]*(constants[70]-constants[67])])
    constants[127] = custom_piecewise([equal(constants[0] , 0.00000), constants[71]+constants[110]*(constants[74]-constants[71]) , equal(constants[0] , 1.00000), constants[72]+constants[110]*(constants[75]-constants[72]) , True, constants[73]+constants[110]*(constants[76]-constants[73])])
    constants[128] = custom_piecewise([equal(constants[0] , 0.00000), constants[77]+constants[110]*(constants[80]-constants[77]) , equal(constants[0] , 1.00000), constants[78]+constants[110]*(constants[81]-constants[78]) , True, constants[79]+constants[110]*(constants[82]-constants[79])])
    constants[129] = custom_piecewise([equal(constants[0] , 0.00000), constants[83]+constants[110]*(constants[86]-constants[83]) , equal(constants[0] , 1.00000), constants[84]+constants[110]*(constants[87]-constants[84]) , True, constants[85]+constants[110]*(constants[88]-constants[85])])
    constants[130] = custom_piecewise([equal(constants[0] , 0.00000), constants[96]+constants[110]*(constants[99]-constants[96]) , equal(constants[0] , 1.00000), constants[97]+constants[110]*(constants[100]-constants[97]) , True, constants[98]+constants[110]*(constants[101]-constants[98])])
    constants[131] = custom_piecewise([equal(constants[0] , 0.00000), constants[103]+constants[110]*(constants[106]-constants[103]) , equal(constants[0] , 1.00000), constants[104]+constants[110]*(constants[107]-constants[104]) , True, constants[105]+constants[110]*(constants[108]-constants[105])])
    constants[132] = (constants[131]*constants[92])/(constants[92]+0.000400000)
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[10] = 1.00000/(1.00000+exp((states[0]+18.6000)/10.1000))
    rates[12] = (algebraic[10]-states[12])/constants[113]
    algebraic[1] = custom_piecewise([equal(constants[0] , 0.00000), power(1.00000/(1.00000+exp(-states[0]/5.46000)), 1.00000/3.00000) , True, power(1.00000/(1.00000+exp(-(states[0]+30.3200)/5.46000)), 1.00000/3.00000)])
    algebraic[14] = custom_piecewise([equal(constants[0] , 0.00000), 0.000624700/(0.832000*exp(-0.335000*(states[0]+56.7000))+0.627000*exp(0.0820000*(states[0]+65.0100)))+4.00000e-05 , True, 0.000624700/(0.832217*exp(-0.335660*(states[0]+56.7062))+0.627400*exp(0.0823000*(states[0]+65.0131)))+4.56900e-05])
    rates[1] = (algebraic[1]-states[1])/algebraic[14]
    algebraic[2] = 1.00000/(1.00000+exp((states[0]+66.1000)/6.40000))
    algebraic[16] = (3.71700e-06*exp(-0.281500*(states[0]+17.1100)))/(1.00000+0.00373200*exp(-0.342600*(states[0]+37.7600)))+0.000597700
    rates[2] = (algebraic[2]-states[2])/algebraic[16]
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+59.3700)/13.1000))
    algebraic[21] = custom_piecewise([equal(constants[0] , 0.00000), 0.0101000+0.0651700/(0.570000*exp(-0.0800000*(states[0]+49.0000)))+2.40000e-05*exp(0.100000*(states[0]+50.9300)) , equal(constants[0] , 1.00000), (0.00100000/3.00000)*(30.3100+195.500/(0.568600*exp(-0.0816100*(states[0]+39.0000+10.0000*constants[110]))+0.717400*exp((0.271900-0.171900*constants[110])*1.00000*(states[0]+40.9300+10.0000*constants[110])))) , True, 0.0101000+0.0651700/(0.568600*exp(-0.0816100*(states[0]+39.0000))+0.717400*exp(0.271900*(states[0]+40.9300)))])
    rates[8] = (algebraic[7]-states[8])/algebraic[21]
    algebraic[8] = 1.00000/(1.00000+exp(-(states[0]-10.9300)/19.7000))
    algebraic[22] = custom_piecewise([equal(constants[0] , 0.00000), 0.00100000*(2.98000+15.5900/(1.03700*exp(0.0900000*(states[0]+30.6100))+0.369000*exp(-0.120000*(states[0]+23.8400)))) , equal(constants[0] , 1.00000), 0.00250000*(1.19100+7.83800/(1.03700*exp(0.0901200*(states[0]+30.6100))+0.369000*exp(-0.119000*(states[0]+23.8400)))) , True, 0.00100000*(2.98000+19.5900/(1.03700*exp(0.0901200*(states[0]+30.6100))+0.369000*exp(-0.119000*(states[0]+23.8400))))])
    rates[9] = (algebraic[8]-states[9])/algebraic[22]
    algebraic[9] = custom_piecewise([constants[0] != 2.00000, 1.00000/(1.00000+exp(-(states[0]+14.2000)/10.6000)) , True, 1.00000/(1.00000+exp(-(states[0]+13.2000)/10.6000))])
    algebraic[23] = custom_piecewise([constants[0] != 2.00000, 1.00000/(37.2000*exp((states[0]-9.00000)/15.9000)+0.960000*exp(-(states[0]-9.00000)/22.5000)) , True, 1.00000/(37.2000*exp((states[0]-10.0000)/15.9000)+0.960000*exp(-(states[0]-10.0000)/22.5000))])
    rates[10] = (algebraic[9]-states[10])/algebraic[23]
    algebraic[11] = 14.0000/(1.00000+exp(-(states[0]-40.0000)/9.00000))
    algebraic[25] = 1.00000*exp(-states[0]/45.0000)
    rates[13] = algebraic[11]*(1.00000-states[13])-algebraic[25]*states[13]
    algebraic[12] = custom_piecewise([equal(constants[0] , 0.00000), 1.00000*exp(-(states[0]+78.9100)/26.6200) , True, 1.00000*exp(-(states[0]+78.9100)/26.6300)])
    algebraic[26] = 1.00000*exp((states[0]+75.1300)/21.2500)
    rates[14] = algebraic[12]*(1.00000-states[14])-algebraic[26]*states[14]
    algebraic[15] = algebraic[2]
    algebraic[28] = (3.18600e-08*exp(-0.621900*(states[0]+18.8000)))/(1.00000+7.18900e-05*exp(-0.668300*(states[0]+34.0700)))+0.00355600
    rates[3] = (algebraic[15]-states[3])/algebraic[28]
    algebraic[24] = algebraic[9]
    algebraic[33] = custom_piecewise([constants[0] != 2.00000, 1.00000/(4.20000*exp((states[0]-9.00000)/17.0000)+0.150000*exp(-(states[0]-9.00000)/21.6000)) , True, 1.00000/(4.20000*exp((states[0]-10.0000)/17.0000)+0.150000*exp(-(states[0]-10.0000)/21.6000))])
    rates[11] = (algebraic[24]-states[11])/algebraic[33]
    algebraic[35] = custom_piecewise([equal(constants[0] , 0.00000), 1.00000/(1.00000+exp(-(states[0]+23.1000)/6.00000)) , equal(constants[0] , 1.00000), 1.00000/(1.00000+exp(-(states[0]+22.3000+0.800000*constants[110])/6.00000)) , True, 1.00000/(1.00000+exp(-(states[0]+22.2000)/6.00000))])
    algebraic[3] = custom_piecewise([equal(constants[0] , 0.00000), (-28.3800*(states[0]+35.0000))/(exp(-(states[0]+35.0000)/2.50000)-1.00000)-(84.9000*states[0])/(exp(-0.208000*states[0])-1.00000) , equal(constants[0] , 1.00000), (-28.3900*(states[0]+35.0000))/(exp(-(states[0]+35.0000)/2.50000)-1.00000)-(84.9000*states[0])/(exp(-0.208000*states[0])-1.00000) , True, (-28.4000*(states[0]+35.0000))/(exp(-(states[0]+35.0000)/2.50000)-1.00000)-(84.9000*states[0])/(exp(-0.208000*states[0])-1.00000)])
    algebraic[17] = custom_piecewise([equal(constants[0] , 1.00000), (11.4300*(states[0]-5.00000))/(exp(0.400000*(states[0]-5.00000))-1.00000) , True, (11.4200*(states[0]-5.00000))/(exp(0.400000*(states[0]-5.00000))-1.00000)])
    algebraic[29] = 2.00000/(algebraic[3]+algebraic[17])
    rates[4] = (algebraic[35]-states[4])/algebraic[29]
    algebraic[36] = 1.00000/(1.00000+exp((states[0]+45.0000)/5.00000))
    algebraic[4] = custom_piecewise([equal(constants[0] , 1.00000), (3.75000*(states[0]+28.0000))/(exp((states[0]+28.0000)/4.00000)-1.00000) , True, (3.12000*(states[0]+28.0000))/(exp((states[0]+28.0000)/4.00000)-1.00000)])
    algebraic[18] = custom_piecewise([equal(constants[0] , 1.00000), 30.0000/(1.00000+exp(-(states[0]+28.0000)/4.00000)) , True, 25.0000/(1.00000+exp(-(states[0]+28.0000)/4.00000))])
    algebraic[30] = custom_piecewise([equal(constants[0] , 1.00000), (1.20000-0.200000*constants[110])/(algebraic[4]+algebraic[18]) , True, 1.00000/(algebraic[4]+algebraic[18])])
    rates[5] = (algebraic[36]-states[5])/algebraic[30]
    algebraic[37] = 1.00000/(1.00000+exp(-(states[0]+37.0000)/6.80000))
    algebraic[5] = 1068.00*exp((states[0]+26.3000)/30.0000)
    algebraic[19] = 1068.00*exp(-(states[0]+26.3000)/30.0000)
    algebraic[31] = 1.00000/(algebraic[5]+algebraic[19])
    rates[6] = (algebraic[37]-states[6])/algebraic[31]
    algebraic[38] = 1.00000/(1.00000+exp((states[0]+71.0000)/9.00000))
    algebraic[6] = custom_piecewise([equal(constants[0] , 1.00000), 15.3000*exp(-(states[0]+71.0000+0.700000*constants[110])/83.3000) , True, 15.3000*exp(-(states[0]+71.7000)/83.3000)])
    algebraic[20] = custom_piecewise([equal(constants[0] , 1.00000), 15.0000*exp((states[0]+71.0000)/15.3800) , True, 15.0000*exp((states[0]+71.7000)/15.3800)])
    algebraic[32] = 1.00000/(algebraic[6]+algebraic[20])
    rates[7] = (algebraic[38]-states[7])/algebraic[32]
    algebraic[0] = custom_piecewise([equal(constants[0] , 0.00000), (0.0952000*exp(-0.0630000*(states[0]+34.4000)))/(1.00000+1.66000*exp(-0.225000*(states[0]+63.7000)))+0.0869000 , True, (0.0951800*exp(-0.0630600*(states[0]+34.4000)))/(1.00000+1.66200*exp(-0.225100*(states[0]+63.7000)))+0.0869300])
    algebraic[13] = (1.00000-algebraic[0])*states[2]+algebraic[0]*states[3]
    algebraic[27] = ((((constants[117]*(power(states[1], 3.00000))*algebraic[13]*constants[14]*(power(constants[5], 2.00000)))/(constants[3]*constants[4]))*(exp(((states[0]-constants[111])*constants[5])/(constants[3]*constants[4]))-1.00000))/(exp((states[0]*constants[5])/(constants[3]*constants[4]))-1.00000))*states[0]
    algebraic[34] = constants[118]*(states[5]*states[4]+0.00600000/(1.00000+exp(-(states[0]+14.1000)/6.00000)))*(states[0]-constants[21])
    algebraic[39] = constants[119]*states[6]*states[7]*(states[0]-constants[28])
    algebraic[40] = constants[120]*states[8]*states[9]*(states[0]-constants[112])
    algebraic[41] = constants[121]*states[9]*(states[0]-constants[112])
    algebraic[42] = 0.600000*states[10]+0.400000*states[11]
    algebraic[43] = constants[122]*algebraic[42]*states[12]*(states[0]-constants[112])
    algebraic[44] = constants[123]*(power(states[13], 2.00000))*(states[0]-constants[114])
    algebraic[45] = custom_piecewise([constants[0] != 2.00000, constants[124]*states[14]*(states[0]-constants[111]) , True, constants[124]*states[14]*(states[0]-77.6000)])
    algebraic[46] = custom_piecewise([constants[0] != 2.00000, constants[125]*states[14]*(states[0]-constants[112]) , True, constants[125]*states[14]*(states[0]+102.000)])
    algebraic[47] = constants[126]*(states[0]-constants[111])
    algebraic[49] = constants[128]*(states[0]-constants[115])
    algebraic[48] = constants[127]*(states[0]-constants[112])
    algebraic[50] = custom_piecewise([equal(constants[0] , 0.00000), (constants[129]*((power(constants[91], 3.00000))*constants[93]*exp(0.0374300*states[0]*constants[90])-(power(constants[14], 3.00000))*constants[92]*exp(0.0374000*states[0]*(constants[90]-1.00000))))/(1.00000+constants[89]*(constants[92]*(power(constants[14], 3.00000))+constants[93]*(power(constants[91], 3.00000)))) , True, (constants[129]*((power(constants[91], 3.00000))*constants[93]*exp(0.0374300*states[0]*constants[90])-(power(constants[14], 3.00000))*constants[92]*exp(0.0374300*states[0]*(constants[90]-1.00000))))/(1.00000+constants[89]*(constants[92]*(power(constants[14], 3.00000))+constants[93]*(power(constants[91], 3.00000))))])
    algebraic[51] = (constants[130]*(power(constants[91]/(constants[94]+constants[91]), 3.00000))*(power(constants[102]/(constants[95]+constants[102]), 2.00000))*1.60000)/(1.50000+exp(-(states[0]+60.0000)/40.0000))
    rates[0] = (-1.00000/constants[116])*(algebraic[27]+algebraic[34]+algebraic[39]+algebraic[40]+algebraic[41]+algebraic[43]+algebraic[44]+algebraic[45]+algebraic[46]+algebraic[47]+algebraic[49]+algebraic[48]+algebraic[50]+algebraic[51]+constants[132])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[10] = 1.00000/(1.00000+exp((states[0]+18.6000)/10.1000))
    algebraic[1] = custom_piecewise([equal(constants[0] , 0.00000), power(1.00000/(1.00000+exp(-states[0]/5.46000)), 1.00000/3.00000) , True, power(1.00000/(1.00000+exp(-(states[0]+30.3200)/5.46000)), 1.00000/3.00000)])
    algebraic[14] = custom_piecewise([equal(constants[0] , 0.00000), 0.000624700/(0.832000*exp(-0.335000*(states[0]+56.7000))+0.627000*exp(0.0820000*(states[0]+65.0100)))+4.00000e-05 , True, 0.000624700/(0.832217*exp(-0.335660*(states[0]+56.7062))+0.627400*exp(0.0823000*(states[0]+65.0131)))+4.56900e-05])
    algebraic[2] = 1.00000/(1.00000+exp((states[0]+66.1000)/6.40000))
    algebraic[16] = (3.71700e-06*exp(-0.281500*(states[0]+17.1100)))/(1.00000+0.00373200*exp(-0.342600*(states[0]+37.7600)))+0.000597700
    algebraic[7] = 1.00000/(1.00000+exp((states[0]+59.3700)/13.1000))
    algebraic[21] = custom_piecewise([equal(constants[0] , 0.00000), 0.0101000+0.0651700/(0.570000*exp(-0.0800000*(states[0]+49.0000)))+2.40000e-05*exp(0.100000*(states[0]+50.9300)) , equal(constants[0] , 1.00000), (0.00100000/3.00000)*(30.3100+195.500/(0.568600*exp(-0.0816100*(states[0]+39.0000+10.0000*constants[110]))+0.717400*exp((0.271900-0.171900*constants[110])*1.00000*(states[0]+40.9300+10.0000*constants[110])))) , True, 0.0101000+0.0651700/(0.568600*exp(-0.0816100*(states[0]+39.0000))+0.717400*exp(0.271900*(states[0]+40.9300)))])
    algebraic[8] = 1.00000/(1.00000+exp(-(states[0]-10.9300)/19.7000))
    algebraic[22] = custom_piecewise([equal(constants[0] , 0.00000), 0.00100000*(2.98000+15.5900/(1.03700*exp(0.0900000*(states[0]+30.6100))+0.369000*exp(-0.120000*(states[0]+23.8400)))) , equal(constants[0] , 1.00000), 0.00250000*(1.19100+7.83800/(1.03700*exp(0.0901200*(states[0]+30.6100))+0.369000*exp(-0.119000*(states[0]+23.8400)))) , True, 0.00100000*(2.98000+19.5900/(1.03700*exp(0.0901200*(states[0]+30.6100))+0.369000*exp(-0.119000*(states[0]+23.8400))))])
    algebraic[9] = custom_piecewise([constants[0] != 2.00000, 1.00000/(1.00000+exp(-(states[0]+14.2000)/10.6000)) , True, 1.00000/(1.00000+exp(-(states[0]+13.2000)/10.6000))])
    algebraic[23] = custom_piecewise([constants[0] != 2.00000, 1.00000/(37.2000*exp((states[0]-9.00000)/15.9000)+0.960000*exp(-(states[0]-9.00000)/22.5000)) , True, 1.00000/(37.2000*exp((states[0]-10.0000)/15.9000)+0.960000*exp(-(states[0]-10.0000)/22.5000))])
    algebraic[11] = 14.0000/(1.00000+exp(-(states[0]-40.0000)/9.00000))
    algebraic[25] = 1.00000*exp(-states[0]/45.0000)
    algebraic[12] = custom_piecewise([equal(constants[0] , 0.00000), 1.00000*exp(-(states[0]+78.9100)/26.6200) , True, 1.00000*exp(-(states[0]+78.9100)/26.6300)])
    algebraic[26] = 1.00000*exp((states[0]+75.1300)/21.2500)
    algebraic[15] = algebraic[2]
    algebraic[28] = (3.18600e-08*exp(-0.621900*(states[0]+18.8000)))/(1.00000+7.18900e-05*exp(-0.668300*(states[0]+34.0700)))+0.00355600
    algebraic[24] = algebraic[9]
    algebraic[33] = custom_piecewise([constants[0] != 2.00000, 1.00000/(4.20000*exp((states[0]-9.00000)/17.0000)+0.150000*exp(-(states[0]-9.00000)/21.6000)) , True, 1.00000/(4.20000*exp((states[0]-10.0000)/17.0000)+0.150000*exp(-(states[0]-10.0000)/21.6000))])
    algebraic[35] = custom_piecewise([equal(constants[0] , 0.00000), 1.00000/(1.00000+exp(-(states[0]+23.1000)/6.00000)) , equal(constants[0] , 1.00000), 1.00000/(1.00000+exp(-(states[0]+22.3000+0.800000*constants[110])/6.00000)) , True, 1.00000/(1.00000+exp(-(states[0]+22.2000)/6.00000))])
    algebraic[3] = custom_piecewise([equal(constants[0] , 0.00000), (-28.3800*(states[0]+35.0000))/(exp(-(states[0]+35.0000)/2.50000)-1.00000)-(84.9000*states[0])/(exp(-0.208000*states[0])-1.00000) , equal(constants[0] , 1.00000), (-28.3900*(states[0]+35.0000))/(exp(-(states[0]+35.0000)/2.50000)-1.00000)-(84.9000*states[0])/(exp(-0.208000*states[0])-1.00000) , True, (-28.4000*(states[0]+35.0000))/(exp(-(states[0]+35.0000)/2.50000)-1.00000)-(84.9000*states[0])/(exp(-0.208000*states[0])-1.00000)])
    algebraic[17] = custom_piecewise([equal(constants[0] , 1.00000), (11.4300*(states[0]-5.00000))/(exp(0.400000*(states[0]-5.00000))-1.00000) , True, (11.4200*(states[0]-5.00000))/(exp(0.400000*(states[0]-5.00000))-1.00000)])
    algebraic[29] = 2.00000/(algebraic[3]+algebraic[17])
    algebraic[36] = 1.00000/(1.00000+exp((states[0]+45.0000)/5.00000))
    algebraic[4] = custom_piecewise([equal(constants[0] , 1.00000), (3.75000*(states[0]+28.0000))/(exp((states[0]+28.0000)/4.00000)-1.00000) , True, (3.12000*(states[0]+28.0000))/(exp((states[0]+28.0000)/4.00000)-1.00000)])
    algebraic[18] = custom_piecewise([equal(constants[0] , 1.00000), 30.0000/(1.00000+exp(-(states[0]+28.0000)/4.00000)) , True, 25.0000/(1.00000+exp(-(states[0]+28.0000)/4.00000))])
    algebraic[30] = custom_piecewise([equal(constants[0] , 1.00000), (1.20000-0.200000*constants[110])/(algebraic[4]+algebraic[18]) , True, 1.00000/(algebraic[4]+algebraic[18])])
    algebraic[37] = 1.00000/(1.00000+exp(-(states[0]+37.0000)/6.80000))
    algebraic[5] = 1068.00*exp((states[0]+26.3000)/30.0000)
    algebraic[19] = 1068.00*exp(-(states[0]+26.3000)/30.0000)
    algebraic[31] = 1.00000/(algebraic[5]+algebraic[19])
    algebraic[38] = 1.00000/(1.00000+exp((states[0]+71.0000)/9.00000))
    algebraic[6] = custom_piecewise([equal(constants[0] , 1.00000), 15.3000*exp(-(states[0]+71.0000+0.700000*constants[110])/83.3000) , True, 15.3000*exp(-(states[0]+71.7000)/83.3000)])
    algebraic[20] = custom_piecewise([equal(constants[0] , 1.00000), 15.0000*exp((states[0]+71.0000)/15.3800) , True, 15.0000*exp((states[0]+71.7000)/15.3800)])
    algebraic[32] = 1.00000/(algebraic[6]+algebraic[20])
    algebraic[0] = custom_piecewise([equal(constants[0] , 0.00000), (0.0952000*exp(-0.0630000*(states[0]+34.4000)))/(1.00000+1.66000*exp(-0.225000*(states[0]+63.7000)))+0.0869000 , True, (0.0951800*exp(-0.0630600*(states[0]+34.4000)))/(1.00000+1.66200*exp(-0.225100*(states[0]+63.7000)))+0.0869300])
    algebraic[13] = (1.00000-algebraic[0])*states[2]+algebraic[0]*states[3]
    algebraic[27] = ((((constants[117]*(power(states[1], 3.00000))*algebraic[13]*constants[14]*(power(constants[5], 2.00000)))/(constants[3]*constants[4]))*(exp(((states[0]-constants[111])*constants[5])/(constants[3]*constants[4]))-1.00000))/(exp((states[0]*constants[5])/(constants[3]*constants[4]))-1.00000))*states[0]
    algebraic[34] = constants[118]*(states[5]*states[4]+0.00600000/(1.00000+exp(-(states[0]+14.1000)/6.00000)))*(states[0]-constants[21])
    algebraic[39] = constants[119]*states[6]*states[7]*(states[0]-constants[28])
    algebraic[40] = constants[120]*states[8]*states[9]*(states[0]-constants[112])
    algebraic[41] = constants[121]*states[9]*(states[0]-constants[112])
    algebraic[42] = 0.600000*states[10]+0.400000*states[11]
    algebraic[43] = constants[122]*algebraic[42]*states[12]*(states[0]-constants[112])
    algebraic[44] = constants[123]*(power(states[13], 2.00000))*(states[0]-constants[114])
    algebraic[45] = custom_piecewise([constants[0] != 2.00000, constants[124]*states[14]*(states[0]-constants[111]) , True, constants[124]*states[14]*(states[0]-77.6000)])
    algebraic[46] = custom_piecewise([constants[0] != 2.00000, constants[125]*states[14]*(states[0]-constants[112]) , True, constants[125]*states[14]*(states[0]+102.000)])
    algebraic[47] = constants[126]*(states[0]-constants[111])
    algebraic[49] = constants[128]*(states[0]-constants[115])
    algebraic[48] = constants[127]*(states[0]-constants[112])
    algebraic[50] = custom_piecewise([equal(constants[0] , 0.00000), (constants[129]*((power(constants[91], 3.00000))*constants[93]*exp(0.0374300*states[0]*constants[90])-(power(constants[14], 3.00000))*constants[92]*exp(0.0374000*states[0]*(constants[90]-1.00000))))/(1.00000+constants[89]*(constants[92]*(power(constants[14], 3.00000))+constants[93]*(power(constants[91], 3.00000)))) , True, (constants[129]*((power(constants[91], 3.00000))*constants[93]*exp(0.0374300*states[0]*constants[90])-(power(constants[14], 3.00000))*constants[92]*exp(0.0374300*states[0]*(constants[90]-1.00000))))/(1.00000+constants[89]*(constants[92]*(power(constants[14], 3.00000))+constants[93]*(power(constants[91], 3.00000))))])
    algebraic[51] = (constants[130]*(power(constants[91]/(constants[94]+constants[91]), 3.00000))*(power(constants[102]/(constants[95]+constants[102]), 2.00000))*1.60000)/(1.50000+exp(-(states[0]+60.0000)/40.0000))
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