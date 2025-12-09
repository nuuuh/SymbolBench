# Size of variable arrays:
sizeAlgebraic = 42
sizeStates = 27
sizeConstants = 64
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component Environment (second)"
    legend_constants[0] = "O2_x in component Environment (molar)"
    legend_constants[1] = "X_AtC in component Environment (mole_per_second_per_l_cell)"
    legend_constants[2] = "RT in component fixed_parameters (kilojoule_per_mole)"
    legend_constants[3] = "F in component fixed_parameters (kilojoule_per_mole_per_millivolt)"
    legend_constants[4] = "dG_C1o in component fixed_parameters (kilojoule_per_mole)"
    legend_constants[5] = "dG_C3o in component fixed_parameters (kilojoule_per_mole)"
    legend_constants[6] = "dG_C4o in component fixed_parameters (kilojoule_per_mole)"
    legend_algebraic[13] = "dG_C1op in component fixed_parameters (kilojoule_per_mole)"
    legend_algebraic[14] = "dG_C3op in component fixed_parameters (kilojoule_per_mole)"
    legend_algebraic[15] = "dG_C4op in component fixed_parameters (kilojoule_per_mole)"
    legend_constants[7] = "dG_F1o in component fixed_parameters (kilojoule_per_mole)"
    legend_algebraic[11] = "dG_H in component fixed_parameters (kilojoule_per_mole)"
    legend_constants[50] = "k_dHPI in component fixed_parameters (molar)"
    legend_constants[51] = "k_dHatp in component fixed_parameters (molar)"
    legend_constants[52] = "k_dHadp in component fixed_parameters (molar)"
    legend_constants[8] = "K_DT in component fixed_parameters (molar)"
    legend_constants[9] = "K_DD in component fixed_parameters (molar)"
    legend_constants[10] = "K_AK in component fixed_parameters (dimensionless)"
    legend_constants[11] = "K_CK in component fixed_parameters (per_molar)"
    legend_constants[56] = "W_m in component fixed_parameters (l_water_per_l_mito)"
    legend_constants[59] = "W_c in component fixed_parameters (l_water_per_l_cyto)"
    legend_constants[57] = "W_x in component fixed_parameters (l_water_per_l_mito)"
    legend_constants[58] = "W_i in component fixed_parameters (l_water_per_l_mito)"
    legend_constants[12] = "V_cyto in component fixed_parameters (cyto_per_cell)"
    legend_constants[13] = "V_mito in component fixed_parameters (mito_per_cell)"
    legend_constants[61] = "Rm_cyto in component fixed_parameters (mito_per_cyto)"
    legend_constants[54] = "Rm_cell in component fixed_parameters (mito_per_cell)"
    legend_constants[62] = "Rc_cell in component fixed_parameters (cyto_per_cell)"
    legend_constants[14] = "n_A in component fixed_parameters (dimensionless)"
    legend_constants[15] = "C_tot in component fixed_parameters (molar)"
    legend_constants[16] = "CR_tot in component fixed_parameters (molar)"
    legend_constants[17] = "Q_tot in component fixed_parameters (molar)"
    legend_constants[18] = "NAD_tot in component fixed_parameters (molar)"
    legend_constants[19] = "pH_C in component fixed_parameters (dimensionless)"
    legend_constants[53] = "H_c in component fixed_parameters (molar)"
    legend_constants[20] = "K_c in component fixed_parameters (molar)"
    legend_constants[60] = "H_i in component fixed_parameters (molar)"
    legend_constants[55] = "K_i in component fixed_parameters (molar)"
    legend_algebraic[9] = "H_x in component dH_x_dt (molar)"
    legend_states[0] = "deltaPsi in component ddeltaPsi_dt (millivolt)"
    legend_constants[21] = "k_PI1 in component adjustable_parameters (molar)"
    legend_constants[22] = "k_PI2 in component adjustable_parameters (molar)"
    legend_constants[23] = "k_PI3 in component adjustable_parameters (molar)"
    legend_constants[24] = "k_PI4 in component adjustable_parameters (molar)"
    legend_constants[25] = "k_PIH in component adjustable_parameters (molar)"
    legend_constants[26] = "r in component adjustable_parameters (dimensionless)"
    legend_constants[27] = "X_DH in component adjustable_parameters (mole_per_second_per_l_mito_per_molar)"
    legend_constants[28] = "X_C1 in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_molar)"
    legend_constants[29] = "X_C3 in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_half_molar)"
    legend_constants[30] = "X_C4 in component adjustable_parameters (mole_per_second_per_l_mito_per_molar)"
    legend_constants[31] = "X_F1 in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_molar)"
    legend_constants[32] = "X_ANT in component adjustable_parameters (mole_per_second_per_l_mito)"
    legend_constants[33] = "X_PI1 in component adjustable_parameters (mole_per_second_per_l_mito_per_molar)"
    legend_constants[34] = "X_KH in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_molar)"
    legend_constants[35] = "X_Hle in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_millivolt)"
    legend_constants[36] = "X_K in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_millivolt)"
    legend_constants[37] = "k_mADP in component adjustable_parameters (molar)"
    legend_constants[38] = "X_AK in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_molar)"
    legend_constants[39] = "X_CK in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_molar)"
    legend_constants[40] = "X_MgA in component adjustable_parameters (mole_per_second_per_l_mito_per_molar_per_molar)"
    legend_constants[41] = "k_O2 in component adjustable_parameters (molar)"
    legend_constants[42] = "x_buff in component adjustable_parameters (per_molar)"
    legend_constants[43] = "C_IM in component adjustable_parameters (mole_per_l_mito_per_millivolt)"
    legend_constants[44] = "X_A in component adjustable_parameters (micron_per_second)"
    legend_constants[45] = "X_PI2 in component adjustable_parameters (micron_per_second)"
    legend_constants[46] = "gamma in component adjustable_parameters (per_micron)"
    legend_algebraic[0] = "NAD_x in component dNAD_x_dt (molar)"
    legend_states[1] = "NADH_x in component dNADH_x_dt (molar)"
    legend_algebraic[2] = "QH2_x in component dQH2_x_dt (molar)"
    legend_states[2] = "Q_x in component dQ_x_dt (molar)"
    legend_algebraic[4] = "Cred_x in component dCred_x_dt (molar)"
    legend_states[3] = "Cox_x in component dCox_x_dt (molar)"
    legend_algebraic[6] = "ATP_fx in component dATP_fx_dt (molar)"
    legend_states[4] = "ATP_x in component dATP_x_dt (molar)"
    legend_states[5] = "ATP_mx in component dATP_mx_dt (molar)"
    legend_algebraic[8] = "ADP_fx in component dADP_fx_dt (molar)"
    legend_states[6] = "ADP_x in component dADP_x_dt (molar)"
    legend_states[7] = "ADP_mx in component dADP_mx_dt (molar)"
    legend_states[8] = "pH_x in component dH_x_dt (dimensionless)"
    legend_algebraic[20] = "J_DH in component J_DH (mole_per_second_per_l_mito)"
    legend_algebraic[22] = "J_C1 in component J_C1 (mole_per_second_per_l_mito)"
    legend_algebraic[23] = "J_C3 in component J_C3 (mole_per_second_per_l_mito)"
    legend_algebraic[24] = "J_C4 in component J_C4 (mole_per_second_per_l_mito)"
    legend_algebraic[25] = "J_F1 in component J_F1 (mole_per_second_per_l_mito)"
    legend_algebraic[33] = "J_PI1 in component J_PI1 (mole_per_second_per_l_mito)"
    legend_algebraic[39] = "J_KH in component J_KH (mole_per_second_per_l_mito)"
    legend_algebraic[35] = "J_Hle in component J_Hle (mole_per_second_per_l_mito)"
    legend_states[9] = "K_x in component dK_x_dt (molar)"
    legend_algebraic[41] = "J_K in component J_K (mole_per_second_per_l_mito)"
    legend_states[10] = "Mg_x in component dMg_x_dt (molar)"
    legend_algebraic[10] = "J_MgATP_x in component J_MgATP_x (mole_per_second_per_l_mito)"
    legend_algebraic[12] = "J_MgADP_x in component J_MgADP_x (mole_per_second_per_l_mito)"
    legend_algebraic[28] = "J_ANT in component J_ANT (mole_per_second_per_l_mito)"
    legend_states[11] = "PI_x in component dPI_x_dt (molar)"
    legend_algebraic[16] = "ATP_fi in component dATP_fi_dt (molar)"
    legend_states[12] = "ATP_i in component dATP_i_dt (molar)"
    legend_states[13] = "ATP_mi in component dATP_mi_dt (molar)"
    legend_algebraic[17] = "ADP_fi in component dADP_fi_dt (molar)"
    legend_states[14] = "ADP_i in component dADP_i_dt (molar)"
    legend_states[15] = "ADP_mi in component dADP_mi_dt (molar)"
    legend_algebraic[34] = "J_ATP in component J_ATP (mole_per_second_per_l_mito)"
    legend_algebraic[36] = "J_AKi in component J_AKi (mole_per_second_per_l_mito)"
    legend_algebraic[32] = "J_ADP in component J_ADP (mole_per_second_per_l_mito)"
    legend_states[16] = "AMP_i in component dAMP_i_dt (molar)"
    legend_algebraic[30] = "J_AMP in component J_AMP (mole_per_second_per_l_mito)"
    legend_algebraic[19] = "J_MgATP_i in component J_MgATP_i (mole_per_second_per_l_mito)"
    legend_algebraic[21] = "J_MgADP_i in component J_MgADP_i (mole_per_second_per_l_mito)"
    legend_states[17] = "PI_i in component dPI_i_dt (molar)"
    legend_algebraic[37] = "J_PI2 in component J_PI2 (mole_per_second_per_l_mito)"
    legend_states[18] = "Mg_i in component dMg_i_dt (molar)"
    legend_states[19] = "ATP_c in component dATP_c_dt (molar)"
    legend_algebraic[38] = "J_AKc in component J_AKc (mole_per_second_per_l_cyto)"
    legend_constants[63] = "J_AtC in component J_AtC (mole_per_second_per_l_cell)"
    legend_algebraic[40] = "J_CKc in component J_CKc (mole_per_second_per_l_cyto)"
    legend_states[20] = "ADP_c in component dADP_c_dt (molar)"
    legend_states[21] = "AMP_c in component dAMP_c_dt (molar)"
    legend_states[22] = "ATP_mc in component dATP_mc_dt (molar)"
    legend_algebraic[5] = "J_MgATP_c in component J_MgATP_c (mole_per_second_per_l_cyto)"
    legend_states[23] = "ADP_mc in component dADP_mc_dt (molar)"
    legend_algebraic[7] = "J_MgADP_c in component J_MgADP_c (mole_per_second_per_l_cyto)"
    legend_states[24] = "PI_c in component dPI_c_dt (molar)"
    legend_states[25] = "Mg_c in component dMg_c_dt (molar)"
    legend_states[26] = "PCr_c in component dPCr_c_dt (molar)"
    legend_algebraic[1] = "ATP_fc in component dATP_fc_dt (molar)"
    legend_algebraic[3] = "ADP_fc in component dADP_fc_dt (molar)"
    legend_algebraic[18] = "Cr_c in component dCr_c_dt (molar)"
    legend_constants[47] = "mincon in component J_ANT (molar)"
    legend_algebraic[26] = "Psi_x in component J_ANT (millivolt)"
    legend_algebraic[27] = "Psi_i in component J_ANT (millivolt)"
    legend_algebraic[29] = "H2PIi in component J_PI1 (molar)"
    legend_algebraic[31] = "H2PIx in component J_PI1 (molar)"
    legend_constants[48] = "mincond in component J_Hle (millivolt)"
    legend_constants[49] = "mincond in component J_K (millivolt)"
    legend_rates[8] = "d/dt pH_x in component dH_x_dt (dimensionless)"
    legend_rates[9] = "d/dt K_x in component dK_x_dt (molar)"
    legend_rates[10] = "d/dt Mg_x in component dMg_x_dt (molar)"
    legend_rates[1] = "d/dt NADH_x in component dNADH_x_dt (molar)"
    legend_rates[2] = "d/dt Q_x in component dQ_x_dt (molar)"
    legend_rates[3] = "d/dt Cox_x in component dCox_x_dt (molar)"
    legend_rates[4] = "d/dt ATP_x in component dATP_x_dt (molar)"
    legend_rates[6] = "d/dt ADP_x in component dADP_x_dt (molar)"
    legend_rates[5] = "d/dt ATP_mx in component dATP_mx_dt (molar)"
    legend_rates[7] = "d/dt ADP_mx in component dADP_mx_dt (molar)"
    legend_rates[11] = "d/dt PI_x in component dPI_x_dt (molar)"
    legend_rates[12] = "d/dt ATP_i in component dATP_i_dt (molar)"
    legend_rates[14] = "d/dt ADP_i in component dADP_i_dt (molar)"
    legend_rates[16] = "d/dt AMP_i in component dAMP_i_dt (molar)"
    legend_rates[13] = "d/dt ATP_mi in component dATP_mi_dt (molar)"
    legend_rates[15] = "d/dt ADP_mi in component dADP_mi_dt (molar)"
    legend_rates[17] = "d/dt PI_i in component dPI_i_dt (molar)"
    legend_rates[18] = "d/dt Mg_i in component dMg_i_dt (molar)"
    legend_rates[19] = "d/dt ATP_c in component dATP_c_dt (molar)"
    legend_rates[20] = "d/dt ADP_c in component dADP_c_dt (molar)"
    legend_rates[21] = "d/dt AMP_c in component dAMP_c_dt (molar)"
    legend_rates[22] = "d/dt ATP_mc in component dATP_mc_dt (molar)"
    legend_rates[23] = "d/dt ADP_mc in component dADP_mc_dt (molar)"
    legend_rates[24] = "d/dt PI_c in component dPI_c_dt (molar)"
    legend_rates[25] = "d/dt Mg_c in component dMg_c_dt (molar)"
    legend_rates[26] = "d/dt PCr_c in component dPCr_c_dt (molar)"
    legend_rates[0] = "d/dt deltaPsi in component ddeltaPsi_dt (millivolt)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 3.48e-5
    constants[1] = 0e-4
    constants[2] = 2.5775
    constants[3] = 0.096484
    constants[4] = -69.37
    constants[5] = -32.53
    constants[6] = -122.94
    constants[7] = 36.03
    constants[8] = 24e-6
    constants[9] = 347e-6
    constants[10] = 0.4331
    constants[11] = 1.66e9
    constants[12] = 0.894
    constants[13] = 0.056
    constants[14] = 3
    constants[15] = 2.7e-3
    constants[16] = 42.7e-3
    constants[17] = 1.35e-3
    constants[18] = 2.97e-3
    constants[19] = 7
    constants[20] = 0.15
    states[0] = 166.67
    constants[21] = 0.1553e-3
    constants[22] = 0.8222e-3
    constants[23] = 0.3601e-3
    constants[24] = 5.651e-3
    constants[25] = 2.542e-4
    constants[26] = 4.559
    constants[27] = 0.0866
    constants[28] = 4.405e3
    constants[29] = 4.887
    constants[30] = 6.766e-5
    constants[31] = 1e3
    constants[32] = 8.123e-3
    constants[33] = 3.85e5
    constants[34] = 5.651e7
    constants[35] = 200
    constants[36] = 0
    constants[37] = 3.5e-6
    constants[38] = 1e7
    constants[39] = 1e7
    constants[40] = 1e7
    constants[41] = 1.2e-4
    constants[42] = 100
    constants[43] = 6.75e-6
    constants[44] = 85
    constants[45] = 327
    constants[46] = 5.99
    states[1] = 0.0015723
    states[2] = 6.75e-4
    states[3] = 1.35e-3
    states[4] = 0.0026657
    states[5] = 0.0026046
    states[6] = 0.0073343
    states[7] = 0.0054765
    states[8] = 7
    states[9] = 0.14661
    states[10] = 1.0229e-3
    states[11] = 1.72e-4
    states[12] = 0.0065339
    states[13] = 0.0063812
    states[14] = 6.5773e-5
    states[15] = 4.8866e-5
    states[16] = 2.8837e-7
    states[17] = 1.72e-4
    states[18] = 0.0010029
    states[19] = 0.0081312
    states[20] = 6.85e-5
    states[21] = 3.0911e-7
    states[22] = 0.0079786
    states[23] = 5.1958e-5
    states[24] = 8.7702e-3
    states[25] = 0.001003
    states[26] = 23.589e-3
    constants[47] = 1e-9
    constants[48] = 1e-6
    constants[49] = 1e-6
    constants[50] = 1.00000*(power(10.0000, -6.75000))
    constants[51] = 1.00000*(power(10.0000, -6.48000))
    constants[52] = 1.00000*(power(10.0000, -6.29000))
    constants[53] = 1.00000*(power(10.0000, -constants[19]))
    constants[54] = constants[13]
    constants[55] = constants[20]
    constants[56] = 0.664000*1.09000
    constants[57] = 0.900000*constants[56]
    constants[58] = 0.100000*constants[56]
    constants[59] = 0.807000*1.04400
    constants[60] = constants[53]
    constants[61] = constants[13]/constants[12]
    constants[62] = constants[12]
    constants[63] = constants[1]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = states[19]-states[22]
    algebraic[5] = constants[40]*1.00000*(algebraic[1]*states[25]-constants[8]*states[22])
    rates[22] = algebraic[5]/constants[59]
    algebraic[3] = states[20]-states[23]
    algebraic[7] = constants[40]*1.00000*(algebraic[3]*states[25]-constants[9]*states[23])
    rates[23] = algebraic[7]/constants[59]
    rates[25] = (-algebraic[5]-algebraic[7])/constants[59]
    algebraic[6] = states[4]-states[5]
    algebraic[10] = constants[40]*(algebraic[6]*states[10]-constants[8]*states[5])
    rates[5] = algebraic[10]/constants[57]
    algebraic[8] = states[6]-states[7]
    algebraic[12] = constants[40]*(algebraic[8]*states[10]-constants[9]*states[7])
    rates[10] = (-algebraic[12]-algebraic[10])/constants[57]
    rates[7] = algebraic[12]/constants[57]
    algebraic[16] = states[12]-states[13]
    algebraic[19] = constants[40]*(algebraic[16]*states[18]-constants[8]*states[13])
    rates[13] = algebraic[19]/constants[58]
    algebraic[17] = states[14]-states[15]
    algebraic[21] = constants[40]*(algebraic[17]*states[18]-constants[9]*states[15])
    rates[15] = algebraic[21]/constants[58]
    rates[18] = (-algebraic[21]-algebraic[19])/constants[58]
    algebraic[0] = constants[18]-states[1]
    algebraic[20] = (constants[27]*(constants[26]*algebraic[0]-states[1])*(1.00000+states[11]/constants[21]))/(1.00000+states[11]/constants[22])
    algebraic[9] = 1.00000*(power(10.0000, -states[8]))
    algebraic[13] = constants[4]-1.00000*constants[2]*log(algebraic[9]/1.00000e-07)
    algebraic[11] = constants[3]*states[0]+constants[2]*log(constants[60]/algebraic[9])
    algebraic[2] = constants[17]-states[2]
    algebraic[22] = constants[28]*(exp(-(algebraic[13]+4.00000*algebraic[11])/constants[2])*states[1]*states[2]-algebraic[0]*algebraic[2])
    rates[1] = (algebraic[20]-algebraic[22])/constants[57]
    algebraic[14] = constants[5]+2.00000*constants[2]*log(algebraic[9]/1.00000e-07)
    algebraic[4] = constants[15]-states[3]
    algebraic[23] = ((constants[29]*(1.00000+states[11]/constants[23]))/(1.00000+states[11]/constants[24]))*(exp(-((algebraic[14]+4.00000*algebraic[11])-2.00000*constants[3]*states[0])/(2.00000*constants[2]))*states[3]*(power(algebraic[2], 1.0/2))-algebraic[4]*(power(states[2], 1.0/2)))
    rates[2] = -(algebraic[22]-algebraic[23])/constants[57]
    algebraic[15] = constants[6]-2.00000*constants[2]*log(algebraic[9]/1.00000e-07)
    algebraic[24] = custom_piecewise([less(constants[0] , 1.00000e-12), 0.00000 , True, ((constants[30]*constants[0])/(constants[0]+constants[41]))*(algebraic[4]/constants[15])*(exp(-(algebraic[15]+2.00000*algebraic[11])/(2.00000*constants[2]))*algebraic[4]*(power(constants[0]/1.00000, 0.250000))-states[3]*exp((constants[3]*states[0])/constants[2]))])
    rates[3] = -(2.00000*algebraic[23]-2.00000*algebraic[24])/constants[57]
    algebraic[25] = constants[31]*(((exp(-(constants[7]-constants[14]*algebraic[11])/constants[2])*constants[9])/constants[8])*states[7]*states[11]-states[5]*1.00000)
    algebraic[26] = -0.650000*states[0]
    algebraic[27] = 0.350000*states[0]
    algebraic[28] = custom_piecewise([greater(algebraic[17] , constants[47]) | greater(algebraic[16] , constants[47]), constants[32]*(algebraic[17]/(algebraic[17]+algebraic[16]*exp((-constants[3]*algebraic[27])/constants[2]))-algebraic[8]/(algebraic[8]+algebraic[6]*exp((-constants[3]*algebraic[26])/constants[2])))*(algebraic[17]/(algebraic[17]+constants[37])) , True, 0.00000])
    rates[4] = (algebraic[25]-algebraic[28])/constants[57]
    rates[6] = (-algebraic[25]+algebraic[28])/constants[57]
    algebraic[29] = (states[17]*constants[60])/(constants[60]+constants[50])
    algebraic[31] = (states[11]*algebraic[9])/(algebraic[9]+constants[50])
    algebraic[33] = (constants[33]*(constants[60]*algebraic[29]-algebraic[9]*algebraic[31]))/(algebraic[29]+constants[25])
    rates[11] = (-algebraic[25]+algebraic[33])/constants[57]
    algebraic[34] = constants[46]*constants[44]*(states[19]-states[12])
    algebraic[36] = constants[38]*(constants[10]*states[14]*states[14]-states[16]*states[12])
    rates[12] = (algebraic[34]+algebraic[28]+algebraic[36])/constants[58]
    algebraic[32] = constants[46]*constants[44]*(states[20]-states[14])
    rates[14] = ((algebraic[32]-algebraic[28])-2.00000*algebraic[36])/constants[58]
    algebraic[30] = constants[46]*constants[44]*(states[21]-states[16])
    rates[16] = (algebraic[30]+algebraic[36])/constants[58]
    algebraic[37] = constants[46]*constants[45]*(states[24]-states[17])
    rates[17] = (-algebraic[33]+algebraic[37])/constants[58]
    algebraic[38] = constants[38]*1.00000*(constants[10]*states[20]*states[20]-states[21]*states[19])
    rates[21] = (-constants[61]*algebraic[30]+algebraic[38])/constants[59]
    rates[24] = (-algebraic[37]*constants[61]+constants[63]/constants[62])/constants[59]
    algebraic[39] = constants[34]*(constants[55]*algebraic[9]-states[9]*constants[60])
    algebraic[35] = custom_piecewise([greater(fabs(states[0]) , constants[48]), (constants[35]*states[0]*(constants[60]*exp((constants[3]*states[0])/constants[2])-algebraic[9]))/(exp((constants[3]*states[0])/constants[2])-1.00000) , True, (constants[35]*constants[2]*(constants[60]-algebraic[9]))/constants[3]])
    rates[8] = ((-1.00000/2.30300)*constants[42]*((((((algebraic[20]-5.00000*algebraic[22])-2.00000*algebraic[23])-4.00000*algebraic[24])+(constants[14]-1.00000)*algebraic[25]+2.00000*algebraic[33])-algebraic[39])+algebraic[35]))/constants[57]
    algebraic[18] = constants[16]-states[26]
    algebraic[40] = constants[39]*1.00000*(constants[11]*states[20]*states[26]*constants[53]-states[19]*algebraic[18])
    rates[19] = (((-constants[61]*algebraic[34]+algebraic[38])-constants[63]/constants[62])+algebraic[40])/constants[59]
    rates[20] = (((-constants[61]*algebraic[32]-2.00000*algebraic[38])+constants[63]/constants[62])-algebraic[40])/constants[59]
    rates[26] = -algebraic[40]/constants[59]
    algebraic[41] = custom_piecewise([greater(fabs(states[0]) , constants[49]), (constants[36]*states[0]*(constants[55]*exp((constants[3]*states[0])/constants[2])-states[9]))/(exp((constants[3]*states[0])/constants[2])-1.00000) , True, (constants[36]*constants[2]*(constants[55]-states[9]))/constants[3]])
    rates[9] = (algebraic[39]+algebraic[41])/constants[57]
    rates[0] = (((((4.00000*algebraic[22]+2.00000*algebraic[23]+4.00000*algebraic[24])-constants[14]*algebraic[25])-algebraic[28])-algebraic[35])-algebraic[41])/constants[43]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = states[19]-states[22]
    algebraic[5] = constants[40]*1.00000*(algebraic[1]*states[25]-constants[8]*states[22])
    algebraic[3] = states[20]-states[23]
    algebraic[7] = constants[40]*1.00000*(algebraic[3]*states[25]-constants[9]*states[23])
    algebraic[6] = states[4]-states[5]
    algebraic[10] = constants[40]*(algebraic[6]*states[10]-constants[8]*states[5])
    algebraic[8] = states[6]-states[7]
    algebraic[12] = constants[40]*(algebraic[8]*states[10]-constants[9]*states[7])
    algebraic[16] = states[12]-states[13]
    algebraic[19] = constants[40]*(algebraic[16]*states[18]-constants[8]*states[13])
    algebraic[17] = states[14]-states[15]
    algebraic[21] = constants[40]*(algebraic[17]*states[18]-constants[9]*states[15])
    algebraic[0] = constants[18]-states[1]
    algebraic[20] = (constants[27]*(constants[26]*algebraic[0]-states[1])*(1.00000+states[11]/constants[21]))/(1.00000+states[11]/constants[22])
    algebraic[9] = 1.00000*(power(10.0000, -states[8]))
    algebraic[13] = constants[4]-1.00000*constants[2]*log(algebraic[9]/1.00000e-07)
    algebraic[11] = constants[3]*states[0]+constants[2]*log(constants[60]/algebraic[9])
    algebraic[2] = constants[17]-states[2]
    algebraic[22] = constants[28]*(exp(-(algebraic[13]+4.00000*algebraic[11])/constants[2])*states[1]*states[2]-algebraic[0]*algebraic[2])
    algebraic[14] = constants[5]+2.00000*constants[2]*log(algebraic[9]/1.00000e-07)
    algebraic[4] = constants[15]-states[3]
    algebraic[23] = ((constants[29]*(1.00000+states[11]/constants[23]))/(1.00000+states[11]/constants[24]))*(exp(-((algebraic[14]+4.00000*algebraic[11])-2.00000*constants[3]*states[0])/(2.00000*constants[2]))*states[3]*(power(algebraic[2], 1.0/2))-algebraic[4]*(power(states[2], 1.0/2)))
    algebraic[15] = constants[6]-2.00000*constants[2]*log(algebraic[9]/1.00000e-07)
    algebraic[24] = custom_piecewise([less(constants[0] , 1.00000e-12), 0.00000 , True, ((constants[30]*constants[0])/(constants[0]+constants[41]))*(algebraic[4]/constants[15])*(exp(-(algebraic[15]+2.00000*algebraic[11])/(2.00000*constants[2]))*algebraic[4]*(power(constants[0]/1.00000, 0.250000))-states[3]*exp((constants[3]*states[0])/constants[2]))])
    algebraic[25] = constants[31]*(((exp(-(constants[7]-constants[14]*algebraic[11])/constants[2])*constants[9])/constants[8])*states[7]*states[11]-states[5]*1.00000)
    algebraic[26] = -0.650000*states[0]
    algebraic[27] = 0.350000*states[0]
    algebraic[28] = custom_piecewise([greater(algebraic[17] , constants[47]) | greater(algebraic[16] , constants[47]), constants[32]*(algebraic[17]/(algebraic[17]+algebraic[16]*exp((-constants[3]*algebraic[27])/constants[2]))-algebraic[8]/(algebraic[8]+algebraic[6]*exp((-constants[3]*algebraic[26])/constants[2])))*(algebraic[17]/(algebraic[17]+constants[37])) , True, 0.00000])
    algebraic[29] = (states[17]*constants[60])/(constants[60]+constants[50])
    algebraic[31] = (states[11]*algebraic[9])/(algebraic[9]+constants[50])
    algebraic[33] = (constants[33]*(constants[60]*algebraic[29]-algebraic[9]*algebraic[31]))/(algebraic[29]+constants[25])
    algebraic[34] = constants[46]*constants[44]*(states[19]-states[12])
    algebraic[36] = constants[38]*(constants[10]*states[14]*states[14]-states[16]*states[12])
    algebraic[32] = constants[46]*constants[44]*(states[20]-states[14])
    algebraic[30] = constants[46]*constants[44]*(states[21]-states[16])
    algebraic[37] = constants[46]*constants[45]*(states[24]-states[17])
    algebraic[38] = constants[38]*1.00000*(constants[10]*states[20]*states[20]-states[21]*states[19])
    algebraic[39] = constants[34]*(constants[55]*algebraic[9]-states[9]*constants[60])
    algebraic[35] = custom_piecewise([greater(fabs(states[0]) , constants[48]), (constants[35]*states[0]*(constants[60]*exp((constants[3]*states[0])/constants[2])-algebraic[9]))/(exp((constants[3]*states[0])/constants[2])-1.00000) , True, (constants[35]*constants[2]*(constants[60]-algebraic[9]))/constants[3]])
    algebraic[18] = constants[16]-states[26]
    algebraic[40] = constants[39]*1.00000*(constants[11]*states[20]*states[26]*constants[53]-states[19]*algebraic[18])
    algebraic[41] = custom_piecewise([greater(fabs(states[0]) , constants[49]), (constants[36]*states[0]*(constants[55]*exp((constants[3]*states[0])/constants[2])-states[9]))/(exp((constants[3]*states[0])/constants[2])-1.00000) , True, (constants[36]*constants[2]*(constants[55]-states[9]))/constants[3]])
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