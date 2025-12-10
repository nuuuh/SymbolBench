# Size of variable arrays:
sizeAlgebraic = 28
sizeStates = 14
sizeConstants = 98
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "actCycA in component actCycA (dimensionless)"
    legend_constants[0] = "k_asa in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[1] = "k_sap in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[2] = "k_sapp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[3] = "k_dia in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[15] = "V_di in component V_di (first_order_rate_constant)"
    legend_algebraic[9] = "V_da in component V_da (first_order_rate_constant)"
    legend_states[1] = "mass in component mass (dimensionless)"
    legend_algebraic[1] = "Tri_A in component Tri_A (dimensionless)"
    legend_algebraic[7] = "TF_E in component TF_E (dimensionless)"
    legend_algebraic[3] = "freeCKI in component freeCKI (dimensionless)"
    legend_states[2] = "actCycB in component actCycB (dimensionless)"
    legend_constants[4] = "k_asb in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[5] = "k_dib in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[18] = "V_sb in component V_sb (first_order_rate_constant)"
    legend_algebraic[19] = "V_db in component V_db (first_order_rate_constant)"
    legend_algebraic[27] = "V_25 in component V_25 (first_order_rate_constant)"
    legend_algebraic[23] = "V_wee in component V_wee (first_order_rate_constant)"
    legend_states[3] = "CycB in component CycB (dimensionless)"
    legend_states[4] = "Tri_B in component Tri_B (dimensionless)"
    legend_states[5] = "preMPF in component preMPF (dimensionless)"
    legend_states[6] = "actCycE in component actCycE (dimensionless)"
    legend_constants[6] = "k_die in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[7] = "k_ase in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[8] = "k_sep in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[9] = "k_sepp in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[8] = "V_de in component V_de (first_order_rate_constant)"
    legend_algebraic[2] = "Tri_E in component Tri_E (dimensionless)"
    legend_states[7] = "CycA in component CycA (dimensionless)"
    legend_states[8] = "CycE in component CycE (dimensionless)"
    legend_algebraic[0] = "CycD in component CycD (dimensionless)"
    legend_constants[10] = "CycD_0 in component CycD (dimensionless)"
    legend_states[9] = "Cdc20_A in component Cdc20_A (dimensionless)"
    legend_constants[11] = "k_a20 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[12] = "k_i20 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[13] = "k_d20 in component kinetic_parameters (first_order_rate_constant)"
    legend_states[10] = "APCP in component APCP (dimensionless)"
    legend_states[11] = "Cdc20_T in component Cdc20_T (dimensionless)"
    legend_constants[14] = "J_a20 in component Cdc20_A (dimensionless)"
    legend_constants[15] = "J_i20 in component Cdc20_A (dimensionless)"
    legend_constants[16] = "k_s20p in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[17] = "k_s20pp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[18] = "n in component Cdc20_T (dimensionless)"
    legend_constants[19] = "J_20 in component Cdc20_T (dimensionless)"
    legend_constants[20] = "k_aAPC in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[21] = "k_iAPC in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[22] = "J_aAPC in component APCP (dimensionless)"
    legend_constants[23] = "J_iAPC in component APCP (dimensionless)"
    legend_states[12] = "CKI in component CKI (dimensionless)"
    legend_states[13] = "Cdh1 in component Cdh1 (dimensionless)"
    legend_constants[24] = "k_ah1p in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[25] = "k_ah1pp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[26] = "k_ih1p in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[27] = "k_ih1pp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[28] = "k_ih1ppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[29] = "k_ih1pppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[30] = "k_ih1ppppp in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[11] = "Cdc14 in component TF_I (dimensionless)"
    legend_constants[31] = "J_ah1 in component Cdh1 (dimensionless)"
    legend_constants[32] = "J_ih1 in component Cdh1 (dimensionless)"
    legend_constants[33] = "mu in component mass (first_order_rate_constant)"
    legend_constants[34] = "maxmass in component mass (dimensionless)"
    legend_algebraic[4] = "V_atf in component V_atf (first_order_rate_constant)"
    legend_constants[35] = "k_atfp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[36] = "k_atfpp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[37] = "k_atfppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[38] = "k_atfpppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[39] = "k_itfp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[40] = "k_itfpp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[41] = "k_itfppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[42] = "J_itf in component TF_E (dimensionless)"
    legend_constants[43] = "J_atf in component TF_E (dimensionless)"
    legend_algebraic[5] = "A1 in component TF_E (first_order_rate_constant)"
    legend_algebraic[6] = "A2 in component TF_E (first_order_rate_constant)"
    legend_constants[88] = "A3 in component TF_E (dimensionless)"
    legend_constants[91] = "A4 in component TF_E (dimensionless)"
    legend_constants[44] = "k_dep in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[45] = "k_depp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[46] = "k_deppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[47] = "k_depppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[48] = "k_dap in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[49] = "k_dapp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[50] = "k_dappp in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[13] = "TF_I in component TF_I (dimensionless)"
    legend_constants[51] = "k_afi in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[52] = "k_ifip in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[53] = "k_ifipp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[54] = "J_ifi in component TF_I (dimensionless)"
    legend_constants[55] = "J_afi in component TF_I (dimensionless)"
    legend_algebraic[12] = "A1 in component TF_I (first_order_rate_constant)"
    legend_algebraic[10] = "A2 in component TF_I (first_order_rate_constant)"
    legend_constants[87] = "A3 in component TF_I (dimensionless)"
    legend_constants[90] = "A4 in component TF_I (dimensionless)"
    legend_algebraic[14] = "V_si in component V_si (first_order_rate_constant)"
    legend_constants[56] = "k_sip in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[57] = "k_sipp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[58] = "k_dip in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[59] = "k_dipp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[60] = "k_dippp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[61] = "k_dipppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[62] = "k_dippppp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[63] = "J_14di in component kinetic_parameters (dimensionless)"
    legend_algebraic[17] = "TF_B in component TF_B (dimensionless)"
    legend_constants[64] = "k_afb in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[65] = "k_ifb in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[66] = "J_ifb in component TF_B (dimensionless)"
    legend_constants[67] = "J_afb in component TF_B (dimensionless)"
    legend_algebraic[16] = "A1 in component TF_B (first_order_rate_constant)"
    legend_constants[89] = "A2 in component TF_B (first_order_rate_constant)"
    legend_constants[92] = "A3 in component TF_B (dimensionless)"
    legend_constants[94] = "A4 in component TF_B (dimensionless)"
    legend_constants[68] = "k_sbp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[69] = "k_sbpp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[70] = "k_dbp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[71] = "k_dbpp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[72] = "k_dbppp in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[22] = "Wee1 in component Wee1 (dimensionless)"
    legend_constants[73] = "k_aweep in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[74] = "k_aweepp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[75] = "k_iwee in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[76] = "J_iwee in component Wee1 (dimensionless)"
    legend_constants[77] = "J_awee in component Wee1 (dimensionless)"
    legend_algebraic[20] = "A1 in component Wee1 (first_order_rate_constant)"
    legend_algebraic[21] = "A2 in component Wee1 (first_order_rate_constant)"
    legend_constants[93] = "A3 in component Wee1 (dimensionless)"
    legend_constants[95] = "A4 in component Wee1 (dimensionless)"
    legend_constants[78] = "k_weep in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[79] = "k_weepp in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[26] = "Cdc25 in component Cdc25 (dimensionless)"
    legend_constants[80] = "k_a25 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[81] = "k_i25p in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[82] = "k_i25pp in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[83] = "J_i25 in component Cdc25 (dimensionless)"
    legend_constants[84] = "J_a25 in component Cdc25 (dimensionless)"
    legend_algebraic[24] = "A1 in component Cdc25 (first_order_rate_constant)"
    legend_algebraic[25] = "A2 in component Cdc25 (first_order_rate_constant)"
    legend_constants[96] = "A3 in component Cdc25 (dimensionless)"
    legend_constants[97] = "A4 in component Cdc25 (dimensionless)"
    legend_constants[85] = "k_25p in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[86] = "k_25pp in component kinetic_parameters (first_order_rate_constant)"
    legend_rates[0] = "d/dt actCycA in component actCycA (dimensionless)"
    legend_rates[2] = "d/dt actCycB in component actCycB (dimensionless)"
    legend_rates[6] = "d/dt actCycE in component actCycE (dimensionless)"
    legend_rates[7] = "d/dt CycA in component CycA (dimensionless)"
    legend_rates[3] = "d/dt CycB in component CycB (dimensionless)"
    legend_rates[8] = "d/dt CycE in component CycE (dimensionless)"
    legend_rates[9] = "d/dt Cdc20_A in component Cdc20_A (dimensionless)"
    legend_rates[11] = "d/dt Cdc20_T in component Cdc20_T (dimensionless)"
    legend_rates[10] = "d/dt APCP in component APCP (dimensionless)"
    legend_rates[13] = "d/dt Cdh1 in component Cdh1 (dimensionless)"
    legend_rates[1] = "d/dt mass in component mass (dimensionless)"
    legend_rates[5] = "d/dt preMPF in component preMPF (dimensionless)"
    legend_rates[4] = "d/dt Tri_B in component Tri_B (dimensionless)"
    legend_rates[12] = "d/dt CKI in component CKI (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.09450308233499527
    constants[0] = 50
    constants[1] = 0.0015
    constants[2] = 0.01
    constants[3] = 0.06
    states[1] = 1.338267803192139
    states[2] = 0.1903585940599442
    constants[4] = 60
    constants[5] = 0.05
    states[3] = 0.09450308233499527
    states[4] = 0.01
    states[5] = 0.01
    states[6] = 0.2092024385929108
    constants[6] = 0
    constants[7] = 0
    constants[8] = 0
    constants[9] = 0.15
    states[7] = 0.09450308233499527
    states[8] = 0.09450308233499527
    constants[10] = 0.108
    states[9] = 0.3572723865509033
    constants[11] = 1
    constants[12] = 0.16
    constants[13] = 0.05
    states[10] = 0.07591402530670166
    states[11] = 0.7702387571334839
    constants[14] = 1
    constants[15] = 1
    constants[16] = 0.001
    constants[17] = 1
    constants[18] = 1
    constants[19] = 10
    constants[20] = 0.1
    constants[21] = 0.15
    constants[22] = 0.1
    constants[23] = 0.1
    states[12] = 0.01
    states[13] = 0.7189393639564514
    constants[24] = 0.02
    constants[25] = 0.8
    constants[26] = 0.001
    constants[27] = 0.35
    constants[28] = 0.1
    constants[29] = 0.06
    constants[30] = 0.005
    constants[31] = 0.03
    constants[32] = 0.03
    constants[33] = 0.005776
    constants[34] = 10000
    constants[35] = 0
    constants[36] = 1.5
    constants[37] = 0.38
    constants[38] = 3
    constants[39] = 0.75
    constants[40] = 8
    constants[41] = 0
    constants[42] = 0.01
    constants[43] = 0.01
    constants[44] = 0.12
    constants[45] = 0
    constants[46] = 0
    constants[47] = 0
    constants[48] = 0.01
    constants[49] = 0.16
    constants[50] = 0
    constants[51] = 6
    constants[52] = 0.008
    constants[53] = 0.05
    constants[54] = 2
    constants[55] = 1
    constants[56] = 0.018
    constants[57] = 0.18
    constants[58] = 0.02
    constants[59] = 0.1
    constants[60] = 0.8
    constants[61] = 0.12
    constants[62] = 0.1
    constants[63] = 12
    constants[64] = 1
    constants[65] = 0.15
    constants[66] = 0.1
    constants[67] = 0.1
    constants[68] = 0.004
    constants[69] = 0.04
    constants[70] = 0.003
    constants[71] = 0.4
    constants[72] = 0.15
    constants[73] = 0.3
    constants[74] = 0
    constants[75] = 1
    constants[76] = 0.05
    constants[77] = 0.05
    constants[78] = 0.02
    constants[79] = 0.2
    constants[80] = 0
    constants[81] = 0.3
    constants[82] = 0
    constants[83] = 0.1
    constants[84] = 0.1
    constants[85] = 0.01
    constants[86] = 5
    constants[87] = constants[55]
    constants[88] = constants[43]
    constants[89] = constants[65]
    constants[90] = constants[54]
    constants[91] = constants[42]
    constants[92] = constants[67]
    constants[93] = constants[77]
    constants[94] = constants[66]
    constants[95] = constants[76]
    constants[96] = constants[84]
    constants[97] = constants[83]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[9] = (constants[11]*states[10]*(states[11]-states[9]))/((constants[14]+states[11])-states[9])-(constants[12]/(constants[15]+states[9])+constants[13])*states[9]
    rates[11] = (constants[16]+constants[17]*(power(states[2], constants[18])))/(power(constants[19], constants[18])+power(states[2], constants[18]))-constants[13]*states[11]
    rates[10] = (constants[20]*states[2]*(1.00000-states[10]))/((constants[22]+1.00000)-states[10])-(constants[21]*states[10])/(constants[23]+states[10])
    rates[1] = constants[33]*states[1]*(1.00000-states[1]/constants[34])
    algebraic[0] = constants[10]*states[1]
    algebraic[4] = constants[35]+constants[36]*states[0]+constants[37]*states[6]+constants[38]*algebraic[0]
    algebraic[5] = algebraic[4]
    algebraic[6] = constants[39]+constants[40]*states[2]+constants[41]*states[0]
    algebraic[7] = (2.00000*constants[91]*algebraic[5])/((algebraic[6]+-algebraic[5]+constants[88]*algebraic[6]+constants[91]*algebraic[5])+power(power(algebraic[6]+-algebraic[5]+constants[88]*algebraic[6]+constants[91]*algebraic[5], 2.00000)-4.00000*(algebraic[6]-algebraic[5])*constants[91]*algebraic[5], 1.0/2))
    algebraic[8] = constants[44]+constants[45]*states[6]+constants[46]*states[0]+constants[47]*states[2]
    rates[8] = (constants[8]+constants[9]*algebraic[7])*states[1]-algebraic[8]*states[8]
    algebraic[9] = constants[48]+constants[49]*states[9]+constants[50]*states[11]
    rates[7] = (constants[1]+constants[2]*algebraic[7])*states[1]-algebraic[9]*states[7]
    algebraic[11] = states[9]
    rates[13] = ((constants[24]+constants[25]*algebraic[11])*(1.00000-states[13]))/((constants[31]+1.00000)-states[13])-((constants[26]+constants[27]*states[0]+constants[28]*states[2]+constants[29]*states[6]+constants[30]*algebraic[0])*states[13])/(constants[32]+states[13])
    algebraic[15] = (constants[58]+constants[59]*states[0]+constants[60]*states[2]+constants[61]*states[6]+constants[62]*algebraic[0])/(1.00000+algebraic[11]/constants[63])
    algebraic[1] = states[7]-states[0]
    algebraic[2] = states[8]-states[6]
    algebraic[3] = states[12]-(states[4]+algebraic[1]+algebraic[2])
    rates[0] = ((constants[1]+constants[2]*algebraic[7])*states[1]+(algebraic[15]+constants[3])*algebraic[1])-(algebraic[9]+constants[0]*algebraic[3])*states[0]
    rates[6] = ((constants[8]+constants[9]*algebraic[7])*states[1]+(algebraic[15]+constants[6])*algebraic[2])-(algebraic[8]+constants[7]*algebraic[3])*states[6]
    algebraic[12] = constants[51]*algebraic[11]
    algebraic[10] = constants[52]+constants[53]*states[2]
    algebraic[13] = (2.00000*constants[90]*algebraic[12])/((algebraic[10]+-algebraic[12]+constants[87]*algebraic[10]+constants[90]*algebraic[12])+power(power(algebraic[10]+-algebraic[12]+constants[87]*algebraic[10]+constants[90]*algebraic[12], 2.00000)-4.00000*(algebraic[10]-algebraic[12])*constants[90]*algebraic[12], 1.0/2))
    algebraic[14] = constants[56]+constants[57]*algebraic[13]
    rates[12] = algebraic[14]-algebraic[15]*states[12]
    algebraic[16] = constants[64]*states[2]
    algebraic[17] = (2.00000*constants[94]*algebraic[16])/((constants[89]+-algebraic[16]+constants[92]*constants[89]+constants[94]*algebraic[16])+power(power(constants[89]+-algebraic[16]+constants[92]*constants[89]+constants[94]*algebraic[16], 2.00000)-4.00000*(constants[89]-algebraic[16])*constants[94]*algebraic[16], 1.0/2))
    algebraic[18] = constants[68]+constants[69]*algebraic[17]
    algebraic[19] = constants[70]+constants[71]*states[13]+constants[72]*states[9]
    rates[3] = algebraic[18]*states[1]-algebraic[19]*states[3]
    rates[4] = constants[4]*(states[3]-states[4])*algebraic[3]-(constants[5]+algebraic[19]+algebraic[15])*states[4]
    algebraic[24] = constants[80]*states[2]
    algebraic[25] = constants[81]+constants[82]*algebraic[11]
    algebraic[26] = (2.00000*constants[97]*algebraic[24])/((algebraic[25]+-algebraic[24]+constants[96]*algebraic[25]+constants[97]*algebraic[24])+power(power(algebraic[25]+-algebraic[24]+constants[96]*algebraic[25]+constants[97]*algebraic[24], 2.00000)-4.00000*(algebraic[25]-algebraic[24])*constants[97]*algebraic[24], 1.0/2))
    algebraic[27] = constants[85]+constants[86]*algebraic[26]
    algebraic[20] = constants[73]+constants[74]*algebraic[11]
    algebraic[21] = constants[75]*states[2]
    algebraic[22] = (2.00000*constants[95]*algebraic[20])/((algebraic[21]+-algebraic[20]+constants[93]*algebraic[21]+constants[95]*algebraic[20])+power(power(algebraic[21]+-algebraic[20]+constants[93]*algebraic[21]+constants[95]*algebraic[20], 2.00000)-4.00000*(algebraic[21]-algebraic[20])*constants[95]*algebraic[20], 1.0/2))
    algebraic[23] = constants[78]+constants[79]*algebraic[22]
    rates[2] = (algebraic[18]*states[1]+algebraic[27]*(states[3]-(states[4]+states[2]))+(constants[5]+algebraic[15])*(states[3]-(states[5]+states[2])))-(algebraic[19]+algebraic[23]+constants[4]*algebraic[3])*states[2]
    rates[5] = algebraic[23]*(states[3]-states[5])-(algebraic[27]+algebraic[19])*states[5]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[10]*states[1]
    algebraic[4] = constants[35]+constants[36]*states[0]+constants[37]*states[6]+constants[38]*algebraic[0]
    algebraic[5] = algebraic[4]
    algebraic[6] = constants[39]+constants[40]*states[2]+constants[41]*states[0]
    algebraic[7] = (2.00000*constants[91]*algebraic[5])/((algebraic[6]+-algebraic[5]+constants[88]*algebraic[6]+constants[91]*algebraic[5])+power(power(algebraic[6]+-algebraic[5]+constants[88]*algebraic[6]+constants[91]*algebraic[5], 2.00000)-4.00000*(algebraic[6]-algebraic[5])*constants[91]*algebraic[5], 1.0/2))
    algebraic[8] = constants[44]+constants[45]*states[6]+constants[46]*states[0]+constants[47]*states[2]
    algebraic[9] = constants[48]+constants[49]*states[9]+constants[50]*states[11]
    algebraic[11] = states[9]
    algebraic[15] = (constants[58]+constants[59]*states[0]+constants[60]*states[2]+constants[61]*states[6]+constants[62]*algebraic[0])/(1.00000+algebraic[11]/constants[63])
    algebraic[1] = states[7]-states[0]
    algebraic[2] = states[8]-states[6]
    algebraic[3] = states[12]-(states[4]+algebraic[1]+algebraic[2])
    algebraic[12] = constants[51]*algebraic[11]
    algebraic[10] = constants[52]+constants[53]*states[2]
    algebraic[13] = (2.00000*constants[90]*algebraic[12])/((algebraic[10]+-algebraic[12]+constants[87]*algebraic[10]+constants[90]*algebraic[12])+power(power(algebraic[10]+-algebraic[12]+constants[87]*algebraic[10]+constants[90]*algebraic[12], 2.00000)-4.00000*(algebraic[10]-algebraic[12])*constants[90]*algebraic[12], 1.0/2))
    algebraic[14] = constants[56]+constants[57]*algebraic[13]
    algebraic[16] = constants[64]*states[2]
    algebraic[17] = (2.00000*constants[94]*algebraic[16])/((constants[89]+-algebraic[16]+constants[92]*constants[89]+constants[94]*algebraic[16])+power(power(constants[89]+-algebraic[16]+constants[92]*constants[89]+constants[94]*algebraic[16], 2.00000)-4.00000*(constants[89]-algebraic[16])*constants[94]*algebraic[16], 1.0/2))
    algebraic[18] = constants[68]+constants[69]*algebraic[17]
    algebraic[19] = constants[70]+constants[71]*states[13]+constants[72]*states[9]
    algebraic[24] = constants[80]*states[2]
    algebraic[25] = constants[81]+constants[82]*algebraic[11]
    algebraic[26] = (2.00000*constants[97]*algebraic[24])/((algebraic[25]+-algebraic[24]+constants[96]*algebraic[25]+constants[97]*algebraic[24])+power(power(algebraic[25]+-algebraic[24]+constants[96]*algebraic[25]+constants[97]*algebraic[24], 2.00000)-4.00000*(algebraic[25]-algebraic[24])*constants[97]*algebraic[24], 1.0/2))
    algebraic[27] = constants[85]+constants[86]*algebraic[26]
    algebraic[20] = constants[73]+constants[74]*algebraic[11]
    algebraic[21] = constants[75]*states[2]
    algebraic[22] = (2.00000*constants[95]*algebraic[20])/((algebraic[21]+-algebraic[20]+constants[93]*algebraic[21]+constants[95]*algebraic[20])+power(power(algebraic[21]+-algebraic[20]+constants[93]*algebraic[21]+constants[95]*algebraic[20], 2.00000)-4.00000*(algebraic[21]-algebraic[20])*constants[95]*algebraic[20], 1.0/2))
    algebraic[23] = constants[78]+constants[79]*algebraic[22]
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