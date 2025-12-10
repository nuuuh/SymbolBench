# Size of variable arrays:
sizeAlgebraic = 42
sizeStates = 31
sizeConstants = 57
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "Casp8 in component Casp8 (micromolar)"
    legend_states[1] = "Casp8Bid in component Casp8 (micromolar)"
    legend_states[2] = "Bid in component Bid (micromolar)"
    legend_algebraic[2] = "J_0 in component Casp8 (flux)"
    legend_algebraic[0] = "J_f0 in component Casp8 (flux)"
    legend_algebraic[1] = "J_Casp8 in component Casp8 (flux)"
    legend_constants[0] = "k_f0 in component model_constant (first_order_rate_constant)"
    legend_constants[1] = "k_10 in component model_constant (second_order_rate_constant)"
    legend_constants[2] = "k_00 in component model_constant (first_order_rate_constant)"
    legend_constants[3] = "u in component model_constant (first_order_rate_constant)"
    legend_states[3] = "Apaf_1 in component Apaf_1 (micromolar)"
    legend_states[4] = "CytcApaf_1 in component CytcApaf_1 (micromolar)"
    legend_states[5] = "Cytc in component Cytc (micromolar)"
    legend_algebraic[3] = "J_Apaf_1 in component Apaf_1 (flux)"
    legend_algebraic[4] = "J_1 in component Apaf_1 (flux)"
    legend_constants[4] = "P_Apaf_1 in component model_constant (flux)"
    legend_constants[5] = "k_11 in component model_constant (second_order_rate_constant)"
    legend_constants[6] = "k_01 in component model_constant (first_order_rate_constant)"
    legend_states[6] = "Apop in component Apop (micromolar)"
    legend_algebraic[6] = "J_1b in component CytcApaf_1 (flux)"
    legend_constants[7] = "k_11b in component model_constant (rate)"
    legend_constants[8] = "k_01b in component model_constant (first_order_rate_constant)"
    legend_constants[9] = "p in component model_constant (dimensionless)"
    legend_algebraic[9] = "J_14 in component Cytc_mito (flux)"
    legend_algebraic[5] = "J_Cytc in component Cytc (flux)"
    legend_states[7] = "Cytc_mito in component Cytc_mito (micromolar)"
    legend_states[8] = "Bax_2 in component Bax_2 (micromolar)"
    legend_algebraic[7] = "J_Cytc_mito in component Cytc_mito (flux)"
    legend_constants[10] = "P_Cytc_mito in component model_constant (flux)"
    legend_constants[11] = "k14 in component model_constant (second_order_rate_constant)"
    legend_algebraic[15] = "J_12b in component tBidBax (flux)"
    legend_algebraic[8] = "J_Bax_2 in component Bax_2 (flux)"
    legend_states[9] = "tBid_mito in component tBid_mito (micromolar)"
    legend_states[10] = "Bax in component Bax (micromolar)"
    legend_algebraic[12] = "J_11 in component tBid (flux)"
    legend_algebraic[11] = "J_12a in component tBid_mito (flux)"
    legend_algebraic[10] = "J_tBid_mito in component tBid_mito (flux)"
    legend_constants[12] = "k12a in component model_constant (second_order_rate_constant)"
    legend_states[11] = "tBid in component tBid (micromolar)"
    legend_algebraic[23] = "J_f8 in component Casp3Bid (flux)"
    legend_algebraic[13] = "J_tBid in component tBid (flux)"
    legend_constants[13] = "k11 in component model_constant (first_order_rate_constant)"
    legend_states[12] = "tBidBax in component tBidBax (micromolar)"
    legend_algebraic[14] = "J_tBidBax in component tBidBax (flux)"
    legend_constants[14] = "k12b in component model_constant (second_order_rate_constant)"
    legend_states[13] = "Bcl_2 in component Bcl_2 (micromolar)"
    legend_algebraic[17] = "J_13 in component Bax (flux)"
    legend_algebraic[16] = "J_Bax in component Bax (flux)"
    legend_constants[55] = "P_Bax in component Bax (flux)"
    legend_constants[15] = "k13 in component model_constant (second_order_rate_constant)"
    legend_constants[16] = "P_oBax in component model_constant (flux)"
    legend_constants[17] = "p53 in component model_constant (micromolar)"
    legend_constants[18] = "p53_thresh in component model_constant (micromolar)"
    legend_constants[19] = "u_Bax in component model_constant (first_order_rate_constant)"
    legend_states[14] = "Casp3 in component Casp3 (micromolar)"
    legend_states[15] = "Casp3Bcl_2 in component Casp3Bcl_2 (micromolar)"
    legend_algebraic[18] = "J_9 in component Bcl_2 (flux)"
    legend_algebraic[20] = "J_Bcl_2 in component Bcl_2 (flux)"
    legend_constants[20] = "k_19 in component model_constant (second_order_rate_constant)"
    legend_constants[21] = "k_09 in component model_constant (first_order_rate_constant)"
    legend_constants[56] = "P_Bcl_2 in component Bcl_2 (flux)"
    legend_constants[22] = "P_oBcl_2 in component model_constant (flux)"
    legend_constants[23] = "u_Bcl_2 in component model_constant (first_order_rate_constant)"
    legend_algebraic[19] = "J_f9 in component Casp3Bcl_2 (flux)"
    legend_constants[24] = "k_f9 in component model_constant (first_order_rate_constant)"
    legend_states[16] = "Casp3Bid in component Casp3Bid (micromolar)"
    legend_algebraic[21] = "J_8 in component Casp3Bid (flux)"
    legend_constants[25] = "k_18 in component model_constant (second_order_rate_constant)"
    legend_constants[26] = "k_08 in component model_constant (first_order_rate_constant)"
    legend_constants[27] = "k_f8 in component model_constant (first_order_rate_constant)"
    legend_algebraic[22] = "J_Bid in component Bid (flux)"
    legend_constants[28] = "P_Bid in component model_constant (flux)"
    legend_states[17] = "IAP in component IAP (micromolar)"
    legend_states[18] = "Casp3IAP in component Casp3 (micromolar)"
    legend_algebraic[32] = "J_f6 in component Casp9Pro3 (flux)"
    legend_algebraic[31] = "J_f6b in component ApopCasp9_2Pro3 (flux)"
    legend_algebraic[25] = "J_7 in component Casp3 (flux)"
    legend_algebraic[24] = "J_Casp3 in component Casp3 (flux)"
    legend_constants[29] = "k_17 in component model_constant (second_order_rate_constant)"
    legend_constants[30] = "k_07 in component model_constant (first_order_rate_constant)"
    legend_algebraic[27] = "J_2 in component Pro9 (flux)"
    legend_algebraic[40] = "J_4b in component ApopCasp9 (flux)"
    legend_states[19] = "Pro9 in component Pro9 (micromolar)"
    legend_states[20] = "ApopPro9 in component Pro9 (micromolar)"
    legend_states[21] = "ApopPro9_2 in component ApopPro9_2 (micromolar)"
    legend_algebraic[28] = "J_3 in component Pro9 (flux)"
    legend_algebraic[26] = "J_Pro9 in component Pro9 (flux)"
    legend_constants[31] = "k_12 in component model_constant (second_order_rate_constant)"
    legend_constants[32] = "k_02 in component model_constant (first_order_rate_constant)"
    legend_constants[33] = "k_13 in component model_constant (second_order_rate_constant)"
    legend_constants[34] = "k_03 in component model_constant (first_order_rate_constant)"
    legend_constants[35] = "P_Pro9 in component model_constant (flux)"
    legend_algebraic[29] = "J_f3 in component ApopPro9_2 (flux)"
    legend_constants[36] = "k_f3 in component model_constant (first_order_rate_constant)"
    legend_states[22] = "ApopCasp9_2Pro3 in component ApopCasp9_2Pro3 (micromolar)"
    legend_states[23] = "ApopCasp9_2 in component ApopCasp9_2 (micromolar)"
    legend_states[24] = "Pro3 in component Pro3 (micromolar)"
    legend_algebraic[30] = "J_6b in component ApopCasp9_2Pro3 (flux)"
    legend_constants[37] = "k_16b in component model_constant (second_order_rate_constant)"
    legend_constants[38] = "k_06b in component model_constant (first_order_rate_constant)"
    legend_constants[39] = "k_f6b in component model_constant (first_order_rate_constant)"
    legend_states[25] = "Casp9Pro3 in component Casp9Pro3 (micromolar)"
    legend_algebraic[33] = "J_6 in component Pro3 (flux)"
    legend_constants[40] = "k_f6 in component model_constant (first_order_rate_constant)"
    legend_states[26] = "Casp9 in component Casp9 (micromolar)"
    legend_algebraic[34] = "J_Pro3 in component Pro3 (flux)"
    legend_constants[41] = "k_16 in component model_constant (second_order_rate_constant)"
    legend_constants[42] = "k_06 in component model_constant (first_order_rate_constant)"
    legend_constants[43] = "P_Pro3 in component model_constant (flux)"
    legend_states[27] = "Casp9IAP in component IAP (micromolar)"
    legend_states[28] = "ApopCasp9IAP in component IAP (micromolar)"
    legend_states[29] = "ApopCasp9_2IAP in component IAP (micromolar)"
    legend_states[30] = "ApopCasp9 in component ApopCasp9 (micromolar)"
    legend_algebraic[35] = "J_5 in component IAP (flux)"
    legend_algebraic[36] = "J_5b in component IAP (flux)"
    legend_algebraic[37] = "J_5c in component IAP (flux)"
    legend_algebraic[39] = "J_IAP in component IAP (flux)"
    legend_constants[44] = "P_IAP in component model_constant (flux)"
    legend_constants[45] = "k_15 in component model_constant (second_order_rate_constant)"
    legend_constants[46] = "k_05 in component model_constant (first_order_rate_constant)"
    legend_constants[47] = "k_15b in component model_constant (second_order_rate_constant)"
    legend_constants[48] = "k_05b in component model_constant (first_order_rate_constant)"
    legend_constants[49] = "k_15c in component model_constant (second_order_rate_constant)"
    legend_constants[50] = "k_05c in component model_constant (first_order_rate_constant)"
    legend_algebraic[38] = "J_4 in component ApopCasp9 (flux)"
    legend_constants[51] = "k_14 in component model_constant (first_order_rate_constant)"
    legend_constants[52] = "k_04 in component model_constant (second_order_rate_constant)"
    legend_constants[53] = "k_14b in component model_constant (first_order_rate_constant)"
    legend_constants[54] = "k_04b in component model_constant (second_order_rate_constant)"
    legend_algebraic[41] = "J_Casp9 in component Casp9 (flux)"
    legend_rates[0] = "d/dt Casp8 in component Casp8 (micromolar)"
    legend_rates[1] = "d/dt Casp8Bid in component Casp8 (micromolar)"
    legend_rates[3] = "d/dt Apaf_1 in component Apaf_1 (micromolar)"
    legend_rates[4] = "d/dt CytcApaf_1 in component CytcApaf_1 (micromolar)"
    legend_rates[5] = "d/dt Cytc in component Cytc (micromolar)"
    legend_rates[7] = "d/dt Cytc_mito in component Cytc_mito (micromolar)"
    legend_rates[8] = "d/dt Bax_2 in component Bax_2 (micromolar)"
    legend_rates[9] = "d/dt tBid_mito in component tBid_mito (micromolar)"
    legend_rates[11] = "d/dt tBid in component tBid (micromolar)"
    legend_rates[12] = "d/dt tBidBax in component tBidBax (micromolar)"
    legend_rates[10] = "d/dt Bax in component Bax (micromolar)"
    legend_rates[13] = "d/dt Bcl_2 in component Bcl_2 (micromolar)"
    legend_rates[15] = "d/dt Casp3Bcl_2 in component Casp3Bcl_2 (micromolar)"
    legend_rates[16] = "d/dt Casp3Bid in component Casp3Bid (micromolar)"
    legend_rates[2] = "d/dt Bid in component Bid (micromolar)"
    legend_rates[18] = "d/dt Casp3IAP in component Casp3 (micromolar)"
    legend_rates[14] = "d/dt Casp3 in component Casp3 (micromolar)"
    legend_rates[6] = "d/dt Apop in component Apop (micromolar)"
    legend_rates[19] = "d/dt Pro9 in component Pro9 (micromolar)"
    legend_rates[20] = "d/dt ApopPro9 in component Pro9 (micromolar)"
    legend_rates[21] = "d/dt ApopPro9_2 in component ApopPro9_2 (micromolar)"
    legend_rates[22] = "d/dt ApopCasp9_2Pro3 in component ApopCasp9_2Pro3 (micromolar)"
    legend_rates[25] = "d/dt Casp9Pro3 in component Casp9Pro3 (micromolar)"
    legend_rates[24] = "d/dt Pro3 in component Pro3 (micromolar)"
    legend_rates[27] = "d/dt Casp9IAP in component IAP (micromolar)"
    legend_rates[28] = "d/dt ApopCasp9IAP in component IAP (micromolar)"
    legend_rates[29] = "d/dt ApopCasp9_2IAP in component IAP (micromolar)"
    legend_rates[17] = "d/dt IAP in component IAP (micromolar)"
    legend_rates[30] = "d/dt ApopCasp9 in component ApopCasp9 (micromolar)"
    legend_rates[26] = "d/dt Casp9 in component Casp9 (micromolar)"
    legend_rates[23] = "d/dt ApopCasp9_2 in component ApopCasp9_2 (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0
    states[1] = 0
    states[2] = 0.004
    constants[0] = 0.1
    constants[1] = 10
    constants[2] = 0.5
    constants[3] = 0.006
    states[3] = 0.004
    states[4] = 0
    states[5] = 0
    constants[4] = 0.0003
    constants[5] = 5
    constants[6] = 0.5
    states[6] = 0
    constants[7] = 50000
    constants[8] = 0.5
    constants[9] = 4
    states[7] = 0.004
    states[8] = 0
    constants[10] = 0.0003
    constants[11] = 10
    states[9] = 0
    states[10] = 0.004
    constants[12] = 10
    states[11] = 0
    constants[13] = 10
    states[12] = 0
    constants[14] = 10
    states[13] = 0.004
    constants[15] = 10
    constants[16] = 0.00003
    constants[17] = 0.0066
    constants[18] = 0.004
    constants[19] = 0.006
    states[14] = 0
    states[15] = 0
    constants[20] = 10
    constants[21] = 0.5
    constants[22] = 0.00008
    constants[23] = 0.006
    constants[24] = 0.1
    states[16] = 0
    constants[25] = 10
    constants[26] = 0.5
    constants[27] = 0.1
    constants[28] = 0.00003
    states[17] = 0.004
    states[18] = 0
    constants[29] = 5
    constants[30] = 0.0035
    states[19] = 0.004
    states[20] = 0
    states[21] = 0
    constants[31] = 10
    constants[32] = 0.5
    constants[33] = 10
    constants[34] = 0.5
    constants[35] = 0.0003
    constants[36] = 0.1
    states[22] = 0
    states[23] = 0
    states[24] = 0.004
    constants[37] = 10
    constants[38] = 0.5
    constants[39] = 0.1
    states[25] = 0
    constants[40] = 0.001
    states[26] = 0
    constants[41] = 10
    constants[42] = 0.5
    constants[43] = 0.0003
    states[27] = 0
    states[28] = 0
    states[29] = 0
    states[30] = 0
    constants[44] = 0.00003
    constants[45] = 5
    constants[46] = 0.0035
    constants[47] = 5
    constants[48] = 0.0035
    constants[49] = 5
    constants[50] = 0.0035
    constants[51] = 5
    constants[52] = 0.5
    constants[53] = 5
    constants[54] = 0.5
    constants[55] = constants[16]*(1.00000+(power(constants[17], 4.00000))/(power(constants[17], 4.00000)+power(constants[18], 4.00000)))
    constants[56] = (constants[22]*(power(constants[18], 4.00000)))/(power(constants[17], 4.00000)+power(constants[18], 4.00000))
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[2] = constants[1]*states[0]*states[2]-constants[2]*states[1]
    algebraic[0] = constants[0]*states[1]
    algebraic[1] = -constants[3]*states[0]
    rates[0] = -algebraic[2]+algebraic[0]+algebraic[1]
    rates[1] = algebraic[2]-algebraic[0]
    algebraic[3] = constants[4]-constants[3]*states[3]
    algebraic[4] = constants[5]*states[5]*states[3]-constants[6]*states[4]
    rates[3] = -algebraic[4]+algebraic[3]
    algebraic[6] = constants[7]*(power(states[4], constants[9]))-constants[8]*states[6]
    rates[4] = algebraic[4]-7.00000*algebraic[6]
    algebraic[9] = constants[11]*states[8]*states[7]
    algebraic[5] = -constants[3]*states[5]
    rates[5] = (algebraic[9]-algebraic[4])+algebraic[5]
    algebraic[7] = constants[10]-constants[3]*states[7]
    rates[7] = -algebraic[9]+algebraic[7]
    algebraic[12] = constants[13]*states[11]
    algebraic[11] = constants[12]*states[9]*states[10]
    algebraic[10] = -constants[3]*states[9]
    rates[9] = (algebraic[12]-algebraic[11])+algebraic[10]
    algebraic[15] = constants[14]*states[12]*states[10]
    algebraic[8] = -constants[3]*states[8]
    rates[8] = algebraic[15]+algebraic[8]
    algebraic[14] = -constants[3]*states[12]
    rates[12] = (algebraic[11]-algebraic[15])+algebraic[14]
    algebraic[17] = constants[15]*states[13]*states[10]
    algebraic[16] = constants[55]-constants[19]*states[10]
    rates[10] = ((-algebraic[11]-algebraic[15])-algebraic[17])+algebraic[16]
    algebraic[18] = constants[20]*states[14]*states[13]-constants[21]*states[15]
    algebraic[20] = constants[56]-constants[23]*states[13]
    rates[13] = (-algebraic[18]-algebraic[17])+algebraic[20]
    algebraic[19] = constants[24]*states[15]
    rates[15] = algebraic[18]-algebraic[19]
    algebraic[23] = constants[27]*states[16]
    algebraic[13] = -constants[3]*states[11]
    rates[11] = ((algebraic[0]+algebraic[23])-algebraic[12])+algebraic[15]+algebraic[13]
    algebraic[21] = constants[25]*states[14]*states[2]-constants[26]*states[16]
    rates[16] = algebraic[21]-algebraic[23]
    algebraic[22] = constants[28]-constants[3]*states[2]
    rates[2] = (-algebraic[2]-algebraic[21])+algebraic[22]
    algebraic[25] = constants[29]*states[14]*states[17]-constants[30]*states[18]
    rates[18] = algebraic[25]
    algebraic[27] = constants[31]*states[6]*states[19]-constants[32]*states[20]
    algebraic[28] = constants[33]*states[20]*states[19]-constants[34]*states[21]
    algebraic[26] = constants[35]-constants[3]*states[19]
    rates[19] = (-algebraic[27]-algebraic[28])+algebraic[26]
    rates[20] = algebraic[27]-algebraic[28]
    algebraic[29] = constants[36]*states[21]
    rates[21] = algebraic[28]-algebraic[29]
    algebraic[31] = constants[39]*states[22]
    algebraic[30] = constants[37]*states[23]*states[24]-constants[38]*states[22]
    rates[22] = algebraic[30]-algebraic[31]
    algebraic[32] = constants[40]*states[25]
    algebraic[24] = -constants[3]*states[14]
    rates[14] = (((((algebraic[32]+algebraic[31])-algebraic[25])-algebraic[21])+algebraic[23])-algebraic[18])+algebraic[19]+algebraic[24]
    algebraic[33] = constants[41]*states[26]*states[24]-constants[42]*states[25]
    rates[25] = algebraic[33]-algebraic[32]
    algebraic[34] = constants[43]-constants[3]*states[24]
    rates[24] = (-algebraic[33]-algebraic[30])+algebraic[34]
    algebraic[35] = constants[45]*states[26]*states[17]-constants[46]*states[27]
    rates[27] = algebraic[35]
    algebraic[36] = constants[47]*states[30]*states[17]-constants[48]*states[28]
    rates[28] = algebraic[36]
    algebraic[37] = constants[49]*states[23]*states[17]-constants[50]*states[29]
    rates[29] = algebraic[37]
    algebraic[39] = constants[44]-constants[3]*states[17]
    rates[17] = (((-algebraic[35]-algebraic[36])-algebraic[37])-algebraic[25])+algebraic[39]
    algebraic[38] = constants[51]*states[23]-constants[52]*states[30]*states[26]
    rates[23] = (((algebraic[29]-algebraic[38])-algebraic[37])-algebraic[30])+algebraic[31]
    algebraic[40] = constants[53]*states[30]-constants[54]*states[6]*states[26]
    rates[6] = (algebraic[6]-algebraic[27])+algebraic[40]
    rates[30] = (algebraic[38]-algebraic[40])-algebraic[36]
    algebraic[41] = -constants[3]*states[26]
    rates[26] = (((algebraic[38]+algebraic[40])-algebraic[35])-algebraic[33])+algebraic[32]+algebraic[41]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[2] = constants[1]*states[0]*states[2]-constants[2]*states[1]
    algebraic[0] = constants[0]*states[1]
    algebraic[1] = -constants[3]*states[0]
    algebraic[3] = constants[4]-constants[3]*states[3]
    algebraic[4] = constants[5]*states[5]*states[3]-constants[6]*states[4]
    algebraic[6] = constants[7]*(power(states[4], constants[9]))-constants[8]*states[6]
    algebraic[9] = constants[11]*states[8]*states[7]
    algebraic[5] = -constants[3]*states[5]
    algebraic[7] = constants[10]-constants[3]*states[7]
    algebraic[12] = constants[13]*states[11]
    algebraic[11] = constants[12]*states[9]*states[10]
    algebraic[10] = -constants[3]*states[9]
    algebraic[15] = constants[14]*states[12]*states[10]
    algebraic[8] = -constants[3]*states[8]
    algebraic[14] = -constants[3]*states[12]
    algebraic[17] = constants[15]*states[13]*states[10]
    algebraic[16] = constants[55]-constants[19]*states[10]
    algebraic[18] = constants[20]*states[14]*states[13]-constants[21]*states[15]
    algebraic[20] = constants[56]-constants[23]*states[13]
    algebraic[19] = constants[24]*states[15]
    algebraic[23] = constants[27]*states[16]
    algebraic[13] = -constants[3]*states[11]
    algebraic[21] = constants[25]*states[14]*states[2]-constants[26]*states[16]
    algebraic[22] = constants[28]-constants[3]*states[2]
    algebraic[25] = constants[29]*states[14]*states[17]-constants[30]*states[18]
    algebraic[27] = constants[31]*states[6]*states[19]-constants[32]*states[20]
    algebraic[28] = constants[33]*states[20]*states[19]-constants[34]*states[21]
    algebraic[26] = constants[35]-constants[3]*states[19]
    algebraic[29] = constants[36]*states[21]
    algebraic[31] = constants[39]*states[22]
    algebraic[30] = constants[37]*states[23]*states[24]-constants[38]*states[22]
    algebraic[32] = constants[40]*states[25]
    algebraic[24] = -constants[3]*states[14]
    algebraic[33] = constants[41]*states[26]*states[24]-constants[42]*states[25]
    algebraic[34] = constants[43]-constants[3]*states[24]
    algebraic[35] = constants[45]*states[26]*states[17]-constants[46]*states[27]
    algebraic[36] = constants[47]*states[30]*states[17]-constants[48]*states[28]
    algebraic[37] = constants[49]*states[23]*states[17]-constants[50]*states[29]
    algebraic[39] = constants[44]-constants[3]*states[17]
    algebraic[38] = constants[51]*states[23]-constants[52]*states[30]*states[26]
    algebraic[40] = constants[53]*states[30]-constants[54]*states[6]*states[26]
    algebraic[41] = -constants[3]*states[26]
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