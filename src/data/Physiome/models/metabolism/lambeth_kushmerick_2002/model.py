# Size of variable arrays:
sizeAlgebraic = 21
sizeStates = 21
sizeConstants = 116
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "GLY in component GLY (millimolar)"
    legend_algebraic[3] = "flux_GP in component flux_GP (flux)"
    legend_states[1] = "G1P in component G1P (millimolar)"
    legend_algebraic[4] = "V_PGLM in component V_PGLM (flux)"
    legend_states[2] = "G6P in component G6P (millimolar)"
    legend_algebraic[5] = "V_PGI in component V_PGI (flux)"
    legend_states[3] = "F6P in component F6P (millimolar)"
    legend_algebraic[9] = "V_PFK in component V_PFK (flux)"
    legend_states[4] = "FBP in component FBP (millimolar)"
    legend_algebraic[10] = "V_ALD in component V_ALD (flux)"
    legend_states[5] = "DHAP in component DHAP (millimolar)"
    legend_algebraic[11] = "V_TPI in component V_TPI (flux)"
    legend_states[6] = "GAP in component GAP (millimolar)"
    legend_algebraic[13] = "V_GAPDH in component V_GAPDH (flux)"
    legend_states[7] = "Thirteen_BPG in component Thirteen_BPG (millimolar)"
    legend_algebraic[14] = "V_PGK in component V_PGK (flux)"
    legend_states[8] = "three_PG in component three_PG (millimolar)"
    legend_algebraic[15] = "V_PGM in component V_PGM (flux)"
    legend_states[9] = "two_PG in component two_PG (millimolar)"
    legend_algebraic[16] = "V_ENOL in component V_ENOL (flux)"
    legend_states[10] = "PEP in component PEP (millimolar)"
    legend_algebraic[17] = "V_PK in component V_PK (flux)"
    legend_states[11] = "PYR in component PYR (millimolar)"
    legend_algebraic[18] = "V_LDH in component V_LDH (flux)"
    legend_states[12] = "LAC in component LAC (millimolar)"
    legend_algebraic[0] = "output in component LAC (flux)"
    legend_states[13] = "Pi in component Pi (millimolar)"
    legend_constants[0] = "V_ATPase in component V_ATPase (flux)"
    legend_states[14] = "ADP in component ADP (millimolar)"
    legend_algebraic[20] = "V_ADK in component V_ADK (flux)"
    legend_algebraic[19] = "V_CK in component V_CK (flux)"
    legend_states[15] = "ATP in component ATP (millimolar)"
    legend_states[16] = "AMP in component AMP (millimolar)"
    legend_states[17] = "PCr in component PCr (millimolar)"
    legend_states[18] = "Cr in component Cr (millimolar)"
    legend_states[19] = "NADH in component NADH (millimolar)"
    legend_states[20] = "NAD in component NAD (millimolar)"
    legend_constants[1] = "frac_a in component flux_GP (dimensionless)"
    legend_constants[2] = "frac_b in component flux_GP (dimensionless)"
    legend_algebraic[1] = "V_GPa in component V_GPa (flux)"
    legend_algebraic[2] = "V_GPb in component V_GPb (flux)"
    legend_constants[3] = "Ki_GLY_f in component V_GPa (millimolar)"
    legend_constants[4] = "Ki_GLY_b in component V_GPa (millimolar)"
    legend_constants[5] = "K_GLY_f in component V_GPa (millimolar)"
    legend_constants[6] = "K_GLY_b in component V_GPa (millimolar)"
    legend_constants[7] = "Keq_GP in component V_GPa (dimensionless)"
    legend_constants[8] = "Ki_Pi in component V_GPa (millimolar)"
    legend_constants[9] = "K_Pi in component V_GPa (millimolar)"
    legend_constants[10] = "Ki_G1P in component V_GPa (millimolar)"
    legend_constants[100] = "V_max_r in component V_GPa (flux)"
    legend_constants[11] = "V_max_f in component V_GPa (flux)"
    legend_constants[12] = "Ki_GLY_f in component V_GPb (millimolar)"
    legend_constants[13] = "Ki_GLY_b in component V_GPb (millimolar)"
    legend_constants[14] = "K_GLY_f in component V_GPb (millimolar)"
    legend_constants[15] = "Keq_GP in component V_GPb (dimensionless)"
    legend_constants[16] = "nH in component V_GPb (dimensionless)"
    legend_constants[17] = "K_AMP in component V_GPb (millimolar_1_75)"
    legend_constants[18] = "K_GLY_b in component V_GPb (millimolar)"
    legend_constants[19] = "Ki_Pi in component V_GPb (millimolar)"
    legend_constants[20] = "K_Pi in component V_GPb (millimolar)"
    legend_constants[21] = "Ki_G1P in component V_GPb (millimolar)"
    legend_constants[22] = "K_G1P in component V_GPb (millimolar)"
    legend_constants[101] = "V_max_r in component V_GPb (flux)"
    legend_constants[23] = "V_max_f in component V_GPb (flux)"
    legend_constants[24] = "K_G1P in component V_PGLM (millimolar)"
    legend_constants[25] = "K_G6P in component V_PGLM (millimolar)"
    legend_constants[26] = "Keq_PGLM in component V_PGLM (dimensionless)"
    legend_constants[102] = "V_max_r in component V_PGLM (flux)"
    legend_constants[27] = "V_max_f in component V_PGLM (flux)"
    legend_constants[28] = "K_G6P in component V_PGI (millimolar)"
    legend_constants[29] = "K_F6P in component V_PGI (millimolar)"
    legend_constants[30] = "Keq_PGI in component V_PGI (dimensionless)"
    legend_constants[31] = "V_max_r in component V_PGI (flux)"
    legend_constants[103] = "V_max_f in component V_PGI (flux)"
    legend_constants[32] = "K_FBP in component V_PFK (millimolar)"
    legend_constants[33] = "K_FBP_ in component V_PFK (millimolar)"
    legend_constants[34] = "K_F6P in component V_PFK (millimolar)"
    legend_constants[35] = "K_ATP in component V_PFK (millimolar)"
    legend_constants[36] = "K_F6P_ in component V_PFK (millimolar)"
    legend_constants[37] = "K_ATP_ in component V_PFK (millimolar)"
    legend_constants[38] = "K_ADP in component V_PFK (millimolar)"
    legend_constants[39] = "K_ADP_ in component V_PFK (millimolar)"
    legend_constants[40] = "Ki_ATP in component V_PFK (millimolar)"
    legend_constants[41] = "Ka_AMP in component V_PFK (millimolar)"
    legend_constants[104] = "V_max_r in component V_PFK (flux)"
    legend_constants[42] = "V_max_f in component V_PFK (flux)"
    legend_constants[105] = "alpha in component V_PFK (dimensionless)"
    legend_algebraic[8] = "L in component V_PFK (dimensionless)"
    legend_constants[43] = "Lo in component V_PFK (dimensionless)"
    legend_algebraic[6] = "delta in component V_PFK (dimensionless)"
    legend_algebraic[7] = "delta_ in component V_PFK (dimensionless)"
    legend_constants[44] = "d in component V_PFK (dimensionless)"
    legend_constants[45] = "var_e in component V_PFK (dimensionless)"
    legend_constants[46] = "K_GAP in component V_ALD (millimolar)"
    legend_constants[47] = "K_FBP in component V_ALD (millimolar)"
    legend_constants[48] = "K_DHAP in component V_ALD (millimolar)"
    legend_constants[49] = "Keq_ALD in component V_ALD (millimolar)"
    legend_constants[106] = "V_max_r in component V_ALD (flux)"
    legend_constants[50] = "V_max_f in component V_ALD (flux)"
    legend_constants[51] = "K_GAP in component V_TPI (millimolar)"
    legend_constants[52] = "K_DHAP in component V_TPI (millimolar)"
    legend_constants[53] = "Keq_TPI in component V_TPI (dimensionless)"
    legend_constants[107] = "V_max_r in component V_TPI (flux)"
    legend_constants[54] = "V_max_f in component V_TPI (flux)"
    legend_algebraic[12] = "D_GAPDH in component V_GAPDH (dimensionless)"
    legend_constants[55] = "K_GAP in component V_GAPDH (millimolar)"
    legend_constants[56] = "K_Pi in component V_GAPDH (millimolar)"
    legend_constants[57] = "K_NAD in component V_GAPDH (millimolar)"
    legend_constants[58] = "K_NADH in component V_GAPDH (millimolar)"
    legend_constants[59] = "K_Thirteen_BPG in component V_GAPDH (millimolar)"
    legend_constants[60] = "Keq_GAPDH in component V_GAPDH (dimensionless)"
    legend_constants[108] = "V_max_r in component V_GAPDH (flux)"
    legend_constants[61] = "V_max_f in component V_GAPDH (flux)"
    legend_constants[62] = "K_three_PG in component V_PGK (millimolar)"
    legend_constants[63] = "K_ADP in component V_PGK (millimolar)"
    legend_constants[64] = "K_ATP in component V_PGK (millimolar)"
    legend_constants[65] = "K_Thirteen_BPG in component V_PGK (millimolar)"
    legend_constants[66] = "Keq_PGK in component V_PGK (dimensionless)"
    legend_constants[67] = "V_max_r in component V_PGK (flux)"
    legend_constants[109] = "V_max_f in component V_PGK (flux)"
    legend_constants[68] = "K_two_PG in component V_PGM (millimolar)"
    legend_constants[69] = "K_three_PG in component V_PGM (millimolar)"
    legend_constants[70] = "Keq_PGM in component V_PGM (dimensionless)"
    legend_constants[110] = "V_max_r in component V_PGM (flux)"
    legend_constants[71] = "V_max_f in component V_PGM (flux)"
    legend_constants[72] = "K_PEP in component V_ENOL (millimolar)"
    legend_constants[73] = "K_two_PG in component V_ENOL (millimolar)"
    legend_constants[74] = "Keq_ENOL in component V_ENOL (dimensionless)"
    legend_constants[111] = "V_max_r in component V_ENOL (flux)"
    legend_constants[75] = "V_max_f in component V_ENOL (flux)"
    legend_constants[76] = "K_PYR in component V_PK (millimolar)"
    legend_constants[77] = "K_ADP in component V_PK (millimolar)"
    legend_constants[78] = "K_ATP in component V_PK (millimolar)"
    legend_constants[79] = "K_PEP in component V_PK (millimolar)"
    legend_constants[80] = "Keq_PK in component V_PK (dimensionless)"
    legend_constants[112] = "V_max_r in component V_PK (flux)"
    legend_constants[81] = "V_max_f in component V_PK (flux)"
    legend_constants[82] = "K_LAC in component V_LDH (millimolar)"
    legend_constants[83] = "K_NADH in component V_LDH (millimolar)"
    legend_constants[84] = "K_NAD in component V_LDH (millimolar)"
    legend_constants[85] = "K_PYR in component V_LDH (millimolar)"
    legend_constants[86] = "Keq_LDH in component V_LDH (dimensionless)"
    legend_constants[113] = "V_max_r in component V_LDH (flux)"
    legend_constants[87] = "V_max_f in component V_LDH (flux)"
    legend_constants[88] = "Ki_ADP in component V_CK (millimolar)"
    legend_constants[89] = "K_Cr in component V_CK (millimolar)"
    legend_constants[90] = "K_PCr in component V_CK (millimolar)"
    legend_constants[91] = "Ki_PCr in component V_CK (millimolar)"
    legend_constants[92] = "Ki_ATP in component V_CK (millimolar)"
    legend_constants[93] = "Keq_CK in component V_CK (dimensionless)"
    legend_constants[94] = "V_max_r in component V_CK (flux)"
    legend_constants[114] = "V_max_f in component V_CK (flux)"
    legend_constants[95] = "K_ADP in component V_ADK (millimolar)"
    legend_constants[96] = "K_AMP in component V_ADK (millimolar)"
    legend_constants[97] = "K_ATP in component V_ADK (millimolar)"
    legend_constants[98] = "Keq_ADK in component V_ADK (dimensionless)"
    legend_constants[115] = "V_max_r in component V_ADK (flux)"
    legend_constants[99] = "V_max_f in component V_ADK (flux)"
    legend_rates[0] = "d/dt GLY in component GLY (millimolar)"
    legend_rates[1] = "d/dt G1P in component G1P (millimolar)"
    legend_rates[2] = "d/dt G6P in component G6P (millimolar)"
    legend_rates[3] = "d/dt F6P in component F6P (millimolar)"
    legend_rates[4] = "d/dt FBP in component FBP (millimolar)"
    legend_rates[5] = "d/dt DHAP in component DHAP (millimolar)"
    legend_rates[6] = "d/dt GAP in component GAP (millimolar)"
    legend_rates[7] = "d/dt Thirteen_BPG in component Thirteen_BPG (millimolar)"
    legend_rates[8] = "d/dt three_PG in component three_PG (millimolar)"
    legend_rates[9] = "d/dt two_PG in component two_PG (millimolar)"
    legend_rates[10] = "d/dt PEP in component PEP (millimolar)"
    legend_rates[11] = "d/dt PYR in component PYR (millimolar)"
    legend_rates[12] = "d/dt LAC in component LAC (millimolar)"
    legend_rates[13] = "d/dt Pi in component Pi (millimolar)"
    legend_rates[14] = "d/dt ADP in component ADP (millimolar)"
    legend_rates[15] = "d/dt ATP in component ATP (millimolar)"
    legend_rates[16] = "d/dt AMP in component AMP (millimolar)"
    legend_rates[17] = "d/dt PCr in component PCr (millimolar)"
    legend_rates[18] = "d/dt Cr in component Cr (millimolar)"
    legend_rates[19] = "d/dt NADH in component NADH (millimolar)"
    legend_rates[20] = "d/dt NAD in component NAD (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 112.0
    states[1] = 0.0589
    states[2] = 0.75
    states[3] = 0.228
    states[4] = 0.0723
    states[5] = 0.0764
    states[6] = 0.0355
    states[7] = 0.065
    states[8] = 0.052
    states[9] = 0.005
    states[10] = 0.0194
    states[11] = 0.0994
    states[12] = 1.3
    states[13] = 4.1
    constants[0] = 600.0
    states[14] = 0.013
    states[15] = 8.2
    states[16] = 0.99
    states[17] = 34.67
    states[18] = 35.0
    states[19] = 0.001
    states[20] = 0.5
    constants[1] = 0.4
    constants[2] = 0.6
    constants[3] = 2.0
    constants[4] = 2.0
    constants[5] = 1.7
    constants[6] = 0.15
    constants[7] = 0.42
    constants[8] = 4.7
    constants[9] = 4.0
    constants[10] = 10.1
    constants[11] = 20.0
    constants[12] = 15.0
    constants[13] = 4.4
    constants[14] = 1.7
    constants[15] = 0.42
    constants[16] = 1.75
    constants[17] = 1.9E-6
    constants[18] = 0.15
    constants[19] = 4.6
    constants[20] = 0.2
    constants[21] = 7.4
    constants[22] = 1.5
    constants[23] = 30.0
    constants[24] = 0.063
    constants[25] = 0.03
    constants[26] = 16.62
    constants[27] = 480.0
    constants[28] = 0.18
    constants[29] = 0.119
    constants[30] = 0.45
    constants[31] = 880.0
    constants[32] = 4.02
    constants[33] = 4.02
    constants[34] = 0.18
    constants[35] = 0.08
    constants[36] = 20.0
    constants[37] = 0.25
    constants[38] = 2.7
    constants[39] = 2.7
    constants[40] = 0.87
    constants[41] = 0.06
    constants[42] = 56.0
    constants[43] = 13.0
    constants[44] = 0.01
    constants[45] = 0.01
    constants[46] = 1.0
    constants[47] = 0.05
    constants[48] = 2.0
    constants[49] = 9.5E-5
    constants[50] = 104.0
    constants[51] = 0.32
    constants[52] = 0.61
    constants[53] = 0.052
    constants[54] = 12000.0
    constants[55] = 0.0025
    constants[56] = 0.29
    constants[57] = 0.09
    constants[58] = 0.0033
    constants[59] = 0.0008
    constants[60] = 0.089
    constants[61] = 1265.0
    constants[62] = 1.2
    constants[63] = 0.008
    constants[64] = 0.35
    constants[65] = 0.002
    constants[66] = 57109.0
    constants[67] = 1120.0
    constants[68] = 0.014
    constants[69] = 0.2
    constants[70] = 0.18
    constants[71] = 1120.0
    constants[72] = 0.37
    constants[73] = 0.1
    constants[74] = 0.49
    constants[75] = 192.0
    constants[76] = 7.05
    constants[77] = 0.3
    constants[78] = 1.13
    constants[79] = 0.08
    constants[80] = 10304.0
    constants[81] = 1440.0
    constants[82] = 17.0
    constants[83] = 0.002
    constants[84] = 0.849
    constants[85] = 0.335
    constants[86] = 16198.0
    constants[87] = 1920.0
    constants[88] = 0.135
    constants[89] = 3.8
    constants[90] = 1.11
    constants[91] = 3.9
    constants[92] = 3.5
    constants[93] = 233.0
    constants[94] = 500.0
    constants[95] = 0.35
    constants[96] = 0.32
    constants[97] = 0.27
    constants[98] = 2.21
    constants[99] = 880.0
    constants[100] = (constants[11]*constants[6]*constants[10])/(constants[3]*constants[9]*constants[7])
    constants[101] = (constants[23]*constants[13]*constants[22])/(constants[14]*constants[19]*constants[15])
    constants[102] = (constants[27]*constants[25])/(constants[24]*constants[26])
    constants[103] = (constants[31]*constants[28]*constants[30])/constants[29]
    constants[104] = (constants[42]*constants[38]*constants[32])/(constants[35]*constants[34])
    constants[105] = (constants[34]*constants[35])/(constants[36]*constants[37])
    constants[106] = (constants[50]*constants[48]*constants[46])/(constants[47]*constants[49])
    constants[107] = (constants[54]*constants[52])/(constants[51]*constants[53])
    constants[108] = (constants[61]*1.00000*constants[59]*constants[58])/(constants[55]*constants[57]*constants[56]*constants[60])
    constants[109] = (constants[67]*constants[65]*constants[63]*constants[66])/(constants[62]*constants[64])
    constants[110] = (constants[71]*constants[68])/(constants[69]*constants[70])
    constants[111] = (constants[75]*constants[72])/(constants[73]*constants[74])
    constants[112] = (constants[81]*constants[78]*constants[76])/(constants[79]*constants[77]*constants[80])
    constants[113] = (constants[87]*constants[82]*constants[84])/(constants[85]*constants[83]*constants[86])
    constants[114] = (constants[94]*constants[92]*constants[89]*constants[93])/(constants[88]*constants[90])
    constants[115] = (constants[99]*(power(constants[95], 2.00000)))/(constants[97]*constants[96]*constants[98])
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = (constants[11]*((states[0]*states[13])/(constants[3]*constants[9]))-constants[100]*((states[0]*states[1])/(constants[4]*constants[10])))/(1.00000+states[0]/constants[3]+states[13]/constants[8]+states[0]/constants[4]+states[1]/constants[10]+(states[0]*states[13])/(constants[5]*constants[8])+(states[0]*states[1])/(constants[6]*constants[10]))
    algebraic[2] = ((constants[23]*((states[0]*states[13])/(constants[12]*constants[20]))-constants[101]*((states[0]*states[1])/(constants[13]*constants[21])))/(1.00000+states[0]/constants[12]+states[13]/constants[19]+states[0]/constants[13]+states[1]/constants[21]+(states[0]*states[13])/(constants[14]*constants[19])+(states[0]*states[1])/(constants[18]*constants[21])))*(((power(states[16], constants[16]))/constants[17])/(1.00000+(power(states[16], constants[16]))/constants[17]))
    algebraic[3] = constants[1]*algebraic[1]+constants[2]*algebraic[2]
    rates[0] = -algebraic[3]
    algebraic[4] = (constants[27]*(states[1]/constants[24])-constants[102]*(states[2]/constants[25]))/(1.00000+states[1]/constants[24]+states[2]/constants[25])
    rates[1] = algebraic[3]-algebraic[4]
    algebraic[5] = (constants[103]*(states[2]/constants[28])-constants[31]*(states[3]/constants[29]))/(1.00000+states[2]/constants[28]+states[3]/constants[29])
    rates[2] = algebraic[4]-algebraic[5]
    algebraic[8] = constants[43]*(power(((1.00000+states[15]/constants[40])/(1.00000+constants[44]*(states[15]/constants[40])))*((1.00000+constants[45]*(states[16]/constants[41]))/(1.00000+states[16]/constants[41])), 4.00000))
    algebraic[6] = (1.00000+states[3]/constants[34])*(1.00000+states[15]/constants[35])+states[14]/constants[38]+(states[4]/constants[32])*(1.00000+states[14]/constants[38])
    algebraic[7] = (1.00000+states[3]/constants[36])*(1.00000+states[15]/constants[37])+states[14]/constants[39]+(states[4]/constants[33])*(1.00000+states[14]/constants[39])
    algebraic[9] = ((constants[42]*((states[15]*states[3])/(constants[35]*constants[34]))-constants[104]*((states[14]*states[4])/(constants[38]*constants[32])))/algebraic[6])*((1.00000+constants[105]*algebraic[8]*(power(algebraic[7]/algebraic[6], 3.00000)))/(1.00000+algebraic[8]*(power(algebraic[7]/algebraic[6], 4.00000))))
    rates[3] = algebraic[5]-algebraic[9]
    algebraic[10] = (constants[50]*(states[4]/constants[47])-constants[106]*((states[5]*states[6])/(constants[48]*constants[46])))/(1.00000+states[4]/constants[47]+states[5]/constants[48]+states[6]/constants[46])
    rates[4] = algebraic[9]-algebraic[10]
    algebraic[11] = (constants[54]*(states[6]/constants[51])-constants[107]*(states[5]/constants[52]))/(1.00000+states[5]/constants[52]+states[6]/constants[51])
    rates[5] = algebraic[10]+algebraic[11]
    algebraic[12] = 1.00000+states[6]/constants[55]+states[20]/constants[57]+states[13]/constants[56]+(states[6]*states[20])/(constants[55]*constants[57])+(states[6]*states[20]*states[13])/(constants[55]*constants[57]*constants[56])+states[7]/constants[59]+states[19]/constants[58]+(states[7]*states[19])/(constants[59]*constants[58])
    algebraic[13] = (constants[61]*((states[6]*states[20]*states[13])/(constants[55]*constants[57]*constants[56]))-constants[108]*((states[7]*states[19])/(constants[59]*constants[58])))/algebraic[12]
    rates[6] = algebraic[10]-(algebraic[11]+algebraic[13])
    rates[13] = constants[0]-(algebraic[3]+algebraic[13])
    algebraic[14] = (constants[109]*((states[7]*states[14])/(constants[65]*constants[63]))-constants[67]*((states[8]*states[15])/(constants[62]*constants[64])))/(1.00000+states[7]/constants[65]+states[14]/constants[63]+(states[7]*states[14])/(constants[65]*constants[63])+states[8]/constants[62]+states[15]/constants[64]+(states[8]*states[15])/(constants[62]*constants[64]))
    rates[7] = algebraic[13]-algebraic[14]
    algebraic[15] = (constants[71]*(states[8]/constants[69])-constants[110]*(states[9]/constants[68]))/(1.00000+states[8]/constants[69]+states[9]/constants[68])
    rates[8] = algebraic[14]-algebraic[15]
    algebraic[16] = (constants[75]*(states[9]/constants[73])-constants[111]*(states[10]/constants[72]))/(1.00000+states[9]/constants[73]+states[10]/constants[72])
    rates[9] = algebraic[15]-algebraic[16]
    algebraic[17] = (constants[81]*((states[10]*states[14])/(constants[79]*constants[77]))-constants[112]*((states[11]*states[15])/(constants[76]*constants[78])))/(1.00000+states[10]/constants[79]+states[14]/constants[77]+(states[10]*states[14])/(constants[79]*constants[77])+states[11]/constants[76]+states[15]/constants[78]+(states[11]*states[15])/(constants[76]*constants[78]))
    rates[10] = algebraic[16]-algebraic[17]
    algebraic[18] = (constants[87]*((states[11]*states[19])/(constants[85]*constants[83]))-constants[113]*((states[12]*states[20])/(constants[82]*constants[84])))/(1.00000+states[11]/constants[85]+states[19]/constants[83]+(states[11]*states[19])/(constants[85]*constants[83])+states[12]/constants[82]+states[20]/constants[84]+(states[12]*states[20])/(constants[82]*constants[84]))
    rates[11] = algebraic[17]-algebraic[18]
    algebraic[0] = 0.200000*states[12]
    rates[12] = algebraic[18]-algebraic[0]
    algebraic[19] = (constants[114]*((states[15]*states[18])/(constants[92]*constants[89]))-constants[94]*((states[14]*states[17])/(constants[88]*constants[90])))/(1.00000+states[14]/constants[88]+states[17]/constants[91]+(states[14]*states[17])/(constants[88]*constants[90])+states[15]/constants[92]+(states[15]*states[18])/(constants[92]*constants[89]))
    rates[17] = algebraic[19]
    rates[18] = -algebraic[19]
    rates[19] = algebraic[13]-algebraic[18]
    rates[20] = algebraic[18]-algebraic[13]
    algebraic[20] = (constants[99]*((states[15]*states[16])/(constants[97]*constants[96]))-constants[115]*((power(states[14], 2.00000))/(power(constants[95], 2.00000))))/(1.00000+states[15]/constants[97]+states[16]/constants[96]+(states[15]*states[16])/(constants[97]*constants[96])+(2.00000*states[14])/constants[95]+(power(states[14], 2.00000))/(power(constants[95], 2.00000)))
    rates[14] = algebraic[9]+2.00000*algebraic[20]+algebraic[19]+constants[0]+-(algebraic[14]+algebraic[17])
    rates[15] = algebraic[14]+algebraic[17]+-(algebraic[9]+algebraic[20]+algebraic[19]+constants[0])
    rates[16] = -algebraic[20]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = (constants[11]*((states[0]*states[13])/(constants[3]*constants[9]))-constants[100]*((states[0]*states[1])/(constants[4]*constants[10])))/(1.00000+states[0]/constants[3]+states[13]/constants[8]+states[0]/constants[4]+states[1]/constants[10]+(states[0]*states[13])/(constants[5]*constants[8])+(states[0]*states[1])/(constants[6]*constants[10]))
    algebraic[2] = ((constants[23]*((states[0]*states[13])/(constants[12]*constants[20]))-constants[101]*((states[0]*states[1])/(constants[13]*constants[21])))/(1.00000+states[0]/constants[12]+states[13]/constants[19]+states[0]/constants[13]+states[1]/constants[21]+(states[0]*states[13])/(constants[14]*constants[19])+(states[0]*states[1])/(constants[18]*constants[21])))*(((power(states[16], constants[16]))/constants[17])/(1.00000+(power(states[16], constants[16]))/constants[17]))
    algebraic[3] = constants[1]*algebraic[1]+constants[2]*algebraic[2]
    algebraic[4] = (constants[27]*(states[1]/constants[24])-constants[102]*(states[2]/constants[25]))/(1.00000+states[1]/constants[24]+states[2]/constants[25])
    algebraic[5] = (constants[103]*(states[2]/constants[28])-constants[31]*(states[3]/constants[29]))/(1.00000+states[2]/constants[28]+states[3]/constants[29])
    algebraic[8] = constants[43]*(power(((1.00000+states[15]/constants[40])/(1.00000+constants[44]*(states[15]/constants[40])))*((1.00000+constants[45]*(states[16]/constants[41]))/(1.00000+states[16]/constants[41])), 4.00000))
    algebraic[6] = (1.00000+states[3]/constants[34])*(1.00000+states[15]/constants[35])+states[14]/constants[38]+(states[4]/constants[32])*(1.00000+states[14]/constants[38])
    algebraic[7] = (1.00000+states[3]/constants[36])*(1.00000+states[15]/constants[37])+states[14]/constants[39]+(states[4]/constants[33])*(1.00000+states[14]/constants[39])
    algebraic[9] = ((constants[42]*((states[15]*states[3])/(constants[35]*constants[34]))-constants[104]*((states[14]*states[4])/(constants[38]*constants[32])))/algebraic[6])*((1.00000+constants[105]*algebraic[8]*(power(algebraic[7]/algebraic[6], 3.00000)))/(1.00000+algebraic[8]*(power(algebraic[7]/algebraic[6], 4.00000))))
    algebraic[10] = (constants[50]*(states[4]/constants[47])-constants[106]*((states[5]*states[6])/(constants[48]*constants[46])))/(1.00000+states[4]/constants[47]+states[5]/constants[48]+states[6]/constants[46])
    algebraic[11] = (constants[54]*(states[6]/constants[51])-constants[107]*(states[5]/constants[52]))/(1.00000+states[5]/constants[52]+states[6]/constants[51])
    algebraic[12] = 1.00000+states[6]/constants[55]+states[20]/constants[57]+states[13]/constants[56]+(states[6]*states[20])/(constants[55]*constants[57])+(states[6]*states[20]*states[13])/(constants[55]*constants[57]*constants[56])+states[7]/constants[59]+states[19]/constants[58]+(states[7]*states[19])/(constants[59]*constants[58])
    algebraic[13] = (constants[61]*((states[6]*states[20]*states[13])/(constants[55]*constants[57]*constants[56]))-constants[108]*((states[7]*states[19])/(constants[59]*constants[58])))/algebraic[12]
    algebraic[14] = (constants[109]*((states[7]*states[14])/(constants[65]*constants[63]))-constants[67]*((states[8]*states[15])/(constants[62]*constants[64])))/(1.00000+states[7]/constants[65]+states[14]/constants[63]+(states[7]*states[14])/(constants[65]*constants[63])+states[8]/constants[62]+states[15]/constants[64]+(states[8]*states[15])/(constants[62]*constants[64]))
    algebraic[15] = (constants[71]*(states[8]/constants[69])-constants[110]*(states[9]/constants[68]))/(1.00000+states[8]/constants[69]+states[9]/constants[68])
    algebraic[16] = (constants[75]*(states[9]/constants[73])-constants[111]*(states[10]/constants[72]))/(1.00000+states[9]/constants[73]+states[10]/constants[72])
    algebraic[17] = (constants[81]*((states[10]*states[14])/(constants[79]*constants[77]))-constants[112]*((states[11]*states[15])/(constants[76]*constants[78])))/(1.00000+states[10]/constants[79]+states[14]/constants[77]+(states[10]*states[14])/(constants[79]*constants[77])+states[11]/constants[76]+states[15]/constants[78]+(states[11]*states[15])/(constants[76]*constants[78]))
    algebraic[18] = (constants[87]*((states[11]*states[19])/(constants[85]*constants[83]))-constants[113]*((states[12]*states[20])/(constants[82]*constants[84])))/(1.00000+states[11]/constants[85]+states[19]/constants[83]+(states[11]*states[19])/(constants[85]*constants[83])+states[12]/constants[82]+states[20]/constants[84]+(states[12]*states[20])/(constants[82]*constants[84]))
    algebraic[0] = 0.200000*states[12]
    algebraic[19] = (constants[114]*((states[15]*states[18])/(constants[92]*constants[89]))-constants[94]*((states[14]*states[17])/(constants[88]*constants[90])))/(1.00000+states[14]/constants[88]+states[17]/constants[91]+(states[14]*states[17])/(constants[88]*constants[90])+states[15]/constants[92]+(states[15]*states[18])/(constants[92]*constants[89]))
    algebraic[20] = (constants[99]*((states[15]*states[16])/(constants[97]*constants[96]))-constants[115]*((power(states[14], 2.00000))/(power(constants[95], 2.00000))))/(1.00000+states[15]/constants[97]+states[16]/constants[96]+(states[15]*states[16])/(constants[97]*constants[96])+(2.00000*states[14])/constants[95]+(power(states[14], 2.00000))/(power(constants[95], 2.00000)))
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