# Size of variable arrays:
sizeAlgebraic = 30
sizeStates = 8
sizeConstants = 81
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "PYR in component PYR (millimolar)"
    legend_algebraic[4] = "V_GLYC in component V_GLYC (flux)"
    legend_algebraic[6] = "V_LDH in component V_LDH (flux)"
    legend_algebraic[8] = "V_PDH in component V_PDH (flux)"
    legend_algebraic[20] = "V_ALS in component V_ALS (flux)"
    legend_states[1] = "ACP in component ACP (millimolar)"
    legend_algebraic[10] = "V_PTA in component V_PTA (flux)"
    legend_algebraic[14] = "V_ACK in component V_ACK (flux)"
    legend_states[2] = "ACAL in component ACAL (millimolar)"
    legend_algebraic[13] = "V_ACALDH in component V_ACALDH (flux)"
    legend_algebraic[17] = "V_ADH in component V_ADH (flux)"
    legend_states[3] = "ACLAC in component ACLAC (millimolar)"
    legend_algebraic[22] = "V_ALDC in component V_ALDC (flux)"
    legend_algebraic[27] = "V_NEALC in component V_NEALC (flux)"
    legend_states[4] = "ACET in component ACET (millimolar)"
    legend_algebraic[26] = "V_ACETDH in component V_ACETDH (flux)"
    legend_algebraic[24] = "V_ACETEFF in component V_ACETEFF (flux)"
    legend_states[5] = "ATP in component ATP (millimolar)"
    legend_algebraic[18] = "V_ATPase in component V_ATPase (flux)"
    legend_algebraic[0] = "ADP in component ADP (millimolar)"
    legend_constants[0] = "A_tot in component ADP (millimolar)"
    legend_states[6] = "NADH in component NADH (millimolar)"
    legend_algebraic[29] = "V_NOX in component V_NOX (flux)"
    legend_algebraic[1] = "NAD in component NAD (millimolar)"
    legend_constants[1] = "NAD_tot in component NAD (millimolar)"
    legend_states[7] = "ACCOA in component ACCOA (millimolar)"
    legend_algebraic[2] = "COA in component COA (millimolar)"
    legend_constants[2] = "C_tot in component COA (millimolar)"
    legend_constants[3] = "AC in component AC (millimolar)"
    legend_constants[4] = "BUT in component BUT (millimolar)"
    legend_constants[5] = "ETOH in component ETOH (millimolar)"
    legend_constants[6] = "GLC in component GLC (millimolar)"
    legend_constants[7] = "LAC in component LAC (millimolar)"
    legend_constants[8] = "O in component O (millimolar)"
    legend_constants[9] = "P in component P (millimolar)"
    legend_algebraic[3] = "V_GLYC_temp in component V_GLYC (flux)"
    legend_constants[10] = "Km_GLC in component V_GLYC (millimolar)"
    legend_constants[11] = "Km_NAD in component V_GLYC (millimolar)"
    legend_constants[12] = "Km_ADP in component V_GLYC (millimolar)"
    legend_constants[13] = "Km_PYR in component V_GLYC (millimolar)"
    legend_constants[14] = "Km_NADH in component V_GLYC (millimolar)"
    legend_constants[15] = "Km_ATP in component V_GLYC (millimolar)"
    legend_constants[16] = "V_GLYC_max in component V_GLYC (flux)"
    legend_algebraic[5] = "V_LDH_temp in component V_LDH (flux)"
    legend_constants[17] = "Keq in component V_LDH (dimensionless)"
    legend_constants[18] = "Km_LAC in component V_LDH (millimolar)"
    legend_constants[19] = "Km_NAD in component V_LDH (millimolar)"
    legend_constants[20] = "Km_PYR in component V_LDH (millimolar)"
    legend_constants[21] = "Km_NADH in component V_LDH (millimolar)"
    legend_constants[22] = "V_LDH_max in component V_LDH (flux)"
    legend_algebraic[7] = "V_PDH_temp in component V_PDH (flux)"
    legend_constants[23] = "Ki in component V_PDH (dimensionless)"
    legend_constants[24] = "Km_NAD in component V_PDH (millimolar)"
    legend_constants[25] = "Km_COA in component V_PDH (millimolar)"
    legend_constants[26] = "Km_PYR in component V_PDH (millimolar)"
    legend_constants[27] = "Km_NADH in component V_PDH (millimolar)"
    legend_constants[28] = "Km_ACCOA in component V_PDH (millimolar)"
    legend_constants[29] = "V_PDH_max in component V_PDH (flux)"
    legend_algebraic[9] = "V_PTA_temp in component V_PTA (flux)"
    legend_constants[30] = "Keq in component V_PTA (dimensionless)"
    legend_constants[31] = "Km_P in component V_PTA (millimolar)"
    legend_constants[32] = "Ki_P in component V_PTA (millimolar)"
    legend_constants[33] = "Ki_COA in component V_PTA (millimolar)"
    legend_constants[34] = "Km_ACP in component V_PTA (millimolar)"
    legend_constants[35] = "Ki_ACP in component V_PTA (millimolar)"
    legend_constants[36] = "Ki_ACCOA in component V_PTA (millimolar)"
    legend_constants[37] = "V_PTA_max in component V_PTA (flux)"
    legend_algebraic[12] = "V_ACK_temp in component V_ACK (flux)"
    legend_constants[38] = "Keq in component V_ACK (dimensionless)"
    legend_constants[39] = "Km_AC in component V_ACK (millimolar)"
    legend_constants[40] = "Km_ATP in component V_ACK (millimolar)"
    legend_constants[41] = "Km_ADP in component V_ACK (millimolar)"
    legend_constants[42] = "Km_ACP in component V_ACK (millimolar)"
    legend_constants[43] = "V_ACK_max in component V_ACK (flux)"
    legend_algebraic[11] = "V_ACALDH_temp in component V_ACALDH (flux)"
    legend_constants[44] = "Keq in component V_ACALDH (millimolar)"
    legend_constants[45] = "Km_NAD in component V_ACALDH (millimolar)"
    legend_constants[46] = "Km_NADH in component V_ACALDH (millimolar)"
    legend_constants[47] = "Km_COA in component V_ACALDH (millimolar)"
    legend_constants[48] = "Km_ACCOA in component V_ACALDH (millimolar)"
    legend_constants[49] = "Km_ACAL in component V_ACALDH (millimolar)"
    legend_constants[50] = "V_ACALDH_max in component V_ACALDH (flux)"
    legend_algebraic[15] = "V_ADH_temp in component V_ADH (flux)"
    legend_constants[51] = "Keq in component V_ADH (dimensionless)"
    legend_constants[52] = "Km_NAD in component V_ADH (millimolar)"
    legend_constants[53] = "Km_NADH in component V_ADH (millimolar)"
    legend_constants[54] = "Km_ETOH in component V_ADH (millimolar)"
    legend_constants[55] = "Km_ACAL in component V_ADH (millimolar)"
    legend_constants[56] = "V_ADH_max in component V_ADH (flux)"
    legend_algebraic[19] = "V_ALS_temp in component V_ALS (flux)"
    legend_constants[57] = "N in component V_ALS (dimensionless)"
    legend_constants[58] = "Keq in component V_ALS (dimensionless)"
    legend_constants[59] = "Km_ACLAC in component V_ALS (millimolar)"
    legend_constants[60] = "Km_PYR in component V_ALS (millimolar)"
    legend_constants[61] = "V_ALS_max in component V_ALS (flux)"
    legend_algebraic[21] = "V_ALDC_temp in component V_ALDC (flux)"
    legend_constants[62] = "Km_ACLAC in component V_ALDC (millimolar)"
    legend_constants[63] = "Km_ACET in component V_ALDC (millimolar)"
    legend_constants[64] = "V_ALDC_max in component V_ALDC (flux)"
    legend_algebraic[23] = "V_ACETEFF_temp in component V_ACETEFF (flux)"
    legend_constants[65] = "Km_ACET in component V_ACETEFF (millimolar)"
    legend_constants[66] = "V_ACETEFF_max in component V_ACETEFF (flux)"
    legend_algebraic[25] = "V_ACETDH_temp in component V_ACETDH (flux)"
    legend_constants[67] = "Keq in component V_ACETDH (dimensionless)"
    legend_constants[68] = "Km_NAD in component V_ACETDH (millimolar)"
    legend_constants[69] = "Km_NADH in component V_ACETDH (millimolar)"
    legend_constants[70] = "Km_BUT in component V_ACETDH (millimolar)"
    legend_constants[71] = "Km_ACET in component V_ACETDH (millimolar)"
    legend_constants[72] = "V_ACETDH_max in component V_ACETDH (flux)"
    legend_algebraic[16] = "V_ATPase_temp in component V_ATPase (flux)"
    legend_constants[73] = "N in component V_ATPase (dimensionless)"
    legend_constants[74] = "Km_ATP in component V_ATPase (dimensionless)"
    legend_constants[75] = "V_ATPase_max in component V_ATPase (flux)"
    legend_algebraic[28] = "V_NOX_temp in component V_NOX (flux)"
    legend_constants[76] = "Km_NAD in component V_NOX (millimolar)"
    legend_constants[77] = "Km_NADH in component V_NOX (millimolar)"
    legend_constants[78] = "Km_O in component V_NOX (millimolar)"
    legend_constants[79] = "V_NOX_max in component V_NOX (flux)"
    legend_constants[80] = "k in component V_NEALC (first_order_rate_constant)"
    legend_rates[0] = "d/dt PYR in component PYR (millimolar)"
    legend_rates[1] = "d/dt ACP in component ACP (millimolar)"
    legend_rates[2] = "d/dt ACAL in component ACAL (millimolar)"
    legend_rates[3] = "d/dt ACLAC in component ACLAC (millimolar)"
    legend_rates[4] = "d/dt ACET in component ACET (millimolar)"
    legend_rates[5] = "d/dt ATP in component ATP (millimolar)"
    legend_rates[6] = "d/dt NADH in component NADH (millimolar)"
    legend_rates[7] = "d/dt ACCOA in component ACCOA (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 1
    states[1] = 0.03145
    states[2] = 0.11
    states[3] = 1e-5
    states[4] = 1e-5
    states[5] = 0.1
    constants[0] = 5
    states[6] = 3.67
    constants[1] = 10
    states[7] = 0.11
    constants[2] = 1
    constants[3] = 0.01
    constants[4] = 0.01
    constants[5] = 0.1
    constants[6] = 15
    constants[7] = 0.1
    constants[8] = 0.2
    constants[9] = 10
    constants[10] = 0.1
    constants[11] = 0.1412
    constants[12] = 0.04699
    constants[13] = 2.5
    constants[14] = 0.08999
    constants[15] = 0.01867
    constants[16] = 2397
    constants[17] = 21120.69
    constants[18] = 100
    constants[19] = 2.4
    constants[20] = 1.5
    constants[21] = 0.08
    constants[22] = 5118
    constants[23] = 46.4159
    constants[24] = 0.4
    constants[25] = 0.014
    constants[26] = 1
    constants[27] = 0.1
    constants[28] = 0.008
    constants[29] = 259
    constants[30] = 0.0065
    constants[31] = 2.6
    constants[32] = 2.6
    constants[33] = 0.029
    constants[34] = 0.7
    constants[35] = 0.2
    constants[36] = 0.2
    constants[37] = 42
    constants[38] = 174.217
    constants[39] = 7
    constants[40] = 0.07
    constants[41] = 0.5
    constants[42] = 0.16
    constants[43] = 2700
    constants[44] = 1
    constants[45] = 0.08
    constants[46] = 0.025
    constants[47] = 0.008
    constants[48] = 0.007
    constants[49] = 10
    constants[50] = 97
    constants[51] = 12354.9
    constants[52] = 0.08
    constants[53] = 0.05
    constants[54] = 1
    constants[55] = 0.03
    constants[56] = 162
    constants[57] = 2.4
    constants[58] = 9e12
    constants[59] = 100
    constants[60] = 50
    constants[61] = 600
    constants[62] = 10
    constants[63] = 100
    constants[64] = 106
    constants[65] = 5
    constants[66] = 200
    constants[67] = 1400
    constants[68] = 0.16
    constants[69] = 0.02
    constants[70] = 2.6
    constants[71] = 0.06
    constants[72] = 105
    constants[73] = 2.58
    constants[74] = 6.196
    constants[75] = 900
    constants[76] = 1
    constants[77] = 0.041
    constants[78] = 0.2
    constants[79] = 118
    constants[80] = 0.0003
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[2] = constants[2]-states[7]
    algebraic[9] = ((constants[37]/(constants[36]*constants[31]))*(states[7]*constants[9]-(states[1]*algebraic[2])/constants[30]))/(1.00000+states[7]/constants[36]+constants[9]/constants[32]+states[1]/constants[35]+algebraic[2]/constants[33]+(states[7]*constants[9])/(constants[36]*constants[31])+(states[1]*algebraic[2])/(constants[34]*constants[33]))
    algebraic[10] = custom_piecewise([greater_equal(algebraic[9] , 0.00000), algebraic[9] , True, algebraic[9]])
    algebraic[0] = constants[0]-states[5]
    algebraic[12] = ((constants[43]/(constants[41]*constants[42]))*(states[1]*algebraic[0]-(constants[3]*states[5])/constants[38]))/((1.00000+states[1]/constants[42]+constants[3]/constants[39])*(1.00000+algebraic[0]/constants[41]+states[5]/constants[40]))
    algebraic[14] = custom_piecewise([greater_equal(algebraic[12] , 0.00000), algebraic[12] , True, algebraic[12]])
    rates[1] = algebraic[10]-algebraic[14]
    algebraic[1] = constants[1]-states[6]
    algebraic[7] = (((((((constants[29]/(1.00000+(constants[23]*states[6])/algebraic[1]))*states[0])/constants[26])*algebraic[1])/constants[24])*algebraic[2])/constants[25])/((1.00000+states[0]/constants[26])*(1.00000+algebraic[1]/constants[24]+states[6]/constants[27])*(1.00000+algebraic[2]/constants[25]+states[7]/constants[28]))
    algebraic[8] = custom_piecewise([greater_equal(algebraic[7] , 0.00000), algebraic[7] , True, algebraic[7]])
    algebraic[11] = ((constants[50]/(constants[48]*constants[46]))*(states[7]*states[6]-(algebraic[1]*algebraic[2]*states[2])/constants[44]))/((1.00000+algebraic[1]/constants[45]+states[6]/constants[46])*(1.00000+states[7]/constants[48]+algebraic[2]/constants[47])*(1.00000+states[2]/constants[49]))
    algebraic[13] = custom_piecewise([greater_equal(algebraic[11] , 0.00000), algebraic[11] , True, algebraic[11]])
    rates[7] = algebraic[8]-(algebraic[13]+algebraic[10])
    algebraic[15] = ((constants[56]/(constants[55]*constants[53]))*(states[2]*states[6]-(constants[5]*algebraic[1])/constants[51]))/((1.00000+algebraic[1]/constants[52]+states[6]/constants[53])*(1.00000+states[2]/constants[55]+constants[5]/constants[54]))
    algebraic[17] = custom_piecewise([greater_equal(algebraic[15] , 0.00000), algebraic[15] , True, algebraic[15]])
    rates[2] = algebraic[13]-algebraic[17]
    algebraic[3] = ((((((constants[16]*constants[6])/constants[10])*algebraic[1])/constants[11])*algebraic[0])/constants[12])/((1.00000+constants[6]/constants[10]+states[0]/constants[13])*(1.00000+algebraic[1]/constants[11]+states[6]/constants[14])*(1.00000+algebraic[0]/constants[12]+states[5]/constants[15]))
    algebraic[4] = custom_piecewise([greater_equal(algebraic[3] , 0.00000), algebraic[3] , True, algebraic[3]])
    algebraic[16] = (constants[75]*(power(states[5]/algebraic[0], constants[73])))/(power(constants[74], constants[73])+power(states[5]/algebraic[0], constants[73]))
    algebraic[18] = custom_piecewise([greater_equal(algebraic[16] , 0.00000), algebraic[16] , True, algebraic[16]])
    rates[5] = (algebraic[4]+algebraic[14])-algebraic[18]
    algebraic[5] = ((constants[22]/(constants[20]*constants[21]))*(states[0]*states[6]-(constants[7]*algebraic[1])/constants[17]))/((1.00000+states[0]/constants[20]+constants[7]/constants[18])*(1.00000+states[6]/constants[21]+algebraic[1]/constants[19]))
    algebraic[6] = custom_piecewise([greater_equal(algebraic[5] , 0.00000), algebraic[5] , True, algebraic[5]])
    algebraic[19] = (((constants[61]*states[0])/constants[60])*(1.00000-states[3]/(states[0]*constants[58]))*(power(states[0]/constants[60]+states[3]/constants[59], constants[57]-1.00000)))/(1.00000+power(states[0]/constants[60]+states[3]/constants[59], constants[57]))
    algebraic[20] = custom_piecewise([greater_equal(algebraic[19] , 0.00000), algebraic[19] , True, algebraic[19]])
    rates[0] = algebraic[4]-(algebraic[6]+algebraic[8]+algebraic[20])
    algebraic[21] = ((constants[64]*states[3])/constants[62])/(1.00000+states[3]/constants[62]+states[4]/constants[63])
    algebraic[22] = custom_piecewise([greater_equal(algebraic[21] , 0.00000), algebraic[21] , True, algebraic[21]])
    algebraic[27] = constants[80]*states[3]
    rates[3] = 0.500000*algebraic[20]-(algebraic[22]+algebraic[27])
    algebraic[25] = ((constants[72]/(constants[71]*constants[69]))*(states[4]*states[6]-(constants[4]*algebraic[1])/constants[67]))/((1.00000+states[4]/constants[71]+constants[4]/constants[70])*(1.00000+states[6]/constants[69]+algebraic[1]/constants[68]))
    algebraic[26] = custom_piecewise([greater_equal(algebraic[25] , 0.00000), algebraic[25] , True, algebraic[25]])
    algebraic[23] = ((constants[66]*states[4])/constants[65])/(1.00000+states[4]/constants[65])
    algebraic[24] = custom_piecewise([greater_equal(algebraic[23] , 0.00000), algebraic[23] , True, algebraic[23]])
    rates[4] = (algebraic[22]+algebraic[27])-(algebraic[26]+algebraic[24])
    algebraic[28] = ((constants[79]*states[6]*constants[8])/(constants[77]*constants[78]))/((1.00000+states[6]/constants[77]+algebraic[1]/constants[76])*(1.00000+constants[8]/constants[78]))
    algebraic[29] = custom_piecewise([greater_equal(algebraic[28] , 0.00000), algebraic[28] , True, algebraic[28]])
    rates[6] = (algebraic[4]+algebraic[8])-(algebraic[6]+algebraic[13]+algebraic[17]+algebraic[26]+algebraic[29])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[2] = constants[2]-states[7]
    algebraic[9] = ((constants[37]/(constants[36]*constants[31]))*(states[7]*constants[9]-(states[1]*algebraic[2])/constants[30]))/(1.00000+states[7]/constants[36]+constants[9]/constants[32]+states[1]/constants[35]+algebraic[2]/constants[33]+(states[7]*constants[9])/(constants[36]*constants[31])+(states[1]*algebraic[2])/(constants[34]*constants[33]))
    algebraic[10] = custom_piecewise([greater_equal(algebraic[9] , 0.00000), algebraic[9] , True, algebraic[9]])
    algebraic[0] = constants[0]-states[5]
    algebraic[12] = ((constants[43]/(constants[41]*constants[42]))*(states[1]*algebraic[0]-(constants[3]*states[5])/constants[38]))/((1.00000+states[1]/constants[42]+constants[3]/constants[39])*(1.00000+algebraic[0]/constants[41]+states[5]/constants[40]))
    algebraic[14] = custom_piecewise([greater_equal(algebraic[12] , 0.00000), algebraic[12] , True, algebraic[12]])
    algebraic[1] = constants[1]-states[6]
    algebraic[7] = (((((((constants[29]/(1.00000+(constants[23]*states[6])/algebraic[1]))*states[0])/constants[26])*algebraic[1])/constants[24])*algebraic[2])/constants[25])/((1.00000+states[0]/constants[26])*(1.00000+algebraic[1]/constants[24]+states[6]/constants[27])*(1.00000+algebraic[2]/constants[25]+states[7]/constants[28]))
    algebraic[8] = custom_piecewise([greater_equal(algebraic[7] , 0.00000), algebraic[7] , True, algebraic[7]])
    algebraic[11] = ((constants[50]/(constants[48]*constants[46]))*(states[7]*states[6]-(algebraic[1]*algebraic[2]*states[2])/constants[44]))/((1.00000+algebraic[1]/constants[45]+states[6]/constants[46])*(1.00000+states[7]/constants[48]+algebraic[2]/constants[47])*(1.00000+states[2]/constants[49]))
    algebraic[13] = custom_piecewise([greater_equal(algebraic[11] , 0.00000), algebraic[11] , True, algebraic[11]])
    algebraic[15] = ((constants[56]/(constants[55]*constants[53]))*(states[2]*states[6]-(constants[5]*algebraic[1])/constants[51]))/((1.00000+algebraic[1]/constants[52]+states[6]/constants[53])*(1.00000+states[2]/constants[55]+constants[5]/constants[54]))
    algebraic[17] = custom_piecewise([greater_equal(algebraic[15] , 0.00000), algebraic[15] , True, algebraic[15]])
    algebraic[3] = ((((((constants[16]*constants[6])/constants[10])*algebraic[1])/constants[11])*algebraic[0])/constants[12])/((1.00000+constants[6]/constants[10]+states[0]/constants[13])*(1.00000+algebraic[1]/constants[11]+states[6]/constants[14])*(1.00000+algebraic[0]/constants[12]+states[5]/constants[15]))
    algebraic[4] = custom_piecewise([greater_equal(algebraic[3] , 0.00000), algebraic[3] , True, algebraic[3]])
    algebraic[16] = (constants[75]*(power(states[5]/algebraic[0], constants[73])))/(power(constants[74], constants[73])+power(states[5]/algebraic[0], constants[73]))
    algebraic[18] = custom_piecewise([greater_equal(algebraic[16] , 0.00000), algebraic[16] , True, algebraic[16]])
    algebraic[5] = ((constants[22]/(constants[20]*constants[21]))*(states[0]*states[6]-(constants[7]*algebraic[1])/constants[17]))/((1.00000+states[0]/constants[20]+constants[7]/constants[18])*(1.00000+states[6]/constants[21]+algebraic[1]/constants[19]))
    algebraic[6] = custom_piecewise([greater_equal(algebraic[5] , 0.00000), algebraic[5] , True, algebraic[5]])
    algebraic[19] = (((constants[61]*states[0])/constants[60])*(1.00000-states[3]/(states[0]*constants[58]))*(power(states[0]/constants[60]+states[3]/constants[59], constants[57]-1.00000)))/(1.00000+power(states[0]/constants[60]+states[3]/constants[59], constants[57]))
    algebraic[20] = custom_piecewise([greater_equal(algebraic[19] , 0.00000), algebraic[19] , True, algebraic[19]])
    algebraic[21] = ((constants[64]*states[3])/constants[62])/(1.00000+states[3]/constants[62]+states[4]/constants[63])
    algebraic[22] = custom_piecewise([greater_equal(algebraic[21] , 0.00000), algebraic[21] , True, algebraic[21]])
    algebraic[27] = constants[80]*states[3]
    algebraic[25] = ((constants[72]/(constants[71]*constants[69]))*(states[4]*states[6]-(constants[4]*algebraic[1])/constants[67]))/((1.00000+states[4]/constants[71]+constants[4]/constants[70])*(1.00000+states[6]/constants[69]+algebraic[1]/constants[68]))
    algebraic[26] = custom_piecewise([greater_equal(algebraic[25] , 0.00000), algebraic[25] , True, algebraic[25]])
    algebraic[23] = ((constants[66]*states[4])/constants[65])/(1.00000+states[4]/constants[65])
    algebraic[24] = custom_piecewise([greater_equal(algebraic[23] , 0.00000), algebraic[23] , True, algebraic[23]])
    algebraic[28] = ((constants[79]*states[6]*constants[8])/(constants[77]*constants[78]))/((1.00000+states[6]/constants[77]+algebraic[1]/constants[76])*(1.00000+constants[8]/constants[78]))
    algebraic[29] = custom_piecewise([greater_equal(algebraic[28] , 0.00000), algebraic[28] , True, algebraic[28]])
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