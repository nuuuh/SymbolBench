# Size of variable arrays:
sizeAlgebraic = 22
sizeStates = 13
sizeConstants = 67
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "GlcI in component GlcI (millimolar)"
    legend_algebraic[5] = "vHK in component vHK (flux)"
    legend_algebraic[4] = "vGlcTr in component vGlcTr (flux)"
    legend_states[1] = "hexose_P in component hexose_P (millimolar)"
    legend_algebraic[6] = "vPFK in component vPFK (flux)"
    legend_states[2] = "Fru16BP in component Fru16BP (millimolar)"
    legend_algebraic[14] = "vALD in component vALD (flux)"
    legend_states[3] = "triose_P in component triose_P (millimolar)"
    legend_algebraic[16] = "vGAPDH in component vGAPDH (flux)"
    legend_algebraic[17] = "vGDH in component vGDH (flux)"
    legend_algebraic[8] = "vGPO in component vGPO (flux)"
    legend_states[4] = "BPGA13 in component BPGA13 (millimolar)"
    legend_algebraic[20] = "vPGK in component vPGK (flux)"
    legend_states[5] = "N in component N (millimolar)"
    legend_algebraic[21] = "vPK in component vPK (flux)"
    legend_states[6] = "Pyr in component Pyr (millimolar)"
    legend_algebraic[9] = "vPyrTr in component vPyrTr (flux)"
    legend_states[7] = "NADH in component NADH (millimolar)"
    legend_states[8] = "NAD in component NAD (millimolar)"
    legend_states[9] = "Gly3P in component Gly3P (millimolar)"
    legend_algebraic[11] = "vGlyK in component vGlyK (flux)"
    legend_states[10] = "Gly in component Gly (millimolar)"
    legend_states[11] = "P in component P (millimolar)"
    legend_algebraic[10] = "vATPase in component vATPase (flux)"
    legend_algebraic[0] = "ATP in component ATP (millimolar)"
    legend_constants[0] = "sumA in component ATP (millimolar)"
    legend_constants[1] = "Keq_AK in component ATP (dimensionless)"
    legend_algebraic[1] = "ADP in component ADP (millimolar)"
    legend_algebraic[12] = "DHAP in component DHAP (millimolar)"
    legend_algebraic[3] = "Fru6P in component Fru6P (millimolar)"
    legend_algebraic[13] = "GAP in component GAP (millimolar)"
    legend_states[12] = "Glc6P in component Glc6P (millimolar)"
    legend_constants[2] = "sumc4 in component DHAP (millimolar)"
    legend_constants[3] = "sumc5 in component DHAP (millimolar)"
    legend_algebraic[2] = "GlcE in component GlcE (millimolar)"
    legend_algebraic[7] = "vPGI in component vPGI (flux)"
    legend_algebraic[18] = "PGA3 in component PGA3 (millimolar)"
    legend_algebraic[19] = "PEP in component PEP (millimolar)"
    legend_constants[4] = "Keq_ENO in component PEP (dimensionless)"
    legend_constants[5] = "Keq_PGM in component PEP (dimensionless)"
    legend_constants[6] = "K_Glc in component vGlcTr (millimolar)"
    legend_constants[7] = "alpha in component vGlcTr (dimensionless)"
    legend_constants[8] = "vGlcTr_max in component vGlcTr (flux)"
    legend_constants[9] = "K_GlcI in component vHK (millimolar)"
    legend_constants[10] = "K_Glc6P in component vHK (millimolar)"
    legend_constants[11] = "K_ATP in component vHK (millimolar)"
    legend_constants[12] = "K_ADP in component vHK (millimolar)"
    legend_constants[13] = "vHK_max in component vHK (flux)"
    legend_constants[14] = "K_Glc6P in component vPGI (millimolar)"
    legend_constants[15] = "K_Fru6P in component vPGI (millimolar)"
    legend_constants[16] = "vPGI_max in component vPGI (flux)"
    legend_constants[17] = "Ki_1 in component vPFK (millimolar)"
    legend_constants[18] = "Ki_2 in component vPFK (millimolar)"
    legend_constants[19] = "KM_Fru6P in component vPFK (millimolar)"
    legend_constants[20] = "KM_ATP in component vPFK (millimolar)"
    legend_constants[21] = "vPFK_max in component vPFK (flux)"
    legend_constants[22] = "sumA in component vALD (millimolar)"
    legend_constants[23] = "KM_GAP in component vALD (millimolar)"
    legend_constants[24] = "Ki_GAP in component vALD (millimolar)"
    legend_constants[25] = "KM_DHAP in component vALD (millimolar)"
    legend_constants[26] = "vALD_max_forward in component vALD (flux)"
    legend_constants[27] = "vALD_max_reverse in component vALD (flux)"
    legend_algebraic[15] = "vTPI in component vTPI (flux)"
    legend_constants[28] = "K_DHAP in component vTPI (millimolar)"
    legend_constants[29] = "K_GAP in component vTPI (millimolar)"
    legend_constants[30] = "vTPI_max in component vTPI (flux)"
    legend_constants[31] = "K_NAD in component vGAPDH (millimolar)"
    legend_constants[32] = "K_GAP in component vGAPDH (millimolar)"
    legend_constants[33] = "K_BPGA13 in component vGAPDH (millimolar)"
    legend_constants[34] = "K_NADH in component vGAPDH (millimolar)"
    legend_constants[35] = "vGAPDH_max_forward in component vGAPDH (flux)"
    legend_constants[36] = "vGAPDH_max_reverse in component vGAPDH (flux)"
    legend_constants[37] = "vGAPDH_max in component vGAPDH (dimensionless)"
    legend_constants[38] = "K_NADH in component vGDH (millimolar)"
    legend_constants[39] = "K_Gly3P in component vGDH (millimolar)"
    legend_constants[40] = "K_DHAP in component vGDH (millimolar)"
    legend_constants[41] = "K_NAD in component vGDH (millimolar)"
    legend_constants[42] = "vGDH_max_forward in component vGDH (flux)"
    legend_constants[43] = "vGDH_max_reverse in component vGDH (flux)"
    legend_constants[44] = "vGDH_max in component vGDH (dimensionless)"
    legend_constants[45] = "K_Gly3P in component vGPO (millimolar)"
    legend_constants[46] = "vGPO_max in component vGPO (flux)"
    legend_constants[47] = "K_pyruvate in component vPyrTr (millimolar)"
    legend_constants[48] = "vPyrTr_max in component vPyrTr (flux)"
    legend_constants[49] = "K_ADP in component vPGK (millimolar)"
    legend_constants[50] = "K_BPGA13 in component vPGK (millimolar)"
    legend_constants[51] = "K_PGA3 in component vPGK (millimolar)"
    legend_constants[52] = "K_ATP in component vPGK (millimolar)"
    legend_constants[53] = "vPGK_max_forward in component vPGK (flux)"
    legend_constants[54] = "vPGK_max_reverse in component vPGK (flux)"
    legend_constants[55] = "vPGK_max in component vPGK (dimensionless)"
    legend_constants[56] = "KM_ADP in component vPK (millimolar)"
    legend_constants[57] = "n in component vPK (dimensionless)"
    legend_constants[58] = "vPK_max in component vPK (flux)"
    legend_constants[59] = "k in component vATPase (flux)"
    legend_constants[60] = "K_ADP in component vGlyK (millimolar)"
    legend_constants[61] = "K_Gly3P in component vGlyK (millimolar)"
    legend_constants[62] = "K_Gly in component vGlyK (millimolar)"
    legend_constants[63] = "K_ATP in component vGlyK (millimolar)"
    legend_constants[64] = "vGlyK_max_forward in component vGlyK (flux)"
    legend_constants[65] = "vGlyK_max_reverse in component vGlyK (flux)"
    legend_constants[66] = "vGlyK_max in component vGlyK (dimensionless)"
    legend_rates[0] = "d/dt GlcI in component GlcI (millimolar)"
    legend_rates[1] = "d/dt hexose_P in component hexose_P (millimolar)"
    legend_rates[2] = "d/dt Fru16BP in component Fru16BP (millimolar)"
    legend_rates[3] = "d/dt triose_P in component triose_P (millimolar)"
    legend_rates[4] = "d/dt BPGA13 in component BPGA13 (millimolar)"
    legend_rates[5] = "d/dt N in component N (millimolar)"
    legend_rates[6] = "d/dt Pyr in component Pyr (millimolar)"
    legend_rates[7] = "d/dt NADH in component NADH (millimolar)"
    legend_rates[8] = "d/dt NAD in component NAD (millimolar)"
    legend_rates[9] = "d/dt Gly3P in component Gly3P (millimolar)"
    legend_rates[10] = "d/dt Gly in component Gly (millimolar)"
    legend_rates[11] = "d/dt P in component P (millimolar)"
    legend_rates[12] = "d/dt Glc6P in component Glc6P (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.0340009
    states[1] = 2.583763
    states[2] = 16.5371
    states[3] = 3.9391429
    states[4] = 0.0326745
    states[5] = 1.59603
    states[6] = 4.77413
    states[7] = 0.0448639
    states[8] = 0.0448639
    states[9] = 0.0
    states[10] = 0.0
    states[11] = 7.63936
    constants[0] = 3.9
    constants[1] = 0.442
    states[12] = 2.07199
    constants[2] = 45.0
    constants[3] = 5.0
    constants[4] = 6.7
    constants[5] = 0.187
    constants[6] = 2.0
    constants[7] = 0.75
    constants[8] = 106.2
    constants[9] = 0.1
    constants[10] = 12.0
    constants[11] = 0.116
    constants[12] = 0.126
    constants[13] = 625.0
    constants[14] = 0.4
    constants[15] = 0.12
    constants[16] = 848.0
    constants[17] = 15.8
    constants[18] = 10.7
    constants[19] = 0.82
    constants[20] = 0.026
    constants[21] = 780.0
    constants[22] = 6.0
    constants[23] = 0.067
    constants[24] = 0.098
    constants[25] = 0.015
    constants[26] = 184.5
    constants[27] = 219.555
    constants[28] = 1.2
    constants[29] = 0.25
    constants[30] = 842.0
    constants[31] = 0.45
    constants[32] = 0.15
    constants[33] = 0.1
    constants[34] = 0.02
    constants[35] = 1470.0
    constants[36] = 984.9
    constants[37] = 1.0
    constants[38] = 0.01
    constants[39] = 2.0
    constants[40] = 0.1
    constants[41] = 0.4
    constants[42] = 533.0
    constants[43] = 149.24
    constants[44] = 1.0
    constants[45] = 1.7
    constants[46] = 368.0
    constants[47] = 1.96
    constants[48] = 200.0
    constants[49] = 0.1
    constants[50] = 0.05
    constants[51] = 1.62
    constants[52] = 0.29
    constants[53] = 640.0
    constants[54] = 18.56
    constants[55] = 1.0
    constants[56] = 0.114
    constants[57] = 2.5
    constants[58] = 2600
    constants[59] = 50
    constants[60] = 0.12
    constants[61] = 5.1
    constants[62] = 0.12
    constants[63] = 0.19
    constants[64] = 220.0
    constants[65] = 334000.0
    constants[66] = 1.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = ((states[11]*(1.00000-4.00000*constants[1])-constants[0])+power(power(constants[0]-(1.00000-4.00000*constants[1])*states[11], 2.00000)+4.00000*(1.00000-4.00000*constants[1])*(-constants[1]*(power(states[11], 2.00000))), 0.500000))/(2.00000*(1.00000-4.00000*constants[1]))
    algebraic[1] = states[11]-2.00000*algebraic[0]
    algebraic[5] = (constants[13]*states[0]*algebraic[0])/(constants[11]*constants[9]*(1.00000+states[12]/constants[10]+states[0]/constants[9])*(1.00000+algebraic[0]/constants[11]+algebraic[1]/constants[12]))
    algebraic[2] = custom_piecewise([greater_equal(voi , 60.0000) & less(voi , 61.0000), 5.00000 , True, 0.0500000])
    algebraic[4] = constants[8]*((algebraic[2]-states[0])/(constants[6]+algebraic[2]+states[0]+constants[7]*algebraic[2]*(states[0]/constants[6])))
    rates[0] = algebraic[4]-algebraic[5]
    algebraic[3] = states[1]-states[12]
    algebraic[6] = (constants[17]*constants[21]*algebraic[3]*algebraic[0])/(constants[20]*constants[19]*(states[2]+constants[17])*(1.00000+states[2]/constants[18]+algebraic[3]/constants[19])*(1.00000+algebraic[0]/constants[20]))
    rates[1] = algebraic[5]-algebraic[6]
    algebraic[7] = (constants[16]*(states[12]/constants[14]-algebraic[3]/constants[15]))/(1.00000+states[12]/constants[14]+algebraic[3]/constants[15])
    rates[12] = algebraic[5]-algebraic[7]
    algebraic[11] = (constants[66]*((constants[64]*algebraic[1]*states[9])/(constants[60]*constants[61])-(constants[65]*algebraic[0]*states[10])/(constants[63]*constants[62])))/((1.00000+states[9]/constants[61]+states[10]/constants[62])*(1.00000+algebraic[1]/constants[60]+algebraic[0]/constants[63]))
    rates[10] = algebraic[11]
    rootfind_0(voi, constants, rates, states, algebraic)
    algebraic[14] = ((constants[26]*states[2])/(0.00900000*(1.00000+algebraic[0]/0.680000+algebraic[1]/1.51000+(constants[22]-(algebraic[0]+algebraic[1]))/3.65000))-(constants[27]*algebraic[13]*algebraic[12])/(constants[25]*constants[23]))/(1.00000+algebraic[13]/constants[23]+algebraic[12]/constants[25]+(algebraic[13]*algebraic[12])/(constants[25]*constants[23])+states[2]/(0.00900000*(1.00000+algebraic[0]/0.680000+algebraic[1]/1.51000+(constants[22]-(algebraic[0]+algebraic[1]))/3.65000))+(states[2]*algebraic[13])/(constants[24]*0.00900000*(1.00000+algebraic[0]/0.680000+algebraic[1]/1.51000+(constants[22]-(algebraic[0]+algebraic[1]))/3.65000)))
    rates[2] = algebraic[6]-algebraic[14]
    algebraic[16] = constants[37]*((constants[35]*(algebraic[13]*((states[8]/constants[32])/constants[31])-(constants[36]/constants[35])*(((states[4]*states[7])/constants[33])/constants[34])))/((1.00000+algebraic[13]/constants[32]+states[4]/constants[33])*(1.00000+states[8]/constants[31]+states[7]/constants[34])))
    algebraic[17] = (constants[44]*constants[42]*((states[7]*algebraic[12])/(constants[38]*constants[40])-(constants[43]*states[8]*states[9])/(constants[39]*constants[41]*constants[42])))/((1.00000+states[8]/constants[41]+states[7]/constants[38])*(1.00000+algebraic[12]/constants[40]+states[9]/constants[39]))
    algebraic[8] = (constants[46]*states[9])/(states[9]+constants[45])
    rates[3] = (2.00000*algebraic[14]+algebraic[8])-(algebraic[16]+algebraic[17])
    rates[7] = algebraic[16]-algebraic[17]
    rates[8] = algebraic[17]-algebraic[16]
    rates[9] = algebraic[17]-(algebraic[11]+algebraic[8])
    rootfind_1(voi, constants, rates, states, algebraic)
    algebraic[20] = (constants[55]*constants[53]*((-constants[54]*algebraic[18]*algebraic[0])/(constants[52]*constants[51]*constants[53])+(states[4]*algebraic[1])/(constants[50]*constants[49])))/((1.00000+states[4]/constants[50]+algebraic[18]/constants[51])*(1.00000+algebraic[1]/constants[49]+algebraic[0]/constants[52]))
    rates[4] = algebraic[16]-algebraic[20]
    algebraic[21] = ((constants[58]*(power(algebraic[19]/(0.340000*(1.00000+algebraic[0]/0.570000+algebraic[1]/0.640000)), constants[57]))*algebraic[1])/constants[56])/((1.00000+power(algebraic[19]/(0.340000*(1.00000+algebraic[0]/0.570000+algebraic[1]/0.640000)), constants[57]))*(1.00000+algebraic[1]/constants[56]))
    rates[5] = algebraic[20]-algebraic[21]
    algebraic[9] = ((constants[48]*states[6])/constants[47])/(1.00000+states[6]/constants[47])
    rates[6] = algebraic[21]-algebraic[9]
    algebraic[10] = (constants[59]*algebraic[0])/algebraic[1]
    rates[11] = (algebraic[20]+algebraic[11]+algebraic[21])-(algebraic[5]+algebraic[6]+algebraic[10])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = ((states[11]*(1.00000-4.00000*constants[1])-constants[0])+power(power(constants[0]-(1.00000-4.00000*constants[1])*states[11], 2.00000)+4.00000*(1.00000-4.00000*constants[1])*(-constants[1]*(power(states[11], 2.00000))), 0.500000))/(2.00000*(1.00000-4.00000*constants[1]))
    algebraic[1] = states[11]-2.00000*algebraic[0]
    algebraic[5] = (constants[13]*states[0]*algebraic[0])/(constants[11]*constants[9]*(1.00000+states[12]/constants[10]+states[0]/constants[9])*(1.00000+algebraic[0]/constants[11]+algebraic[1]/constants[12]))
    algebraic[2] = custom_piecewise([greater_equal(voi , 60.0000) & less(voi , 61.0000), 5.00000 , True, 0.0500000])
    algebraic[4] = constants[8]*((algebraic[2]-states[0])/(constants[6]+algebraic[2]+states[0]+constants[7]*algebraic[2]*(states[0]/constants[6])))
    algebraic[3] = states[1]-states[12]
    algebraic[6] = (constants[17]*constants[21]*algebraic[3]*algebraic[0])/(constants[20]*constants[19]*(states[2]+constants[17])*(1.00000+states[2]/constants[18]+algebraic[3]/constants[19])*(1.00000+algebraic[0]/constants[20]))
    algebraic[7] = (constants[16]*(states[12]/constants[14]-algebraic[3]/constants[15]))/(1.00000+states[12]/constants[14]+algebraic[3]/constants[15])
    algebraic[11] = (constants[66]*((constants[64]*algebraic[1]*states[9])/(constants[60]*constants[61])-(constants[65]*algebraic[0]*states[10])/(constants[63]*constants[62])))/((1.00000+states[9]/constants[61]+states[10]/constants[62])*(1.00000+algebraic[1]/constants[60]+algebraic[0]/constants[63]))
    algebraic[14] = ((constants[26]*states[2])/(0.00900000*(1.00000+algebraic[0]/0.680000+algebraic[1]/1.51000+(constants[22]-(algebraic[0]+algebraic[1]))/3.65000))-(constants[27]*algebraic[13]*algebraic[12])/(constants[25]*constants[23]))/(1.00000+algebraic[13]/constants[23]+algebraic[12]/constants[25]+(algebraic[13]*algebraic[12])/(constants[25]*constants[23])+states[2]/(0.00900000*(1.00000+algebraic[0]/0.680000+algebraic[1]/1.51000+(constants[22]-(algebraic[0]+algebraic[1]))/3.65000))+(states[2]*algebraic[13])/(constants[24]*0.00900000*(1.00000+algebraic[0]/0.680000+algebraic[1]/1.51000+(constants[22]-(algebraic[0]+algebraic[1]))/3.65000)))
    algebraic[16] = constants[37]*((constants[35]*(algebraic[13]*((states[8]/constants[32])/constants[31])-(constants[36]/constants[35])*(((states[4]*states[7])/constants[33])/constants[34])))/((1.00000+algebraic[13]/constants[32]+states[4]/constants[33])*(1.00000+states[8]/constants[31]+states[7]/constants[34])))
    algebraic[17] = (constants[44]*constants[42]*((states[7]*algebraic[12])/(constants[38]*constants[40])-(constants[43]*states[8]*states[9])/(constants[39]*constants[41]*constants[42])))/((1.00000+states[8]/constants[41]+states[7]/constants[38])*(1.00000+algebraic[12]/constants[40]+states[9]/constants[39]))
    algebraic[8] = (constants[46]*states[9])/(states[9]+constants[45])
    algebraic[20] = (constants[55]*constants[53]*((-constants[54]*algebraic[18]*algebraic[0])/(constants[52]*constants[51]*constants[53])+(states[4]*algebraic[1])/(constants[50]*constants[49])))/((1.00000+states[4]/constants[50]+algebraic[18]/constants[51])*(1.00000+algebraic[1]/constants[49]+algebraic[0]/constants[52]))
    algebraic[21] = ((constants[58]*(power(algebraic[19]/(0.340000*(1.00000+algebraic[0]/0.570000+algebraic[1]/0.640000)), constants[57]))*algebraic[1])/constants[56])/((1.00000+power(algebraic[19]/(0.340000*(1.00000+algebraic[0]/0.570000+algebraic[1]/0.640000)), constants[57]))*(1.00000+algebraic[1]/constants[56]))
    algebraic[9] = ((constants[48]*states[6])/constants[47])/(1.00000+states[6]/constants[47])
    algebraic[10] = (constants[59]*algebraic[0])/algebraic[1]
    algebraic[15] = (constants[30]*(algebraic[12]/constants[28]-(5.70000*algebraic[13])/constants[29]))/(1.00000+algebraic[13]/constants[29]+algebraic[12]/constants[28])
    return algebraic

initialGuess0 = None
def rootfind_0(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess0
    if initialGuess0 is None: initialGuess0 = ones(2)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_0, initialGuess0, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess0 = soln
        algebraic[12] = soln[0]
        algebraic[13] = soln[1]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess0 = soln
            algebraic[12][i] = soln[0]
            algebraic[13][i] = soln[1]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 2)
    algebraic[12] = algebraicCandidate[0]
    algebraic[13] = algebraicCandidate[1]
    resid[0] = (algebraic[12]-(constants[3]*algebraic[12])/((constants[2]+constants[3])-(states[4]+2.00000*states[2]+algebraic[3]+algebraic[13]+states[12]+states[11])))
    resid[1] = (algebraic[13]-(states[3]-algebraic[12]))
    return resid

initialGuess1 = None
def rootfind_1(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess1
    if initialGuess1 is None: initialGuess1 = ones(2)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_1, initialGuess1, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess1 = soln
        algebraic[18] = soln[0]
        algebraic[19] = soln[1]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_1, initialGuess1, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess1 = soln
            algebraic[18][i] = soln[0]
            algebraic[19][i] = soln[1]

def residualSN_1(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 2)
    algebraic[18] = algebraicCandidate[0]
    algebraic[19] = algebraicCandidate[1]
    resid[0] = (algebraic[18]-(states[5]-algebraic[19]))
    resid[1] = (algebraic[19]-constants[4]*constants[5]*algebraic[18])
    return resid

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