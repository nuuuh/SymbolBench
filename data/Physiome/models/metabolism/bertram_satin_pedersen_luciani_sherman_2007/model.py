# Size of variable arrays:
sizeAlgebraic = 79
sizeStates = 11
sizeConstants = 75
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "Vm in component membrane (millivolt)"
    legend_constants[0] = "cm in component membrane (femtofarad)"
    legend_algebraic[3] = "Ica in component Ica (femtoampere)"
    legend_algebraic[1] = "Ik in component Ik (femtoampere)"
    legend_algebraic[4] = "Ikca in component Ikca (femtoampere)"
    legend_algebraic[43] = "Ikatp in component Ikatp (femtoampere)"
    legend_constants[1] = "gK in component Ik (picosiemens)"
    legend_constants[2] = "VK in component model_parameters (millivolt)"
    legend_states[1] = "n in component n (dimensionless)"
    legend_algebraic[0] = "n_infinity in component n (dimensionless)"
    legend_constants[3] = "tau_n in component n (millisecond)"
    legend_constants[4] = "vn in component n (millivolt)"
    legend_constants[5] = "sn in component n (millivolt)"
    legend_constants[6] = "gCa in component Ica (picosiemens)"
    legend_constants[7] = "VCa in component model_parameters (millivolt)"
    legend_algebraic[2] = "m_infinity in component m (dimensionless)"
    legend_constants[8] = "v in component m (millivolt)"
    legend_constants[9] = "sm in component m (millivolt)"
    legend_constants[10] = "gkCa in component Ikca (picosiemens)"
    legend_constants[11] = "kd in component Ikca (micromolar)"
    legend_states[2] = "c in component c (micromolar)"
    legend_constants[12] = "gkATP_ in component Ikatp (picosiemens)"
    legend_algebraic[41] = "katpo in component Ikatp (dimensionless)"
    legend_algebraic[6] = "topo in component Ikatp (dimensionless)"
    legend_algebraic[39] = "bottomo in component Ikatp (dimensionless)"
    legend_algebraic[5] = "mgadp in component Ikatp (micromolar)"
    legend_algebraic[7] = "adp3m in component Ikatp (micromolar)"
    legend_algebraic[36] = "atp4m in component Ikatp (micromolar)"
    legend_algebraic[35] = "atp in component atp (micromolar)"
    legend_states[3] = "adp in component adp (micromolar)"
    legend_algebraic[8] = "JGPDH in component JGPDH (micromolar_millisecond)"
    legend_constants[13] = "kGPDH in component JGPDH (micromolar_millisecond)"
    legend_states[4] = "FBP in component FBP (micromolar)"
    legend_states[5] = "G6P in component G6P (micromolar)"
    legend_constants[74] = "JGK_ms in component JGK (micromolar_millisecond)"
    legend_algebraic[78] = "JPFK_ms in component JPFK (micromolar_millisecond)"
    legend_algebraic[9] = "F6P in component F6P (micromolar)"
    legend_algebraic[77] = "JPFK in component JPFK (micromolar_second)"
    legend_constants[14] = "bottom1 in component JPFK (dimensionless)"
    legend_constants[15] = "topa1 in component JPFK (dimensionless)"
    legend_constants[16] = "k1 in component JPFK (micromolar)"
    legend_constants[17] = "k2 in component JPFK (micromolar)"
    legend_constants[18] = "k3 in component JPFK (micromolar)"
    legend_constants[19] = "k4 in component JPFK (micromolar)"
    legend_constants[20] = "VmaxPFK in component JPFK (micromolar_millisecond)"
    legend_algebraic[37] = "weight2 in component JPFK (dimensionless)"
    legend_constants[71] = "topa2 in component JPFK (dimensionless)"
    legend_algebraic[40] = "bottom2 in component JPFK (dimensionless)"
    legend_algebraic[13] = "topa3 in component JPFK (dimensionless)"
    legend_algebraic[11] = "weight3 in component JPFK (dimensionless)"
    legend_algebraic[42] = "bottom3 in component JPFK (dimensionless)"
    legend_constants[21] = "f13 in component JPFK (dimensionless)"
    legend_constants[22] = "f43 in component JPFK (dimensionless)"
    legend_constants[23] = "f23 in component JPFK (dimensionless)"
    legend_constants[24] = "f42 in component JPFK (dimensionless)"
    legend_constants[25] = "f41 in component JPFK (dimensionless)"
    legend_algebraic[44] = "weight4 in component JPFK (dimensionless)"
    legend_algebraic[45] = "topa4 in component JPFK (dimensionless)"
    legend_algebraic[46] = "bottom4 in component JPFK (dimensionless)"
    legend_algebraic[15] = "weight5 in component JPFK (dimensionless)"
    legend_algebraic[47] = "topa5 in component JPFK (dimensionless)"
    legend_algebraic[48] = "bottom5 in component JPFK (dimensionless)"
    legend_algebraic[49] = "weight6 in component JPFK (dimensionless)"
    legend_algebraic[50] = "topa6 in component JPFK (dimensionless)"
    legend_algebraic[51] = "bottom6 in component JPFK (dimensionless)"
    legend_algebraic[16] = "weight7 in component JPFK (dimensionless)"
    legend_algebraic[52] = "topa7 in component JPFK (dimensionless)"
    legend_algebraic[53] = "bottom7 in component JPFK (dimensionless)"
    legend_algebraic[54] = "weight8 in component JPFK (dimensionless)"
    legend_algebraic[55] = "topa8 in component JPFK (dimensionless)"
    legend_algebraic[56] = "bottom8 in component JPFK (dimensionless)"
    legend_constants[73] = "weight9 in component JPFK (dimensionless)"
    legend_algebraic[57] = "topa9 in component JPFK (dimensionless)"
    legend_algebraic[58] = "bottom9 in component JPFK (dimensionless)"
    legend_algebraic[59] = "weight10 in component JPFK (dimensionless)"
    legend_algebraic[60] = "topa10 in component JPFK (dimensionless)"
    legend_algebraic[61] = "bottom10 in component JPFK (dimensionless)"
    legend_algebraic[17] = "weight11 in component JPFK (dimensionless)"
    legend_algebraic[62] = "topa11 in component JPFK (dimensionless)"
    legend_algebraic[63] = "bottom11 in component JPFK (dimensionless)"
    legend_algebraic[64] = "weight12 in component JPFK (dimensionless)"
    legend_algebraic[65] = "topa12 in component JPFK (dimensionless)"
    legend_algebraic[66] = "bottom12 in component JPFK (dimensionless)"
    legend_algebraic[18] = "weight13 in component JPFK (dimensionless)"
    legend_algebraic[67] = "topa13 in component JPFK (dimensionless)"
    legend_algebraic[68] = "bottom13 in component JPFK (dimensionless)"
    legend_algebraic[69] = "weight14 in component JPFK (dimensionless)"
    legend_algebraic[70] = "topa14 in component JPFK (dimensionless)"
    legend_algebraic[71] = "bottom14 in component JPFK (dimensionless)"
    legend_algebraic[19] = "weight15 in component JPFK (dimensionless)"
    legend_algebraic[72] = "topa15 in component JPFK (dimensionless)"
    legend_algebraic[73] = "bottom15 in component JPFK (dimensionless)"
    legend_algebraic[74] = "weight16 in component JPFK (dimensionless)"
    legend_algebraic[75] = "topa16 in component JPFK (dimensionless)"
    legend_algebraic[76] = "bottom16 in component JPFK (dimensionless)"
    legend_algebraic[20] = "topb in component JPFK (dimensionless)"
    legend_constants[26] = "AMP in component model_parameters (micromolar)"
    legend_constants[27] = "lambda in component JPFK (dimensionless)"
    legend_constants[28] = "kappa in component model_parameters (dimensionless)"
    legend_algebraic[14] = "JPDH in component JPDH (micromolar_millisecond)"
    legend_constants[29] = "p1 in component JPDH (dimensionless)"
    legend_constants[30] = "p2 in component JPDH (dimensionless)"
    legend_constants[31] = "p3 in component JPDH (micromolar)"
    legend_constants[32] = "JGPDHbas in component JPDH (micromolar_millisecond)"
    legend_algebraic[12] = "NADm in component NADm (millimolar)"
    legend_states[6] = "Cam in component Cam (micromolar)"
    legend_states[7] = "NADHm in component NADHm (millimolar)"
    legend_constants[33] = "gamma in component model_parameters (dimensionless)"
    legend_algebraic[10] = "JO in component JO (micromolar_millisecond)"
    legend_constants[34] = "p4 in component JO (micromolar_millisecond)"
    legend_constants[35] = "p5 in component JO (millimolar)"
    legend_constants[36] = "p6 in component JO (millivolt)"
    legend_constants[37] = "p7 in component JO (millivolt)"
    legend_states[8] = "delta_psi in component delta_psi (millivolt)"
    legend_constants[38] = "NADmtot in component NADm (millimolar)"
    legend_constants[39] = "Cmito in component delta_psi (micromolar_millivolt)"
    legend_algebraic[21] = "JHres in component JHres (micromolar_millisecond)"
    legend_algebraic[29] = "JHatp in component JHatp (micromolar_millisecond)"
    legend_algebraic[33] = "JANT in component JANT (micromolar_millisecond)"
    legend_algebraic[22] = "JHleak in component JHleak (micromolar_millisecond)"
    legend_algebraic[24] = "JNaCa in component JNaCa (micromolar_millisecond)"
    legend_algebraic[23] = "Juni in component Juni (micromolar_millisecond)"
    legend_constants[40] = "p8 in component JHres (micromolar_millisecond)"
    legend_constants[41] = "p9 in component JHres (millimolar)"
    legend_constants[42] = "p10 in component JHres (millivolt)"
    legend_constants[43] = "p11 in component JHres (millivolt)"
    legend_algebraic[27] = "JF1F0 in component JF1F0 (micromolar_millisecond)"
    legend_constants[44] = "p13 in component JF1F0 (millimolar)"
    legend_constants[45] = "p14 in component JF1F0 (millivolt)"
    legend_constants[46] = "p15 in component JF1F0 (millivolt)"
    legend_constants[47] = "p16 in component JF1F0 (micromolar_millisecond)"
    legend_algebraic[25] = "ATPm in component ATPm (millimolar)"
    legend_constants[48] = "JGK in component JGK (micromolar_second)"
    legend_constants[49] = "p17 in component JHleak (micromolar_millisecond_millivolt)"
    legend_constants[50] = "p18 in component JHleak (micromolar_millisecond)"
    legend_constants[51] = "p19 in component JANT (micromolar_millisecond)"
    legend_constants[52] = "p20 in component JANT (dimensionless)"
    legend_constants[53] = "FRT in component JANT (per_millivolt)"
    legend_algebraic[31] = "RATm in component RATm (dimensionless)"
    legend_states[9] = "ADPm in component ADPm (millimolar)"
    legend_constants[54] = "p21 in component Juni (per_micromolar_millisecond_millivolt)"
    legend_constants[55] = "p22 in component Juni (second_order_rate_constant)"
    legend_constants[56] = "p23 in component JNaCa (micromolar_millisecond)"
    legend_constants[57] = "p24 in component JNaCa (per_millivolt)"
    legend_constants[58] = "fmito in component Cam (dimensionless)"
    legend_algebraic[26] = "Jmito in component Jmito (micromolar_millisecond)"
    legend_constants[59] = "Amtot in component ATPm (millimolar)"
    legend_constants[72] = "delta in component model_parameters (dimensionless)"
    legend_algebraic[38] = "Jhyd in component Jhyd (micromolar_millisecond)"
    legend_constants[60] = "khyd in component Jhyd (second_order_rate_constant)"
    legend_constants[61] = "khydbas in component Jhyd (first_order_rate_constant)"
    legend_constants[62] = "atot in component atp (micromolar)"
    legend_constants[63] = "fcyt in component c (dimensionless)"
    legend_algebraic[34] = "Jer in component Jer (micromolar_millisecond)"
    legend_algebraic[28] = "Jmem in component Jmem (micromolar_millisecond)"
    legend_constants[64] = "kPMCA in component Jmem (first_order_rate_constant)"
    legend_constants[65] = "alpha in component Jmem (micromolar_millisecond_femtoampere)"
    legend_constants[66] = "Cbas in component Jmem (micromolar)"
    legend_algebraic[30] = "Jleak in component Jleak (micromolar_millisecond)"
    legend_constants[67] = "pleak in component Jleak (first_order_rate_constant)"
    legend_states[10] = "Caer in component Caer (micromolar)"
    legend_algebraic[32] = "JSERCA in component JSERCA (micromolar_millisecond)"
    legend_constants[68] = "kSERCA in component JSERCA (first_order_rate_constant)"
    legend_constants[69] = "fer in component Caer (dimensionless)"
    legend_constants[70] = "Vc_Ver in component Caer (dimensionless)"
    legend_rates[0] = "d/dt Vm in component membrane (millivolt)"
    legend_rates[1] = "d/dt n in component n (dimensionless)"
    legend_rates[5] = "d/dt G6P in component G6P (micromolar)"
    legend_rates[4] = "d/dt FBP in component FBP (micromolar)"
    legend_rates[7] = "d/dt NADHm in component NADHm (millimolar)"
    legend_rates[8] = "d/dt delta_psi in component delta_psi (millivolt)"
    legend_rates[6] = "d/dt Cam in component Cam (micromolar)"
    legend_rates[9] = "d/dt ADPm in component ADPm (millimolar)"
    legend_rates[3] = "d/dt adp in component adp (micromolar)"
    legend_rates[2] = "d/dt c in component c (micromolar)"
    legend_rates[10] = "d/dt Caer in component Caer (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -60.0
    constants[0] = 5300.0
    constants[1] = 2700.0
    constants[2] = -75.0
    states[1] = 0.0
    constants[3] = 20.0
    constants[4] = -16.0
    constants[5] = 5.0
    constants[6] = 1000.0
    constants[7] = 25.0
    constants[8] = -20.0
    constants[9] = 12.0
    constants[10] = 300.0
    constants[11] = 0.5
    states[2] = 0.17
    constants[12] = 16000.0
    states[3] = 1137.0
    constants[13] = 0.0005
    states[4] = 2.16
    states[5] = 301.0
    constants[14] = 1
    constants[15] = 0
    constants[16] = 30
    constants[17] = 1
    constants[18] = 50000
    constants[19] = 1000
    constants[20] = 5.0
    constants[21] = 0.02
    constants[22] = 20
    constants[23] = 0.2
    constants[24] = 20
    constants[25] = 20
    constants[26] = 500.0
    constants[27] = 0.06
    constants[28] = 0.001
    constants[29] = 400.0
    constants[30] = 1.0
    constants[31] = 0.01
    constants[32] = 0.0005
    states[6] = 0.2
    states[7] = 0.4
    constants[33] = 0.001
    constants[34] = 0.6
    constants[35] = 0.1
    constants[36] = 177.0
    constants[37] = 5.0
    states[8] = 164.0
    constants[38] = 10.0
    constants[39] = 1.8
    constants[40] = 7.0
    constants[41] = 0.1
    constants[42] = 177.0
    constants[43] = 5.0
    constants[44] = 10.0
    constants[45] = 190.0
    constants[46] = 8.5
    constants[47] = 35.0
    constants[48] = 0.4
    constants[49] = 0.002
    constants[50] = -0.03
    constants[51] = 0.35
    constants[52] = 2.0
    constants[53] = 0.037410133
    states[9] = 11.1
    constants[54] = 0.04
    constants[55] = 1.1
    constants[56] = 0.01
    constants[57] = 0.016
    constants[58] = 0.01
    constants[59] = 15.0
    constants[60] = 0.00005
    constants[61] = 0.00005
    constants[62] = 2500.0
    constants[63] = 0.01
    constants[64] = 0.1
    constants[65] = 4.5E-6
    constants[66] = 0.05
    constants[67] = 0.0002
    states[10] = 345.0
    constants[68] = 0.4
    constants[69] = 0.01
    constants[70] = 31.0
    constants[71] = constants[15]
    constants[72] = 3.90000/53.2000
    constants[73] = constants[26]/constants[16]
    constants[74] = constants[28]*constants[48]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = 1.00000/(1.00000+exp((constants[4]-states[0])/constants[5]))
    rates[1] = (algebraic[0]-states[1])/constants[3]
    algebraic[8] = constants[13]*(power(fabs(states[4]/1.00000), 1.0/2))
    algebraic[12] = constants[38]-states[7]
    algebraic[14] = (constants[29]/(constants[30]+states[7]/algebraic[12]))*(states[6]/(constants[31]+states[6]))*(algebraic[8]+constants[32])
    algebraic[10] = (constants[34]*(states[7]/(constants[35]+states[7])))/(1.00000+exp((states[8]-constants[36])/constants[37]))
    rates[7] = constants[33]*(algebraic[14]-algebraic[10])
    algebraic[24] = (constants[56]/states[2])*states[6]*exp(constants[57]*states[8])
    algebraic[23] = (constants[54]*states[8]-constants[55])*(power(states[2], 2.00000))
    algebraic[26] = algebraic[24]-algebraic[23]
    rates[6] = -constants[58]*algebraic[26]
    algebraic[21] = (constants[40]*(states[7]/(constants[41]+states[7])))/(1.00000+exp((states[8]-constants[42])/constants[43]))
    algebraic[25] = constants[59]-states[9]
    algebraic[27] = ((constants[47]*constants[44])/(constants[44]+algebraic[25]))/(1.00000+exp((constants[45]-states[8])/constants[46]))
    algebraic[29] = 3.00000*algebraic[27]
    algebraic[31] = algebraic[25]/states[9]
    algebraic[33] = constants[51]*((algebraic[31]/(algebraic[31]+constants[52]))/exp(-0.500000*constants[53]*states[8]))
    algebraic[22] = constants[49]*states[8]+constants[50]
    rates[8] = (algebraic[21]-(algebraic[29]+algebraic[33]+algebraic[22]+algebraic[24]+2.00000*algebraic[23]))/constants[39]
    rates[9] = constants[33]*(algebraic[33]-algebraic[27])
    algebraic[30] = constants[67]*(states[10]-states[2])
    algebraic[32] = constants[68]*states[2]
    algebraic[34] = algebraic[30]-algebraic[32]
    algebraic[2] = 1.00000/(1.00000+exp((constants[8]-states[0])/constants[9]))
    algebraic[3] = constants[6]*algebraic[2]*(states[0]-constants[7])
    algebraic[28] = -(constants[65]*algebraic[3]+constants[64]*(states[2]-constants[66]))
    rates[2] = constants[63]*(algebraic[28]+algebraic[34]+constants[72]*algebraic[26])
    rates[10] = -constants[69]*constants[70]*algebraic[34]
    algebraic[35] = constants[62]-states[3]
    algebraic[38] = (constants[60]*states[2]+constants[61])*algebraic[35]
    rates[3] = -constants[72]*algebraic[33]+algebraic[38]
    algebraic[1] = constants[1]*states[1]*(states[0]-constants[2])
    algebraic[4] = (constants[10]/(1.00000+power(constants[11]/states[2], 2.00000)))*(states[0]-constants[2])
    algebraic[5] = 0.165000*states[3]
    algebraic[6] = 0.0800000*(1.00000+(2.00000*algebraic[5])/17.0000)+0.890000*(power(algebraic[5]/17.0000, 2.00000))
    algebraic[7] = 0.135000*states[3]
    algebraic[36] = 0.0500000*algebraic[35]
    algebraic[39] = (power(1.00000+algebraic[5]/17.0000, 2.00000))*(1.00000+algebraic[7]/26.0000+algebraic[36]/1.00000)
    algebraic[41] = algebraic[6]/algebraic[39]
    algebraic[43] = constants[12]*algebraic[41]*(states[0]-constants[2])
    rates[0] = -(algebraic[1]+algebraic[3]+algebraic[4]+algebraic[43])/constants[0]
    algebraic[9] = 0.300000*states[5]
    algebraic[11] = (power(algebraic[9], 2.00000))/(constants[18]*1.00000)
    algebraic[13] = constants[71]+algebraic[11]
    algebraic[44] = (power(algebraic[9]*algebraic[35], 2.00000))/(constants[22]*constants[18]*constants[19]*(power(1.00000, 2.00000)))
    algebraic[45] = algebraic[13]+algebraic[44]
    algebraic[47] = algebraic[45]
    algebraic[50] = algebraic[47]
    algebraic[16] = (states[4]*(power(algebraic[9], 2.00000)))/(constants[17]*constants[18]*constants[23]*1.00000)
    algebraic[52] = algebraic[50]+algebraic[16]
    algebraic[54] = (states[4]*(power(algebraic[9], 2.00000))*(power(algebraic[35], 2.00000)))/(constants[17]*constants[18]*constants[19]*constants[23]*constants[24]*constants[22]*(power(1.00000, 2.00000)))
    algebraic[55] = algebraic[52]+algebraic[54]
    algebraic[57] = algebraic[55]
    algebraic[60] = algebraic[57]
    algebraic[17] = (constants[26]*(power(algebraic[9], 2.00000)))/(constants[16]*constants[18]*constants[21]*1.00000)
    algebraic[62] = algebraic[60]+algebraic[17]
    algebraic[64] = (constants[26]*(power(algebraic[9], 2.00000))*(power(algebraic[35], 2.00000)))/(constants[16]*constants[18]*constants[19]*constants[21]*constants[25]*constants[22]*(power(1.00000, 2.00000)))
    algebraic[65] = algebraic[62]+algebraic[64]
    algebraic[67] = algebraic[65]
    algebraic[70] = algebraic[67]
    algebraic[72] = algebraic[70]
    algebraic[74] = (constants[26]*states[4]*(power(algebraic[9], 2.00000))*(power(algebraic[35], 2.00000)))/(constants[16]*constants[17]*constants[18]*constants[19]*constants[23]*constants[21]*constants[24]*constants[25]*constants[22]*(power(1.00000, 2.00000)))
    algebraic[75] = algebraic[72]+algebraic[74]
    algebraic[37] = (power(algebraic[35], 2.00000))/(constants[19]*1.00000)
    algebraic[40] = constants[14]+algebraic[37]
    algebraic[42] = algebraic[40]+algebraic[11]
    algebraic[46] = algebraic[42]+algebraic[44]
    algebraic[15] = states[4]/constants[17]
    algebraic[48] = algebraic[46]+algebraic[15]
    algebraic[49] = (states[4]*(power(algebraic[35], 2.00000)))/(constants[17]*constants[19]*constants[24]*1.00000)
    algebraic[51] = algebraic[48]+algebraic[49]
    algebraic[53] = algebraic[51]+algebraic[16]
    algebraic[56] = algebraic[53]+algebraic[54]
    algebraic[58] = algebraic[56]+constants[73]
    algebraic[59] = (constants[26]*(power(algebraic[35], 2.00000)))/(constants[16]*constants[19]*constants[25]*1.00000)
    algebraic[61] = algebraic[58]+algebraic[59]
    algebraic[63] = algebraic[61]+algebraic[17]
    algebraic[66] = algebraic[63]+algebraic[64]
    algebraic[18] = (constants[26]*states[4])/(constants[16]*constants[17])
    algebraic[68] = algebraic[66]+algebraic[18]
    algebraic[69] = (constants[26]*states[4]*(power(algebraic[35], 2.00000)))/(constants[16]*constants[17]*constants[19]*constants[24]*constants[25]*1.00000)
    algebraic[71] = algebraic[68]+algebraic[69]
    algebraic[19] = (constants[26]*states[4]*(power(algebraic[9], 2.00000)))/(constants[16]*constants[17]*constants[18]*constants[23]*constants[21]*1.00000)
    algebraic[73] = algebraic[71]+algebraic[19]
    algebraic[76] = algebraic[73]+algebraic[74]
    algebraic[20] = algebraic[19]
    algebraic[77] = (constants[27]*constants[20]*algebraic[75]+constants[20]*algebraic[20])/algebraic[76]
    algebraic[78] = constants[28]*algebraic[77]
    rates[5] = constants[74]-algebraic[78]
    rates[4] = algebraic[78]-0.500000*algebraic[8]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = 1.00000/(1.00000+exp((constants[4]-states[0])/constants[5]))
    algebraic[8] = constants[13]*(power(fabs(states[4]/1.00000), 1.0/2))
    algebraic[12] = constants[38]-states[7]
    algebraic[14] = (constants[29]/(constants[30]+states[7]/algebraic[12]))*(states[6]/(constants[31]+states[6]))*(algebraic[8]+constants[32])
    algebraic[10] = (constants[34]*(states[7]/(constants[35]+states[7])))/(1.00000+exp((states[8]-constants[36])/constants[37]))
    algebraic[24] = (constants[56]/states[2])*states[6]*exp(constants[57]*states[8])
    algebraic[23] = (constants[54]*states[8]-constants[55])*(power(states[2], 2.00000))
    algebraic[26] = algebraic[24]-algebraic[23]
    algebraic[21] = (constants[40]*(states[7]/(constants[41]+states[7])))/(1.00000+exp((states[8]-constants[42])/constants[43]))
    algebraic[25] = constants[59]-states[9]
    algebraic[27] = ((constants[47]*constants[44])/(constants[44]+algebraic[25]))/(1.00000+exp((constants[45]-states[8])/constants[46]))
    algebraic[29] = 3.00000*algebraic[27]
    algebraic[31] = algebraic[25]/states[9]
    algebraic[33] = constants[51]*((algebraic[31]/(algebraic[31]+constants[52]))/exp(-0.500000*constants[53]*states[8]))
    algebraic[22] = constants[49]*states[8]+constants[50]
    algebraic[30] = constants[67]*(states[10]-states[2])
    algebraic[32] = constants[68]*states[2]
    algebraic[34] = algebraic[30]-algebraic[32]
    algebraic[2] = 1.00000/(1.00000+exp((constants[8]-states[0])/constants[9]))
    algebraic[3] = constants[6]*algebraic[2]*(states[0]-constants[7])
    algebraic[28] = -(constants[65]*algebraic[3]+constants[64]*(states[2]-constants[66]))
    algebraic[35] = constants[62]-states[3]
    algebraic[38] = (constants[60]*states[2]+constants[61])*algebraic[35]
    algebraic[1] = constants[1]*states[1]*(states[0]-constants[2])
    algebraic[4] = (constants[10]/(1.00000+power(constants[11]/states[2], 2.00000)))*(states[0]-constants[2])
    algebraic[5] = 0.165000*states[3]
    algebraic[6] = 0.0800000*(1.00000+(2.00000*algebraic[5])/17.0000)+0.890000*(power(algebraic[5]/17.0000, 2.00000))
    algebraic[7] = 0.135000*states[3]
    algebraic[36] = 0.0500000*algebraic[35]
    algebraic[39] = (power(1.00000+algebraic[5]/17.0000, 2.00000))*(1.00000+algebraic[7]/26.0000+algebraic[36]/1.00000)
    algebraic[41] = algebraic[6]/algebraic[39]
    algebraic[43] = constants[12]*algebraic[41]*(states[0]-constants[2])
    algebraic[9] = 0.300000*states[5]
    algebraic[11] = (power(algebraic[9], 2.00000))/(constants[18]*1.00000)
    algebraic[13] = constants[71]+algebraic[11]
    algebraic[44] = (power(algebraic[9]*algebraic[35], 2.00000))/(constants[22]*constants[18]*constants[19]*(power(1.00000, 2.00000)))
    algebraic[45] = algebraic[13]+algebraic[44]
    algebraic[47] = algebraic[45]
    algebraic[50] = algebraic[47]
    algebraic[16] = (states[4]*(power(algebraic[9], 2.00000)))/(constants[17]*constants[18]*constants[23]*1.00000)
    algebraic[52] = algebraic[50]+algebraic[16]
    algebraic[54] = (states[4]*(power(algebraic[9], 2.00000))*(power(algebraic[35], 2.00000)))/(constants[17]*constants[18]*constants[19]*constants[23]*constants[24]*constants[22]*(power(1.00000, 2.00000)))
    algebraic[55] = algebraic[52]+algebraic[54]
    algebraic[57] = algebraic[55]
    algebraic[60] = algebraic[57]
    algebraic[17] = (constants[26]*(power(algebraic[9], 2.00000)))/(constants[16]*constants[18]*constants[21]*1.00000)
    algebraic[62] = algebraic[60]+algebraic[17]
    algebraic[64] = (constants[26]*(power(algebraic[9], 2.00000))*(power(algebraic[35], 2.00000)))/(constants[16]*constants[18]*constants[19]*constants[21]*constants[25]*constants[22]*(power(1.00000, 2.00000)))
    algebraic[65] = algebraic[62]+algebraic[64]
    algebraic[67] = algebraic[65]
    algebraic[70] = algebraic[67]
    algebraic[72] = algebraic[70]
    algebraic[74] = (constants[26]*states[4]*(power(algebraic[9], 2.00000))*(power(algebraic[35], 2.00000)))/(constants[16]*constants[17]*constants[18]*constants[19]*constants[23]*constants[21]*constants[24]*constants[25]*constants[22]*(power(1.00000, 2.00000)))
    algebraic[75] = algebraic[72]+algebraic[74]
    algebraic[37] = (power(algebraic[35], 2.00000))/(constants[19]*1.00000)
    algebraic[40] = constants[14]+algebraic[37]
    algebraic[42] = algebraic[40]+algebraic[11]
    algebraic[46] = algebraic[42]+algebraic[44]
    algebraic[15] = states[4]/constants[17]
    algebraic[48] = algebraic[46]+algebraic[15]
    algebraic[49] = (states[4]*(power(algebraic[35], 2.00000)))/(constants[17]*constants[19]*constants[24]*1.00000)
    algebraic[51] = algebraic[48]+algebraic[49]
    algebraic[53] = algebraic[51]+algebraic[16]
    algebraic[56] = algebraic[53]+algebraic[54]
    algebraic[58] = algebraic[56]+constants[73]
    algebraic[59] = (constants[26]*(power(algebraic[35], 2.00000)))/(constants[16]*constants[19]*constants[25]*1.00000)
    algebraic[61] = algebraic[58]+algebraic[59]
    algebraic[63] = algebraic[61]+algebraic[17]
    algebraic[66] = algebraic[63]+algebraic[64]
    algebraic[18] = (constants[26]*states[4])/(constants[16]*constants[17])
    algebraic[68] = algebraic[66]+algebraic[18]
    algebraic[69] = (constants[26]*states[4]*(power(algebraic[35], 2.00000)))/(constants[16]*constants[17]*constants[19]*constants[24]*constants[25]*1.00000)
    algebraic[71] = algebraic[68]+algebraic[69]
    algebraic[19] = (constants[26]*states[4]*(power(algebraic[9], 2.00000)))/(constants[16]*constants[17]*constants[18]*constants[23]*constants[21]*1.00000)
    algebraic[73] = algebraic[71]+algebraic[19]
    algebraic[76] = algebraic[73]+algebraic[74]
    algebraic[20] = algebraic[19]
    algebraic[77] = (constants[27]*constants[20]*algebraic[75]+constants[20]*algebraic[20])/algebraic[76]
    algebraic[78] = constants[28]*algebraic[77]
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