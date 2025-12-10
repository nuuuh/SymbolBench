# Size of variable arrays:
sizeAlgebraic = 34
sizeStates = 15
sizeConstants = 59
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_algebraic[0] = "L_iso in component beta_1_adrenergic_parameters (uM)"
    legend_constants[0] = "K_H in component beta_1_adrenergic_parameters (uM)"
    legend_constants[1] = "K_L in component beta_1_adrenergic_parameters (uM)"
    legend_constants[2] = "K_C in component beta_1_adrenergic_parameters (uM)"
    legend_algebraic[1] = "L_ach in component muscarinic_parameters (uM)"
    legend_constants[3] = "K_H in component muscarinic_parameters (uM)"
    legend_constants[4] = "K_L in component muscarinic_parameters (uM)"
    legend_constants[5] = "K_C in component muscarinic_parameters (uM)"
    legend_constants[6] = "k_PDE2 in component PDE_parameters (per_sec)"
    legend_constants[7] = "Km_PDE2 in component PDE_parameters (uM)"
    legend_constants[8] = "k_PDE3 in component PDE_parameters (per_sec)"
    legend_constants[9] = "Km_PDE3 in component PDE_parameters (uM)"
    legend_constants[10] = "k_PDE4 in component PDE_parameters (per_sec)"
    legend_constants[11] = "Km_PDE4 in component PDE_parameters (uM)"
    legend_constants[12] = "k_act1 in component G_s_parameters (per_sec)"
    legend_constants[13] = "k_act2 in component G_s_parameters (per_sec)"
    legend_constants[14] = "k_hydr in component G_s_parameters (per_sec)"
    legend_constants[15] = "k_reas in component G_s_parameters (per_uM_per_sec)"
    legend_constants[16] = "k_act1 in component G_i_parameters (per_sec)"
    legend_constants[17] = "k_act2 in component G_i_parameters (per_sec)"
    legend_constants[18] = "k_hydr in component G_i_parameters (per_sec)"
    legend_constants[19] = "k_reas in component G_i_parameters (per_uM_per_sec)"
    legend_algebraic[12] = "R in component caveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[13] = "LR in component caveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[14] = "LRG in component caveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[15] = "RG in component caveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_constants[20] = "R_Total in component caveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[5] = "Gs_alpha_beta_gamma in component caveolar_G_s_protein_activation_module (uM)"
    legend_algebraic[16] = "R in component caveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[17] = "LR in component caveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[18] = "LRG in component caveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[19] = "RG in component caveolar_muscarinic_receptor_module (uM)"
    legend_constants[21] = "R_Total in component caveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[6] = "Gi_alpha_beta_gamma in component caveolar_G_i_protein_activation_module (uM)"
    legend_states[0] = "Gs_alpha_GTP in component caveolar_G_s_protein_activation_module (uM)"
    legend_states[1] = "Gs_beta_gamma in component caveolar_G_s_protein_activation_module (uM)"
    legend_states[2] = "Gs_alpha_GDP in component caveolar_G_s_protein_activation_module (uM)"
    legend_constants[22] = "Gs_Total in component caveolar_G_s_protein_activation_module (uM)"
    legend_states[3] = "Gi_alpha_GTP in component caveolar_G_i_protein_activation_module (uM)"
    legend_states[4] = "Gi_beta_gamma in component caveolar_G_i_protein_activation_module (uM)"
    legend_states[5] = "Gi_alpha_GDP in component caveolar_G_i_protein_activation_module (uM)"
    legend_constants[23] = "Gi_Total in component caveolar_G_i_protein_activation_module (uM)"
    legend_algebraic[20] = "R in component extracaveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[21] = "LR in component extracaveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[22] = "LRG in component extracaveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[23] = "RG in component extracaveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_constants[24] = "R_Total in component extracaveolar_beta_1_adrenergic_receptor_module (uM)"
    legend_algebraic[7] = "Gs_alpha_beta_gamma in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_algebraic[24] = "R in component extracaveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[25] = "LR in component extracaveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[26] = "LRG in component extracaveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[27] = "RG in component extracaveolar_muscarinic_receptor_module (uM)"
    legend_constants[25] = "R_Total in component extracaveolar_muscarinic_receptor_module (uM)"
    legend_algebraic[8] = "Gi_alpha_beta_gamma in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_states[6] = "Gs_alpha_GTP in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_states[7] = "Gs_beta_gamma in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_states[8] = "Gs_alpha_GDP in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_constants[26] = "Gs_Total in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_states[9] = "Gi_alpha_GTP in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_states[10] = "Gi_beta_gamma in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_states[11] = "Gi_alpha_GDP in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_constants[27] = "Gi_Total in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_algebraic[9] = "dcAMP_AC_56_dt in component AC56_module (uM_per_sec)"
    legend_algebraic[2] = "k_AC56 in component AC56_module (per_sec)"
    legend_constants[28] = "AC_56 in component AC56_module (uM)"
    legend_constants[29] = "AF56 in component AC56_module (dimensionless)"
    legend_constants[30] = "MW_AC56 in component AC56_module (kDa)"
    legend_constants[31] = "ATP in component AC56_module (uM)"
    legend_constants[32] = "Km_ATP in component AC56_module (uM)"
    legend_algebraic[10] = "dcAMP_AC_47_ecav_dt in component AC47_ecav_module (uM_per_sec)"
    legend_algebraic[3] = "k_AC47_ecav in component AC47_ecav_module (per_sec)"
    legend_constants[33] = "AC_47_ecav in component AC47_ecav_module (uM)"
    legend_constants[34] = "AF47 in component AC47_ecav_module (dimensionless)"
    legend_constants[35] = "MW_AC47 in component AC47_ecav_module (kDa)"
    legend_constants[36] = "ATP in component AC47_ecav_module (uM)"
    legend_constants[37] = "Km_ATP in component AC47_ecav_module (uM)"
    legend_constants[55] = "dcAMP_AC_47_cyt_dt in component AC47_cyt_module (uM_per_sec)"
    legend_constants[38] = "k_AC47_cyt in component AC47_cyt_module (per_sec)"
    legend_constants[39] = "AC_47_cyt in component AC47_cyt_module (uM)"
    legend_constants[40] = "AF47 in component AC47_cyt_module (dimensionless)"
    legend_constants[41] = "ATP in component AC47_cyt_module (uM)"
    legend_constants[42] = "Km_ATP in component AC47_cyt_module (uM)"
    legend_algebraic[28] = "dcAMP_cav_PDE2_dt in component caveolar_PDE_module (uM_per_sec)"
    legend_algebraic[31] = "dcAMP_cav_PDE3_dt in component caveolar_PDE_module (uM_per_sec)"
    legend_algebraic[33] = "dcAMP_cav_PDE4_dt in component caveolar_PDE_module (uM_per_sec)"
    legend_states[12] = "cAMP_cav in component cAMP_flux_module (uM)"
    legend_constants[43] = "PDE2 in component caveolar_PDE_module (uM)"
    legend_constants[44] = "PDE3 in component caveolar_PDE_module (uM)"
    legend_constants[45] = "PDE4 in component caveolar_PDE_module (uM)"
    legend_algebraic[29] = "dcAMP_ecav_PDE2_dt in component extracaveolar_PDE_module (uM_per_sec)"
    legend_algebraic[32] = "dcAMP_ecav_PDE4_dt in component extracaveolar_PDE_module (uM_per_sec)"
    legend_states[13] = "cAMP_ecav in component cAMP_flux_module (uM)"
    legend_constants[46] = "PDE2 in component extracaveolar_PDE_module (uM)"
    legend_constants[47] = "PDE4 in component extracaveolar_PDE_module (uM)"
    legend_algebraic[4] = "dcAMP_cyt_PDE2_dt in component bulk_cytoplasmic_PDE_module (uM_per_sec)"
    legend_algebraic[11] = "dcAMP_cyt_PDE3_dt in component bulk_cytoplasmic_PDE_module (uM_per_sec)"
    legend_algebraic[30] = "dcAMP_cyt_PDE4_dt in component bulk_cytoplasmic_PDE_module (uM_per_sec)"
    legend_states[14] = "cAMP_cyt in component cAMP_flux_module (uM)"
    legend_constants[48] = "PDE2 in component bulk_cytoplasmic_PDE_module (uM)"
    legend_constants[49] = "PDE3 in component bulk_cytoplasmic_PDE_module (uM)"
    legend_constants[50] = "PDE4 in component bulk_cytoplasmic_PDE_module (uM)"
    legend_constants[56] = "V_cav in component cAMP_flux_module (liter)"
    legend_constants[57] = "V_ecav in component cAMP_flux_module (liter)"
    legend_constants[58] = "V_cyt in component cAMP_flux_module (liter)"
    legend_constants[51] = "V_cell in component cAMP_flux_module (liter)"
    legend_constants[52] = "J_cav_ecav in component cAMP_flux_module (liters_per_second)"
    legend_constants[53] = "J_cav_cyt in component cAMP_flux_module (liters_per_second)"
    legend_constants[54] = "J_ecav_cyt in component cAMP_flux_module (liters_per_second)"
    legend_rates[0] = "d/dt Gs_alpha_GTP in component caveolar_G_s_protein_activation_module (uM)"
    legend_rates[1] = "d/dt Gs_beta_gamma in component caveolar_G_s_protein_activation_module (uM)"
    legend_rates[2] = "d/dt Gs_alpha_GDP in component caveolar_G_s_protein_activation_module (uM)"
    legend_rates[3] = "d/dt Gi_alpha_GTP in component caveolar_G_i_protein_activation_module (uM)"
    legend_rates[4] = "d/dt Gi_beta_gamma in component caveolar_G_i_protein_activation_module (uM)"
    legend_rates[5] = "d/dt Gi_alpha_GDP in component caveolar_G_i_protein_activation_module (uM)"
    legend_rates[6] = "d/dt Gs_alpha_GTP in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_rates[7] = "d/dt Gs_beta_gamma in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_rates[8] = "d/dt Gs_alpha_GDP in component extracaveolar_G_s_protein_activation_module (uM)"
    legend_rates[9] = "d/dt Gi_alpha_GTP in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_rates[10] = "d/dt Gi_beta_gamma in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_rates[11] = "d/dt Gi_alpha_GDP in component extracaveolar_G_i_protein_activation_module (uM)"
    legend_rates[12] = "d/dt cAMP_cav in component cAMP_flux_module (uM)"
    legend_rates[13] = "d/dt cAMP_ecav in component cAMP_flux_module (uM)"
    legend_rates[14] = "d/dt cAMP_cyt in component cAMP_flux_module (uM)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 0.035
    constants[1] = 0.386
    constants[2] = 8.809
    constants[3] = 0.16
    constants[4] = 11
    constants[5] = 30
    constants[6] = 20
    constants[7] = 50
    constants[8] = 1.25
    constants[9] = 0.08
    constants[10] = 2.5
    constants[11] = 2.2
    constants[12] = 5
    constants[13] = 0.1
    constants[14] = 0.8
    constants[15] = 1.21e3
    constants[16] = 2.5
    constants[17] = 0.05
    constants[18] = 0.8
    constants[19] = 1.21e3
    constants[20] = 0.633
    constants[21] = 0.633
    states[0] = 0.041983438
    states[1] = 0.042634499
    states[2] = 0.000651061
    constants[22] = 10
    states[3] = 0.012644961
    states[4] = 0.013274751
    states[5] = 0.00062979
    constants[23] = 20
    constants[24] = 0.633
    constants[25] = 0.633
    states[6] = 0.083866891
    states[7] = 0.084522918
    states[8] = 0.000656025
    constants[26] = 10
    states[9] = 0.001018705
    states[10] = 0.001475253
    states[11] = 0.000456548
    constants[27] = 1
    constants[28] = 3.379
    constants[29] = 500
    constants[30] = 130
    constants[31] = 5000
    constants[32] = 315
    constants[33] = 0.2
    constants[34] = 130
    constants[35] = 130
    constants[36] = 5000
    constants[37] = 315
    constants[38] = 1.08e-3
    constants[39] = 0.136
    constants[40] = 130
    constants[41] = 5000
    constants[42] = 315
    states[12] = 0.11750433
    constants[43] = 4.5
    constants[44] = 5.6
    constants[45] = 2
    states[13] = 1.092200547
    constants[46] = 0.02
    constants[47] = 0.16
    states[14] = 0.992583576
    constants[48] = 5e-3
    constants[49] = 7.5e-3
    constants[50] = 5e-3
    constants[51] = 38e-12
    constants[52] = 7.5e-15
    constants[53] = 7.5e-14
    constants[54] = 1.5e-17
    constants[55] = (constants[38]*constants[39]*constants[40]*constants[41])/(constants[42]+constants[41])
    constants[56] = 0.0100000*constants[51]
    constants[57] = 0.0200000*constants[51]
    constants[58] = 0.500000*constants[51]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[2] = states[0]*constants[14]-states[2]*states[1]*constants[15]
    rates[5] = states[3]*constants[18]-states[5]*states[4]*constants[19]
    rates[8] = states[6]*constants[14]-states[8]*states[7]*constants[15]
    rates[11] = states[9]*constants[18]-states[11]*states[10]*constants[19]
    algebraic[0] = custom_piecewise([greater(voi , 120.000) & less_equal(voi , 720.000), 1.00000 , True, 1.00000])
    algebraic[5] = (constants[22]-states[0])-states[2]
    rootfind_0(voi, constants, rates, states, algebraic)
    rates[0] = (algebraic[15]*constants[13]+algebraic[14]*constants[12])-states[0]*constants[14]
    rates[1] = (algebraic[15]*constants[13]+algebraic[14]*constants[12])-states[2]*states[1]*constants[15]
    algebraic[1] = custom_piecewise([greater(voi , 240.000) & less_equal(voi , 540.000), 0.00000 , True, 0.00000])
    algebraic[6] = (constants[23]-states[3])-states[5]
    rootfind_1(voi, constants, rates, states, algebraic)
    rates[3] = (algebraic[19]*constants[17]+algebraic[18]*constants[16])-states[3]*constants[18]
    rates[4] = (algebraic[19]*constants[17]+algebraic[18]*constants[16])-states[5]*states[4]*constants[19]
    algebraic[7] = (constants[26]-states[6])-states[8]
    rootfind_2(voi, constants, rates, states, algebraic)
    rates[6] = (algebraic[23]*constants[13]+algebraic[22]*constants[12])-states[6]*constants[14]
    rates[7] = (algebraic[23]*constants[13]+algebraic[22]*constants[12])-states[8]*states[7]*constants[15]
    algebraic[8] = (constants[27]-states[9])-states[11]
    rootfind_3(voi, constants, rates, states, algebraic)
    rates[9] = (algebraic[27]*constants[17]+algebraic[26]*constants[16])-states[9]*constants[18]
    rates[10] = (algebraic[27]*constants[17]+algebraic[26]*constants[16])-states[11]*states[10]*constants[19]
    algebraic[4] = (constants[6]*constants[48]*states[14])/(constants[7]+states[14])
    algebraic[11] = (constants[8]*constants[49]*states[14])/(constants[9]+states[14])
    algebraic[30] = (constants[10]*constants[50]*states[14])/(constants[11]+states[14])
    rates[14] = (constants[55]-(algebraic[4]+algebraic[11]+algebraic[30]))+(constants[53]*(states[12]-states[14]))/constants[58]+(constants[54]*(states[13]-states[14]))/constants[58]
    algebraic[3] = (((0.0630000+(2.01000*(power(states[6]*1000.00, 1.00430)))/(31.5440+power(states[6]*1000.00, 1.00430)))*(1.00000+((1.00000/3.01000)*49.1000*(power(states[10]*1000.00, 0.892100)))/(25.4400+power(states[10]*1000.00, 0.892100)))*constants[35])/60.0000)*0.00100000
    algebraic[10] = (algebraic[3]*constants[33]*constants[34]*constants[36])/(constants[37]+constants[36])
    algebraic[29] = (constants[6]*constants[46]*states[13])/(constants[7]+states[13])
    algebraic[32] = (constants[10]*constants[47]*states[13])/(constants[11]+states[13])
    rates[13] = ((algebraic[10]-(algebraic[29]+algebraic[32]))+(constants[52]*(states[12]-states[13]))/constants[57])-(constants[54]*(states[13]-states[14]))/constants[57]
    algebraic[2] = (((0.700000+(3.82340*(power(states[0]/1.00000, 0.978700)))/(0.198600+power(states[0]/1.00000, 0.978700)))*(1.00000+((1.00000/1.44320)*-1.00610*(power(states[3]/1.00000, 0.835600)))/(0.191800+power(states[3]/1.00000, 0.835600)))*constants[30])/60.0000)*0.00100000
    algebraic[9] = (algebraic[2]*constants[28]*constants[29]*constants[31])/(constants[32]+constants[31])
    algebraic[28] = (constants[6]*constants[43]*states[12])/(constants[7]+states[12])
    algebraic[31] = (constants[8]*constants[44]*states[12])/(constants[9]+states[12])
    algebraic[33] = (constants[10]*constants[45]*states[12])/(constants[11]+states[12])
    rates[12] = ((algebraic[9]-(algebraic[28]+algebraic[31]+algebraic[33]))-(constants[52]*(states[12]-states[13]))/constants[56])-(constants[53]*(states[12]-states[14]))/constants[56]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater(voi , 120.000) & less_equal(voi , 720.000), 1.00000 , True, 1.00000])
    algebraic[5] = (constants[22]-states[0])-states[2]
    algebraic[1] = custom_piecewise([greater(voi , 240.000) & less_equal(voi , 540.000), 0.00000 , True, 0.00000])
    algebraic[6] = (constants[23]-states[3])-states[5]
    algebraic[7] = (constants[26]-states[6])-states[8]
    algebraic[8] = (constants[27]-states[9])-states[11]
    algebraic[4] = (constants[6]*constants[48]*states[14])/(constants[7]+states[14])
    algebraic[11] = (constants[8]*constants[49]*states[14])/(constants[9]+states[14])
    algebraic[30] = (constants[10]*constants[50]*states[14])/(constants[11]+states[14])
    algebraic[3] = (((0.0630000+(2.01000*(power(states[6]*1000.00, 1.00430)))/(31.5440+power(states[6]*1000.00, 1.00430)))*(1.00000+((1.00000/3.01000)*49.1000*(power(states[10]*1000.00, 0.892100)))/(25.4400+power(states[10]*1000.00, 0.892100)))*constants[35])/60.0000)*0.00100000
    algebraic[10] = (algebraic[3]*constants[33]*constants[34]*constants[36])/(constants[37]+constants[36])
    algebraic[29] = (constants[6]*constants[46]*states[13])/(constants[7]+states[13])
    algebraic[32] = (constants[10]*constants[47]*states[13])/(constants[11]+states[13])
    algebraic[2] = (((0.700000+(3.82340*(power(states[0]/1.00000, 0.978700)))/(0.198600+power(states[0]/1.00000, 0.978700)))*(1.00000+((1.00000/1.44320)*-1.00610*(power(states[3]/1.00000, 0.835600)))/(0.191800+power(states[3]/1.00000, 0.835600)))*constants[30])/60.0000)*0.00100000
    algebraic[9] = (algebraic[2]*constants[28]*constants[29]*constants[31])/(constants[32]+constants[31])
    algebraic[28] = (constants[6]*constants[43]*states[12])/(constants[7]+states[12])
    algebraic[31] = (constants[8]*constants[44]*states[12])/(constants[9]+states[12])
    algebraic[33] = (constants[10]*constants[45]*states[12])/(constants[11]+states[12])
    return algebraic

initialGuess0 = None
def rootfind_0(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess0
    if initialGuess0 is None: initialGuess0 = ones(4)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_0, initialGuess0, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess0 = soln
        algebraic[12] = soln[0]
        algebraic[13] = soln[1]
        algebraic[14] = soln[2]
        algebraic[15] = soln[3]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess0 = soln
            algebraic[12][i] = soln[0]
            algebraic[13][i] = soln[1]
            algebraic[14][i] = soln[2]
            algebraic[15][i] = soln[3]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 4)
    algebraic[12] = algebraicCandidate[0]
    algebraic[13] = algebraicCandidate[1]
    algebraic[14] = algebraicCandidate[2]
    algebraic[15] = algebraicCandidate[3]
    resid[0] = (algebraic[12]-(((constants[20]-algebraic[13])-algebraic[14])-algebraic[15]))
    resid[1] = (algebraic[13]-(algebraic[0]*algebraic[12])/constants[1])
    resid[2] = (algebraic[14]-(algebraic[0]*algebraic[12]*algebraic[5])/(constants[0]*constants[2]))
    resid[3] = (algebraic[15]-(algebraic[12]*algebraic[5])/constants[2])
    return resid

initialGuess1 = None
def rootfind_1(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess1
    if initialGuess1 is None: initialGuess1 = ones(4)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_1, initialGuess1, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess1 = soln
        algebraic[16] = soln[0]
        algebraic[17] = soln[1]
        algebraic[18] = soln[2]
        algebraic[19] = soln[3]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_1, initialGuess1, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess1 = soln
            algebraic[16][i] = soln[0]
            algebraic[17][i] = soln[1]
            algebraic[18][i] = soln[2]
            algebraic[19][i] = soln[3]

def residualSN_1(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 4)
    algebraic[16] = algebraicCandidate[0]
    algebraic[17] = algebraicCandidate[1]
    algebraic[18] = algebraicCandidate[2]
    algebraic[19] = algebraicCandidate[3]
    resid[0] = (algebraic[16]-(((constants[21]-algebraic[17])-algebraic[18])-algebraic[19]))
    resid[1] = (algebraic[17]-(algebraic[1]*algebraic[16])/constants[4])
    resid[2] = (algebraic[18]-(algebraic[1]*algebraic[16]*algebraic[6])/(constants[3]*constants[5]))
    resid[3] = (algebraic[19]-(algebraic[16]*algebraic[6])/constants[5])
    return resid

initialGuess2 = None
def rootfind_2(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess2
    if initialGuess2 is None: initialGuess2 = ones(4)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_2, initialGuess2, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess2 = soln
        algebraic[20] = soln[0]
        algebraic[21] = soln[1]
        algebraic[22] = soln[2]
        algebraic[23] = soln[3]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_2, initialGuess2, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess2 = soln
            algebraic[20][i] = soln[0]
            algebraic[21][i] = soln[1]
            algebraic[22][i] = soln[2]
            algebraic[23][i] = soln[3]

def residualSN_2(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 4)
    algebraic[20] = algebraicCandidate[0]
    algebraic[21] = algebraicCandidate[1]
    algebraic[22] = algebraicCandidate[2]
    algebraic[23] = algebraicCandidate[3]
    resid[0] = (algebraic[20]-(((constants[24]-algebraic[21])-algebraic[22])-algebraic[23]))
    resid[1] = (algebraic[21]-(algebraic[0]*algebraic[20])/constants[1])
    resid[2] = (algebraic[22]-(algebraic[0]*algebraic[20]*algebraic[7])/(constants[0]*constants[2]))
    resid[3] = (algebraic[23]-(algebraic[20]*algebraic[7])/constants[2])
    return resid

initialGuess3 = None
def rootfind_3(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess3
    if initialGuess3 is None: initialGuess3 = ones(4)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_3, initialGuess3, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess3 = soln
        algebraic[24] = soln[0]
        algebraic[25] = soln[1]
        algebraic[26] = soln[2]
        algebraic[27] = soln[3]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_3, initialGuess3, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess3 = soln
            algebraic[24][i] = soln[0]
            algebraic[25][i] = soln[1]
            algebraic[26][i] = soln[2]
            algebraic[27][i] = soln[3]

def residualSN_3(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 4)
    algebraic[24] = algebraicCandidate[0]
    algebraic[25] = algebraicCandidate[1]
    algebraic[26] = algebraicCandidate[2]
    algebraic[27] = algebraicCandidate[3]
    resid[0] = (algebraic[24]-(((constants[25]-algebraic[25])-algebraic[26])-algebraic[27]))
    resid[1] = (algebraic[25]-(algebraic[1]*algebraic[24])/constants[4])
    resid[2] = (algebraic[26]-(algebraic[1]*algebraic[24]*algebraic[8])/(constants[3]*constants[5]))
    resid[3] = (algebraic[27]-(algebraic[24]*algebraic[8])/constants[5])
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