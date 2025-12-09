# Size of variable arrays:
sizeAlgebraic = 40
sizeStates = 7
sizeConstants = 38
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
    legend_constants[0] = "V_T in component membrane (millivolt)"
    legend_constants[1] = "V_S in component membrane (millivolt)"
    legend_constants[2] = "C_m in component membrane (mF_per_cm_squared)"
    legend_constants[3] = "F in component membrane (coulomb_per_mole)"
    legend_constants[4] = "R in component membrane (joule_per_mole_per_kelvin)"
    legend_constants[5] = "T in component membrane (kelvin)"
    legend_algebraic[16] = "I_leak in component I_leak (milliampere_per_cm_squared)"
    legend_algebraic[21] = "I_Na in component I_Na (milliampere_per_cm_squared)"
    legend_algebraic[22] = "I_KD in component I_KD (milliampere_per_cm_squared)"
    legend_algebraic[23] = "I_KM in component I_KM (milliampere_per_cm_squared)"
    legend_algebraic[25] = "I_CaL in component I_CaL (milliampere_per_cm_squared)"
    legend_algebraic[39] = "I_h in component I_h (milliampere_per_cm_squared)"
    legend_algebraic[11] = "I_app in component stimulus_protocol (milliampere_per_cm_squared)"
    legend_constants[6] = "i_stimStart in component stimulus_protocol (second)"
    legend_constants[7] = "i_stimEnd in component stimulus_protocol (second)"
    legend_constants[8] = "i_stimAmplitude in component stimulus_protocol (milliampere_per_cm_squared)"
    legend_algebraic[5] = "tau in component stimulus_protocol (second)"
    legend_constants[9] = "period in component stimulus_protocol (second)"
    legend_constants[10] = "g_leak in component I_leak (millisiemens_per_cm_squared)"
    legend_constants[11] = "E_leak in component I_leak (millivolt)"
    legend_constants[12] = "g_Na in component I_Na (millisiemens_per_cm_squared)"
    legend_constants[13] = "E_Na in component I_Na (millivolt)"
    legend_states[1] = "m in component Na_m_gate (dimensionless)"
    legend_states[2] = "h in component Na_h_gate (dimensionless)"
    legend_algebraic[0] = "alpha in component Na_m_gate (per_second)"
    legend_algebraic[6] = "beta in component Na_m_gate (per_second)"
    legend_algebraic[12] = "tau_m in component Na_m_gate (second)"
    legend_algebraic[17] = "m_inf in component Na_m_gate (dimensionless)"
    legend_algebraic[1] = "alpha_h in component Na_h_gate (per_second)"
    legend_algebraic[7] = "beta_h in component Na_h_gate (per_second)"
    legend_algebraic[18] = "h_inf in component Na_h_gate (dimensionless)"
    legend_algebraic[13] = "tau_h in component Na_h_gate (second)"
    legend_constants[14] = "g_KD in component I_KD (millisiemens_per_cm_squared)"
    legend_constants[15] = "E_K in component I_KD (millivolt)"
    legend_states[3] = "n in component KD_n_gate (dimensionless)"
    legend_algebraic[2] = "alpha_n in component KD_n_gate (per_second)"
    legend_algebraic[8] = "beta_n in component KD_n_gate (per_second)"
    legend_algebraic[14] = "tau_n in component KD_n_gate (second)"
    legend_algebraic[19] = "n_inf in component KD_n_gate (dimensionless)"
    legend_constants[16] = "g_KM in component I_KM (millisiemens_per_cm_squared)"
    legend_states[4] = "p in component KM_p_gate (dimensionless)"
    legend_algebraic[3] = "p_inf in component KM_p_gate (dimensionless)"
    legend_algebraic[9] = "tau_p in component KM_p_gate (second)"
    legend_constants[17] = "tau_max in component KM_p_gate (second)"
    legend_constants[18] = "P_Ca in component I_CaL (cm_per_second)"
    legend_algebraic[24] = "G in component G_nonlin (coulomb_per_cm_cubed)"
    legend_states[5] = "q in component CaL_q_gate (dimensionless)"
    legend_algebraic[4] = "alpha_q in component CaL_q_gate (per_second)"
    legend_algebraic[10] = "beta_q in component CaL_q_gate (per_second)"
    legend_algebraic[15] = "tau_q in component CaL_q_gate (second)"
    legend_algebraic[20] = "q_inf in component CaL_q_gate (dimensionless)"
    legend_constants[19] = "Z in component G_nonlin (dimensionless)"
    legend_constants[20] = "Ca_o in component G_nonlin (mM)"
    legend_states[6] = "Ca_i in component dCa_i_dt (mM)"
    legend_constants[21] = "Ca_inf in component dCa_i_dt (mM)"
    legend_constants[22] = "tau_r in component dCa_i_dt (second)"
    legend_constants[23] = "d in component dCa_i_dt (centimeter)"
    legend_algebraic[26] = "drive_channel in component dCa_i_dt (mM_per_second)"
    legend_constants[24] = "k in component dCa_i_dt (fixer)"
    legend_algebraic[38] = "m in component I_h (dimensionless)"
    legend_constants[25] = "E_h in component I_h (millivolt)"
    legend_constants[26] = "g_hbar in component I_h (millisiemens_per_cm_squared)"
    legend_constants[27] = "cac in component I_h (mM)"
    legend_constants[28] = "V_S in component I_h (millivolt)"
    legend_algebraic[35] = "o_1 in component kinetic (dimensionless)"
    legend_algebraic[36] = "o_2 in component kinetic (dimensionless)"
    legend_constants[29] = "g_inc in component I_h (dimensionless)"
    legend_algebraic[32] = "p_0 in component kinetic (dimensionless)"
    legend_algebraic[33] = "p_1 in component kinetic (dimensionless)"
    legend_algebraic[37] = "c_1 in component kinetic (dimensionless)"
    legend_algebraic[29] = "alpha in component rate_constants (dimensionless)"
    legend_algebraic[30] = "beta in component rate_constants (dimensionless)"
    legend_algebraic[31] = "k_1Ca in component rate_constants (per_second)"
    legend_constants[30] = "k_2 in component rate_constants (per_second)"
    legend_algebraic[34] = "k_3p in component rate_constants (per_second)"
    legend_constants[31] = "k_4 in component rate_constants (per_second)"
    legend_algebraic[27] = "h_inf in component rate_constants (second)"
    legend_algebraic[28] = "tau_s in component rate_constants (second)"
    legend_constants[32] = "P_c in component rate_constants (dimensionless)"
    legend_constants[33] = "n_Ca in component rate_constants (dimensionless)"
    legend_constants[34] = "n_exp in component rate_constants (dimensionless)"
    legend_constants[35] = "p_C in component rate_constants (dimensionless)"
    legend_constants[36] = "Ca_c in component rate_constants (mM)"
    legend_constants[37] = "tau_m in component rate_constants (second)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt m in component Na_m_gate (dimensionless)"
    legend_rates[2] = "d/dt h in component Na_h_gate (dimensionless)"
    legend_rates[3] = "d/dt n in component KD_n_gate (dimensionless)"
    legend_rates[4] = "d/dt p in component KM_p_gate (dimensionless)"
    legend_rates[5] = "d/dt q in component CaL_q_gate (dimensionless)"
    legend_rates[6] = "d/dt Ca_i in component dCa_i_dt (mM)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -70
    constants[0] = -55
    constants[1] = 0
    constants[2] = 1e-3
    constants[3] = 96489
    constants[4] = 8.314
    constants[5] = 296.65
    constants[6] = 5
    constants[7] = 9
    constants[8] = -0.3
    constants[9] = 9
    constants[10] = 1
    constants[11] = -70
    constants[12] = 70
    constants[13] = 50
    states[1] = 0
    states[2] = 0
    constants[14] = 7
    constants[15] = -95
    states[3] = 0
    constants[16] = 0.004
    states[4] = 0
    constants[17] = 4
    constants[18] = 2.76e-4
    states[5] = 0.00247262
    constants[19] = 2
    constants[20] = 2
    states[6] = 100e-6
    constants[21] = 100e-6
    constants[22] = 17e-3
    constants[23] = 1e-4
    constants[24] = 0.1
    constants[25] = -20
    constants[26] = 0.02
    constants[27] = 0.006
    constants[28] = 0
    constants[29] = 2
    constants[30] = 0.1
    constants[31] = 1
    constants[32] = 0.01
    constants[33] = 4
    constants[34] = 1
    constants[35] = 0.01
    constants[36] = 0.006
    constants[37] = 20e-3
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[3] = custom_piecewise([less(-(states[0]+35.0000)/10.0000 , 25.0000) & greater(-(states[0]+35.0000)/10.0000 , -25.0000), 1.00000/(1.00000+exp(-(states[0]+35.0000)/10.0000)) , True, 1.00000])
    algebraic[9] = custom_piecewise([less((states[0]+35.0000)/20.0000 , 25.0000) & greater((states[0]+35.0000)/20.0000 , -25.0000), constants[17]/(3.30000*exp((states[0]+35.0000)/20.0000)+exp(-(states[0]+35.0000)/20.0000)) , True, 1.00000])
    rates[4] = (algebraic[3]-states[4])/algebraic[9]
    algebraic[0] = custom_piecewise([less(fabs(((13.0000+constants[0])-states[0])/4.00000) , 1.00000e-06), 0.320000*4.00000*(1.00000-((13.0000+constants[0])-states[0])/(2.00000*4.00000)) , True, (0.320000*((13.0000+constants[0])-states[0]))/(exp(((13.0000+constants[0])-states[0])/4.00000)-1.00000)])
    algebraic[6] = custom_piecewise([less(fabs(-((states[0]-constants[0])-40.0000)/5.00000) , 1.00000e-06), -0.280000*5.00000*(1.00000--((states[0]-constants[0])-40.0000)/(2.00000*5.00000)) , True, (-0.280000*((states[0]-constants[0])-40.0000))/(exp(-((states[0]-constants[0])-40.0000)/5.00000)-1.00000)])
    algebraic[12] = 1.00000/(algebraic[0]+algebraic[6])
    algebraic[17] = algebraic[0]/(algebraic[0]+algebraic[6])
    rates[1] = (algebraic[17]-states[1])/algebraic[12]
    algebraic[1] = 0.128000*exp(((17.0000+constants[0]+constants[1])-states[0])/18.0000)
    algebraic[7] = 4.00000/(1.00000+exp(((40.0000+constants[1]+constants[0])-states[0])/5.00000))
    algebraic[18] = algebraic[1]/(algebraic[1]+algebraic[7])
    algebraic[13] = 1.00000/(algebraic[1]+algebraic[7])
    rates[2] = (algebraic[18]-states[2])/algebraic[13]
    algebraic[2] = custom_piecewise([less(fabs(((states[0]-constants[0])-15.0000)/5.00000) , 1.00000e-06), -0.0320000*5.00000*(1.00000-((states[0]-constants[0])-15.0000)/(2.00000*5.00000)) , True, (-0.0320000*((states[0]-constants[0])-15.0000))/(exp(((states[0]-constants[0])-15.0000)/5.00000)-1.00000)])
    algebraic[8] = custom_piecewise([less(fabs(-((states[0]-constants[0])-10.0000)/40.0000) , 1.00000e-06), 0.500000*40.0000*(1.00000+((states[0]-constants[0])-10.0000)/(2.00000*40.0000)) , True, (0.500000*-((states[0]-constants[0])-10.0000))/(exp(-((states[0]-constants[0])-10.0000)/40.0000)-1.00000)])
    algebraic[14] = 1.00000/(algebraic[2]+algebraic[8])
    algebraic[19] = algebraic[2]/(algebraic[2]+algebraic[8])
    rates[3] = (algebraic[19]-states[3])/algebraic[14]
    algebraic[4] = 6.32000/(1.00000+exp(-(states[0]-5.00000)/13.8900))
    algebraic[10] = custom_piecewise([less(fabs((1.31000-states[0])/5.36000) , 1.00000e-06), 0.0200000*(5.36000+(1.31000-states[0])/2.00000) , True, (0.0200000*(1.31000-states[0]))/(1.00000-exp((states[0]-1.31000)/5.36000))])
    algebraic[15] = 1.00000/(algebraic[4]+algebraic[10])
    algebraic[20] = 1.00000/(1.00000+exp((states[0]+10.0000)/-10.0000))
    rates[5] = (algebraic[20]-states[5])/algebraic[15]
    algebraic[24] = ((((power(constants[19], 2.00000))*(power(constants[3], 2.00000))*0.00100000*states[0])/(constants[4]*constants[5]))*1.00000e-06*(states[6]-constants[20]*exp((constants[19]*constants[3]*0.00100000*states[0])/(constants[4]*constants[5]))))/(1.00000-exp((0.00100000*constants[19]*constants[3]*states[0])/(constants[4]*constants[5])))
    algebraic[25] = 1000.00*constants[18]*(power(states[5], 2.00000))*algebraic[24]
    algebraic[26] = (constants[24]*algebraic[25])/(2.00000*constants[3]*constants[23])
    rates[6] = custom_piecewise([less_equal(algebraic[26] , 0.00000), (constants[21]-states[6])/constants[22] , True, algebraic[26]+(constants[21]-states[6])/constants[22]])
    algebraic[16] = 1000.00*constants[10]*(states[0]-constants[11])
    algebraic[21] = 1000.00*constants[12]*(power(states[1], 3.00000))*states[2]*(states[0]-constants[13])
    algebraic[22] = 1000.00*constants[14]*(power(states[3], 4.00000))*(states[0]-constants[15])
    algebraic[23] = 1000.00*constants[16]*states[4]*(states[0]-constants[15])
    algebraic[27] = 1.00000/(1.00000+exp(((states[0]+75.0000)-constants[28])/5.50000))
    algebraic[28] = constants[37]+1000.00/(exp(((states[0]+71.5500)-constants[28])/14.2000)+exp(-((states[0]+89.0000)-constants[28])/11.6000))
    algebraic[29] = algebraic[27]/algebraic[28]
    algebraic[30] = (1.00000-algebraic[27])/algebraic[28]
    algebraic[31] = constants[30]*(power(states[6]/constants[36], constants[33]))
    rootfind_0(voi, constants, rates, states, algebraic)
    algebraic[34] = constants[31]*(power(algebraic[33]/constants[35], constants[34]))
    rootfind_1(voi, constants, rates, states, algebraic)
    algebraic[38] = algebraic[35]+constants[29]+algebraic[36]
    algebraic[39] = 1000.00*constants[26]*algebraic[38]*(states[0]-constants[25])
    algebraic[5] = voi-constants[9]*floor(voi/constants[9])
    algebraic[11] = custom_piecewise([greater_equal(algebraic[5] , constants[6]) & less_equal(algebraic[5] , constants[7]), constants[8] , True, 0.00000])
    rates[0] = (0.00100000*((((((algebraic[11]+-algebraic[16])-algebraic[21])-algebraic[22])-algebraic[23])-algebraic[25])-algebraic[39]))/constants[2]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = custom_piecewise([less(-(states[0]+35.0000)/10.0000 , 25.0000) & greater(-(states[0]+35.0000)/10.0000 , -25.0000), 1.00000/(1.00000+exp(-(states[0]+35.0000)/10.0000)) , True, 1.00000])
    algebraic[9] = custom_piecewise([less((states[0]+35.0000)/20.0000 , 25.0000) & greater((states[0]+35.0000)/20.0000 , -25.0000), constants[17]/(3.30000*exp((states[0]+35.0000)/20.0000)+exp(-(states[0]+35.0000)/20.0000)) , True, 1.00000])
    algebraic[0] = custom_piecewise([less(fabs(((13.0000+constants[0])-states[0])/4.00000) , 1.00000e-06), 0.320000*4.00000*(1.00000-((13.0000+constants[0])-states[0])/(2.00000*4.00000)) , True, (0.320000*((13.0000+constants[0])-states[0]))/(exp(((13.0000+constants[0])-states[0])/4.00000)-1.00000)])
    algebraic[6] = custom_piecewise([less(fabs(-((states[0]-constants[0])-40.0000)/5.00000) , 1.00000e-06), -0.280000*5.00000*(1.00000--((states[0]-constants[0])-40.0000)/(2.00000*5.00000)) , True, (-0.280000*((states[0]-constants[0])-40.0000))/(exp(-((states[0]-constants[0])-40.0000)/5.00000)-1.00000)])
    algebraic[12] = 1.00000/(algebraic[0]+algebraic[6])
    algebraic[17] = algebraic[0]/(algebraic[0]+algebraic[6])
    algebraic[1] = 0.128000*exp(((17.0000+constants[0]+constants[1])-states[0])/18.0000)
    algebraic[7] = 4.00000/(1.00000+exp(((40.0000+constants[1]+constants[0])-states[0])/5.00000))
    algebraic[18] = algebraic[1]/(algebraic[1]+algebraic[7])
    algebraic[13] = 1.00000/(algebraic[1]+algebraic[7])
    algebraic[2] = custom_piecewise([less(fabs(((states[0]-constants[0])-15.0000)/5.00000) , 1.00000e-06), -0.0320000*5.00000*(1.00000-((states[0]-constants[0])-15.0000)/(2.00000*5.00000)) , True, (-0.0320000*((states[0]-constants[0])-15.0000))/(exp(((states[0]-constants[0])-15.0000)/5.00000)-1.00000)])
    algebraic[8] = custom_piecewise([less(fabs(-((states[0]-constants[0])-10.0000)/40.0000) , 1.00000e-06), 0.500000*40.0000*(1.00000+((states[0]-constants[0])-10.0000)/(2.00000*40.0000)) , True, (0.500000*-((states[0]-constants[0])-10.0000))/(exp(-((states[0]-constants[0])-10.0000)/40.0000)-1.00000)])
    algebraic[14] = 1.00000/(algebraic[2]+algebraic[8])
    algebraic[19] = algebraic[2]/(algebraic[2]+algebraic[8])
    algebraic[4] = 6.32000/(1.00000+exp(-(states[0]-5.00000)/13.8900))
    algebraic[10] = custom_piecewise([less(fabs((1.31000-states[0])/5.36000) , 1.00000e-06), 0.0200000*(5.36000+(1.31000-states[0])/2.00000) , True, (0.0200000*(1.31000-states[0]))/(1.00000-exp((states[0]-1.31000)/5.36000))])
    algebraic[15] = 1.00000/(algebraic[4]+algebraic[10])
    algebraic[20] = 1.00000/(1.00000+exp((states[0]+10.0000)/-10.0000))
    algebraic[24] = ((((power(constants[19], 2.00000))*(power(constants[3], 2.00000))*0.00100000*states[0])/(constants[4]*constants[5]))*1.00000e-06*(states[6]-constants[20]*exp((constants[19]*constants[3]*0.00100000*states[0])/(constants[4]*constants[5]))))/(1.00000-exp((0.00100000*constants[19]*constants[3]*states[0])/(constants[4]*constants[5])))
    algebraic[25] = 1000.00*constants[18]*(power(states[5], 2.00000))*algebraic[24]
    algebraic[26] = (constants[24]*algebraic[25])/(2.00000*constants[3]*constants[23])
    algebraic[16] = 1000.00*constants[10]*(states[0]-constants[11])
    algebraic[21] = 1000.00*constants[12]*(power(states[1], 3.00000))*states[2]*(states[0]-constants[13])
    algebraic[22] = 1000.00*constants[14]*(power(states[3], 4.00000))*(states[0]-constants[15])
    algebraic[23] = 1000.00*constants[16]*states[4]*(states[0]-constants[15])
    algebraic[27] = 1.00000/(1.00000+exp(((states[0]+75.0000)-constants[28])/5.50000))
    algebraic[28] = constants[37]+1000.00/(exp(((states[0]+71.5500)-constants[28])/14.2000)+exp(-((states[0]+89.0000)-constants[28])/11.6000))
    algebraic[29] = algebraic[27]/algebraic[28]
    algebraic[30] = (1.00000-algebraic[27])/algebraic[28]
    algebraic[31] = constants[30]*(power(states[6]/constants[36], constants[33]))
    algebraic[34] = constants[31]*(power(algebraic[33]/constants[35], constants[34]))
    algebraic[38] = algebraic[35]+constants[29]+algebraic[36]
    algebraic[39] = 1000.00*constants[26]*algebraic[38]*(states[0]-constants[25])
    algebraic[5] = voi-constants[9]*floor(voi/constants[9])
    algebraic[11] = custom_piecewise([greater_equal(algebraic[5] , constants[6]) & less_equal(algebraic[5] , constants[7]), constants[8] , True, 0.00000])
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
        algebraic[32] = soln[0]
        algebraic[33] = soln[1]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess0 = soln
            algebraic[32][i] = soln[0]
            algebraic[33][i] = soln[1]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 2)
    algebraic[32] = algebraicCandidate[0]
    algebraic[33] = algebraicCandidate[1]
    resid[0] = (algebraic[32]-(algebraic[33]*constants[30])/algebraic[31])
    resid[1] = (algebraic[33]-(1.00000-algebraic[32]))
    return resid

initialGuess1 = None
def rootfind_1(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess1
    if initialGuess1 is None: initialGuess1 = ones(3)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_1, initialGuess1, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess1 = soln
        algebraic[35] = soln[0]
        algebraic[36] = soln[1]
        algebraic[37] = soln[2]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_1, initialGuess1, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess1 = soln
            algebraic[35][i] = soln[0]
            algebraic[36][i] = soln[1]
            algebraic[37][i] = soln[2]

def residualSN_1(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 3)
    algebraic[35] = algebraicCandidate[0]
    algebraic[36] = algebraicCandidate[1]
    algebraic[37] = algebraicCandidate[2]
    resid[0] = (algebraic[37]-(algebraic[30]/algebraic[29])*algebraic[35])
    resid[1] = (algebraic[35]-(constants[31]/algebraic[34])*algebraic[36])
    resid[2] = (algebraic[36]-((1.00000-algebraic[37])-algebraic[35]))
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