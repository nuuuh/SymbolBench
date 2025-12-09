# Size of variable arrays:
sizeAlgebraic = 31
sizeStates = 16
sizeConstants = 72
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "C in component membrane (uF)"
    legend_algebraic[0] = "i_Na in component sodium_current (nanoA)"
    legend_algebraic[15] = "i_NaP in component persistent_sodium_current (nanoA)"
    legend_algebraic[22] = "i_K in component delayed_rectifier_current (nanoA)"
    legend_algebraic[23] = "i_leak in component leak_current (nanoA)"
    legend_algebraic[24] = "i_T in component LVA_calcium_current (nanoA)"
    legend_algebraic[25] = "i_N in component N_HVA_calcium_current (nanoA)"
    legend_algebraic[26] = "i_P in component P_HVA_calcium_current (nanoA)"
    legend_algebraic[27] = "i_SK in component calcium_dependent_potassium_current (nanoA)"
    legend_algebraic[28] = "i_A in component fast_transient_potassium_current (nanoA)"
    legend_algebraic[29] = "i_H in component hyperpolarization_activated_current (nanoA)"
    legend_algebraic[30] = "i_app in component stimulus_protocol (nanoA)"
    legend_constants[1] = "g_Na in component sodium_current (uS)"
    legend_constants[2] = "E_Na in component sodium_current (millivolt)"
    legend_states[1] = "m in component sodium_current_m_gate (dimensionless)"
    legend_states[2] = "h in component sodium_current_h_gate (dimensionless)"
    legend_constants[3] = "theta_h in component sodium_current_h_gate (millivolt)"
    legend_constants[4] = "sigma_h in component sodium_current_h_gate (millivolt)"
    legend_constants[5] = "theta_1 in component sodium_current_h_gate (millivolt)"
    legend_constants[6] = "sigma_1 in component sodium_current_h_gate (millivolt)"
    legend_constants[7] = "sigma_2 in component sodium_current_h_gate (millivolt)"
    legend_algebraic[1] = "tau_h in component sodium_current_h_gate (millisecond)"
    legend_algebraic[16] = "h_infinity in component sodium_current_h_gate (dimensionless)"
    legend_constants[8] = "theta_m in component sodium_current_m_gate (millivolt)"
    legend_constants[9] = "sigma_m in component sodium_current_m_gate (millivolt)"
    legend_constants[10] = "tau_m in component sodium_current_m_gate (millisecond)"
    legend_algebraic[2] = "m_infinity in component sodium_current_m_gate (dimensionless)"
    legend_constants[11] = "g_NaP in component persistent_sodium_current (uS)"
    legend_states[3] = "m in component persistent_sodium_current_m_gate (dimensionless)"
    legend_states[4] = "h in component persistent_sodium_current_h_gate (dimensionless)"
    legend_constants[12] = "theta_h in component persistent_sodium_current_h_gate (millivolt)"
    legend_constants[13] = "sigma_h in component persistent_sodium_current_h_gate (millivolt)"
    legend_constants[14] = "tau_h in component persistent_sodium_current_h_gate (millisecond)"
    legend_algebraic[3] = "h_infinity in component persistent_sodium_current_h_gate (dimensionless)"
    legend_constants[15] = "theta_m in component persistent_sodium_current_m_gate (millivolt)"
    legend_constants[16] = "sigma_m in component persistent_sodium_current_m_gate (millivolt)"
    legend_constants[17] = "tau_m in component persistent_sodium_current_m_gate (millisecond)"
    legend_algebraic[4] = "m_infinity in component persistent_sodium_current_m_gate (dimensionless)"
    legend_constants[18] = "g_K in component delayed_rectifier_current (uS)"
    legend_constants[19] = "E_K in component delayed_rectifier_current (millivolt)"
    legend_states[5] = "n in component delayed_rectifier_current_n_gate (dimensionless)"
    legend_constants[20] = "theta_n in component delayed_rectifier_current_n_gate (millivolt)"
    legend_constants[21] = "sigma_n in component delayed_rectifier_current_n_gate (millivolt)"
    legend_constants[22] = "theta_1 in component delayed_rectifier_current_n_gate (millivolt)"
    legend_constants[23] = "sigma_1 in component delayed_rectifier_current_n_gate (millivolt)"
    legend_constants[24] = "sigma_2 in component delayed_rectifier_current_n_gate (millivolt)"
    legend_algebraic[5] = "tau_n in component delayed_rectifier_current_n_gate (millisecond)"
    legend_algebraic[17] = "n_infinity in component delayed_rectifier_current_n_gate (dimensionless)"
    legend_constants[25] = "g_leak in component leak_current (uS)"
    legend_constants[26] = "E_leak in component leak_current (millivolt)"
    legend_constants[27] = "g_T in component LVA_calcium_current (uS)"
    legend_constants[28] = "E_Ca in component LVA_calcium_current (millivolt)"
    legend_states[6] = "m in component LVA_calcium_current_m_gate (dimensionless)"
    legend_states[7] = "h in component LVA_calcium_current_h_gate (dimensionless)"
    legend_constants[29] = "theta_m in component LVA_calcium_current_m_gate (millivolt)"
    legend_constants[30] = "sigma_m in component LVA_calcium_current_m_gate (millivolt)"
    legend_constants[31] = "theta_1 in component LVA_calcium_current_m_gate (millivolt)"
    legend_constants[32] = "sigma_1 in component LVA_calcium_current_m_gate (millivolt)"
    legend_constants[33] = "sigma_2 in component LVA_calcium_current_m_gate (millivolt)"
    legend_algebraic[6] = "tau_m in component LVA_calcium_current_m_gate (millisecond)"
    legend_algebraic[18] = "m_infinity in component LVA_calcium_current_m_gate (dimensionless)"
    legend_constants[34] = "theta_h in component LVA_calcium_current_h_gate (millivolt)"
    legend_constants[35] = "sigma_h in component LVA_calcium_current_h_gate (millivolt)"
    legend_constants[36] = "theta_1 in component LVA_calcium_current_h_gate (millivolt)"
    legend_constants[37] = "sigma_1 in component LVA_calcium_current_h_gate (millivolt)"
    legend_algebraic[7] = "tau_h in component LVA_calcium_current_h_gate (millisecond)"
    legend_algebraic[19] = "h_infinity in component LVA_calcium_current_h_gate (dimensionless)"
    legend_constants[38] = "g_N in component N_HVA_calcium_current (uS)"
    legend_states[8] = "m in component N_HVA_calcium_current_m_gate (dimensionless)"
    legend_states[9] = "h in component N_HVA_calcium_current_h_gate (dimensionless)"
    legend_constants[39] = "theta_m in component N_HVA_calcium_current_m_gate (millivolt)"
    legend_constants[40] = "sigma_m in component N_HVA_calcium_current_m_gate (millivolt)"
    legend_constants[41] = "tau_m in component N_HVA_calcium_current_m_gate (millisecond)"
    legend_algebraic[8] = "m_infinity in component N_HVA_calcium_current_m_gate (dimensionless)"
    legend_constants[42] = "theta_h in component N_HVA_calcium_current_h_gate (millivolt)"
    legend_constants[43] = "sigma_h in component N_HVA_calcium_current_h_gate (millivolt)"
    legend_constants[44] = "tau_h in component N_HVA_calcium_current_h_gate (millisecond)"
    legend_algebraic[9] = "h_infinity in component N_HVA_calcium_current_h_gate (dimensionless)"
    legend_constants[45] = "g_P in component P_HVA_calcium_current (uS)"
    legend_states[10] = "m in component P_HVA_calcium_current_m_gate (dimensionless)"
    legend_constants[46] = "theta_m in component P_HVA_calcium_current_m_gate (millivolt)"
    legend_constants[47] = "sigma_m in component P_HVA_calcium_current_m_gate (millivolt)"
    legend_constants[48] = "tau_m in component P_HVA_calcium_current_m_gate (millisecond)"
    legend_algebraic[10] = "m_infinity in component P_HVA_calcium_current_m_gate (dimensionless)"
    legend_constants[49] = "g_SK in component calcium_dependent_potassium_current (uS)"
    legend_states[11] = "z in component calcium_dependent_potassium_current_z_gate (dimensionless)"
    legend_constants[50] = "K1 in component calcium_dependent_potassium_current_z_gate (uM_per_nanocoulomb)"
    legend_constants[51] = "K2 in component calcium_dependent_potassium_current_z_gate (per_ms)"
    legend_states[12] = "Ca_conc in component calcium_dependent_potassium_current_z_gate (uM)"
    legend_constants[52] = "tau_z in component calcium_dependent_potassium_current_z_gate (millisecond)"
    legend_algebraic[11] = "z_infinity in component calcium_dependent_potassium_current_z_gate (dimensionless)"
    legend_constants[53] = "g_A in component fast_transient_potassium_current (uS)"
    legend_states[13] = "m in component fast_transient_potassium_current_m_gate (dimensionless)"
    legend_states[14] = "h in component fast_transient_potassium_current_h_gate (dimensionless)"
    legend_constants[54] = "theta_m in component fast_transient_potassium_current_m_gate (millivolt)"
    legend_constants[55] = "sigma_m in component fast_transient_potassium_current_m_gate (millivolt)"
    legend_constants[56] = "theta_1 in component fast_transient_potassium_current_m_gate (millivolt)"
    legend_constants[57] = "theta_2 in component fast_transient_potassium_current_m_gate (millivolt)"
    legend_constants[58] = "sigma_1 in component fast_transient_potassium_current_m_gate (millivolt)"
    legend_constants[59] = "sigma_2 in component fast_transient_potassium_current_m_gate (millivolt)"
    legend_algebraic[12] = "tau_m in component fast_transient_potassium_current_m_gate (millisecond)"
    legend_algebraic[20] = "m_infinity in component fast_transient_potassium_current_m_gate (dimensionless)"
    legend_constants[60] = "theta_h in component fast_transient_potassium_current_h_gate (millivolt)"
    legend_constants[61] = "sigma_h in component fast_transient_potassium_current_h_gate (millivolt)"
    legend_constants[62] = "tau_h in component fast_transient_potassium_current_h_gate (millisecond)"
    legend_algebraic[13] = "h_infinity in component fast_transient_potassium_current_h_gate (dimensionless)"
    legend_constants[63] = "g_H in component hyperpolarization_activated_current (uS)"
    legend_constants[64] = "E_H in component hyperpolarization_activated_current (millivolt)"
    legend_states[15] = "m in component hyperpolarization_activated_current_m_gate (dimensionless)"
    legend_constants[65] = "theta_m in component hyperpolarization_activated_current_m_gate (millivolt)"
    legend_constants[66] = "sigma_m in component hyperpolarization_activated_current_m_gate (millivolt)"
    legend_constants[67] = "theta_1 in component hyperpolarization_activated_current_m_gate (millivolt)"
    legend_constants[68] = "sigma_1 in component hyperpolarization_activated_current_m_gate (millivolt)"
    legend_algebraic[14] = "tau_m in component hyperpolarization_activated_current_m_gate (millisecond)"
    legend_algebraic[21] = "m_infinity in component hyperpolarization_activated_current_m_gate (dimensionless)"
    legend_constants[69] = "i_stimStart in component stimulus_protocol (millisecond)"
    legend_constants[70] = "i_stimEnd in component stimulus_protocol (millisecond)"
    legend_constants[71] = "i_stimAmplitude in component stimulus_protocol (nanoA)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[2] = "d/dt h in component sodium_current_h_gate (dimensionless)"
    legend_rates[1] = "d/dt m in component sodium_current_m_gate (dimensionless)"
    legend_rates[4] = "d/dt h in component persistent_sodium_current_h_gate (dimensionless)"
    legend_rates[3] = "d/dt m in component persistent_sodium_current_m_gate (dimensionless)"
    legend_rates[5] = "d/dt n in component delayed_rectifier_current_n_gate (dimensionless)"
    legend_rates[6] = "d/dt m in component LVA_calcium_current_m_gate (dimensionless)"
    legend_rates[7] = "d/dt h in component LVA_calcium_current_h_gate (dimensionless)"
    legend_rates[8] = "d/dt m in component N_HVA_calcium_current_m_gate (dimensionless)"
    legend_rates[9] = "d/dt h in component N_HVA_calcium_current_h_gate (dimensionless)"
    legend_rates[10] = "d/dt m in component P_HVA_calcium_current_m_gate (dimensionless)"
    legend_rates[12] = "d/dt Ca_conc in component calcium_dependent_potassium_current_z_gate (uM)"
    legend_rates[11] = "d/dt z in component calcium_dependent_potassium_current_z_gate (dimensionless)"
    legend_rates[13] = "d/dt m in component fast_transient_potassium_current_m_gate (dimensionless)"
    legend_rates[14] = "d/dt h in component fast_transient_potassium_current_h_gate (dimensionless)"
    legend_rates[15] = "d/dt m in component hyperpolarization_activated_current_m_gate (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -71.847
    constants[0] = 0.04
    constants[1] = 0.7
    constants[2] = 60
    states[1] = 0.015
    states[2] = 0.981
    constants[3] = 44.1
    constants[4] = 7
    constants[5] = 35
    constants[6] = 4
    constants[7] = 25
    constants[8] = 36
    constants[9] = 8.5
    constants[10] = 0.1
    constants[11] = 0.05
    states[3] = 0.002
    states[4] = 0.797
    constants[12] = 65
    constants[13] = 5
    constants[14] = 150
    constants[15] = 47.1
    constants[16] = 4.1
    constants[17] = 0.1
    constants[18] = 1.3
    constants[19] = -80
    states[5] = 0.158
    constants[20] = 30
    constants[21] = 25
    constants[22] = 30
    constants[23] = 40
    constants[24] = 50
    constants[25] = 0.005
    constants[26] = -50
    constants[27] = 0.1
    constants[28] = 40
    states[6] = 0.001
    states[7] = 0.562
    constants[29] = 38
    constants[30] = 5
    constants[31] = 28
    constants[32] = 25
    constants[33] = 70
    constants[34] = 70.1
    constants[35] = 7
    constants[36] = 70
    constants[37] = 65
    constants[38] = 0.05
    states[8] = 0.001
    states[9] = 0.649
    constants[39] = 30
    constants[40] = 6
    constants[41] = 5
    constants[42] = 70
    constants[43] = 3
    constants[44] = 25
    constants[45] = 0.05
    states[10] = 0
    constants[46] = 17
    constants[47] = 3
    constants[48] = 10
    constants[49] = 0.3
    states[11] = 0
    constants[50] = -500
    constants[51] = 0.04
    states[12] = 0.0604
    constants[52] = 1
    constants[53] = 1
    states[13] = 0.057
    states[14] = 0.287
    constants[54] = 27
    constants[55] = 16
    constants[56] = 40
    constants[57] = 74
    constants[58] = 5
    constants[59] = 7.5
    constants[60] = 80
    constants[61] = 11
    constants[62] = 20
    constants[63] = 0.005
    constants[64] = -38.8
    states[15] = 0.182
    constants[65] = 79.8
    constants[66] = 5.3
    constants[67] = 70
    constants[68] = 11
    constants[69] = 10
    constants[70] = 11
    constants[71] = 10
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[2] = 1.00000/(1.00000+exp((states[0]+constants[8])/-constants[9]))
    rates[1] = (algebraic[2]-states[1])/constants[10]
    algebraic[3] = 1.00000/(1.00000+exp((states[0]+constants[12])/constants[13]))
    rates[4] = (algebraic[3]-states[4])/constants[14]
    algebraic[4] = 1.00000/(1.00000+exp(-(states[0]+constants[15])/constants[16]))
    rates[3] = (algebraic[4]-states[3])/constants[17]
    algebraic[8] = 1.00000/(1.00000+exp(-(states[0]+constants[39])/constants[40]))
    rates[8] = (algebraic[8]-states[8])/constants[41]
    algebraic[9] = 1.00000/(1.00000+exp((states[0]+constants[42])/constants[43]))
    rates[9] = (algebraic[9]-states[9])/constants[44]
    algebraic[10] = 1.00000/(1.00000+exp(-(states[0]+constants[46])/constants[47]))
    rates[10] = (algebraic[10]-states[10])/constants[48]
    algebraic[11] = 1.00000/(1.00000+power(0.00300000/states[12], 2.00000))
    rates[11] = (algebraic[11]-states[11])/constants[52]
    algebraic[13] = 1.00000/(1.00000+exp((states[0]+constants[60])/constants[61]))
    rates[14] = (algebraic[13]-states[14])/constants[62]
    algebraic[1] = 3.50000/(exp((states[0]+constants[5])/constants[6])+exp(-(states[0]+constants[5])/constants[7]))+1.00000
    algebraic[16] = 1.00000/(1.00000+exp((states[0]+constants[3])/constants[4]))
    rates[2] = (algebraic[16]-states[2])/algebraic[1]
    algebraic[5] = 2.50000/(exp((states[0]+constants[22])/constants[23])+exp(-(states[0]+constants[22])/constants[24]))+0.0100000
    algebraic[17] = 1.00000/(1.00000+exp(-(states[0]+constants[20])/constants[21]))
    rates[5] = (algebraic[17]-states[5])/algebraic[5]
    algebraic[6] = 5.00000/(exp((states[0]+constants[31])/constants[32])+exp(-(states[0]+constants[31])/constants[33]))+2.00000
    algebraic[18] = 1.00000/(1.00000+exp(-(states[0]+constants[29])/constants[30]))
    rates[6] = (algebraic[18]-states[6])/algebraic[6]
    algebraic[7] = 20.0000/(exp((states[0]+constants[36])/constants[37])+exp(-(states[0]+constants[36])/constants[37]))+1.00000
    algebraic[19] = 1.00000/(1.00000+exp((states[0]+constants[34])/constants[35]))
    rates[7] = (algebraic[19]-states[7])/algebraic[7]
    algebraic[12] = 1.00000/(exp((states[0]+constants[56])/constants[58])+exp(-(states[0]+constants[57])/constants[59]))+0.370000
    algebraic[20] = 1.00000/(1.00000+exp(-(states[0]+constants[54])/constants[55]))
    rates[13] = (algebraic[20]-states[13])/algebraic[12]
    algebraic[14] = 1.00000/(exp((states[0]+constants[67])/constants[68])+exp(-(states[0]+constants[67])/constants[68]))+50.0000
    algebraic[21] = 1.00000/(1.00000+exp((states[0]+constants[65])/constants[66]))
    rates[15] = (algebraic[21]-states[15])/algebraic[14]
    algebraic[24] = constants[27]*states[6]*states[7]*(states[0]-constants[28])
    algebraic[25] = constants[38]*states[8]*states[9]*(states[0]-constants[28])
    algebraic[26] = constants[45]*states[10]*(states[0]-constants[28])
    rates[12] = (1.00000/1000.00)*constants[50]*(algebraic[24]+algebraic[25]+algebraic[26])-constants[51]*states[12]
    algebraic[0] = constants[1]*(power(states[1], 3.00000))*states[2]*(states[0]-constants[2])
    algebraic[15] = constants[11]*states[3]*states[4]*(states[0]-constants[2])
    algebraic[22] = constants[18]*(power(states[5], 4.00000))*(states[0]-constants[19])
    algebraic[23] = constants[25]*(states[0]-constants[26])
    algebraic[27] = constants[49]*(power(states[11], 2.00000))*(states[0]-constants[19])
    algebraic[28] = constants[53]*states[13]*states[14]*(states[0]-constants[19])
    algebraic[29] = constants[63]*states[15]*(states[0]-constants[64])
    algebraic[30] = custom_piecewise([greater_equal(voi , constants[69]) & less_equal(voi , constants[70]), constants[71] , True, 0.00000])
    rates[0] = (-(algebraic[0]+algebraic[15]+algebraic[22]+algebraic[23]+algebraic[24]+algebraic[25]+algebraic[26]+algebraic[27]+algebraic[28]+algebraic[29])+algebraic[30])/constants[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[2] = 1.00000/(1.00000+exp((states[0]+constants[8])/-constants[9]))
    algebraic[3] = 1.00000/(1.00000+exp((states[0]+constants[12])/constants[13]))
    algebraic[4] = 1.00000/(1.00000+exp(-(states[0]+constants[15])/constants[16]))
    algebraic[8] = 1.00000/(1.00000+exp(-(states[0]+constants[39])/constants[40]))
    algebraic[9] = 1.00000/(1.00000+exp((states[0]+constants[42])/constants[43]))
    algebraic[10] = 1.00000/(1.00000+exp(-(states[0]+constants[46])/constants[47]))
    algebraic[11] = 1.00000/(1.00000+power(0.00300000/states[12], 2.00000))
    algebraic[13] = 1.00000/(1.00000+exp((states[0]+constants[60])/constants[61]))
    algebraic[1] = 3.50000/(exp((states[0]+constants[5])/constants[6])+exp(-(states[0]+constants[5])/constants[7]))+1.00000
    algebraic[16] = 1.00000/(1.00000+exp((states[0]+constants[3])/constants[4]))
    algebraic[5] = 2.50000/(exp((states[0]+constants[22])/constants[23])+exp(-(states[0]+constants[22])/constants[24]))+0.0100000
    algebraic[17] = 1.00000/(1.00000+exp(-(states[0]+constants[20])/constants[21]))
    algebraic[6] = 5.00000/(exp((states[0]+constants[31])/constants[32])+exp(-(states[0]+constants[31])/constants[33]))+2.00000
    algebraic[18] = 1.00000/(1.00000+exp(-(states[0]+constants[29])/constants[30]))
    algebraic[7] = 20.0000/(exp((states[0]+constants[36])/constants[37])+exp(-(states[0]+constants[36])/constants[37]))+1.00000
    algebraic[19] = 1.00000/(1.00000+exp((states[0]+constants[34])/constants[35]))
    algebraic[12] = 1.00000/(exp((states[0]+constants[56])/constants[58])+exp(-(states[0]+constants[57])/constants[59]))+0.370000
    algebraic[20] = 1.00000/(1.00000+exp(-(states[0]+constants[54])/constants[55]))
    algebraic[14] = 1.00000/(exp((states[0]+constants[67])/constants[68])+exp(-(states[0]+constants[67])/constants[68]))+50.0000
    algebraic[21] = 1.00000/(1.00000+exp((states[0]+constants[65])/constants[66]))
    algebraic[24] = constants[27]*states[6]*states[7]*(states[0]-constants[28])
    algebraic[25] = constants[38]*states[8]*states[9]*(states[0]-constants[28])
    algebraic[26] = constants[45]*states[10]*(states[0]-constants[28])
    algebraic[0] = constants[1]*(power(states[1], 3.00000))*states[2]*(states[0]-constants[2])
    algebraic[15] = constants[11]*states[3]*states[4]*(states[0]-constants[2])
    algebraic[22] = constants[18]*(power(states[5], 4.00000))*(states[0]-constants[19])
    algebraic[23] = constants[25]*(states[0]-constants[26])
    algebraic[27] = constants[49]*(power(states[11], 2.00000))*(states[0]-constants[19])
    algebraic[28] = constants[53]*states[13]*states[14]*(states[0]-constants[19])
    algebraic[29] = constants[63]*states[15]*(states[0]-constants[64])
    algebraic[30] = custom_piecewise([greater_equal(voi , constants[69]) & less_equal(voi , constants[70]), constants[71] , True, 0.00000])
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