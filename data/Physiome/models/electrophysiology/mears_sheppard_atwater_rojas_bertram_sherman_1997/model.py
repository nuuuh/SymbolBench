# Size of variable arrays:
sizeAlgebraic = 24
sizeStates = 5
sizeConstants = 37
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
    legend_constants[0] = "Cm in component membrane (femtoF)"
    legend_algebraic[1] = "i_K in component K_current (picoA)"
    legend_algebraic[11] = "i_K_Ca in component K_Ca_current (picoA)"
    legend_algebraic[4] = "i_K_ATP in component K_ATP_current (picoA)"
    legend_algebraic[15] = "i_CRAC in component CRAC_current (picoA)"
    legend_algebraic[10] = "i_Ca in component Ca_current_total (picoA)"
    legend_algebraic[17] = "i_leak in component leak_current (picoA)"
    legend_constants[1] = "V_K in component K_current (millivolt)"
    legend_constants[2] = "g_K in component K_current (picoS)"
    legend_states[1] = "n in component K_channel_n_gate (dimensionless)"
    legend_algebraic[0] = "n_infinity in component K_channel_n_gate (dimensionless)"
    legend_algebraic[3] = "tau_n in component K_channel_n_gate (millisecond)"
    legend_constants[3] = "Vn in component K_channel_n_gate (millivolt)"
    legend_constants[4] = "Sn in component K_channel_n_gate (millivolt)"
    legend_constants[5] = "lambda_n in component K_channel_n_gate (dimensionless)"
    legend_constants[6] = "g_K_ATP in component K_ATP_current (picoS)"
    legend_algebraic[7] = "i_Ca_f in component fast_Ca_current (picoA)"
    legend_constants[7] = "V_Ca in component fast_Ca_current (millivolt)"
    legend_constants[8] = "g_Ca_f in component fast_Ca_current (picoS)"
    legend_algebraic[6] = "m_f_infinity in component fast_Ca_channel_m_gate (dimensionless)"
    legend_constants[9] = "Vm_f in component fast_Ca_channel_m_gate (millivolt)"
    legend_constants[10] = "Sm_f in component fast_Ca_channel_m_gate (millivolt)"
    legend_algebraic[9] = "i_Ca_s in component slow_Ca_current (picoA)"
    legend_constants[11] = "g_Ca_s in component slow_Ca_current (picoS)"
    legend_algebraic[8] = "m_s_infinity in component slow_Ca_channel_m_gate (dimensionless)"
    legend_states[2] = "jm in component slow_Ca_channel_j_gate (dimensionless)"
    legend_constants[12] = "Vm_s in component slow_Ca_channel_m_gate (millivolt)"
    legend_constants[13] = "Sm_s in component slow_Ca_channel_m_gate (millivolt)"
    legend_algebraic[2] = "jm_infinity in component slow_Ca_channel_j_gate (dimensionless)"
    legend_constants[14] = "Vj in component slow_Ca_channel_j_gate (millivolt)"
    legend_algebraic[5] = "tau_j in component slow_Ca_channel_j_gate (millisecond)"
    legend_constants[15] = "Sj in component slow_Ca_channel_j_gate (millivolt)"
    legend_constants[16] = "g_K_Ca in component K_Ca_current (picoS)"
    legend_states[3] = "Ca_i in component Ca_equations (micromolar)"
    legend_constants[17] = "kdkca in component K_Ca_current (micromolar)"
    legend_constants[18] = "g_CRAC in component CRAC_current (picoS)"
    legend_constants[19] = "V_CRAC in component CRAC_current (millivolt)"
    legend_states[4] = "Ca_er in component Ca_equations (micromolar)"
    legend_algebraic[13] = "r_infinity in component CRAC_r_gate (dimensionless)"
    legend_constants[20] = "Ca_er_bar in component CRAC_r_gate (micromolar)"
    legend_constants[21] = "sloper in component CRAC_r_gate (micromolar)"
    legend_constants[22] = "g_leak in component leak_current (picoS)"
    legend_algebraic[12] = "J_er_p in component ER_parameters (micromolar_per_millisecond)"
    legend_constants[23] = "IP3 in component ER_parameters (micromolar)"
    legend_constants[24] = "kerp in component ER_parameters (micromolar)"
    legend_constants[25] = "verp in component ER_parameters (micromolar_per_millisecond)"
    legend_constants[26] = "dact in component ER_parameters (micromolar)"
    legend_constants[27] = "dinh in component ER_parameters (micromolar)"
    legend_constants[28] = "dip3 in component ER_parameters (micromolar)"
    legend_algebraic[14] = "a_infinity in component ER_parameters (dimensionless)"
    legend_constants[36] = "b_infinity in component ER_parameters (dimensionless)"
    legend_algebraic[16] = "h_infinity in component ER_parameters (dimensionless)"
    legend_algebraic[18] = "O in component ER_parameters (per_millisecond)"
    legend_algebraic[21] = "J_er_tot in component Ca_equations (micromolar_per_millisecond)"
    legend_algebraic[20] = "J_er_IP3 in component Ca_equations (micromolar_per_millisecond)"
    legend_algebraic[19] = "J_er_leak in component Ca_equations (micromolar_per_millisecond)"
    legend_algebraic[23] = "J_mem_tot in component Ca_membrane_flux (micromolar_per_millisecond)"
    legend_constants[29] = "perl in component Ca_equations (per_millisecond)"
    legend_constants[30] = "lambda_er in component Ca_equations (dimensionless)"
    legend_constants[31] = "sigma_er in component Ca_equations (dimensionless)"
    legend_constants[32] = "kmp in component Ca_membrane_flux (micromolar)"
    legend_constants[33] = "vmp in component Ca_membrane_flux (micromolar)"
    legend_constants[34] = "gamma in component Ca_membrane_flux (micromolar_per_picoA)"
    legend_algebraic[22] = "Jmp in component Ca_membrane_flux (micromolar)"
    legend_constants[35] = "f in component Ca_membrane_flux (per_millisecond)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt n in component K_channel_n_gate (dimensionless)"
    legend_rates[2] = "d/dt jm in component slow_Ca_channel_j_gate (dimensionless)"
    legend_rates[4] = "d/dt Ca_er in component Ca_equations (micromolar)"
    legend_rates[3] = "d/dt Ca_i in component Ca_equations (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -61
    constants[0] = 6158
    constants[1] = -70
    constants[2] = 3900
    states[1] = 0.0005
    constants[3] = -15
    constants[4] = 6
    constants[5] = 1.85
    constants[6] = 150
    constants[7] = 100
    constants[8] = 810
    constants[9] = -20
    constants[10] = 7.5
    constants[11] = 510
    states[2] = 0.12
    constants[12] = -16
    constants[13] = 10
    constants[14] = -53
    constants[15] = 2
    constants[16] = 1200
    states[3] = 0.11
    constants[17] = 0.55
    constants[18] = 75
    constants[19] = 0
    states[4] = 60
    constants[20] = 40
    constants[21] = 3
    constants[22] = 0
    constants[23] = 0
    constants[24] = 0.09
    constants[25] = 0.24
    constants[26] = 0.35
    constants[27] = 0.4
    constants[28] = 0.2
    constants[29] = 0.003
    constants[30] = 250
    constants[31] = 1
    constants[32] = 0.35
    constants[33] = 0.08
    constants[34] = 0.000003607
    constants[35] = 0.01
    constants[36] = constants[23]/(constants[23]+constants[28])
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = 1.00000/(1.00000+exp((constants[3]-states[0])/constants[4]))
    algebraic[3] = 9.09000/(1.00000+exp((states[0]-constants[3])/constants[4]))
    rates[1] = (constants[5]*(algebraic[0]-states[1]))/algebraic[3]
    algebraic[2] = 1.00000-1.00000/(1.00000+exp((states[0]-constants[14])/constants[15]))
    algebraic[5] = 50000.0/(exp((states[0]-constants[14])/4.00000)+exp((constants[14]-states[0])/4.00000))+1500.00
    rates[2] = (algebraic[2]-states[2])/algebraic[5]
    algebraic[1] = constants[2]*states[1]*(states[0]-constants[1])
    algebraic[11] = ((constants[16]*(power(states[3], 5.00000)))/(power(states[3], 5.00000)+power(constants[17], 5.00000)))*(states[0]-constants[1])
    algebraic[4] = constants[6]*(states[0]-constants[1])
    algebraic[13] = 1.00000/(1.00000+exp((states[4]-constants[20])/constants[21]))
    algebraic[15] = constants[18]*algebraic[13]*(states[0]-constants[19])
    algebraic[6] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[7] = constants[8]*algebraic[6]*(states[0]-constants[7])
    algebraic[8] = 1.00000/(1.00000+exp((constants[12]-states[0])/constants[13]))
    algebraic[9] = constants[11]*algebraic[8]*(1.00000-states[2])*(states[0]-constants[7])
    algebraic[10] = algebraic[7]+algebraic[9]
    algebraic[17] = constants[22]*(states[0]-constants[19])
    rates[0] = -(algebraic[10]+algebraic[1]+algebraic[4]+algebraic[11]+algebraic[15]+algebraic[17])/constants[0]
    algebraic[12] = (constants[25]*(power(states[3], 2.00000)))/(power(states[3], 2.00000)+power(constants[24], 2.00000))
    algebraic[14] = 1.00000/(1.00000+constants[26]/states[3])
    algebraic[16] = 1.00000/(1.00000+states[3]/constants[27])
    algebraic[18] = (power(algebraic[14], 3.00000))*(power(constants[36], 3.00000))*(power(algebraic[16], 3.00000))*1.00000
    algebraic[20] = algebraic[18]*(states[4]-states[3])
    algebraic[19] = constants[29]*(states[4]-states[3])
    algebraic[21] = (algebraic[19]+algebraic[20])-algebraic[12]
    rates[4] = -algebraic[21]/(constants[30]*constants[31])
    algebraic[22] = (constants[33]*(power(states[3], 2.00000)))/(power(states[3], 2.00000)+power(constants[32], 2.00000))
    algebraic[23] = -constants[35]*(constants[34]*algebraic[10]+algebraic[22])
    rates[3] = algebraic[21]/constants[30]+algebraic[23]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = 1.00000/(1.00000+exp((constants[3]-states[0])/constants[4]))
    algebraic[3] = 9.09000/(1.00000+exp((states[0]-constants[3])/constants[4]))
    algebraic[2] = 1.00000-1.00000/(1.00000+exp((states[0]-constants[14])/constants[15]))
    algebraic[5] = 50000.0/(exp((states[0]-constants[14])/4.00000)+exp((constants[14]-states[0])/4.00000))+1500.00
    algebraic[1] = constants[2]*states[1]*(states[0]-constants[1])
    algebraic[11] = ((constants[16]*(power(states[3], 5.00000)))/(power(states[3], 5.00000)+power(constants[17], 5.00000)))*(states[0]-constants[1])
    algebraic[4] = constants[6]*(states[0]-constants[1])
    algebraic[13] = 1.00000/(1.00000+exp((states[4]-constants[20])/constants[21]))
    algebraic[15] = constants[18]*algebraic[13]*(states[0]-constants[19])
    algebraic[6] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[7] = constants[8]*algebraic[6]*(states[0]-constants[7])
    algebraic[8] = 1.00000/(1.00000+exp((constants[12]-states[0])/constants[13]))
    algebraic[9] = constants[11]*algebraic[8]*(1.00000-states[2])*(states[0]-constants[7])
    algebraic[10] = algebraic[7]+algebraic[9]
    algebraic[17] = constants[22]*(states[0]-constants[19])
    algebraic[12] = (constants[25]*(power(states[3], 2.00000)))/(power(states[3], 2.00000)+power(constants[24], 2.00000))
    algebraic[14] = 1.00000/(1.00000+constants[26]/states[3])
    algebraic[16] = 1.00000/(1.00000+states[3]/constants[27])
    algebraic[18] = (power(algebraic[14], 3.00000))*(power(constants[36], 3.00000))*(power(algebraic[16], 3.00000))*1.00000
    algebraic[20] = algebraic[18]*(states[4]-states[3])
    algebraic[19] = constants[29]*(states[4]-states[3])
    algebraic[21] = (algebraic[19]+algebraic[20])-algebraic[12]
    algebraic[22] = (constants[33]*(power(states[3], 2.00000)))/(power(states[3], 2.00000)+power(constants[32], 2.00000))
    algebraic[23] = -constants[35]*(constants[34]*algebraic[10]+algebraic[22])
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