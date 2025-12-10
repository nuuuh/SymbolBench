# Size of variable arrays:
sizeAlgebraic = 9
sizeStates = 5
sizeConstants = 28
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
    legend_algebraic[5] = "i_Ca in component calcium_current (picoA)"
    legend_algebraic[0] = "i_K in component rapidly_activating_K_current (picoA)"
    legend_algebraic[6] = "i_slow in component slow_K_current (picoA)"
    legend_algebraic[8] = "i_Na_Ca in component Na_Ca_exchanger_current (picoA)"
    legend_constants[1] = "V_K in component rapidly_activating_K_current (millivolt)"
    legend_constants[2] = "g_K in component rapidly_activating_K_current (picoS)"
    legend_states[1] = "n in component rapidly_activating_K_current_n_gate (dimensionless)"
    legend_algebraic[1] = "n_infinity in component rapidly_activating_K_current_n_gate (dimensionless)"
    legend_constants[3] = "lamda in component rapidly_activating_K_current_n_gate (dimensionless)"
    legend_constants[4] = "tau_n in component rapidly_activating_K_current_n_gate (millisecond)"
    legend_constants[5] = "V_n in component rapidly_activating_K_current_n_gate (millivolt)"
    legend_constants[6] = "S_n in component rapidly_activating_K_current_n_gate (millivolt)"
    legend_constants[7] = "V_Ca in component calcium_current (millivolt)"
    legend_constants[8] = "g_Ca in component calcium_current (picoS)"
    legend_algebraic[3] = "m_infinity in component calcium_current_m_gate (dimensionless)"
    legend_constants[9] = "V_m in component calcium_current_m_gate (millivolt)"
    legend_constants[10] = "S_m in component calcium_current_m_gate (millivolt)"
    legend_constants[11] = "g_s in component slow_K_current (picoS)"
    legend_states[2] = "s in component slow_K_current_s_gate (dimensionless)"
    legend_algebraic[2] = "s_infinity in component slow_K_current_s_gate (dimensionless)"
    legend_constants[12] = "tau_s in component slow_K_current_s_gate (millisecond)"
    legend_constants[13] = "V_s in component slow_K_current_s_gate (millivolt)"
    legend_constants[14] = "S_s in component slow_K_current_s_gate (millivolt)"
    legend_constants[15] = "R_s in component slow_K_current_s_gate (dimensionless)"
    legend_algebraic[4] = "S_V_R_s in component slow_K_current_s_gate (dimensionless)"
    legend_constants[16] = "g_Na_Ca in component Na_Ca_exchanger_current (picoS)"
    legend_constants[17] = "K_1_2 in component Na_Ca_exchanger_current (micromolar)"
    legend_algebraic[7] = "V_Na_Ca in component Na_Ca_exchanger_current (millivolt)"
    legend_constants[18] = "RT_F in component Na_Ca_exchanger_current (millivolt)"
    legend_constants[19] = "nH in component Na_Ca_exchanger_current (dimensionless)"
    legend_states[3] = "Ca_i in component ionic_concentrations (micromolar)"
    legend_constants[20] = "Ca_o in component ionic_concentrations (micromolar)"
    legend_constants[21] = "Na_i in component ionic_concentrations (millimolar)"
    legend_constants[22] = "Na_o in component ionic_concentrations (millimolar)"
    legend_states[4] = "Ca_ret in component ionic_concentrations (micromolar)"
    legend_constants[23] = "f in component ionic_concentrations (dimensionless)"
    legend_constants[24] = "k_Ca in component ionic_concentrations (per_millisecond)"
    legend_constants[25] = "k_rel in component ionic_concentrations (per_millisecond)"
    legend_constants[26] = "k_pump in component ionic_concentrations (per_millisecond)"
    legend_constants[27] = "alpha in component ionic_concentrations (mole_per_microlitre_coulomb)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt n in component rapidly_activating_K_current_n_gate (dimensionless)"
    legend_rates[2] = "d/dt s in component slow_K_current_s_gate (dimensionless)"
    legend_rates[3] = "d/dt Ca_i in component ionic_concentrations (micromolar)"
    legend_rates[4] = "d/dt Ca_ret in component ionic_concentrations (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -76.0
    constants[0] = 5310.0
    constants[1] = -75.0
    constants[2] = 2700.0
    states[1] = 0.1
    constants[3] = 1.0
    constants[4] = 20.0
    constants[5] = -16.0
    constants[6] = 5.6
    constants[7] = 25.0
    constants[8] = 1000.0
    constants[9] = -20.0
    constants[10] = 12.0
    constants[11] = 200.0
    states[2] = 0.1
    constants[12] = 12000.0
    constants[13] = -52.0
    constants[14] = 10.0
    constants[15] = 0.58
    constants[16] = 350.0
    constants[17] = 1.5
    constants[18] = 26.54
    constants[19] = 5.0
    states[3] = 0.52
    constants[20] = 2600.0
    constants[21] = 10.0
    constants[22] = 140.0
    states[4] = 0.7
    constants[23] = 0.02
    constants[24] = 0.64
    constants[25] = 0.0006
    constants[26] = 0.2
    constants[27] = 0.00006
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[4] = -constants[25]*(states[4]-states[3])+constants[26]*states[3]
    algebraic[1] = 1.00000/(1.00000+exp((constants[5]-states[0])/constants[6]))
    rates[1] = constants[3]*((algebraic[1]-states[1])/constants[4])
    algebraic[2] = 1.00000/(1.00000+exp((constants[13]-states[0])/constants[14]))
    algebraic[4] = algebraic[2]+constants[15]
    rates[2] = (algebraic[4]-states[2])/constants[12]
    algebraic[3] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[5] = constants[8]*algebraic[3]*(states[0]-constants[7])
    algebraic[0] = constants[2]*states[1]*(states[0]-constants[1])
    algebraic[6] = constants[11]*states[2]*(states[0]-constants[1])
    algebraic[7] = constants[18]*(3.00000*log(constants[22]/constants[21]-log(constants[20]/states[3])))
    algebraic[8] = constants[16]*((power(states[3], constants[19]))/(power(constants[17], constants[19])+power(states[3], constants[19])))*(states[0]-algebraic[7])
    rates[0] = -(algebraic[0]+algebraic[5]+algebraic[6]+algebraic[8])/constants[0]
    rates[3] = (constants[23]*(-constants[27]*(algebraic[5]-2.00000*algebraic[8])-constants[24]*states[3])+constants[25]*(states[4]-states[3]))-constants[26]*states[3]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = 1.00000/(1.00000+exp((constants[5]-states[0])/constants[6]))
    algebraic[2] = 1.00000/(1.00000+exp((constants[13]-states[0])/constants[14]))
    algebraic[4] = algebraic[2]+constants[15]
    algebraic[3] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[5] = constants[8]*algebraic[3]*(states[0]-constants[7])
    algebraic[0] = constants[2]*states[1]*(states[0]-constants[1])
    algebraic[6] = constants[11]*states[2]*(states[0]-constants[1])
    algebraic[7] = constants[18]*(3.00000*log(constants[22]/constants[21]-log(constants[20]/states[3])))
    algebraic[8] = constants[16]*((power(states[3], constants[19]))/(power(constants[17], constants[19])+power(states[3], constants[19])))*(states[0]-algebraic[7])
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