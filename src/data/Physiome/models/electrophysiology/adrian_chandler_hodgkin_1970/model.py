# Size of variable arrays:
sizeAlgebraic = 12
sizeStates = 5
sizeConstants = 31
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "t in component interface (ms)"
    legend_constants[0] = "Cm in component interface (uFpmmsq)"
    legend_constants[1] = "Am in component interface (pmm)"
    legend_algebraic[0] = "Istim in component interface (uApmmcu)"
    legend_states[0] = "Vm in component membrane (mV)"
    legend_states[1] = "Vt in component Ttubular_current_Vt_var (mV)"
    legend_states[2] = "m in component sodium_current_m_gate (dimensionless)"
    legend_states[3] = "h in component sodium_current_h_gate (dimensionless)"
    legend_states[4] = "n in component potassium_current_n_gate (dimensionless)"
    legend_algebraic[4] = "INa in component sodium_current (uApmmsq)"
    legend_algebraic[9] = "IK in component potassium_current (uApmmsq)"
    legend_algebraic[10] = "IL in component leak_current (uApmmsq)"
    legend_algebraic[11] = "IT in component Ttubular_current (uApmmsq)"
    legend_algebraic[5] = "IStimC in component interface (uApmmcu)"
    legend_constants[30] = "AmC in component interface (pmm)"
    legend_constants[2] = "IstimStart in component interface (ms)"
    legend_constants[3] = "IstimEnd in component interface (ms)"
    legend_constants[4] = "IstimAmplitude in component interface (uApmmcu)"
    legend_constants[5] = "IstimPeriod in component interface (ms)"
    legend_constants[6] = "IstimPulseDuration in component interface (ms)"
    legend_constants[7] = "gNa_max in component sodium_current (mSpmmsq)"
    legend_constants[8] = "ENa in component sodium_current (mV)"
    legend_algebraic[1] = "alpha_m in component sodium_current_m_gate (pms)"
    legend_algebraic[6] = "beta_m in component sodium_current_m_gate (pms)"
    legend_constants[9] = "alpha_m_max in component sodium_current_m_gate (pms)"
    legend_constants[10] = "beta_m_max in component sodium_current_m_gate (pms)"
    legend_constants[11] = "Em in component sodium_current_m_gate (mV)"
    legend_constants[12] = "v_alpha_m in component sodium_current_m_gate (dimensionless)"
    legend_constants[13] = "v_beta_m in component sodium_current_m_gate (mV)"
    legend_algebraic[2] = "alpha_h in component sodium_current_h_gate (pms)"
    legend_algebraic[7] = "beta_h in component sodium_current_h_gate (pms)"
    legend_constants[14] = "alpha_h_max in component sodium_current_h_gate (pms)"
    legend_constants[15] = "beta_h_max in component sodium_current_h_gate (pms)"
    legend_constants[16] = "Eh in component sodium_current_h_gate (mV)"
    legend_constants[17] = "v_alpha_h in component sodium_current_h_gate (mV)"
    legend_constants[18] = "v_beta_h in component sodium_current_h_gate (mV)"
    legend_constants[19] = "gK_max in component potassium_current (mSpmmsq)"
    legend_constants[20] = "EK in component potassium_current (mV)"
    legend_algebraic[3] = "alpha_n in component potassium_current_n_gate (pms)"
    legend_algebraic[8] = "beta_n in component potassium_current_n_gate (pms)"
    legend_constants[21] = "alpha_n_max in component potassium_current_n_gate (pms)"
    legend_constants[22] = "beta_n_max in component potassium_current_n_gate (pms)"
    legend_constants[23] = "En in component potassium_current_n_gate (mV)"
    legend_constants[24] = "v_alpha_n in component potassium_current_n_gate (dimensionless)"
    legend_constants[25] = "v_beta_n in component potassium_current_n_gate (mV)"
    legend_constants[26] = "EL in component leak_current (mV)"
    legend_constants[27] = "gL_max in component leak_current (mSpmmsq)"
    legend_constants[28] = "Rs in component Ttubular_current (mmsqpmS)"
    legend_constants[29] = "Ct in component Ttubular_current_Vt_var (uFpmmsq)"
    legend_rates[0] = "d/dt Vm in component membrane (mV)"
    legend_rates[2] = "d/dt m in component sodium_current_m_gate (dimensionless)"
    legend_rates[3] = "d/dt h in component sodium_current_h_gate (dimensionless)"
    legend_rates[4] = "d/dt n in component potassium_current_n_gate (dimensionless)"
    legend_rates[1] = "d/dt Vt in component Ttubular_current_Vt_var (mV)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 0.009
    constants[1] = 200.0
    states[0] = -95.0
    states[1] = -95.0
    states[2] = 0.0
    states[3] = 1.0
    states[4] = 0.0
    constants[2] = 10
    constants[3] = 50000
    constants[4] = 0.5
    constants[5] = 1000
    constants[6] = 1
    constants[7] = 1.8
    constants[8] = 50.0
    constants[9] = 0.208
    constants[10] = 2.081
    constants[11] = -42.0
    constants[12] = 10.0
    constants[13] = 18.0
    constants[14] = 0.0156
    constants[15] = 3.382
    constants[16] = -41.0
    constants[17] = 14.7
    constants[18] = 7.6
    constants[19] = 0.415
    constants[20] = -70.0
    constants[21] = 0.0229
    constants[22] = 0.09616
    constants[23] = -40.0
    constants[24] = 7.0
    constants[25] = 40.0
    constants[26] = -95.0
    constants[27] = 0.0024
    constants[28] = 15.0
    constants[29] = 0.04
    constants[30] = constants[1]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = (states[0]-states[1])/(constants[28]*constants[29])
    algebraic[1] = (constants[9]*(states[0]-constants[11]))/(1.00000-exp((constants[11]-states[0])/constants[12]))
    algebraic[6] = constants[10]*exp((constants[11]-states[0])/constants[13])
    rates[2] = algebraic[1]*(1.00000-states[2])-algebraic[6]*states[2]
    algebraic[2] = constants[14]*exp((constants[16]-states[0])/constants[17])
    algebraic[7] = constants[15]/(1.00000+exp((constants[16]-states[0])/constants[18]))
    rates[3] = algebraic[2]*(1.00000-states[3])-algebraic[7]*states[3]
    algebraic[3] = (constants[21]*(states[0]-constants[23]))/(1.00000-exp((constants[23]-states[0])/constants[24]))
    algebraic[8] = constants[22]*exp((constants[23]-states[0])/constants[25])
    rates[4] = algebraic[3]*(1.00000-states[4])-algebraic[8]*states[4]
    algebraic[0] = custom_piecewise([greater_equal(voi , constants[2]) & less_equal(voi , constants[3]) & less_equal((voi-constants[2])-floor((voi-constants[2])/constants[5])*constants[5] , constants[6]), constants[4] , True, 0.00000])
    algebraic[4] = constants[7]*states[2]*states[2]*states[2]*states[3]*(states[0]-constants[8])
    algebraic[9] = constants[19]*states[4]*states[4]*states[4]*states[4]*(states[0]-constants[20])
    algebraic[10] = constants[27]*(states[0]-constants[26])
    algebraic[11] = (states[0]-states[1])/constants[28]
    rates[0] = (algebraic[0]-(algebraic[4]+algebraic[9]+algebraic[10]+algebraic[11]))/constants[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = (constants[9]*(states[0]-constants[11]))/(1.00000-exp((constants[11]-states[0])/constants[12]))
    algebraic[6] = constants[10]*exp((constants[11]-states[0])/constants[13])
    algebraic[2] = constants[14]*exp((constants[16]-states[0])/constants[17])
    algebraic[7] = constants[15]/(1.00000+exp((constants[16]-states[0])/constants[18]))
    algebraic[3] = (constants[21]*(states[0]-constants[23]))/(1.00000-exp((constants[23]-states[0])/constants[24]))
    algebraic[8] = constants[22]*exp((constants[23]-states[0])/constants[25])
    algebraic[0] = custom_piecewise([greater_equal(voi , constants[2]) & less_equal(voi , constants[3]) & less_equal((voi-constants[2])-floor((voi-constants[2])/constants[5])*constants[5] , constants[6]), constants[4] , True, 0.00000])
    algebraic[4] = constants[7]*states[2]*states[2]*states[2]*states[3]*(states[0]-constants[8])
    algebraic[9] = constants[19]*states[4]*states[4]*states[4]*states[4]*(states[0]-constants[20])
    algebraic[10] = constants[27]*(states[0]-constants[26])
    algebraic[11] = (states[0]-states[1])/constants[28]
    algebraic[5] = algebraic[0]
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