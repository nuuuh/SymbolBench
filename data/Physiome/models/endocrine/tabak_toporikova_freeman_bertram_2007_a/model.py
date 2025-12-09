# Size of variable arrays:
sizeAlgebraic = 11
sizeStates = 4
sizeConstants = 23
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
    legend_constants[0] = "C in component membrane (picofarad)"
    legend_algebraic[5] = "I_Ca in component I_Ca (picoampere)"
    legend_algebraic[0] = "I_K in component I_K (picoampere)"
    legend_algebraic[7] = "I_SK in component I_SK (picoampere)"
    legend_algebraic[10] = "I_DA in component I_DA (picoampere)"
    legend_constants[1] = "gK in component I_K (nanosiemens)"
    legend_constants[2] = "VK in component model_parameters (millivolt)"
    legend_states[1] = "n in component n (dimensionless)"
    legend_algebraic[1] = "n_infinity in component n (dimensionless)"
    legend_constants[3] = "lambda in component n (dimensionless)"
    legend_constants[4] = "tau_n in component n (millisecond)"
    legend_constants[5] = "vn in component n (millivolt)"
    legend_constants[6] = "sn in component n (millivolt)"
    legend_constants[7] = "gCa in component I_Ca (nanosiemens)"
    legend_constants[8] = "VCa in component model_parameters (millivolt)"
    legend_algebraic[4] = "m_infinity in component m (dimensionless)"
    legend_constants[9] = "vm in component m (millivolt)"
    legend_constants[10] = "sm in component m (millivolt)"
    legend_constants[11] = "gSK in component I_SK (nanosiemens)"
    legend_algebraic[6] = "s_infinity in component I_SK (dimensionless)"
    legend_constants[12] = "ks in component I_SK (micromolar)"
    legend_states[2] = "Ca in component Ca (micromolar)"
    legend_algebraic[9] = "I_BK in component I_DA (picoampere)"
    legend_constants[13] = "gBK in component I_DA (nanosiemens)"
    legend_algebraic[8] = "f_infinity in component f (dimensionless)"
    legend_constants[14] = "vf in component f (millivolt)"
    legend_constants[15] = "sf in component f (millivolt)"
    legend_states[3] = "h in component h (dimensionless)"
    legend_algebraic[2] = "h_infinity in component h (dimensionless)"
    legend_constants[16] = "tau_h in component h (millisecond)"
    legend_constants[17] = "vh in component h (millivolt)"
    legend_constants[18] = "sh in component h (millivolt)"
    legend_constants[19] = "fc in component Ca (dimensionless)"
    legend_constants[20] = "alpha in component Ca (micromolar_femtocoulomb)"
    legend_constants[21] = "kc in component Ca (first_order_rate_constant)"
    legend_algebraic[3] = "PRL in component PRL (dimensionless)"
    legend_constants[22] = "kPRL in component PRL (micromolar_4)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt n in component n (dimensionless)"
    legend_rates[3] = "d/dt h in component h (dimensionless)"
    legend_rates[2] = "d/dt Ca in component Ca (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -60
    constants[0] = 10
    constants[1] = 4
    constants[2] = -75
    states[1] = 0.1
    constants[3] = 0.7
    constants[4] = 30
    constants[5] = -5
    constants[6] = 10
    constants[7] = 2
    constants[8] = 50
    constants[9] = -20
    constants[10] = 12
    constants[11] = 1.7
    constants[12] = 0.5
    states[2] = 0.1
    constants[13] = 0.2
    constants[14] = -20
    constants[15] = 5.6
    states[3] = 0.1
    constants[16] = 20
    constants[17] = -60
    constants[18] = 5
    constants[19] = 0.01
    constants[20] = 0.0015
    constants[21] = 0.16
    constants[22] = 1
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = 1.00000/(1.00000+exp((constants[5]-states[0])/constants[6]))
    rates[1] = (constants[3]*(algebraic[1]-states[1]))/constants[4]
    algebraic[2] = 1.00000/(1.00000+exp((states[0]-constants[17])/constants[18]))
    rates[3] = (algebraic[2]-states[3])/constants[16]
    algebraic[4] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[5] = constants[7]*algebraic[4]*(states[0]-constants[8])
    rates[2] = -constants[19]*(constants[20]*algebraic[5]+constants[21]*states[2])
    algebraic[0] = constants[1]*states[1]*(states[0]-constants[2])
    algebraic[6] = (power(states[2], 2.00000))/(power(states[2], 2.00000)+power(constants[12], 2.00000))
    algebraic[7] = constants[11]*algebraic[6]*(states[0]-constants[2])
    algebraic[8] = 1.00000/(1.00000+exp((constants[14]-states[0])/constants[15]))
    algebraic[9] = constants[13]*algebraic[8]*(states[0]-constants[2])
    algebraic[10] = algebraic[9]
    rates[0] = -(algebraic[5]+algebraic[0]+algebraic[7]+algebraic[10])/constants[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = 1.00000/(1.00000+exp((constants[5]-states[0])/constants[6]))
    algebraic[2] = 1.00000/(1.00000+exp((states[0]-constants[17])/constants[18]))
    algebraic[4] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[5] = constants[7]*algebraic[4]*(states[0]-constants[8])
    algebraic[0] = constants[1]*states[1]*(states[0]-constants[2])
    algebraic[6] = (power(states[2], 2.00000))/(power(states[2], 2.00000)+power(constants[12], 2.00000))
    algebraic[7] = constants[11]*algebraic[6]*(states[0]-constants[2])
    algebraic[8] = 1.00000/(1.00000+exp((constants[14]-states[0])/constants[15]))
    algebraic[9] = constants[13]*algebraic[8]*(states[0]-constants[2])
    algebraic[10] = algebraic[9]
    algebraic[3] = constants[22]*(power(states[2], 4.00000))
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