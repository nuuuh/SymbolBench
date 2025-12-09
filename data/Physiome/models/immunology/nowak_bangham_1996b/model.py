# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 4
sizeConstants = 9
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (day)"
    legend_constants[0] = "lambda in component uninfected (first_order_rate_constant)"
    legend_constants[1] = "d in component uninfected (first_order_rate_constant)"
    legend_constants[2] = "beta in component infected (first_order_rate_constant)"
    legend_states[0] = "v in component virus (dimensionless)"
    legend_states[1] = "x in component uninfected (dimensionless)"
    legend_constants[3] = "a in component infected (first_order_rate_constant)"
    legend_constants[4] = "p in component infected (first_order_rate_constant)"
    legend_states[2] = "z in component CTL (dimensionless)"
    legend_states[3] = "y in component infected (dimensionless)"
    legend_constants[5] = "k in component virus (first_order_rate_constant)"
    legend_constants[6] = "u in component virus (first_order_rate_constant)"
    legend_constants[7] = "c in component CTL (first_order_rate_constant)"
    legend_constants[8] = "b in component CTL (first_order_rate_constant)"
    legend_rates[1] = "d/dt x in component uninfected (dimensionless)"
    legend_rates[3] = "d/dt y in component infected (dimensionless)"
    legend_rates[0] = "d/dt v in component virus (dimensionless)"
    legend_rates[2] = "d/dt z in component CTL (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1
    constants[1] = 0.01
    constants[2] = 0.02
    states[0] = 10
    states[1] = 100
    constants[3] = 0.5
    constants[4] = 1
    states[2] = 1
    states[3] = 0
    constants[5] = 1
    constants[6] = 1
    constants[7] = 1
    constants[8] = 0.05
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = (constants[0]-constants[1]*states[1])-constants[2]*states[1]*states[0]
    rates[3] = (constants[2]*states[1]*states[0]-constants[3]*states[3])-constants[4]*states[3]*states[2]
    rates[0] = constants[5]*states[3]-constants[6]*states[0]
    rates[2] = constants[7]*states[3]*states[2]-constants[8]*states[2]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
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