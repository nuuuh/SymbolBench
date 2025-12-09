# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 3
sizeConstants = 11
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "X1 in component X1 (ng_ml)"
    legend_constants[0] = "a0 in component X1 (ng_ml_min)"
    legend_constants[1] = "a1 in component X1 (ng_ml_min)"
    legend_constants[2] = "a2 in component X1 (dl_microg)"
    legend_constants[3] = "b1 in component X1 (first_order_rate_constant)"
    legend_states[1] = "X3 in component X3 (microg_dl)"
    legend_states[2] = "X2 in component X2 (pg_ml)"
    legend_constants[4] = "a3 in component X2 (pg_ml_min)"
    legend_constants[5] = "a4 in component X2 (first_order_rate_constant)"
    legend_constants[6] = "a5 in component X2 (dl_microg)"
    legend_constants[7] = "b2 in component X2 (first_order_rate_constant)"
    legend_constants[8] = "a6 in component X3 (microg_dl_min)"
    legend_constants[9] = "a7 in component X3 (first_order_rate_constant)"
    legend_constants[10] = "b3 in component X3 (first_order_rate_constant)"
    legend_rates[0] = "d/dt X1 in component X1 (ng_ml)"
    legend_rates[2] = "d/dt X2 in component X2 (pg_ml)"
    legend_rates[1] = "d/dt X3 in component X3 (microg_dl)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 10.0
    constants[0] = 0.0014
    constants[1] = 0.000517
    constants[2] = 0.0164
    constants[3] = 0.0598
    states[1] = 0.0
    states[2] = 0.0
    constants[4] = 1.38
    constants[5] = 0.60
    constants[6] = 0.00498
    constants[7] = 0.053
    constants[8] = 0.0084
    constants[9] = 0.0081
    constants[10] = 0.0138
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (constants[0]+constants[1]/(1.00000+constants[2]*states[1]))-constants[3]*states[0]
    rates[2] = (constants[4]+constants[5]*states[0])/(1.00000+constants[6]*states[1])-constants[7]*states[2]
    rates[1] = (constants[8]+constants[9]*states[2])-constants[10]*states[1]
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