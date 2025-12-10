# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 2
sizeConstants = 3
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "t in component Main (millisecond)"
    legend_states[0] = "v in component Main (dimensionless)"
    legend_states[1] = "w in component Main (dimensionless)"
    legend_constants[0] = "alpha in component Main (dimensionless)"
    legend_constants[1] = "gamma in component Main (dimensionless)"
    legend_constants[2] = "epsilon in component Main (dimensionless)"
    legend_algebraic[0] = "I in component Main (dimensionless)"
    legend_rates[0] = "d/dt v in component Main (dimensionless)"
    legend_rates[1] = "d/dt w in component Main (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0
    states[1] = 0
    constants[0] = -0.08
    constants[1] = 3
    constants[2] = 0.005
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = 1.00000*constants[2]*(states[0]-constants[1]*states[1])
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less_equal(voi , 0.500000), -80.0000 , True, 0.00000])
    rates[0] = 1.00000*((states[0]*(states[0]-constants[0])*(1.00000-states[0])-states[1])+algebraic[0])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less_equal(voi , 0.500000), -80.0000 , True, 0.00000])
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