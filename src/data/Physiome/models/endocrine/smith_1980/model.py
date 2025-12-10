# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 3
sizeConstants = 7
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_states[0] = "R in component R (ng_ml)"
    legend_constants[0] = "h in component R (first_order_rate_constant)"
    legend_constants[1] = "c in component R (ng_ml_hr)"
    legend_constants[2] = "b1 in component R (first_order_rate_constant)"
    legend_algebraic[1] = "H in component R (dimensionless)"
    legend_algebraic[0] = "x in component R (ng_ml)"
    legend_states[1] = "T in component T (ng_ml)"
    legend_states[2] = "L in component L (ng_ml)"
    legend_constants[3] = "g1 in component L (first_order_rate_constant)"
    legend_constants[4] = "b2 in component L (first_order_rate_constant)"
    legend_constants[5] = "g2 in component T (first_order_rate_constant)"
    legend_constants[6] = "b3 in component T (first_order_rate_constant)"
    legend_rates[0] = "d/dt R in component R (ng_ml)"
    legend_rates[2] = "d/dt L in component L (ng_ml)"
    legend_rates[1] = "d/dt T in component T (ng_ml)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.5
    constants[0] = 12
    constants[1] = 100
    constants[2] = 1.29
    states[1] = 15
    states[2] = 22
    constants[3] = 10
    constants[4] = 0.97
    constants[5] = 0.7
    constants[6] = 1.39
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[2] = constants[3]*states[0]-constants[4]*states[2]
    rates[1] = constants[5]*states[2]-constants[6]*states[1]
    algebraic[0] = states[1]-constants[1]/constants[0]
    algebraic[1] = custom_piecewise([greater(algebraic[0] , 0.00000), 1.00000 , True, 0.00000])
    rates[0] = (constants[1]-constants[0]*states[1])*(1.00000-algebraic[1])-constants[2]*states[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[1]-constants[1]/constants[0]
    algebraic[1] = custom_piecewise([greater(algebraic[0] , 0.00000), 1.00000 , True, 0.00000])
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