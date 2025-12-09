# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 4
sizeConstants = 7
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "C1 in component C1 (dimensionless)"
    legend_constants[0] = "O2C1 in component reaction_constants (second_order_rate_constant)"
    legend_constants[1] = "C1O2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[2] = "C1C2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[3] = "C2C1 in component reaction_constants (first_order_rate_constant)"
    legend_states[1] = "C2 in component C2 (dimensionless)"
    legend_states[2] = "O2 in component O2 (dimensionless)"
    legend_constants[4] = "Ca in component reaction_constants (micromolar)"
    legend_states[3] = "O1 in component O1 (dimensionless)"
    legend_constants[5] = "O1O2 in component reaction_constants (second_order_rate_constant)"
    legend_constants[6] = "O2O1 in component reaction_constants (first_order_rate_constant)"
    legend_rates[0] = "d/dt C1 in component C1 (dimensionless)"
    legend_rates[1] = "d/dt C2 in component C2 (dimensionless)"
    legend_rates[3] = "d/dt O1 in component O1 (dimensionless)"
    legend_rates[2] = "d/dt O2 in component O2 (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.25
    constants[0] = 2.62
    constants[1] = 2277.0
    constants[2] = 60.8
    constants[3] = 198.0
    states[1] = 0.25
    states[2] = 0.25
    constants[4] = 50.0
    states[3] = 0.25
    constants[5] = 2.41
    constants[6] = 85.1
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (constants[0]*constants[4]*states[2]+constants[3]*states[1])-(constants[1]*states[0]+constants[2]*states[0])
    rates[1] = constants[2]*states[0]-constants[3]*states[1]
    rates[3] = constants[6]*states[2]-constants[5]*constants[4]*states[3]
    rates[2] = (constants[5]*constants[4]*states[3]+constants[1]*states[0])-(constants[6]*states[2]+constants[0]*constants[4]*states[2])
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