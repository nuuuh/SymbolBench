# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 8
sizeConstants = 17
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
    legend_constants[0] = "C1C2 in component reaction_constants (second_order_rate_constant)"
    legend_constants[1] = "C2C1 in component reaction_constants (first_order_rate_constant)"
    legend_states[1] = "C2 in component C2 (dimensionless)"
    legend_constants[2] = "Ca in component reaction_constants (micromolar)"
    legend_constants[3] = "C2C3 in component reaction_constants (second_order_rate_constant)"
    legend_constants[4] = "C3C2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[5] = "C2C5 in component reaction_constants (first_order_rate_constant)"
    legend_constants[6] = "C5C2 in component reaction_constants (first_order_rate_constant)"
    legend_states[2] = "C3 in component C3 (dimensionless)"
    legend_states[3] = "C5 in component C5 (dimensionless)"
    legend_constants[7] = "O1C3 in component reaction_constants (first_order_rate_constant)"
    legend_constants[8] = "C3O1 in component reaction_constants (first_order_rate_constant)"
    legend_constants[9] = "O2C3 in component reaction_constants (first_order_rate_constant)"
    legend_constants[10] = "C3O2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[11] = "O3C3 in component reaction_constants (first_order_rate_constant)"
    legend_constants[12] = "C3O3 in component reaction_constants (first_order_rate_constant)"
    legend_states[4] = "O1 in component O1 (dimensionless)"
    legend_states[5] = "O2 in component O2 (dimensionless)"
    legend_states[6] = "O3 in component O3 (dimensionless)"
    legend_states[7] = "C4 in component C4 (dimensionless)"
    legend_constants[13] = "O2C4 in component reaction_constants (first_order_rate_constant)"
    legend_constants[14] = "C4O2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[15] = "O3C4 in component reaction_constants (first_order_rate_constant)"
    legend_constants[16] = "C4O3 in component reaction_constants (first_order_rate_constant)"
    legend_rates[0] = "d/dt C1 in component C1 (dimensionless)"
    legend_rates[1] = "d/dt C2 in component C2 (dimensionless)"
    legend_rates[2] = "d/dt C3 in component C3 (dimensionless)"
    legend_rates[7] = "d/dt C4 in component C4 (dimensionless)"
    legend_rates[3] = "d/dt C5 in component C5 (dimensionless)"
    legend_rates[4] = "d/dt O1 in component O1 (dimensionless)"
    legend_rates[5] = "d/dt O2 in component O2 (dimensionless)"
    legend_rates[6] = "d/dt O3 in component O3 (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.125
    constants[0] = 2.5
    constants[1] = 13.3
    states[1] = 0.125
    constants[2] = 20.0
    constants[3] = 68.0
    constants[4] = 8000.0
    constants[5] = 0.13
    constants[6] = 3.6
    states[2] = 0.125
    states[3] = 0.125
    constants[7] = 3400.0
    constants[8] = 1100.0
    constants[9] = 92.0
    constants[10] = 17.0
    constants[11] = 138.0
    constants[12] = 14.0
    states[4] = 0.125
    states[5] = 0.125
    states[6] = 0.125
    states[7] = 0.125
    constants[13] = 1900.0
    constants[14] = 520.0
    constants[15] = 300.0
    constants[16] = 46.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[1]*states[1]-constants[0]*constants[2]*states[1]
    rates[1] = (constants[0]*constants[2]*states[0]+constants[4]*states[2]+constants[6]*states[3])-(constants[1]*states[1]+constants[3]*constants[2]*states[1]+constants[5]*states[1])
    rates[2] = (constants[3]*constants[2]*states[1]+constants[7]*states[4]+constants[11]*states[6]+constants[9]*states[5])-(constants[4]*states[2]+constants[8]*states[2]+constants[10]*states[2]+constants[12]*states[2])
    rates[7] = (constants[13]*states[5]+constants[15]*states[6])-(constants[14]*states[7]+constants[16]*states[7])
    rates[3] = constants[5]*states[1]-constants[6]*states[3]
    rates[4] = constants[8]*states[2]-constants[7]*states[4]
    rates[5] = (constants[10]*states[2]+constants[14]*states[7])-(constants[9]*states[5]+constants[13]*states[5])
    rates[6] = (constants[12]*states[2]+constants[16]*states[7])-(constants[11]*states[6]+constants[15]*states[6])
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