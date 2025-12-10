# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 9
sizeConstants = 10
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "x1 in component x1 (nanomolar)"
    legend_states[1] = "x2 in component x2 (nanomolar)"
    legend_states[2] = "x3 in component x3 (nanomolar)"
    legend_states[3] = "x4 in component x4 (nanomolar)"
    legend_states[4] = "x5 in component x5 (nanomolar)"
    legend_states[5] = "x6 in component x6 (nanomolar)"
    legend_states[6] = "x7 in component x7 (nanomolar)"
    legend_states[7] = "x8 in component x8 (nanomolar)"
    legend_states[8] = "x9 in component x9 (nanomolar)"
    legend_constants[0] = "k1 in component rate_variables (second_order_rate_constant)"
    legend_constants[1] = "k1_ in component rate_variables (first_order_rate_constant)"
    legend_constants[2] = "k2 in component rate_variables (second_order_rate_constant)"
    legend_constants[3] = "k2_ in component rate_variables (first_order_rate_constant)"
    legend_constants[4] = "k4 in component rate_variables (second_order_rate_constant)"
    legend_constants[5] = "k4_ in component rate_variables (first_order_rate_constant)"
    legend_constants[6] = "k5 in component rate_variables (second_order_rate_constant)"
    legend_constants[7] = "k5_ in component rate_variables (first_order_rate_constant)"
    legend_algebraic[0] = "scatchard in component x1 (dimensionless)"
    legend_constants[8] = "k3 in component rate_variables (second_order_rate_constant)"
    legend_constants[9] = "k3_ in component rate_variables (first_order_rate_constant)"
    legend_rates[0] = "d/dt x1 in component x1 (nanomolar)"
    legend_rates[1] = "d/dt x2 in component x2 (nanomolar)"
    legend_rates[2] = "d/dt x3 in component x3 (nanomolar)"
    legend_rates[3] = "d/dt x4 in component x4 (nanomolar)"
    legend_rates[4] = "d/dt x5 in component x5 (nanomolar)"
    legend_rates[5] = "d/dt x6 in component x6 (nanomolar)"
    legend_rates[6] = "d/dt x7 in component x7 (nanomolar)"
    legend_rates[7] = "d/dt x8 in component x8 (nanomolar)"
    legend_rates[8] = "d/dt x9 in component x9 (nanomolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 10
    states[1] = 0.1
    states[2] = 1.0
    states[3] = 1.0
    states[4] = 1.0
    states[5] = 1.0
    states[6] = 1.0
    states[7] = 1.0
    states[8] = 1.0
    constants[0] = 1000000
    constants[1] = 0.0004
    constants[2] = 1000000
    constants[3] = 0.04
    constants[4] = 1000000
    constants[5] = 0.0004
    constants[6] = 10000000
    constants[7] = 0.004
    constants[8] = 1000000
    constants[9] = 0.0004
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = ((((((constants[1]*states[2]-constants[0]*states[0]*states[1])+constants[3]*states[3])-constants[2]*states[0]*states[2])+constants[5]*(states[5]+states[6]))-constants[4]*states[0]*(states[4]+states[8]))+constants[7]*(states[6]+states[7]+states[8]))-constants[6]*states[0]*(states[4]+states[5]+states[6])
    rates[1] = ((constants[1]*states[2]-constants[0]*states[0]*states[1])+constants[9]*(states[4]+states[8]))-constants[8]*states[1]*(states[2]+states[3])
    rates[2] = (((((constants[0]*states[0]*states[1]-constants[1]*states[2])+constants[3]*states[3])-constants[2]*states[0]*states[2])+constants[9]*states[4])-constants[8]*states[1]*states[2])
    rates[3] = ((constants[2]*states[0]*states[2]-constants[3]*states[3])+constants[9]*states[8])-constants[2]*states[1]*states[3]
    rates[4] = ((((constants[8]*states[1]*states[2]-constants[9]*states[4])+constants[5]*states[5])-constants[4]*states[0]*states[4])+constants[7]*states[8])-constants[6]*states[0]*states[4]
    rates[5] = ((constants[4]*states[0]*states[4]-constants[5]*states[5])+constants[7]*states[6])-constants[6]*states[0]*states[5]
    rates[6] = (constants[4]*states[0]*states[8]-constants[5]*states[6])+constants[7]*(states[7]-states[6])+constants[6]*states[0]*(states[5]-states[6])
    rates[7] = constants[6]*states[0]*states[6]-constants[7]*states[7]
    rates[8] = ((((constants[8]*states[1]*states[3]-constants[9]*states[8])+constants[5]*states[6])-constants[4]*states[0]*states[8])+constants[6]*states[0]*states[4])-constants[7]*states[8]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = (states[2]+states[3]+states[6]+states[7]+states[8])/states[0]
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