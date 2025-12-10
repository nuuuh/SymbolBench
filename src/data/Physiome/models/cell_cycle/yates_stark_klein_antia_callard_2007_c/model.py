# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 4
sizeConstants = 10
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (day)"
    legend_states[0] = "x in component x (dimensionless)"
    legend_constants[0] = "d1 in component model_parameters (first_order_rate_constant)"
    legend_constants[1] = "a in component model_parameters (first_order_rate_constant)"
    legend_constants[8] = "alpha in component model_parameters (first_order_rate_constant)"
    legend_constants[2] = "mu in component model_parameters (first_order_rate_constant)"
    legend_constants[3] = "r in component model_parameters (first_order_rate_constant)"
    legend_states[1] = "y in component y (dimensionless)"
    legend_states[2] = "w in component w (dimensionless)"
    legend_constants[4] = "kappa in component model_parameters (dimensionless)"
    legend_constants[9] = "d2 in component y (first_order_rate_constant)"
    legend_constants[5] = "p in component model_parameters (first_order_rate_constant)"
    legend_states[3] = "z in component z (dimensionless)"
    legend_constants[6] = "f in component w (dimensionless)"
    legend_constants[7] = "v in component z (first_order_rate_constant)"
    legend_rates[0] = "d/dt x in component x (dimensionless)"
    legend_rates[1] = "d/dt y in component y (dimensionless)"
    legend_rates[2] = "d/dt w in component w (dimensionless)"
    legend_rates[3] = "d/dt z in component z (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 10E-1
    constants[0] = 0.005
    constants[1] = 0.03333
    constants[2] = 0.023
    constants[3] = 1.0
    states[1] = 0.0
    states[2] = 0.0
    constants[4] = 1.0
    constants[5] = 200.0
    states[3] = 0.01
    constants[6] = 100.0
    constants[7] = 0.5
    constants[8] = 0.0100000*constants[1]
    constants[9] = (-(99.0000*constants[1]*constants[0])+constants[1]*constants[3]+constants[0]*constants[3])/(constants[1]-constants[0])
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (2.00000*constants[3]*states[1]+0.100000*constants[2]*states[2])-(constants[1]*states[0]*(1.00000-states[0]/constants[4])+constants[0]*states[0]*(states[0]/constants[4])+constants[8]*states[0])
    rates[1] = constants[1]*states[0]*(1.00000-states[0]/constants[4])-((constants[3]+constants[9])*states[1]+constants[5]*states[3]*states[1])
    rates[2] = constants[6]*constants[8]*states[0]-(constants[2]*states[2]+constants[5]*states[2]*states[3])
    rates[3] = constants[5]*states[3]*(states[1]+states[2])-constants[7]*states[3]
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