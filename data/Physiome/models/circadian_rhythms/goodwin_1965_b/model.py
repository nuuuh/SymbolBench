# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 4
sizeConstants = 14
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "X1 in component X1 (nanomolar)"
    legend_constants[0] = "a1 in component X1 (flux)"
    legend_constants[1] = "b1 in component X1 (flux)"
    legend_constants[2] = "A1 in component X1 (dimensionless)"
    legend_constants[3] = "k11 in component X1 (per_nanomolar)"
    legend_constants[4] = "k12 in component X1 (per_nanomolar)"
    legend_states[1] = "Y1 in component Y1 (nanomolar)"
    legend_states[2] = "Y2 in component Y2 (nanomolar)"
    legend_constants[5] = "beta_1 in component Y1 (flux)"
    legend_constants[6] = "alpha_1 in component Y1 (first_order_rate_constant)"
    legend_states[3] = "X2 in component X2 (nanomolar)"
    legend_constants[7] = "a2 in component X2 (flux)"
    legend_constants[8] = "b2 in component X2 (flux)"
    legend_constants[9] = "A2 in component X2 (dimensionless)"
    legend_constants[10] = "k21 in component X2 (per_nanomolar)"
    legend_constants[11] = "k22 in component X2 (per_nanomolar)"
    legend_constants[12] = "beta_2 in component Y2 (flux)"
    legend_constants[13] = "alpha_2 in component Y2 (first_order_rate_constant)"
    legend_rates[0] = "d/dt X1 in component X1 (nanomolar)"
    legend_rates[1] = "d/dt Y1 in component Y1 (nanomolar)"
    legend_rates[3] = "d/dt X2 in component X2 (nanomolar)"
    legend_rates[2] = "d/dt Y2 in component Y2 (nanomolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 7
    constants[0] = 360
    constants[1] = 10
    constants[2] = 36
    constants[3] = 1
    constants[4] = 0
    states[1] = -10
    states[2] = -10
    constants[5] = 0
    constants[6] = 0.5
    states[3] = 7
    constants[7] = 360
    constants[8] = 10
    constants[9] = 43
    constants[10] = 0
    constants[11] = 1
    constants[12] = 0
    constants[13] = 0.6
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[0]/(constants[2]+constants[3]*states[1]+constants[4]*states[2])-constants[1]
    rates[1] = constants[6]*states[0]-constants[5]
    rates[3] = constants[7]/(constants[9]+constants[10]*states[1]+constants[11]*states[2])-constants[8]
    rates[2] = constants[13]*states[3]-constants[12]
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