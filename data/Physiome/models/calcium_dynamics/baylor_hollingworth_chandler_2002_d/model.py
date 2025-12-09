# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 5
sizeConstants = 8
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "Ca in component Ca (micromolar)"
    legend_constants[0] = "k1 in component reaction_constants (second_order_rate_constant)"
    legend_constants[1] = "k1_ in component reaction_constants (first_order_rate_constant)"
    legend_constants[2] = "k3 in component reaction_constants (second_order_rate_constant)"
    legend_constants[3] = "k3_ in component reaction_constants (first_order_rate_constant)"
    legend_constants[4] = "k4 in component reaction_constants (first_order_rate_constant)"
    legend_constants[5] = "k4_ in component reaction_constants (second_order_rate_constant)"
    legend_states[1] = "E in component E (micromolar)"
    legend_states[2] = "CaE in component CaE (micromolar)"
    legend_states[3] = "Ca2E_ in component Ca2E_ (micromolar)"
    legend_states[4] = "CaE_ in component CaE_ (micromolar)"
    legend_constants[6] = "k2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[7] = "k2_ in component reaction_constants (first_order_rate_constant)"
    legend_rates[0] = "d/dt Ca in component Ca (micromolar)"
    legend_rates[1] = "d/dt E in component E (micromolar)"
    legend_rates[2] = "d/dt CaE in component CaE (micromolar)"
    legend_rates[3] = "d/dt Ca2E_ in component Ca2E_ (micromolar)"
    legend_rates[4] = "d/dt CaE_ in component CaE_ (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.05
    constants[0] = 1.0E14
    constants[1] = 1000.0
    constants[2] = 1.0E14
    constants[3] = 10.0
    constants[4] = 20.0
    constants[5] = 0.0
    states[1] = 240.0
    states[2] = 0.01
    states[3] = 1.0
    states[4] = 0.01
    constants[6] = 500.0
    constants[7] = 1200.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (constants[1]*states[2]+constants[3]*states[3]+constants[4]*states[3])-(constants[0]*states[0]*states[1]+constants[2]*states[0]*states[4]+constants[5]*2.00000*states[0]*states[1])
    rates[1] = (constants[1]*states[2]+constants[4]*states[3])-(constants[0]*states[0]*states[1]+constants[5]*2.00000*states[0]*states[1])
    rates[2] = (constants[0]*states[0]*states[1]+constants[7]*states[4])-(constants[1]*states[2]+constants[6]*states[2])
    rates[3] = (constants[2]*states[4]*states[0]+constants[5]*states[1]*2.00000*states[0])-(constants[3]*states[3]+constants[4]*states[3])
    rates[4] = (constants[3]*states[3]+constants[6]*states[2])-(constants[7]*states[4]+constants[2]*states[0]*states[4])
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