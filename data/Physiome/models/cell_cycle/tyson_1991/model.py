# Size of variable arrays:
sizeAlgebraic = 7
sizeStates = 6
sizeConstants = 10
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "C2 in component C2 (millimolar)"
    legend_constants[0] = "k6 in component reaction_constants (first_order_rate_constant)"
    legend_constants[8] = "k8 in component reaction_constants (second_order_rate_constant)"
    legend_constants[1] = "k9 in component reaction_constants (first_order_rate_constant)"
    legend_states[1] = "M in component M (millimolar)"
    legend_constants[2] = "P in component reaction_constants (millimolar)"
    legend_states[2] = "CP in component CP (millimolar)"
    legend_algebraic[4] = "k3 in component reaction_constants (second_order_rate_constant)"
    legend_states[3] = "Y in component Y (millimolar)"
    legend_states[4] = "pM in component pM (millimolar)"
    legend_constants[9] = "k5 in component reaction_constants (second_order_rate_constant)"
    legend_algebraic[6] = "F in component reaction_constants (first_order_rate_constant)"
    legend_algebraic[1] = "k1 in component reaction_constants (first_order_rate_constant)"
    legend_constants[3] = "k2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[4] = "aa in component reaction_constants (millimolar)"
    legend_states[5] = "YP in component YP (millimolar)"
    legend_constants[5] = "k7 in component reaction_constants (first_order_rate_constant)"
    legend_constants[6] = "k4 in component reaction_constants (first_order_rate_constant)"
    legend_constants[7] = "k4_ in component reaction_constants (first_order_rate_constant)"
    legend_algebraic[0] = "CT in component reaction_constants (millimolar)"
    legend_algebraic[2] = "YT in component reaction_constants (millimolar)"
    legend_algebraic[5] = "YT_CT in component reaction_constants (dimensionless)"
    legend_algebraic[3] = "M_CT in component reaction_constants (dimensionless)"
    legend_rates[0] = "d/dt C2 in component C2 (millimolar)"
    legend_rates[2] = "d/dt CP in component CP (millimolar)"
    legend_rates[4] = "d/dt pM in component pM (millimolar)"
    legend_rates[1] = "d/dt M in component M (millimolar)"
    legend_rates[3] = "d/dt Y in component Y (millimolar)"
    legend_rates[5] = "d/dt YP in component YP (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.001
    constants[0] = 1.0
    constants[1] = 100.0
    states[1] = 0.001
    constants[2] = 1.0
    states[2] = 0.001
    states[3] = 0.001
    states[4] = 0.001
    constants[3] = 0.0
    constants[4] = 1.0
    states[5] = 0.001
    constants[5] = 0.6
    constants[6] = 180.0
    constants[7] = 0.018
    constants[8] = constants[1]/constants[2]
    constants[9] = 0.00000/constants[2]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (constants[0]*states[1]+constants[1]*states[2])-constants[8]*constants[2]*states[0]
    rates[5] = constants[0]*states[1]-constants[5]*states[5]
    algebraic[0] = states[0]+states[2]+states[4]+states[1]
    algebraic[4] = 200.000/algebraic[0]
    rates[2] = constants[8]*constants[2]*states[0]-(constants[1]*states[2]+algebraic[4]*states[2]*states[3])
    algebraic[1] = (0.0150000*algebraic[0])/constants[4]
    rates[3] = constants[4]*algebraic[1]-(algebraic[4]*states[2]*states[3]+constants[3]*states[3])
    algebraic[6] = constants[7]+constants[6]*(power(states[1]/algebraic[0], 2.00000))
    rates[4] = (algebraic[4]*states[2]*states[3]+constants[9]*constants[2]*states[1])-states[4]*algebraic[6]
    rates[1] = states[4]*algebraic[6]-(constants[9]*constants[2]*states[1]+constants[0]*states[1])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[0]+states[2]+states[4]+states[1]
    algebraic[4] = 200.000/algebraic[0]
    algebraic[1] = (0.0150000*algebraic[0])/constants[4]
    algebraic[6] = constants[7]+constants[6]*(power(states[1]/algebraic[0], 2.00000))
    algebraic[2] = states[3]+states[5]+states[4]+states[1]
    algebraic[3] = states[1]/algebraic[0]
    algebraic[5] = algebraic[2]/algebraic[0]
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