# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 2
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
    legend_states[0] = "GSH in component GSH (millimolar)"
    legend_constants[0] = "V in component GSH (flux)"
    legend_constants[1] = "Kms in component GSH (millimolar)"
    legend_constants[2] = "Kmp in component GSH (millimolar)"
    legend_constants[3] = "Kmq in component GSH (millimolar)"
    legend_constants[4] = "Keq in component GSH (millimolar)"
    legend_states[1] = "SDLGSH in component SDLGSH (millimolar)"
    legend_constants[5] = "D_lactate in component D_lactate (millimolar)"
    legend_constants[6] = "HTA in component HTA (millimolar)"
    legend_constants[7] = "V in component SDLGSH (flux)"
    legend_constants[8] = "Kms in component SDLGSH (millimolar)"
    legend_constants[9] = "Kmp in component SDLGSH (millimolar)"
    legend_constants[10] = "Keq in component SDLGSH (dimensionless)"
    legend_rates[0] = "d/dt GSH in component GSH (millimolar)"
    legend_rates[1] = "d/dt SDLGSH in component SDLGSH (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 1.0
    constants[0] = 3.44E-3
    constants[1] = 0.49
    constants[2] = 0.49
    constants[3] = 0.49
    constants[4] = 0.49
    states[1] = 1.0
    constants[5] = 0.0
    constants[6] = 1.0
    constants[7] = 8.12E-2
    constants[8] = 0.61
    constants[9] = 0.61
    constants[10] = 0.61
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = ((constants[0]/constants[1])*(states[1]-(states[0]*constants[5])/constants[4]))/(1.00000+states[1]/constants[1]+states[0]/constants[2]+constants[5]/constants[3])
    rates[1] = ((constants[7]/constants[8])*(constants[6]-states[1]/constants[10]))/(1.00000+constants[6]/constants[8]+states[1]/constants[9])
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