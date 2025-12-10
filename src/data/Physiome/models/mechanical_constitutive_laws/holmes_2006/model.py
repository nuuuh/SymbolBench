# Size of variable arrays:
sizeAlgebraic = 4
sizeStates = 1
sizeConstants = 5
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "a in component contraction (mNpermmsq)"
    legend_constants[1] = "b in component contraction (pms)"
    legend_constants[2] = "Po in component contraction (mNpermmsq)"
    legend_constants[3] = "alpha in component contraction (mNpermmsq)"
    legend_constants[4] = "L_se_o in component contraction (dimensionless)"
    legend_algebraic[0] = "L in component contraction (dimensionless)"
    legend_algebraic[3] = "v in component contraction (pms)"
    legend_algebraic[1] = "L_se in component contraction (dimensionless)"
    legend_states[0] = "L_ce in component contraction (dimensionless)"
    legend_algebraic[2] = "P in component contraction (mNpermmsq)"
    legend_rates[0] = "d/dt L_ce in component contraction (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 37.24
    constants[1] = 0.325
    constants[2] = 144.9
    constants[3] = 1449.027
    constants[4] = 0.3
    states[0] = 0.7
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = custom_piecewise([less_equal(voi , 1.00000), 1.00000 , greater(voi , 1.00000) & less(voi , 5.00000), 0.920000 , True, 0.900000])
    algebraic[1] = algebraic[0]-states[0]
    algebraic[2] = constants[3]*(algebraic[1]-constants[4])
    algebraic[3] = (-constants[1]*(constants[2]-algebraic[2]))/(algebraic[2]+constants[0])
    rates[0] = algebraic[3]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less_equal(voi , 1.00000), 1.00000 , greater(voi , 1.00000) & less(voi , 5.00000), 0.920000 , True, 0.900000])
    algebraic[1] = algebraic[0]-states[0]
    algebraic[2] = constants[3]*(algebraic[1]-constants[4])
    algebraic[3] = (-constants[1]*(constants[2]-algebraic[2]))/(algebraic[2]+constants[0])
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