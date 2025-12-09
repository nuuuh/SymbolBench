# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 5
sizeConstants = 6
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_algebraic[0] = "GnRH in component GnRH (nanomolar)"
    legend_states[0] = "F in component F (dimensionless)"
    legend_constants[0] = "kfb in component model_parameters (second_order_rate_constant)"
    legend_constants[1] = "kdf in component model_parameters (first_order_rate_constant)"
    legend_states[1] = "D in component D (dimensionless)"
    legend_constants[2] = "kbd in component model_parameters (first_order_rate_constant)"
    legend_states[2] = "B in component B (dimensionless)"
    legend_states[3] = "R in component R (dimensionless)"
    legend_constants[3] = "s in component model_parameters (first_order_rate_constant)"
    legend_constants[4] = "a1 in component model_parameters (first_order_rate_constant)"
    legend_constants[5] = "a2 in component model_parameters (first_order_rate_constant)"
    legend_states[4] = "C in component C (dimensionless)"
    legend_rates[0] = "d/dt F in component F (dimensionless)"
    legend_rates[1] = "d/dt D in component D (dimensionless)"
    legend_rates[2] = "d/dt B in component B (dimensionless)"
    legend_rates[3] = "d/dt R in component R (dimensionless)"
    legend_rates[4] = "d/dt C in component C (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 1.0
    constants[0] = 19.35
    constants[1] = 2.52
    states[1] = 0.0
    constants[2] = 9.91
    states[2] = 0.0
    states[3] = 2115.0
    constants[3] = 6.80
    constants[4] = 0.0328
    constants[5] = 0.0830
    states[4] = 0.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = constants[2]*states[2]-constants[1]*states[1]
    rates[3] = constants[3]-(constants[4]+constants[5]*states[2])*states[3]
    rates[4] = (constants[4]+constants[5]*states[2])*states[3]
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 0.0666667), 0.500000 , greater_equal(voi , 0.0666667) & less(voi , 0.400000), 0.00000 , greater_equal(voi , 0.400000) & less(voi , 0.466667), 0.500000 , greater_equal(voi , 0.466667) & less(voi , 2.46667), 0.00000 , greater_equal(voi , 2.46667) & less(voi , 2.53333), 0.500000 , greater_equal(voi , 2.53333) & less(voi , 2.61667), 0.00000 , greater_equal(voi , 2.61667) & less(voi , 2.68333), 0.500000 , greater_equal(voi , 2.68333) & less(voi , 4.68333), 0.00000 , greater_equal(voi , 4.68333) & less(voi , 4.75000), 0.500000 , greater_equal(voi , 4.75000) & less(voi , 4.91667), 0.00000 , greater_equal(voi , 4.91667) & less(voi , 4.98333), 0.500000 , greater_equal(voi , 4.98333) & less(voi , 6.98333), 0.00000 , greater_equal(voi , 6.98333) & less(voi , 7.06667), 0.500000 , greater_equal(voi , 7.06667) & less(voi , 7.73333), 0.00000 , greater_equal(voi , 7.73333) & less(voi , 7.80000), 0.500000 , greater_equal(voi , 7.80000) & less(voi , 9.80000), 0.00000 , True, float('nan')])
    rates[0] = constants[1]*states[1]-constants[0]*states[0]*algebraic[0]
    rates[2] = constants[0]*states[0]*algebraic[0]-constants[2]*states[2]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 0.0666667), 0.500000 , greater_equal(voi , 0.0666667) & less(voi , 0.400000), 0.00000 , greater_equal(voi , 0.400000) & less(voi , 0.466667), 0.500000 , greater_equal(voi , 0.466667) & less(voi , 2.46667), 0.00000 , greater_equal(voi , 2.46667) & less(voi , 2.53333), 0.500000 , greater_equal(voi , 2.53333) & less(voi , 2.61667), 0.00000 , greater_equal(voi , 2.61667) & less(voi , 2.68333), 0.500000 , greater_equal(voi , 2.68333) & less(voi , 4.68333), 0.00000 , greater_equal(voi , 4.68333) & less(voi , 4.75000), 0.500000 , greater_equal(voi , 4.75000) & less(voi , 4.91667), 0.00000 , greater_equal(voi , 4.91667) & less(voi , 4.98333), 0.500000 , greater_equal(voi , 4.98333) & less(voi , 6.98333), 0.00000 , greater_equal(voi , 6.98333) & less(voi , 7.06667), 0.500000 , greater_equal(voi , 7.06667) & less(voi , 7.73333), 0.00000 , greater_equal(voi , 7.73333) & less(voi , 7.80000), 0.500000 , greater_equal(voi , 7.80000) & less(voi , 9.80000), 0.00000 , True, float('nan')])
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