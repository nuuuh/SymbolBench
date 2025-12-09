# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 6
sizeConstants = 19
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "y1 in component y1 (micromolar)"
    legend_constants[0] = "a1 in component y1 (flux)"
    legend_algebraic[0] = "g1 in component y1 (micromolar)"
    legend_constants[1] = "b1 in component y1 (micromolar)"
    legend_constants[2] = "d1 in component y1 (first_order_rate_constant)"
    legend_states[1] = "y2 in component y2 (micromolar)"
    legend_constants[3] = "a2 in component y2 (flux)"
    legend_algebraic[1] = "g2 in component y2 (micromolar)"
    legend_constants[4] = "b2 in component y2 (micromolar)"
    legend_constants[5] = "d2 in component y2 (first_order_rate_constant)"
    legend_states[2] = "y3 in component y3 (micromolar)"
    legend_constants[6] = "f53 in component y3 (second_order_rate_constant)"
    legend_constants[7] = "f13 in component y3 (second_order_rate_constant)"
    legend_constants[8] = "h36 in component y3 (second_order_rate_constant)"
    legend_constants[9] = "d3 in component y3 (first_order_rate_constant)"
    legend_constants[10] = "E in component y3 (micromolar)"
    legend_states[3] = "y5 in component y5 (micromolar)"
    legend_states[4] = "y6 in component y6 (micromolar)"
    legend_states[5] = "y4 in component y4 (micromolar)"
    legend_constants[11] = "f14 in component y4 (first_order_rate_constant)"
    legend_constants[12] = "f24 in component y4 (first_order_rate_constant)"
    legend_constants[13] = "d4 in component y4 (first_order_rate_constant)"
    legend_constants[14] = "f35 in component y5 (first_order_rate_constant)"
    legend_constants[15] = "f45 in component y5 (first_order_rate_constant)"
    legend_constants[16] = "d5 in component y5 (first_order_rate_constant)"
    legend_constants[17] = "h36 in component y6 (second_order_rate_constant)"
    legend_constants[18] = "d6 in component y6 (first_order_rate_constant)"
    legend_rates[0] = "d/dt y1 in component y1 (micromolar)"
    legend_rates[1] = "d/dt y2 in component y2 (micromolar)"
    legend_rates[2] = "d/dt y3 in component y3 (micromolar)"
    legend_rates[5] = "d/dt y4 in component y4 (micromolar)"
    legend_rates[3] = "d/dt y5 in component y5 (micromolar)"
    legend_rates[4] = "d/dt y6 in component y6 (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.0
    constants[0] = 10.0
    constants[1] = 10.0
    constants[2] = 0.2
    states[1] = 0.0
    constants[3] = 10.0
    constants[4] = 10.0
    constants[5] = 0.1
    states[2] = 0.0
    constants[6] = 1.5
    constants[7] = 0.6
    constants[8] = 0.1
    constants[9] = 1.0
    constants[10] = 10.0
    states[3] = 0.0
    states[4] = 0.0
    states[5] = 0.0
    constants[11] = 0.1
    constants[12] = 0.8
    constants[13] = 1.1
    constants[14] = 0.3
    constants[15] = 0.1
    constants[16] = 1.0
    constants[17] = 0.1
    constants[18] = 0.001
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[2] = (constants[7]*(constants[10]-(states[2]+states[4]))*states[0]+constants[6]*(constants[10]-(states[2]+states[4]))*states[3])-(constants[8]*states[1]*states[2]+constants[9]*states[2])
    rates[5] = (constants[11]*states[0]+constants[12]*states[1])-constants[13]*states[5]
    rates[3] = (constants[14]*states[2]+constants[15]*states[5])-constants[16]*states[3]
    rates[4] = constants[17]*states[1]*states[2]-constants[18]*states[4]
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 5.00000), 0.00000 , greater_equal(voi , 5.00000) & less_equal(voi , 10.0000), 1.00000 , True, 0.00000])
    rates[0] = constants[0]*(algebraic[0]/(constants[1]+algebraic[0]))-constants[2]*states[0]
    algebraic[1] = custom_piecewise([greater_equal(voi , 0.00000) & less_equal(voi , 5.00000), 1.00000 , True, 0.00000])
    rates[1] = constants[3]*(algebraic[1]/(constants[4]+algebraic[1]))-constants[5]*states[1]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 5.00000), 0.00000 , greater_equal(voi , 5.00000) & less_equal(voi , 10.0000), 1.00000 , True, 0.00000])
    algebraic[1] = custom_piecewise([greater_equal(voi , 0.00000) & less_equal(voi , 5.00000), 1.00000 , True, 0.00000])
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