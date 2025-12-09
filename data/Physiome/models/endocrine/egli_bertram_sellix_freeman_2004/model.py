# Size of variable arrays:
sizeAlgebraic = 3
sizeStates = 3
sizeConstants = 11
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_states[0] = "PRL in component PRL (nanog_ml)"
    legend_constants[0] = "kD in component PRL (nanog_ml)"
    legend_constants[1] = "kO in component PRL (picog_ml)"
    legend_constants[2] = "rP in component PRL (nanog_ml_hr)"
    legend_constants[3] = "qP in component PRL (first_order_rate_constant)"
    legend_states[1] = "OT in component OT (picog_ml)"
    legend_states[2] = "DA in component DA (nanog_ml)"
    legend_algebraic[0] = "vD in component DA (nanog_ml_hr)"
    legend_constants[4] = "vDbar in component DA (nanog_ml_hr)"
    legend_constants[5] = "DA_infinity in component DA (nanog_ml)"
    legend_constants[6] = "qD in component DA (first_order_rate_constant)"
    legend_constants[7] = "kx in component OT (picog_ml)"
    legend_algebraic[1] = "vO in component OT (picog_ml_hr)"
    legend_constants[8] = "vObar in component OT (picog_ml_hr)"
    legend_constants[9] = "rO in component OT (picog_ml_hr)"
    legend_constants[10] = "qO in component OT (first_order_rate_constant)"
    legend_algebraic[2] = "x in component x (picog_ml)"
    legend_rates[0] = "d/dt PRL in component PRL (nanog_ml)"
    legend_rates[2] = "d/dt DA in component DA (nanog_ml)"
    legend_rates[1] = "d/dt OT in component OT (picog_ml)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 20.0
    constants[0] = 300.0
    constants[1] = 9.0
    constants[2] = 300000.0
    constants[3] = 0.5
    states[1] = 25.0
    states[2] = 20000.0
    constants[4] = 10000.0
    constants[5] = 20000.0
    constants[6] = 0.2
    constants[7] = 50.0
    constants[8] = 500.0
    constants[9] = 1100.0
    constants[10] = 1.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[2]*(1.00000/(constants[0]+states[2]))*(power(states[1]/(constants[1]+states[1]), 2.00000))-constants[3]*states[0]
    algebraic[0] = custom_piecewise([greater_equal(voi , 2.00000) & less_equal(voi , 4.00000), constants[4] , True, 0.00000])
    rates[2] = constants[6]*(constants[5]-states[2])-algebraic[0]
    algebraic[1] = custom_piecewise([greater_equal(voi , 2.00000) & less_equal(voi , 4.00000), constants[8] , True, 0.00000])
    algebraic[2] = custom_piecewise([greater_equal(voi , 2.00000) & less_equal(voi , 4.00000), 51.0000 , greater_equal(voi , 16.0000) & less_equal(voi , 18.0000), 51.0000 , True, 1.00000])
    rates[1] = constants[9]*(algebraic[2]/(constants[7]+algebraic[2]))-(constants[10]*states[1]+algebraic[1])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 2.00000) & less_equal(voi , 4.00000), constants[4] , True, 0.00000])
    algebraic[1] = custom_piecewise([greater_equal(voi , 2.00000) & less_equal(voi , 4.00000), constants[8] , True, 0.00000])
    algebraic[2] = custom_piecewise([greater_equal(voi , 2.00000) & less_equal(voi , 4.00000), 51.0000 , greater_equal(voi , 16.0000) & less_equal(voi , 18.0000), 51.0000 , True, 1.00000])
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