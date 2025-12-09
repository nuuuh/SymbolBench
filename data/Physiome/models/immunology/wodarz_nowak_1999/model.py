# Size of variable arrays:
sizeAlgebraic = 3
sizeStates = 4
sizeConstants = 9
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (day)"
    legend_states[0] = "x in component x (per_mm3)"
    legend_constants[0] = "lamda in component x (per_mm3_per_day)"
    legend_constants[1] = "d in component x (per_day)"
    legend_constants[2] = "beta in component kinetic_parameters (mm3_per_day)"
    legend_algebraic[0] = "s in component kinetic_parameters (dimensionless)"
    legend_states[1] = "y in component y (per_mm3)"
    legend_constants[3] = "a in component y (per_day)"
    legend_constants[4] = "p in component kinetic_parameters (mm3_per_day)"
    legend_states[2] = "z in component z (per_mm3)"
    legend_algebraic[1] = "log_y in component y (dimensionless)"
    legend_states[3] = "w in component w (per_mm3)"
    legend_constants[5] = "b in component w (per_day)"
    legend_constants[6] = "c in component kinetic_parameters (mm3_mm3_per_day)"
    legend_constants[7] = "q in component kinetic_parameters (per_mm3)"
    legend_algebraic[2] = "log_w in component w (dimensionless)"
    legend_constants[8] = "h in component z (per_day)"
    legend_rates[0] = "d/dt x in component x (per_mm3)"
    legend_rates[1] = "d/dt y in component y (per_mm3)"
    legend_rates[3] = "d/dt w in component w (per_mm3)"
    legend_rates[2] = "d/dt z in component z (per_mm3)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 10
    constants[0] = 1
    constants[1] = 0.1
    constants[2] = 0.5
    states[1] = 0.1
    constants[3] = 0.2
    constants[4] = 1
    states[2] = 0
    states[3] = 0.001
    constants[5] = 0.01
    constants[6] = 0.1
    constants[7] = 0.5
    constants[8] = 0.1
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[3] = constants[6]*states[0]*states[1]*states[3]-(constants[6]*constants[7]*states[1]*states[3]+constants[5]*states[3])
    rates[2] = constants[6]*constants[7]*states[1]*states[3]-constants[8]*states[2]
    algebraic[0] = custom_piecewise([less_equal(voi , 15.0000), 1.00000 , greater_equal(voi , 40.0000), 1.00000 , True, 0.00420000])
    rates[0] = constants[0]-(constants[1]*states[0]+algebraic[0]*constants[2]*states[0]*states[1])
    rates[1] = algebraic[0]*constants[2]*states[0]*states[1]-(constants[3]*states[1]+constants[4]*states[1]*states[2])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less_equal(voi , 15.0000), 1.00000 , greater_equal(voi , 40.0000), 1.00000 , True, 0.00420000])
    algebraic[1] = log(states[1]*1.00000, 10)
    algebraic[2] = log(states[3]*1.00000, 10)
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