# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 7
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
    legend_states[0] = "Ca_t in component equations (uM_per_kg)"
    legend_states[1] = "TnCa_t in component equations (uM_per_kg)"
    legend_states[2] = "CB_on_t in component equations (uM_per_kg)"
    legend_states[3] = "Ca_released in component equations (uM_per_kg)"
    legend_states[4] = "Ca_sequestered in component equations (uM_per_kg)"
    legend_states[5] = "cumCB_on_t in component equations (uM_per_kg)"
    legend_states[6] = "cumCB_off_t in component equations (uM_per_kg)"
    legend_algebraic[0] = "Ca_release_rate in component equations (uM_per_kg_per_second)"
    legend_algebraic[1] = "dTnCa_t_dt in component equations (uM_per_kg_per_second)"
    legend_constants[0] = "Ca_tot_released in component equations (uM_per_kg)"
    legend_constants[1] = "total_Tn in component equations (uM_per_kg)"
    legend_constants[2] = "total_CB in component equations (uM_per_kg)"
    legend_constants[3] = "k_1 in component equations (kg_per_uM_per_second)"
    legend_constants[4] = "k_2 in component equations (per_second)"
    legend_constants[5] = "k_3 in component equations (per_second)"
    legend_constants[6] = "f in component equations (kg_per_uM_per_second)"
    legend_constants[7] = "g in component equations (per_second)"
    legend_rates[0] = "d/dt Ca_t in component equations (uM_per_kg)"
    legend_rates[1] = "d/dt TnCa_t in component equations (uM_per_kg)"
    legend_rates[2] = "d/dt CB_on_t in component equations (uM_per_kg)"
    legend_rates[3] = "d/dt Ca_released in component equations (uM_per_kg)"
    legend_rates[4] = "d/dt Ca_sequestered in component equations (uM_per_kg)"
    legend_rates[5] = "d/dt cumCB_on_t in component equations (uM_per_kg)"
    legend_rates[6] = "d/dt cumCB_off_t in component equations (uM_per_kg)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0
    states[1] = 0
    states[2] = 0
    states[3] = 0
    states[4] = 0
    states[5] = 0
    states[6] = 0
    constants[0] = 35
    constants[1] = 70
    constants[2] = 150
    constants[3] = 5e6
    constants[4] = 10
    constants[5] = 1000
    constants[6] = 0.4e6
    constants[7] = 10
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = constants[3]*states[0]*(constants[1]-states[1])-constants[4]*states[1]
    rates[2] = constants[6]*states[1]*(constants[2]-states[2])-constants[7]*states[2]
    rates[4] = constants[5]*states[0]
    rates[5] = constants[6]*states[1]*(constants[2]-states[2])
    rates[6] = constants[7]*states[2]
    algebraic[0] = custom_piecewise([greater(voi , 0.100000), 0.00000 , True, 20.0000*constants[0]*(1.00000-10.0000*voi)])
    rates[3] = algebraic[0]
    algebraic[1] = constants[3]*states[0]*(constants[1]-states[1])-constants[4]*states[1]
    rates[0] = (algebraic[0]-constants[5]*states[0])-algebraic[1]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater(voi , 0.100000), 0.00000 , True, 20.0000*constants[0]*(1.00000-10.0000*voi)])
    algebraic[1] = constants[3]*states[0]*(constants[1]-states[1])-constants[4]*states[1]
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