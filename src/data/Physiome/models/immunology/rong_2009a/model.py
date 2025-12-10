# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 4
sizeConstants = 10
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (day)"
    legend_constants[0] = "lambda in component uninfected (per_ml_day)"
    legend_constants[1] = "d_T in component uninfected (per_day)"
    legend_constants[2] = "efficacy in component drug_efficacy (dimensionless)"
    legend_constants[3] = "k in component uninfected (ml_per_day)"
    legend_states[0] = "V in component viral_load (per_ml)"
    legend_states[1] = "T in component uninfected (per_ml)"
    legend_constants[4] = "d_0 in component latently_infected (per_day)"
    legend_constants[5] = "a_L in component latently_infected (per_day)"
    legend_constants[6] = "eta in component latently_infected (dimensionless)"
    legend_states[2] = "L in component latently_infected (per_ml)"
    legend_constants[7] = "delta in component productively_infected (per_day)"
    legend_states[3] = "T_star in component productively_infected (per_ml)"
    legend_constants[8] = "N in component viral_load (dimensionless)"
    legend_constants[9] = "c in component viral_load (per_day)"
    legend_rates[1] = "d/dt T in component uninfected (per_ml)"
    legend_rates[2] = "d/dt L in component latently_infected (per_ml)"
    legend_rates[3] = "d/dt T_star in component productively_infected (per_ml)"
    legend_rates[0] = "d/dt V in component viral_load (per_ml)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1E4
    constants[1] = 0.01
    constants[2] = 0.4
    constants[3] = 2.4E-8
    states[0] = 50
    states[1] = 600000
    constants[4] = 0.001
    constants[5] = 0.1
    constants[6] = 0.001
    states[2] = 2
    constants[7] = 1
    states[3] = 0.3
    constants[8] = 2000
    constants[9] = 23
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = (constants[0]-constants[1]*states[1])-(1.00000-constants[2])*constants[3]*states[0]*states[1]
    rates[2] = (constants[6]*(1.00000-constants[2])*constants[3]*states[0]*states[1]-constants[4]*states[2])-constants[5]*states[2]
    rates[3] = ((1.00000-constants[6])*(1.00000-constants[2])*constants[3]*states[0]*states[1]-constants[7]*states[3])+constants[5]*states[2]
    rates[0] = constants[8]*constants[7]*states[3]-constants[9]*states[0]
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