# Size of variable arrays:
sizeAlgebraic = 3
sizeStates = 2
sizeConstants = 14
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (day)"
    legend_constants[0] = "s in component uninfected (per_day_mm3)"
    legend_constants[1] = "p in component uninfected (per_day)"
    legend_constants[2] = "gamma in component uninfected (per_day)"
    legend_constants[13] = "beta in component uninfected (dimensionless)"
    legend_constants[3] = "N in component free_virus_particle (dimensionless)"
    legend_constants[4] = "k_1 in component latently_infected (mm3_per_day)"
    legend_constants[5] = "k_2 in component actively_infected (per_day)"
    legend_constants[6] = "k_3 in component latently_infected (per_day)"
    legend_constants[7] = "mu_V in component free_virus_particle (per_day)"
    legend_states[0] = "T_1 in component latently_infected (per_mm3)"
    legend_constants[8] = "mu_b in component actively_infected (per_day)"
    legend_states[1] = "T in component uninfected (per_mm3)"
    legend_constants[9] = "k_4 in component latently_infected (per_day)"
    legend_constants[10] = "T_0 in component latently_infected (per_mm3)"
    legend_constants[11] = "V_0 in component latently_infected (per_mm3)"
    legend_constants[12] = "t_min in component latently_infected (day)"
    legend_algebraic[0] = "T_1_t in component latently_infected (per_mm3)"
    legend_algebraic[1] = "T_2 in component actively_infected (per_mm3)"
    legend_algebraic[2] = "V in component free_virus_particle (per_mm3)"
    legend_rates[1] = "d/dt T in component uninfected (per_mm3)"
    legend_rates[0] = "d/dt T_1 in component latently_infected (per_mm3)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 10
    constants[1] = 0.01
    constants[2] = 2E-5
    constants[3] = 1000
    constants[4] = 2.4E-5
    constants[5] = 3E-3
    constants[6] = 0.023
    constants[7] = 2.4
    states[0] = 0
    constants[8] = 0.24
    states[1] = 1000
    constants[9] = 2.424
    constants[10] = 1000
    constants[11] = 1E-3
    constants[12] = 2
    constants[13] = (constants[2]/constants[6])*(1.00000+constants[5]/constants[8])
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = ((constants[0]+constants[1]*states[1])-constants[2]*(power(states[1], 2.00000)))-(constants[6]*constants[13]+(constants[3]*constants[4]*constants[5])/(constants[4]*states[1]+constants[7]))*states[1]*states[0]
    algebraic[0] = ((constants[4]*constants[10]*constants[11])/(constants[9]-constants[6]))*(exp(-constants[6]*voi)-exp(-constants[9]*voi))
    rates[0] = custom_piecewise([less_equal(voi , constants[12]), algebraic[0] , True, ((constants[3]*constants[4]*constants[5])/(constants[4]*states[1]+constants[7]))*states[1]*states[0]-constants[6]*states[0]])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = ((constants[4]*constants[10]*constants[11])/(constants[9]-constants[6]))*(exp(-constants[6]*voi)-exp(-constants[9]*voi))
    algebraic[1] = (constants[5]*states[0])/constants[8]
    algebraic[2] = (constants[3]*constants[5]*states[0])/(constants[4]*states[1]+constants[7])
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