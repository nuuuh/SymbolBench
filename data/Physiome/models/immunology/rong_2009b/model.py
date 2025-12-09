# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 5
sizeConstants = 15
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
    legend_constants[4] = "eta in component latently_infected (dimensionless)"
    legend_constants[5] = "d_0 in component latently_infected (per_day)"
    legend_algebraic[0] = "f in component activated (dimensionless)"
    legend_constants[6] = "a in component latently_infected (per_day)"
    legend_constants[7] = "rho in component latently_infected (per_day)"
    legend_states[2] = "L_a in component activated (per_ml)"
    legend_states[3] = "L_0 in component latently_infected (per_ml)"
    legend_constants[8] = "p in component activated (per_day)"
    legend_constants[9] = "a_L in component activated (per_day)"
    legend_constants[10] = "T_off in component activated (day)"
    legend_constants[11] = "T_on in component activated (day)"
    legend_constants[12] = "delta in component productively_infected (per_day)"
    legend_states[4] = "T_star in component productively_infected (per_ml)"
    legend_constants[13] = "p_v in component viral_load (per_day)"
    legend_constants[14] = "c in component viral_load (per_day)"
    legend_rates[1] = "d/dt T in component uninfected (per_ml)"
    legend_rates[3] = "d/dt L_0 in component latently_infected (per_ml)"
    legend_rates[2] = "d/dt L_a in component activated (per_ml)"
    legend_rates[4] = "d/dt T_star in component productively_infected (per_ml)"
    legend_rates[0] = "d/dt V in component viral_load (per_ml)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1E4
    constants[1] = 0.01
    constants[2] = 0.85
    constants[3] = 2.4E-8
    states[0] = 50
    states[1] = 600000
    constants[4] = 0.001
    constants[5] = 0.001
    constants[6] = 0.03
    constants[7] = 0.01
    states[2] = 0
    states[3] = 2
    constants[8] = 1.4
    constants[9] = 0.1
    constants[10] = 54
    constants[11] = 50
    constants[12] = 1
    states[4] = 0.3
    constants[13] = 2000
    constants[14] = 23
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = (constants[0]-constants[1]*states[1])-(1.00000-constants[2])*constants[3]*states[0]*states[1]
    rates[4] = ((1.00000-constants[4])*(1.00000-constants[2])*constants[3]*states[0]*states[1]-1.00000*constants[12]*states[4]*states[4])+constants[9]*states[2]
    rates[0] = constants[13]*states[4]-constants[14]*states[0]
    algebraic[0] = custom_piecewise([less(constants[10] , voi), 0.00000 , less(voi , constants[11]), 0.00000 , True, 1.00000])
    rates[3] = ((constants[4]*(1.00000-constants[2])*constants[3]*states[0]*states[1]-constants[5]*states[3])-algebraic[0]*constants[6]*states[3])+(1.00000-algebraic[0])*constants[7]*states[2]
    rates[2] = (algebraic[0]*(constants[6]*states[3]+constants[8]*states[2])-(1.00000-algebraic[0])*(constants[6]+constants[7])*states[2])-constants[9]*states[2]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less(constants[10] , voi), 0.00000 , less(voi , constants[11]), 0.00000 , True, 1.00000])
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