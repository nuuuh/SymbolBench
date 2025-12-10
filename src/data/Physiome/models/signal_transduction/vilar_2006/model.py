# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 6
sizeConstants = 8
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_algebraic[0] = "l in component l (dimensionless)"
    legend_states[0] = "l_RI_RII in component l_RI_RII (dimensionless)"
    legend_states[1] = "RI in component RI (dimensionless)"
    legend_states[2] = "RII in component RII (dimensionless)"
    legend_constants[0] = "ka in component model_parameters (per_minute)"
    legend_constants[1] = "ki in component model_parameters (per_minute)"
    legend_constants[2] = "kcd in component model_parameters (per_minute)"
    legend_constants[3] = "klid in component model_parameters (per_minute)"
    legend_states[3] = "RI_endo in component RI_endo (dimensionless)"
    legend_states[4] = "l_RI_RII_endo in component l_RI_RII_endo (dimensionless)"
    legend_constants[4] = "kr in component model_parameters (per_minute)"
    legend_constants[5] = "p_RI in component model_parameters (per_minute)"
    legend_constants[6] = "alpha in component model_parameters (dimensionless)"
    legend_states[5] = "RII_endo in component RII_endo (dimensionless)"
    legend_constants[7] = "p_RII in component model_parameters (per_minute)"
    legend_rates[0] = "d/dt l_RI_RII in component l_RI_RII (dimensionless)"
    legend_rates[1] = "d/dt RI in component RI (dimensionless)"
    legend_rates[2] = "d/dt RII in component RII (dimensionless)"
    legend_rates[4] = "d/dt l_RI_RII_endo in component l_RI_RII_endo (dimensionless)"
    legend_rates[3] = "d/dt RI_endo in component RI_endo (dimensionless)"
    legend_rates[5] = "d/dt RII_endo in component RII_endo (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.0
    states[1] = 20.0
    states[2] = 20.0
    constants[0] = 1
    constants[1] = 0.333333333
    constants[2] = 0.0277777778
    constants[3] = 0.25
    states[3] = 0.0
    states[4] = 40.0
    constants[4] = 0.033333333
    constants[5] = 8
    constants[6] = 1
    states[5] = 0.0
    constants[7] = 4
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[4] = constants[1]*states[0]-constants[4]*states[4]
    rates[3] = constants[1]*states[1]-constants[4]*states[3]
    rates[5] = constants[1]*states[2]-constants[4]*states[5]
    algebraic[0] = custom_piecewise([greater_equal(voi , 2500.00), 0.0100000 , True, 3.00000e-05])
    rates[0] = constants[0]*algebraic[0]*states[1]*states[2]-(constants[2]+constants[3]+constants[1])*states[0]
    rates[1] = (constants[5]+constants[4]*states[3]+constants[6]*constants[4]*states[4])-(constants[0]*algebraic[0]*states[1]*states[2]+(constants[2]+constants[1])*states[1])
    rates[2] = (constants[7]+constants[4]*states[5]+constants[6]*constants[4]*states[4])-(constants[0]*algebraic[0]*states[1]*states[2]+(constants[2]+constants[1])*states[2])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 2500.00), 0.0100000 , True, 3.00000e-05])
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