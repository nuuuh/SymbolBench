# Size of variable arrays:
sizeAlgebraic = 5
sizeStates = 2
sizeConstants = 17
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "x1 in component x1 (picomole)"
    legend_constants[14] = "k in component k (flux)"
    legend_algebraic[4] = "lambda_Ca in component lambda_Ca (per_minute)"
    legend_constants[0] = "lambda1 in component model_parameters (per_minute)"
    legend_constants[15] = "A in component A (per_minute)"
    legend_constants[16] = "B in component B (per_minute)"
    legend_algebraic[0] = "Ca in component Ca (millimolar)"
    legend_algebraic[3] = "S in component S (millimolar)"
    legend_algebraic[2] = "m_Ca in component m_Ca (dimensionless)"
    legend_constants[1] = "m1 in component model_parameters (dimensionless)"
    legend_constants[2] = "m2 in component model_parameters (dimensionless)"
    legend_constants[3] = "beta in component model_parameters (litre_per_millimole)"
    legend_constants[4] = "R in component model_parameters (millimolar)"
    legend_states[1] = "x2 in component x2 (picomole)"
    legend_algebraic[1] = "PTH in component x2 (picomole)"
    legend_constants[5] = "lambda2 in component model_parameters (per_minute)"
    legend_constants[6] = "Ca_0 in component model_parameters (millimolar)"
    legend_constants[7] = "Ca_1 in component model_parameters (millimolar)"
    legend_constants[8] = "alpha in component model_parameters (per_minute)"
    legend_constants[9] = "t0 in component model_parameters (minute)"
    legend_constants[10] = "x1_n in component model_parameters (picomole)"
    legend_constants[11] = "x2_n in component model_parameters (picomole)"
    legend_constants[12] = "x2_max in component model_parameters (picomole)"
    legend_constants[13] = "x2_min in component model_parameters (picomole)"
    legend_rates[0] = "d/dt x1 in component x1 (picomole)"
    legend_rates[1] = "d/dt x2 in component x2 (picomole)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.00
    constants[0] = 0.0125
    constants[1] = 112.5200
    constants[2] = 15.00
    constants[3] = 1e6
    constants[4] = 1.2162
    states[1] = 0.00
    constants[5] = 0.5595
    constants[6] = 1.255
    constants[7] = 0.1817
    constants[8] = 0.0442
    constants[9] = 575.0
    constants[10] = 490.7800
    constants[11] = 6.6290
    constants[12] = 14.0430
    constants[13] = 0.6697
    constants[14] = constants[5]*constants[11]+constants[0]*constants[10]
    constants[15] = (constants[0]*constants[5]*constants[12])/(constants[14]-constants[5]*constants[12])
    constants[16] = (constants[0]*constants[5]*constants[13])/(constants[14]-constants[5]*constants[13])
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = custom_piecewise([less(voi , constants[9]), constants[6] , True, constants[6]-constants[7]*(1.00000-exp(-constants[8]*(voi-constants[9])))])
    algebraic[2] = constants[1]/(1.00000+exp(-constants[3]*(constants[4]-algebraic[0])))+constants[2]
    algebraic[3] = constants[6]*(power(-((constants[10]*constants[16]-constants[5]*constants[11])/(constants[10]*constants[15]-constants[5]*constants[11])), 1.00000/algebraic[2]))
    algebraic[4] = (constants[15]-constants[16])/(1.00000+power(algebraic[0]/algebraic[3], algebraic[2]))+constants[16]
    rates[0] = constants[14]-(algebraic[4]*states[0]+constants[0]*states[0])
    rates[1] = algebraic[4]*states[0]-constants[5]*states[1]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less(voi , constants[9]), constants[6] , True, constants[6]-constants[7]*(1.00000-exp(-constants[8]*(voi-constants[9])))])
    algebraic[2] = constants[1]/(1.00000+exp(-constants[3]*(constants[4]-algebraic[0])))+constants[2]
    algebraic[3] = constants[6]*(power(-((constants[10]*constants[16]-constants[5]*constants[11])/(constants[10]*constants[15]-constants[5]*constants[11])), 1.00000/algebraic[2]))
    algebraic[4] = (constants[15]-constants[16])/(1.00000+power(algebraic[0]/algebraic[3], algebraic[2]))+constants[16]
    algebraic[1] = states[1]/2.75000
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