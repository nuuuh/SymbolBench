# Size of variable arrays:
sizeAlgebraic = 13
sizeStates = 14
sizeConstants = 18
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "O2 in component O2 (molar)"
    legend_algebraic[0] = "v1 in component v1 (flux)"
    legend_algebraic[4] = "v5 in component v5 (flux)"
    legend_algebraic[6] = "v7 in component v7 (flux)"
    legend_algebraic[10] = "v11 in component v11 (flux)"
    legend_algebraic[12] = "v13_back in component v13_back (flux)"
    legend_constants[17] = "v13 in component v13 (flux)"
    legend_states[1] = "O2_radical in component O2_radical (molar)"
    legend_algebraic[5] = "v6 in component v6 (flux)"
    legend_states[2] = "NADH in component NADH (molar)"
    legend_constants[16] = "v12 in component v12 (flux)"
    legend_algebraic[11] = "v14 in component v14 (flux)"
    legend_states[3] = "NAD in component NAD (molar)"
    legend_algebraic[7] = "v8 in component v8 (flux)"
    legend_algebraic[9] = "v10 in component v10 (flux)"
    legend_states[4] = "NAD_radical in component NAD_radical (molar)"
    legend_algebraic[8] = "v9 in component v9 (flux)"
    legend_states[5] = "NAD2 in component NAD2 (molar)"
    legend_states[6] = "H2O2 in component H2O2 (molar)"
    legend_algebraic[1] = "v2 in component v2 (flux)"
    legend_states[7] = "Per3 in component Per3 (molar)"
    legend_algebraic[3] = "v4 in component v4 (flux)"
    legend_states[8] = "Per2 in component Per2 (molar)"
    legend_states[9] = "coI in component coI (molar)"
    legend_algebraic[2] = "v3 in component v3 (flux)"
    legend_states[10] = "coII in component coII (molar)"
    legend_states[11] = "coIII in component coIII (molar)"
    legend_states[12] = "Ar_radical in component Ar_radical (molar)"
    legend_states[13] = "ArH in component ArH (molar)"
    legend_constants[0] = "k1 in component v1 (second_order_rate_constant)"
    legend_constants[1] = "k2 in component v2 (second_order_rate_constant)"
    legend_constants[2] = "k3 in component v3 (second_order_rate_constant)"
    legend_constants[3] = "k4 in component v4 (second_order_rate_constant)"
    legend_constants[4] = "k5 in component v5 (second_order_rate_constant)"
    legend_constants[5] = "k6 in component v6 (second_order_rate_constant)"
    legend_constants[6] = "k7 in component v7 (second_order_rate_constant)"
    legend_constants[7] = "k8 in component v8 (second_order_rate_constant)"
    legend_constants[8] = "k9 in component v9 (second_order_rate_constant)"
    legend_constants[9] = "k10 in component v10 (second_order_rate_constant)"
    legend_constants[10] = "k11 in component v11 (second_order_rate_constant)"
    legend_constants[11] = "k12 in component v12 (flux)"
    legend_constants[12] = "O2eq in component v13 (molar)"
    legend_constants[13] = "k13 in component v13 (first_order_rate_constant)"
    legend_constants[14] = "k13_ in component v13_back (first_order_rate_constant)"
    legend_constants[15] = "k14 in component v14 (second_order_rate_constant)"
    legend_rates[0] = "d/dt O2 in component O2 (molar)"
    legend_rates[1] = "d/dt O2_radical in component O2_radical (molar)"
    legend_rates[2] = "d/dt NADH in component NADH (molar)"
    legend_rates[3] = "d/dt NAD in component NAD (molar)"
    legend_rates[4] = "d/dt NAD_radical in component NAD_radical (molar)"
    legend_rates[5] = "d/dt NAD2 in component NAD2 (molar)"
    legend_rates[6] = "d/dt H2O2 in component H2O2 (molar)"
    legend_rates[7] = "d/dt Per3 in component Per3 (molar)"
    legend_rates[8] = "d/dt Per2 in component Per2 (molar)"
    legend_rates[9] = "d/dt coI in component coI (molar)"
    legend_rates[10] = "d/dt coII in component coII (molar)"
    legend_rates[11] = "d/dt coIII in component coIII (molar)"
    legend_rates[12] = "d/dt Ar_radical in component Ar_radical (molar)"
    legend_rates[13] = "d/dt ArH in component ArH (molar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 5.33e-6
    states[1] = 0.026e-6
    states[2] = 0.95e-6
    states[3] = 558e-6
    states[4] = 5e-10
    states[5] = 0.12e-6
    states[6] = 0.013e-6
    states[7] = 0.059e-6
    states[8] = 1e-10
    states[9] = 8.9e-10
    states[10] = 0.026e-6
    states[11] = 1.31e-6
    states[12] = 0.12e-6
    states[13] = 299.88e-6
    constants[0] = 3
    constants[1] = 1.8e7
    constants[2] = 1.5e5
    constants[3] = 5.2e3
    constants[4] = 2e7
    constants[5] = 1.7e7
    constants[6] = 2e7
    constants[7] = 4e7
    constants[8] = 6e7
    constants[9] = 1.8e6
    constants[10] = 1e5
    constants[11] = 0.08e-6
    constants[12] = 1.2e-5
    constants[13] = 6e-3
    constants[14] = 6e-3
    constants[15] = 7e5
    constants[16] = constants[11]
    constants[17] = constants[13]*constants[12]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[3] = constants[3]*states[10]*states[13]
    algebraic[2] = constants[2]*states[9]*states[13]
    rates[10] = algebraic[2]-algebraic[3]
    algebraic[4] = constants[4]*states[4]*states[0]
    algebraic[6] = constants[6]*(power(states[1], 2.00000))
    algebraic[5] = constants[5]*states[1]*states[7]
    rates[1] = (algebraic[4]-algebraic[5])-2.00000*algebraic[6]
    algebraic[0] = constants[0]*states[2]*states[0]
    algebraic[1] = constants[1]*states[6]*states[7]
    rates[6] = (algebraic[0]-algebraic[1])-algebraic[6]
    algebraic[7] = constants[7]*states[11]*states[4]
    rates[9] = (algebraic[1]-algebraic[2])+algebraic[7]
    algebraic[8] = constants[8]*(power(states[4], 2.00000))
    rates[5] = algebraic[8]
    algebraic[9] = constants[9]*states[7]*states[4]
    rates[3] = algebraic[0]+algebraic[4]+algebraic[7]+algebraic[9]
    rates[7] = ((-algebraic[1]+algebraic[3])-algebraic[5])-algebraic[9]
    algebraic[11] = constants[15]*states[12]*states[2]
    rates[2] = (-algebraic[0]+constants[16])-algebraic[11]
    rates[4] = (((-2.00000*algebraic[8]-algebraic[4])-algebraic[7])-algebraic[9])+algebraic[11]
    algebraic[10] = constants[10]*states[8]*states[0]
    rates[8] = -algebraic[10]+algebraic[9]
    rates[11] = (algebraic[5]-algebraic[7])+algebraic[10]
    rates[12] = (algebraic[2]+algebraic[3])-algebraic[11]
    rates[13] = (-algebraic[2]-algebraic[3])+algebraic[11]
    algebraic[12] = constants[14]*states[0]
    rates[0] = ((((-algebraic[0]-algebraic[4])+algebraic[6])-algebraic[10])+constants[17])-algebraic[12]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = constants[3]*states[10]*states[13]
    algebraic[2] = constants[2]*states[9]*states[13]
    algebraic[4] = constants[4]*states[4]*states[0]
    algebraic[6] = constants[6]*(power(states[1], 2.00000))
    algebraic[5] = constants[5]*states[1]*states[7]
    algebraic[0] = constants[0]*states[2]*states[0]
    algebraic[1] = constants[1]*states[6]*states[7]
    algebraic[7] = constants[7]*states[11]*states[4]
    algebraic[8] = constants[8]*(power(states[4], 2.00000))
    algebraic[9] = constants[9]*states[7]*states[4]
    algebraic[11] = constants[15]*states[12]*states[2]
    algebraic[10] = constants[10]*states[8]*states[0]
    algebraic[12] = constants[14]*states[0]
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