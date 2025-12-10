# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 5
sizeConstants = 36
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "x1 in component x1 (microg_l)"
    legend_constants[0] = "lambda_1 in component model_parameters (first_order_rate_constant)"
    legend_states[1] = "x3 in component x3 (microg_l)"
    legend_constants[1] = "a1 in component model_parameters (flux)"
    legend_constants[2] = "a2 in component model_parameters (flux)"
    legend_constants[3] = "a3 in component model_parameters (first_order_rate_constant)"
    legend_constants[4] = "a4 in component model_parameters (second_order_rate_constant)"
    legend_constants[5] = "a5 in component model_parameters (per_microg_l)"
    legend_constants[6] = "a6 in component model_parameters (per_microg_l2)"
    legend_constants[7] = "a7 in component model_parameters (per_microg_l)"
    legend_constants[8] = "a8 in component model_parameters (per_microg_l2)"
    legend_states[2] = "x2 in component x2 (microg_l)"
    legend_constants[9] = "lambda_2 in component model_parameters (first_order_rate_constant)"
    legend_constants[10] = "a9 in component model_parameters (flux)"
    legend_constants[11] = "a10 in component model_parameters (first_order_rate_constant)"
    legend_constants[12] = "a11 in component model_parameters (second_order_rate_constant)"
    legend_constants[13] = "a12 in component model_parameters (per_microg_l)"
    legend_constants[14] = "a13 in component model_parameters (per_microg_l2)"
    legend_constants[15] = "a14 in component model_parameters (per_microg_l)"
    legend_constants[16] = "a15 in component model_parameters (per_microg_l2)"
    legend_constants[33] = "lambda_3_ in component model_parameters (first_order_rate_constant)"
    legend_states[3] = "x4 in component x4 (microg_l)"
    legend_states[4] = "x5 in component x5 (microg_l)"
    legend_constants[17] = "a16 in component model_parameters (flux)"
    legend_constants[18] = "a17 in component model_parameters (first_order_rate_constant)"
    legend_constants[19] = "a18 in component model_parameters (second_order_rate_constant)"
    legend_constants[20] = "a19 in component model_parameters (first_order_rate_constant)"
    legend_constants[21] = "a20 in component model_parameters (second_order_rate_constant)"
    legend_constants[22] = "a21 in component model_parameters (per_microg_l)"
    legend_constants[23] = "a22 in component model_parameters (per_microg_l2)"
    legend_constants[24] = "a23 in component model_parameters (per_microg_l)"
    legend_constants[25] = "a24 in component model_parameters (per_microg_l2)"
    legend_constants[26] = "a25 in component model_parameters (first_order_rate_constant)"
    legend_constants[27] = "a26 in component model_parameters (first_order_rate_constant)"
    legend_constants[34] = "lambda_4_ in component model_parameters (first_order_rate_constant)"
    legend_constants[28] = "a27 in component model_parameters (first_order_rate_constant)"
    legend_constants[35] = "lambda_5_ in component model_parameters (first_order_rate_constant)"
    legend_constants[29] = "a28 in component model_parameters (first_order_rate_constant)"
    legend_constants[30] = "lambda_3 in component model_parameters (first_order_rate_constant)"
    legend_constants[31] = "lambda_4 in component model_parameters (first_order_rate_constant)"
    legend_constants[32] = "lambda_5 in component model_parameters (first_order_rate_constant)"
    legend_rates[0] = "d/dt x1 in component x1 (microg_l)"
    legend_rates[2] = "d/dt x2 in component x2 (microg_l)"
    legend_rates[1] = "d/dt x3 in component x3 (microg_l)"
    legend_rates[3] = "d/dt x4 in component x4 (microg_l)"
    legend_rates[4] = "d/dt x5 in component x5 (microg_l)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.01067
    constants[0] = 0.059
    states[1] = 6.51
    constants[1] = 0.000017
    constants[2] = 0.0023
    constants[3] = 0.6
    constants[4] = 45
    constants[5] = 36
    constants[6] = 216
    constants[7] = 0.28
    constants[8] = 0.36
    states[2] = 0.04665
    constants[9] = 0.028
    constants[10] = 0.0003
    constants[11] = 0.18
    constants[12] = 150
    constants[13] = 18
    constants[14] = 460
    constants[15] = 0.46
    constants[16] = 0.1
    states[3] = 60.61
    states[4] = 12.61
    constants[17] = 0.04
    constants[18] = 150
    constants[19] = 3800
    constants[20] = 57
    constants[21] = 2600
    constants[22] = 200
    constants[23] = 9400
    constants[24] = 10
    constants[25] = 320
    constants[26] = 0.04
    constants[27] = 0.00097
    constants[28] = 0.57
    constants[29] = 0.0017
    constants[30] = 0.0986
    constants[31] = 0.024
    constants[32] = 3e-5
    constants[33] = constants[30]+constants[28]+constants[29]
    constants[34] = constants[31]+constants[26]
    constants[35] = constants[32]+constants[27]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (constants[1]+(constants[2]+constants[3]*states[0]+constants[4]*(power(states[0], 2.00000)))/(1.00000+constants[5]*states[0]+constants[6]*(power(states[0], 2.00000))+constants[7]*states[1]+constants[8]*(power(states[1], 2.00000))))-constants[0]*states[0]
    rates[2] = (constants[10]+constants[11]*states[0]+constants[12]*(power(states[0], 2.00000)))/(1.00000+constants[13]*states[0]+constants[14]*(power(states[0], 2.00000))+constants[15]*states[1]+constants[16]*(power(states[1], 2.00000)))-constants[9]*states[2]
    rates[1] = (constants[17]+(constants[18]*states[0]+constants[19]*(power(states[0], 2.00000))+constants[20]*states[2]+constants[21]*(power(states[2], 2.00000)))/(1.00000+constants[22]*states[0]+constants[23]*(power(states[0], 2.00000))+constants[24]*states[2]+constants[25]*(power(states[2], 2.00000)))+constants[26]*states[3]+constants[27]*states[4])-constants[33]*states[1]
    rates[3] = constants[28]*states[1]-constants[34]*states[3]
    rates[4] = constants[29]*states[1]-constants[35]*states[4]
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