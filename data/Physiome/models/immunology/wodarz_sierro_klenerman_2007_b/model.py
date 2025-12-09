# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 15
sizeConstants = 16
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_states[0] = "x in component x (dimensionless)"
    legend_constants[0] = "d in component model_parameters (first_order_rate_constant)"
    legend_constants[1] = "gamma in component model_parameters (first_order_rate_constant)"
    legend_constants[2] = "lambda in component model_parameters (first_order_rate_constant)"
    legend_constants[3] = "beta in component model_parameters (first_order_rate_constant)"
    legend_states[1] = "v in component v (dimensionless)"
    legend_states[2] = "y0 in component y0 (dimensionless)"
    legend_constants[4] = "a0 in component model_parameters (first_order_rate_constant)"
    legend_constants[5] = "eta in component model_parameters (first_order_rate_constant)"
    legend_constants[6] = "phi in component model_parameters (first_order_rate_constant)"
    legend_states[3] = "L in component L (dimensionless)"
    legend_constants[7] = "pa in component model_parameters (first_order_rate_constant)"
    legend_states[4] = "za in component za (dimensionless)"
    legend_states[5] = "y1 in component y1 (dimensionless)"
    legend_constants[8] = "a1 in component model_parameters (first_order_rate_constant)"
    legend_constants[9] = "k in component model_parameters (first_order_rate_constant)"
    legend_constants[10] = "u in component model_parameters (first_order_rate_constant)"
    legend_states[6] = "m_8 in component m_8 (dimensionless)"
    legend_constants[11] = "r in component model_parameters (first_order_rate_constant)"
    legend_states[7] = "m_7 in component m_7 (dimensionless)"
    legend_states[8] = "m_6 in component m_6 (dimensionless)"
    legend_states[9] = "m_5 in component m_5 (dimensionless)"
    legend_states[10] = "m_4 in component m_4 (dimensionless)"
    legend_states[11] = "m_3 in component m_3 (dimensionless)"
    legend_states[12] = "m_2 in component m_2 (dimensionless)"
    legend_states[13] = "m_1 in component m_1 (dimensionless)"
    legend_states[14] = "m_0 in component m_0 (dimensionless)"
    legend_constants[12] = "alpha in component model_parameters (first_order_rate_constant)"
    legend_constants[13] = "ca in component model_parameters (first_order_rate_constant)"
    legend_constants[14] = "ba in component model_parameters (first_order_rate_constant)"
    legend_algebraic[0] = "log_za in component za (dimensionless)"
    legend_constants[15] = "R0 in component R0 (dimensionless)"
    legend_rates[0] = "d/dt x in component x (dimensionless)"
    legend_rates[2] = "d/dt y0 in component y0 (dimensionless)"
    legend_rates[5] = "d/dt y1 in component y1 (dimensionless)"
    legend_rates[3] = "d/dt L in component L (dimensionless)"
    legend_rates[1] = "d/dt v in component v (dimensionless)"
    legend_rates[6] = "d/dt m_8 in component m_8 (dimensionless)"
    legend_rates[7] = "d/dt m_7 in component m_7 (dimensionless)"
    legend_rates[8] = "d/dt m_6 in component m_6 (dimensionless)"
    legend_rates[9] = "d/dt m_5 in component m_5 (dimensionless)"
    legend_rates[10] = "d/dt m_4 in component m_4 (dimensionless)"
    legend_rates[11] = "d/dt m_3 in component m_3 (dimensionless)"
    legend_rates[12] = "d/dt m_2 in component m_2 (dimensionless)"
    legend_rates[13] = "d/dt m_1 in component m_1 (dimensionless)"
    legend_rates[14] = "d/dt m_0 in component m_0 (dimensionless)"
    legend_rates[4] = "d/dt za in component za (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 1
    constants[0] = 0.1
    constants[1] = 0.5
    constants[2] = 10
    constants[3] = 0.1
    states[1] = 1
    states[2] = 0
    constants[4] = 0.1
    constants[5] = 0.01
    constants[6] = 0.1
    states[3] = 0
    constants[7] = 0.000001
    states[4] = 1
    states[5] = 0
    constants[8] = 0.2
    constants[9] = 1
    constants[10] = 1
    states[6] = 0
    constants[11] = 1
    states[7] = 0
    states[8] = 0
    states[9] = 0
    states[10] = 0
    states[11] = 0
    states[12] = 0
    states[13] = 0
    states[14] = 1
    constants[12] = 0.2
    constants[13] = 15.5
    constants[14] = 0.1
    constants[15] = ((constants[2]*constants[5])/(constants[0]*constants[8]*(constants[4]+constants[5])))*(constants[3]+(constants[1]*constants[6])/(constants[6]+constants[0]))
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[2]-(constants[0]*states[0]+constants[3]*states[0]*states[1]+constants[1]*states[0]*states[1])
    rates[2] = (constants[3]*states[0]*states[1]-(constants[4]*states[2]+constants[5]*states[2]+constants[7]*states[2]*states[4]))+constants[6]*states[3]
    rates[5] = constants[5]*states[2]-(constants[8]*states[5]+constants[7]*states[5]*states[4])
    rates[3] = constants[1]*states[0]*states[1]-(constants[6]*states[3]+constants[0]*states[3])
    rates[1] = constants[9]*states[5]-constants[10]*states[1]
    rates[6] = 2.00000*constants[11]*states[7]-constants[11]*states[6]
    rates[7] = 2.00000*constants[11]*states[8]-constants[11]*states[7]
    rates[8] = 2.00000*constants[11]*states[9]-constants[11]*states[8]
    rates[9] = 2.00000*constants[11]*states[10]-constants[11]*states[9]
    rates[10] = 2.00000*constants[11]*states[11]-constants[11]*states[10]
    rates[11] = 2.00000*constants[11]*states[12]-constants[11]*states[11]
    rates[12] = 2.00000*constants[11]*states[13]-constants[11]*states[12]
    rates[13] = 2.00000*constants[11]*states[14]-constants[11]*states[13]
    rates[14] = -constants[11]*states[14]
    rates[4] = (constants[12]*states[6]+constants[13]*(states[2]+states[5])*states[4])-constants[14]*states[4]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = log(states[4], 10)
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