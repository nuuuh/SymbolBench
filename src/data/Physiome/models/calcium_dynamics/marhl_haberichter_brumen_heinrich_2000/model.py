# Size of variable arrays:
sizeAlgebraic = 7
sizeStates = 3
sizeConstants = 17
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_algebraic[0] = "Ca_Pr in component total_calcium (micromolar)"
    legend_constants[0] = "Ca_tot in component total_calcium (micromolar)"
    legend_constants[1] = "rho_ER in component ER_calcium (dimensionless)"
    legend_constants[2] = "beta_ER in component ER_calcium (dimensionless)"
    legend_constants[3] = "rho_m in component mitochondrial_calcium (dimensionless)"
    legend_constants[4] = "beta_m in component mitochondrial_calcium (dimensionless)"
    legend_states[0] = "Ca_cyt in component cytosolic_calcium (micromolar)"
    legend_states[1] = "Ca_ER in component ER_calcium (micromolar)"
    legend_states[2] = "Ca_m in component mitochondrial_calcium (micromolar)"
    legend_algebraic[1] = "Pr in component total_protein (micromolar)"
    legend_constants[5] = "Pr_tot in component total_protein (micromolar)"
    legend_constants[6] = "k_plus in component cytosolic_calcium (second_order_rate_constant)"
    legend_constants[7] = "k_minus in component cytosolic_calcium (first_order_rate_constant)"
    legend_algebraic[3] = "J_ch in component Ca_efflux_from_the_ER (flux)"
    legend_algebraic[4] = "J_leak in component Ca_leak_flux_from_the_ER (flux)"
    legend_algebraic[2] = "J_pump in component ATP_dependent_Ca_uptake_into_the_ER (flux)"
    legend_algebraic[6] = "J_out in component mitochondrial_Ca_release (flux)"
    legend_algebraic[5] = "J_in in component mitochondrial_Ca_uptake (flux)"
    legend_constants[8] = "k_pump in component ATP_dependent_Ca_uptake_into_the_ER (first_order_rate_constant)"
    legend_constants[9] = "k_ch in component Ca_efflux_from_the_ER (first_order_rate_constant)"
    legend_constants[10] = "K1 in component Ca_efflux_from_the_ER (micromolar)"
    legend_constants[11] = "k_leak in component Ca_leak_flux_from_the_ER (first_order_rate_constant)"
    legend_constants[12] = "k_in in component mitochondrial_Ca_uptake (flux)"
    legend_constants[13] = "K2 in component mitochondrial_Ca_uptake (micromolar)"
    legend_constants[14] = "k_out in component mitochondrial_Ca_release (first_order_rate_constant)"
    legend_constants[15] = "k_m in component mitochondrial_Ca_release (first_order_rate_constant)"
    legend_constants[16] = "K3 in component mitochondrial_Ca_release (micromolar)"
    legend_rates[0] = "d/dt Ca_cyt in component cytosolic_calcium (micromolar)"
    legend_rates[1] = "d/dt Ca_ER in component ER_calcium (micromolar)"
    legend_rates[2] = "d/dt Ca_m in component mitochondrial_calcium (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 90.0
    constants[1] = 0.01
    constants[2] = 0.0025
    constants[3] = 0.01
    constants[4] = 0.0025
    states[0] = 0.05
    states[1] = 1.0
    states[2] = 0.4
    constants[5] = 120.0
    constants[6] = 0.1
    constants[7] = 0.01
    constants[8] = 20.0
    constants[9] = 4100.0
    constants[10] = 5.0
    constants[11] = 0.05
    constants[12] = 300.0
    constants[13] = 0.8
    constants[14] = 125.0
    constants[15] = 0.00625
    constants[16] = 5.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[3] = constants[9]*((power(states[0], 2.00000))/(power(constants[10], 2.00000)+power(states[0], 2.00000)))*(states[1]-states[0])
    algebraic[4] = constants[11]*(states[1]-states[0])
    algebraic[2] = constants[8]*states[0]
    rates[1] = (constants[2]/constants[1])*(algebraic[2]-(algebraic[3]+algebraic[4]))
    algebraic[0] = constants[0]-(states[0]+(constants[1]/constants[2])*states[1]+(constants[3]/constants[4])*states[2])
    algebraic[1] = constants[5]-algebraic[0]
    algebraic[6] = (constants[14]*((power(states[0], 2.00000))/(power(constants[16], 2.00000)+power(states[0], 2.00000)))+constants[15])*states[2]
    algebraic[5] = constants[12]*((power(states[0], 8.00000))/(power(constants[13], 8.00000)+power(states[0], 8.00000)))
    rates[0] = (algebraic[3]+algebraic[4]+algebraic[6]+constants[7]*algebraic[0])-(algebraic[2]+algebraic[5]+constants[6]*states[0]*algebraic[1])
    rates[2] = (constants[4]/constants[3])*(algebraic[5]-algebraic[6])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = constants[9]*((power(states[0], 2.00000))/(power(constants[10], 2.00000)+power(states[0], 2.00000)))*(states[1]-states[0])
    algebraic[4] = constants[11]*(states[1]-states[0])
    algebraic[2] = constants[8]*states[0]
    algebraic[0] = constants[0]-(states[0]+(constants[1]/constants[2])*states[1]+(constants[3]/constants[4])*states[2])
    algebraic[1] = constants[5]-algebraic[0]
    algebraic[6] = (constants[14]*((power(states[0], 2.00000))/(power(constants[16], 2.00000)+power(states[0], 2.00000)))+constants[15])*states[2]
    algebraic[5] = constants[12]*((power(states[0], 8.00000))/(power(constants[13], 8.00000)+power(states[0], 8.00000)))
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