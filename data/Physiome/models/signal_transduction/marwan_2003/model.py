# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 9
sizeConstants = 11
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_states[0] = "Pfr in component Pfr (micromolar)"
    legend_states[1] = "Pr in component Pr (micromolar)"
    legend_constants[0] = "Ifr_sigma_fr_phi_fr in component model_parameters (dimensionless)"
    legend_constants[1] = "Ir_sigma_r_phi_r in component model_parameters (dimensionless)"
    legend_constants[2] = "kd in component Pr (first_order_rate_constant)"
    legend_states[2] = "Xi in component Xi (micromolar)"
    legend_states[3] = "Xa in component Xa (micromolar)"
    legend_constants[3] = "kia in component model_parameters (second_order_rate_constant)"
    legend_constants[4] = "kai in component model_parameters (first_order_rate_constant)"
    legend_states[4] = "prepreS in component prepreS (micromolar)"
    legend_constants[5] = "kx in component model_parameters (second_order_rate_constant)"
    legend_states[5] = "preS in component preS (micromolar)"
    legend_states[6] = "Ya in component Ya (micromolar)"
    legend_constants[6] = "ky in component model_parameters (second_order_rate_constant)"
    legend_states[7] = "S in component S (micromolar)"
    legend_constants[7] = "alpha1 in component S (micromolar)"
    legend_states[8] = "V in component V (micromolar)"
    legend_constants[8] = "alpha2 in component V (micromolar)"
    legend_constants[9] = "kG in component Ya (second_order_rate_constant)"
    legend_constants[10] = "glucose in component model_parameters (micromolar)"
    legend_rates[0] = "d/dt Pfr in component Pfr (micromolar)"
    legend_rates[1] = "d/dt Pr in component Pr (micromolar)"
    legend_rates[2] = "d/dt Xi in component Xi (micromolar)"
    legend_rates[3] = "d/dt Xa in component Xa (micromolar)"
    legend_rates[4] = "d/dt prepreS in component prepreS (micromolar)"
    legend_rates[5] = "d/dt preS in component preS (micromolar)"
    legend_rates[7] = "d/dt S in component S (micromolar)"
    legend_rates[8] = "d/dt V in component V (micromolar)"
    legend_rates[6] = "d/dt Ya in component Ya (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 10.0
    states[1] = 0.0
    constants[0] = 0.1
    constants[1] = 0.0
    constants[2] = 0.1
    states[2] = 6.0
    states[3] = 0.0
    constants[3] = 0.1
    constants[4] = 0.8
    states[4] = 200.0
    constants[5] = 0.2
    states[5] = 0.0
    states[6] = 0.9
    constants[6] = 1.0
    states[7] = 0.0
    constants[7] = 30.0
    states[8] = 50.0
    constants[8] = 50.0
    constants[9] = 0.1
    constants[10] = 1.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[1]*states[1]-constants[0]*states[0]
    rates[1] = constants[0]*states[1]-(constants[1]*states[1]+constants[2]*states[1])
    rates[2] = constants[4]*states[3]-constants[3]*states[1]*states[2]
    rates[3] = constants[3]*states[1]*states[2]-constants[4]*states[3]
    rates[4] = -(constants[5]*states[3]*states[4])
    rates[5] = constants[5]*states[3]*states[4]-constants[6]*states[6]*states[5]
    rates[7] = (constants[6]*states[6]*states[5]+constants[7]/(1.00000+power(states[8], 3.00000)))-states[7]
    rates[8] = constants[8]/(1.00000+power(states[7], 3.00000))-states[8]
    rates[6] = -(constants[9]*constants[10]*states[6])
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