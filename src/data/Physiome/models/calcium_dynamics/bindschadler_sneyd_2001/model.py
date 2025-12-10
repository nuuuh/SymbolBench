# Size of variable arrays:
sizeAlgebraic = 13
sizeStates = 4
sizeConstants = 14
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_algebraic[4] = "phi3_c1 in component phi (per_second)"
    legend_states[0] = "h1 in component h1 (dimensionless)"
    legend_algebraic[0] = "phi1_c1 in component phi (second_order_rate)"
    legend_algebraic[2] = "phi2_c1 in component phi (per_second)"
    legend_constants[0] = "p in component model_parameters (micro_molar)"
    legend_algebraic[1] = "phi_1_c1 in component phi (per_second)"
    legend_algebraic[8] = "phi3_c2 in component phi (per_second)"
    legend_states[1] = "h2 in component h2 (dimensionless)"
    legend_algebraic[3] = "phi1_c2 in component phi (second_order_rate)"
    legend_algebraic[6] = "phi2_c2 in component phi (per_second)"
    legend_algebraic[5] = "phi_1_c2 in component phi (per_second)"
    legend_constants[1] = "r2 in component model_parameters (second_order_rate)"
    legend_constants[2] = "R1 in component model_parameters (micro_molar)"
    legend_constants[3] = "k1 in component model_parameters (micro_molar_per_second)"
    legend_constants[4] = "R3 in component model_parameters (micro_molar)"
    legend_constants[5] = "k2 in component model_parameters (micro_molar_per_second)"
    legend_constants[6] = "r4 in component model_parameters (per_second)"
    legend_constants[7] = "k3 in component model_parameters (micro_molar_per_second)"
    legend_constants[8] = "R5 in component model_parameters (micro_molar)"
    legend_states[2] = "c1 in component c1 (micro_molar)"
    legend_states[3] = "c2 in component c2 (micro_molar)"
    legend_constants[9] = "Vp in component model_parameters (micro_molar_per_second)"
    legend_constants[10] = "Kp in component model_parameters (micro_molar)"
    legend_algebraic[7] = "j_pump_c1 in component j_pump (micro_molar_per_second)"
    legend_algebraic[9] = "j_pump_c2 in component j_pump (micro_molar_per_second)"
    legend_constants[11] = "kf in component model_parameters (micro_molar_per_second)"
    legend_algebraic[10] = "j_receptor_c1 in component j_receptor (micro_molar_per_second)"
    legend_algebraic[11] = "j_receptor_c2 in component j_receptor (micro_molar_per_second)"
    legend_algebraic[12] = "j_diffusion in component j_diffusion (micro_molar_per_second)"
    legend_constants[12] = "D in component model_parameters (per_second)"
    legend_constants[13] = "j_leak in component model_parameters (micro_molar_per_second)"
    legend_rates[0] = "d/dt h1 in component h1 (dimensionless)"
    legend_rates[1] = "d/dt h2 in component h2 (dimensionless)"
    legend_rates[2] = "d/dt c1 in component c1 (micro_molar)"
    legend_rates[3] = "d/dt c2 in component c2 (micro_molar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.8
    constants[0] = 0.2778
    states[1] = 0.1
    constants[1] = 100
    constants[2] = 6
    constants[3] = 44
    constants[4] = 50
    constants[5] = 26.5
    constants[6] = 20
    constants[7] = 1.6
    constants[8] = 1.6
    states[2] = 0.3
    states[3] = 0.1
    constants[9] = 1.2
    constants[10] = 0.18
    constants[11] = 28
    constants[12] = 0.01
    constants[13] = 0.2
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[4] = constants[7]/(constants[8]+states[2])
    algebraic[0] = (constants[1]*states[2])/(constants[2]+states[2])
    algebraic[2] = (constants[5]+constants[6]*states[2])/(constants[4]+states[2])
    algebraic[1] = constants[3]/(constants[4]+states[2])
    rates[0] = algebraic[4]*(1.00000-states[0])-(algebraic[0]*algebraic[2]*states[0]*constants[0])/(algebraic[0]*constants[0]+algebraic[1])
    algebraic[8] = constants[7]/(constants[8]+states[3])
    algebraic[3] = (constants[1]*states[3])/(constants[2]+states[3])
    algebraic[6] = (constants[5]+constants[6]*states[3])/(constants[4]+states[3])
    algebraic[5] = constants[3]/(constants[4]+states[3])
    rates[1] = algebraic[8]*(1.00000-states[1])-(algebraic[3]*algebraic[6]*states[1]*constants[0])/(algebraic[3]*constants[0]+algebraic[5])
    algebraic[7] = (constants[9]*(power(states[2], 2.00000)))/(power(constants[10], 2.00000)+power(states[2], 2.00000))
    algebraic[10] = constants[11]*(power((constants[0]*states[0]*algebraic[0])/(algebraic[0]*constants[0]+algebraic[1]), 4.00000))
    algebraic[12] = constants[12]*(states[3]-states[2])
    rates[2] = (algebraic[10]-algebraic[7])+constants[13]+algebraic[12]
    algebraic[9] = (constants[9]*(power(states[3], 2.00000)))/(power(constants[10], 2.00000)+power(states[3], 2.00000))
    algebraic[11] = constants[11]*(power((constants[0]*states[1]*algebraic[3])/(algebraic[3]*constants[0]+algebraic[5]), 4.00000))
    rates[3] = (algebraic[11]-algebraic[9])+constants[13]+algebraic[12]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[4] = constants[7]/(constants[8]+states[2])
    algebraic[0] = (constants[1]*states[2])/(constants[2]+states[2])
    algebraic[2] = (constants[5]+constants[6]*states[2])/(constants[4]+states[2])
    algebraic[1] = constants[3]/(constants[4]+states[2])
    algebraic[8] = constants[7]/(constants[8]+states[3])
    algebraic[3] = (constants[1]*states[3])/(constants[2]+states[3])
    algebraic[6] = (constants[5]+constants[6]*states[3])/(constants[4]+states[3])
    algebraic[5] = constants[3]/(constants[4]+states[3])
    algebraic[7] = (constants[9]*(power(states[2], 2.00000)))/(power(constants[10], 2.00000)+power(states[2], 2.00000))
    algebraic[10] = constants[11]*(power((constants[0]*states[0]*algebraic[0])/(algebraic[0]*constants[0]+algebraic[1]), 4.00000))
    algebraic[12] = constants[12]*(states[3]-states[2])
    algebraic[9] = (constants[9]*(power(states[3], 2.00000)))/(power(constants[10], 2.00000)+power(states[3], 2.00000))
    algebraic[11] = constants[11]*(power((constants[0]*states[1]*algebraic[3])/(algebraic[3]*constants[0]+algebraic[5]), 4.00000))
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