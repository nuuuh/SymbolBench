# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 10
sizeConstants = 38
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_constants[0] = "v_sP in component nucleus (nanomolar_hour)"
    legend_constants[1] = "v_mP in component nucleus (nanomolar_hour)"
    legend_constants[2] = "K_IP in component nucleus (nanomolar)"
    legend_constants[3] = "K_mP in component nucleus (nanomolar)"
    legend_constants[4] = "v_sT in component nucleus (nanomolar_hour)"
    legend_constants[5] = "v_mT in component nucleus (nanomolar_hour)"
    legend_constants[6] = "K_IT in component nucleus (nanomolar)"
    legend_constants[7] = "K_mT in component nucleus (nanomolar)"
    legend_constants[8] = "k_d in component cytosol (per_hour)"
    legend_constants[9] = "n in component nucleus (dimensionless)"
    legend_constants[10] = "k_1 in component cytosol (per_hour)"
    legend_constants[11] = "k_2 in component cytosol (per_hour)"
    legend_constants[12] = "k_dN in component nucleus (per_hour)"
    legend_states[0] = "C in component cytosol (nanomolar)"
    legend_states[1] = "M_P in component nucleus (nanomolar)"
    legend_states[2] = "M_T in component nucleus (nanomolar)"
    legend_states[3] = "C_N in component nucleus (nanomolar)"
    legend_constants[13] = "k_3 in component cytosol (per_nanomolar_hour)"
    legend_constants[14] = "k_4 in component cytosol (per_hour)"
    legend_constants[15] = "k_dC in component cytosol (per_hour)"
    legend_states[4] = "P_0 in component PER (nanomolar)"
    legend_states[5] = "P_1 in component PER (nanomolar)"
    legend_states[6] = "P_2 in component PER (nanomolar)"
    legend_states[7] = "T_0 in component TIM (nanomolar)"
    legend_states[8] = "T_1 in component TIM (nanomolar)"
    legend_states[9] = "T_2 in component TIM (nanomolar)"
    legend_constants[16] = "V_1P in component PER (nanomolar_hour)"
    legend_constants[17] = "V_2P in component PER (nanomolar_hour)"
    legend_constants[18] = "V_3P in component PER (nanomolar_hour)"
    legend_constants[19] = "V_4P in component PER (nanomolar_hour)"
    legend_constants[20] = "K_1P in component PER (nanomolar)"
    legend_constants[21] = "K_2P in component PER (nanomolar)"
    legend_constants[22] = "K_3P in component PER (nanomolar)"
    legend_constants[23] = "K_4P in component PER (nanomolar)"
    legend_constants[24] = "K_dP in component PER (nanomolar)"
    legend_constants[25] = "v_dP in component PER (nanomolar_hour)"
    legend_constants[26] = "k_sP in component PER (per_hour)"
    legend_constants[27] = "V_1T in component TIM (nanomolar_hour)"
    legend_constants[28] = "V_2T in component TIM (nanomolar_hour)"
    legend_constants[29] = "V_3T in component TIM (nanomolar_hour)"
    legend_constants[30] = "V_4T in component TIM (nanomolar_hour)"
    legend_constants[31] = "K_1T in component TIM (nanomolar)"
    legend_constants[32] = "K_2T in component TIM (nanomolar)"
    legend_constants[33] = "K_3T in component TIM (nanomolar)"
    legend_constants[34] = "K_4T in component TIM (nanomolar)"
    legend_constants[35] = "K_dT in component TIM (nanomolar)"
    legend_constants[36] = "v_dT in component TIM (nanomolar_hour)"
    legend_constants[37] = "k_sT in component TIM (per_hour)"
    legend_algebraic[0] = "P_t in component PER_total (nanomolar)"
    legend_algebraic[1] = "T_t in component TIM_total (nanomolar)"
    legend_rates[1] = "d/dt M_P in component nucleus (nanomolar)"
    legend_rates[2] = "d/dt M_T in component nucleus (nanomolar)"
    legend_rates[3] = "d/dt C_N in component nucleus (nanomolar)"
    legend_rates[0] = "d/dt C in component cytosol (nanomolar)"
    legend_rates[4] = "d/dt P_0 in component PER (nanomolar)"
    legend_rates[5] = "d/dt P_1 in component PER (nanomolar)"
    legend_rates[6] = "d/dt P_2 in component PER (nanomolar)"
    legend_rates[7] = "d/dt T_0 in component TIM (nanomolar)"
    legend_rates[8] = "d/dt T_1 in component TIM (nanomolar)"
    legend_rates[9] = "d/dt T_2 in component TIM (nanomolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1
    constants[1] = 0.7
    constants[2] = 1
    constants[3] = 0.2
    constants[4] = 1
    constants[5] = 0.7
    constants[6] = 1
    constants[7] = 0.2
    constants[8] = 0.01
    constants[9] = 4
    constants[10] = 0.6
    constants[11] = 0.2
    constants[12] = 0.01
    states[0] = 0.344
    states[1] = 0.031
    states[2] = 0.031
    states[3] = 1.77
    constants[13] = 1.2
    constants[14] = 0.6
    constants[15] = 0.01
    states[4] = 0.0114
    states[5] = 0.0178
    states[6] = 0.0322
    states[7] = 0.0114
    states[8] = 0.0178
    states[9] = 0.0324
    constants[16] = 8
    constants[17] = 1
    constants[18] = 8
    constants[19] = 1
    constants[20] = 2
    constants[21] = 2
    constants[22] = 2
    constants[23] = 2
    constants[24] = 0.2
    constants[25] = 2
    constants[26] = 0.9
    constants[27] = 8
    constants[28] = 1
    constants[29] = 8
    constants[30] = 1
    constants[31] = 2
    constants[32] = 2
    constants[33] = 2
    constants[34] = 2
    constants[35] = 0.2
    constants[36] = 2
    constants[37] = 0.9
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = (constants[0]*((power(constants[2], constants[9]))/(power(constants[2], constants[9])+power(states[3], constants[9])))-constants[1]*(states[1]/(constants[3]+states[1])))-constants[8]*states[1]
    rates[2] = (constants[4]*((power(constants[6], constants[9]))/(power(constants[6], constants[9])+power(states[3], constants[9])))-constants[5]*(states[2]/(constants[7]+states[2])))-constants[8]*states[2]
    rates[3] = (constants[10]*states[0]-constants[11]*states[3])-constants[12]*states[3]
    rates[0] = (((constants[13]*states[6]*states[9]-constants[14]*states[0])-constants[10]*states[0])+constants[11]*states[3])-constants[15]*states[0]
    rates[4] = ((constants[26]*states[1]-constants[16]*(states[4]/(constants[20]+states[4])))+constants[17]*(states[5]/(constants[21]+states[5])))-constants[8]*states[4]
    rates[5] = (((constants[16]*(states[4]/(constants[20]+states[4]))-constants[17]*(states[5]/(constants[21]+states[5])))-constants[18]*(states[5]/(constants[22]+states[5])))+constants[19]*(states[6]/(constants[23]+states[6])))-constants[8]*states[5]
    rates[6] = ((((constants[18]*(states[5]/(constants[22]+states[5]))-constants[19]*(states[6]/(constants[23]+states[6])))-constants[13]*states[6]*states[9])+constants[14]*states[0])-constants[25]*(states[6]/(constants[24]+states[6])))-constants[8]*states[6]
    rates[7] = ((constants[37]*states[2]-constants[27]*(states[7]/(constants[31]+states[7])))+constants[28]*(states[8]/(constants[32]+states[8])))-constants[8]*states[7]
    rates[8] = (((constants[27]*(states[7]/(constants[31]+states[7]))-constants[28]*(states[8]/(constants[32]+states[8])))-constants[29]*(states[8]/(constants[33]+states[8])))+constants[30]*(states[9]/(constants[34]+states[9])))-constants[8]*states[8]
    rates[9] = ((((constants[29]*(states[8]/(constants[33]+states[8]))-constants[30]*(states[9]/(constants[34]+states[9])))-constants[13]*states[6]*states[9])+constants[14]*states[0])-constants[36]*(states[9]/(constants[35]+states[9])))-constants[8]*states[9]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[4]+states[5]+states[6]+states[0]+states[3]
    algebraic[1] = states[7]+states[8]+states[9]+states[0]+states[3]
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