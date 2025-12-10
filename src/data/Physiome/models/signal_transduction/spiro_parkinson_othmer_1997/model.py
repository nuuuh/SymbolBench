# Size of variable arrays:
sizeAlgebraic = 9
sizeStates = 22
sizeConstants = 45
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "T2 in component T2 (micromolar)"
    legend_algebraic[0] = "Vmax in component T2 (flux)"
    legend_constants[37] = "KR in component T2 (molar)"
    legend_states[1] = "LT2 in component LT2 (micromolar)"
    legend_states[2] = "T2R in component T2R (micromolar)"
    legend_constants[0] = "R in component R (micromolar)"
    legend_states[3] = "T3 in component T3 (micromolar)"
    legend_constants[1] = "L in component L (micromolar)"
    legend_states[4] = "Bp in component Bp (micromolar)"
    legend_states[5] = "T2p in component T2p (micromolar)"
    legend_constants[2] = "Y0 in component Y0 (micromolar)"
    legend_states[6] = "Yp in component Yp (micromolar)"
    legend_constants[3] = "B0 in component B0 (micromolar)"
    legend_constants[4] = "k5 in component model_parameters (second_order_rate_constant)"
    legend_constants[5] = "k_5 in component model_parameters (first_order_rate_constant)"
    legend_constants[6] = "k8 in component model_parameters (first_order_rate_constant)"
    legend_constants[7] = "k_1 in component model_parameters (second_order_rate_constant)"
    legend_constants[8] = "ky in component model_parameters (second_order_rate_constant)"
    legend_constants[9] = "kb in component model_parameters (second_order_rate_constant)"
    legend_constants[10] = "k1a in component model_parameters (first_order_rate_constant)"
    legend_constants[11] = "k1b in component model_parameters (second_order_rate_constant)"
    legend_constants[12] = "k1c in component model_parameters (first_order_rate_constant)"
    legend_algebraic[1] = "Vmax in component T3 (flux)"
    legend_constants[38] = "KR in component T3 (molar)"
    legend_states[7] = "LT3 in component LT3 (micromolar)"
    legend_states[8] = "T3R in component T3R (micromolar)"
    legend_states[9] = "T4 in component T4 (micromolar)"
    legend_states[10] = "T3p in component T3p (micromolar)"
    legend_constants[13] = "k6 in component model_parameters (second_order_rate_constant)"
    legend_constants[14] = "k_6 in component model_parameters (first_order_rate_constant)"
    legend_constants[15] = "k9 in component model_parameters (first_order_rate_constant)"
    legend_constants[16] = "k_2 in component model_parameters (second_order_rate_constant)"
    legend_constants[17] = "k2a in component model_parameters (first_order_rate_constant)"
    legend_constants[18] = "k2b in component model_parameters (second_order_rate_constant)"
    legend_constants[19] = "k2c in component model_parameters (first_order_rate_constant)"
    legend_states[11] = "LT4 in component LT4 (micromolar)"
    legend_states[12] = "T4p in component T4p (micromolar)"
    legend_constants[20] = "k7 in component model_parameters (second_order_rate_constant)"
    legend_constants[21] = "k_7 in component model_parameters (first_order_rate_constant)"
    legend_constants[22] = "k10 in component model_parameters (first_order_rate_constant)"
    legend_algebraic[2] = "Vmax in component T2p (flux)"
    legend_constants[39] = "KR in component T2p (molar)"
    legend_states[13] = "LT2p in component LT2p (micromolar)"
    legend_states[14] = "T2pR in component T2pR (micromolar)"
    legend_algebraic[3] = "Vmax in component T3p (flux)"
    legend_constants[40] = "KR in component T3p (molar)"
    legend_states[15] = "LT3p in component LT3p (micromolar)"
    legend_states[16] = "T3pR in component T3pR (micromolar)"
    legend_states[17] = "LT4p in component LT4p (micromolar)"
    legend_algebraic[4] = "Vmax in component LT2 (flux)"
    legend_constants[41] = "KR in component LT2 (molar)"
    legend_states[18] = "LT2R in component LT2R (micromolar)"
    legend_constants[23] = "k11 in component model_parameters (first_order_rate_constant)"
    legend_constants[24] = "k_3 in component model_parameters (second_order_rate_constant)"
    legend_constants[25] = "k3a in component model_parameters (first_order_rate_constant)"
    legend_constants[26] = "k3b in component model_parameters (second_order_rate_constant)"
    legend_constants[27] = "k3c in component model_parameters (first_order_rate_constant)"
    legend_algebraic[5] = "Vmax in component LT3 (flux)"
    legend_constants[42] = "KR in component LT3 (molar)"
    legend_states[19] = "LT3R in component LT3R (micromolar)"
    legend_constants[28] = "k12 in component model_parameters (first_order_rate_constant)"
    legend_constants[29] = "k_4 in component model_parameters (second_order_rate_constant)"
    legend_constants[30] = "k4a in component model_parameters (first_order_rate_constant)"
    legend_constants[31] = "k4b in component model_parameters (second_order_rate_constant)"
    legend_constants[32] = "k4c in component model_parameters (first_order_rate_constant)"
    legend_constants[33] = "k13 in component model_parameters (first_order_rate_constant)"
    legend_algebraic[6] = "Vmax in component LT2p (flux)"
    legend_constants[43] = "KR in component LT2p (molar)"
    legend_states[20] = "LT2pR in component LT2pR (micromolar)"
    legend_algebraic[7] = "Vmax in component LT3p (flux)"
    legend_constants[44] = "KR in component LT3p (molar)"
    legend_states[21] = "LT3pR in component LT3pR (micromolar)"
    legend_algebraic[8] = "P in component P (micromolar)"
    legend_constants[34] = "Z in component Z (micromolar)"
    legend_constants[35] = "k_y in component model_parameters (second_order_rate_constant)"
    legend_constants[36] = "k_b in component model_parameters (first_order_rate_constant)"
    legend_rates[0] = "d/dt T2 in component T2 (micromolar)"
    legend_rates[3] = "d/dt T3 in component T3 (micromolar)"
    legend_rates[9] = "d/dt T4 in component T4 (micromolar)"
    legend_rates[5] = "d/dt T2p in component T2p (micromolar)"
    legend_rates[10] = "d/dt T3p in component T3p (micromolar)"
    legend_rates[12] = "d/dt T4p in component T4p (micromolar)"
    legend_rates[1] = "d/dt LT2 in component LT2 (micromolar)"
    legend_rates[7] = "d/dt LT3 in component LT3 (micromolar)"
    legend_rates[11] = "d/dt LT4 in component LT4 (micromolar)"
    legend_rates[13] = "d/dt LT2p in component LT2p (micromolar)"
    legend_rates[15] = "d/dt LT3p in component LT3p (micromolar)"
    legend_rates[17] = "d/dt LT4p in component LT4p (micromolar)"
    legend_rates[6] = "d/dt Yp in component Yp (micromolar)"
    legend_rates[4] = "d/dt Bp in component Bp (micromolar)"
    legend_rates[2] = "d/dt T2R in component T2R (micromolar)"
    legend_rates[8] = "d/dt T3R in component T3R (micromolar)"
    legend_rates[14] = "d/dt T2pR in component T2pR (micromolar)"
    legend_rates[16] = "d/dt T3pR in component T3pR (micromolar)"
    legend_rates[18] = "d/dt LT2R in component LT2R (micromolar)"
    legend_rates[19] = "d/dt LT3R in component LT3R (micromolar)"
    legend_rates[20] = "d/dt LT2pR in component LT2pR (micromolar)"
    legend_rates[21] = "d/dt LT3pR in component LT3pR (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.01
    states[1] = 0.01
    states[2] = 0.01
    constants[0] = 0.3
    states[3] = 0.01
    constants[1] = 3.0
    states[4] = 0.01
    states[5] = 0.01
    constants[2] = 20.0
    states[6] = 0.01
    constants[3] = 1.7
    constants[4] = 7E7
    constants[5] = 70.0
    constants[6] = 15.0
    constants[7] = 4E5
    constants[8] = 3E7
    constants[9] = 8E5
    constants[10] = 10
    constants[11] = 17
    constants[12] = 0.17
    states[7] = 0.01
    states[8] = 0.01
    states[9] = 0.01
    states[10] = 0.01
    constants[13] = 7E7
    constants[14] = 70.0
    constants[15] = 45.0
    constants[16] = 3E4
    constants[17] = 10
    constants[18] = 17
    constants[19] = 0.017
    states[11] = 0.01
    states[12] = 0.01
    constants[20] = 7E7
    constants[21] = 70.0
    constants[22] = 48.0
    states[13] = 0.01
    states[14] = 0.01
    states[15] = 0.01
    states[16] = 0.01
    states[17] = 0.01
    states[18] = 0.01
    constants[23] = 0.0
    constants[24] = 4E5
    constants[25] = 10
    constants[26] = 17
    constants[27] = 5.1
    states[19] = 0.01
    constants[28] = 16.5
    constants[29] = 3E4
    constants[30] = 10
    constants[31] = 17
    constants[32] = 0.51
    constants[33] = 34.56
    states[20] = 0.01
    states[21] = 0.01
    constants[34] = 40.0
    constants[35] = 5E5
    constants[36] = 0.35
    constants[37] = (constants[10]+constants[12])/constants[11]
    constants[38] = (constants[17]+constants[19])/constants[18]
    constants[39] = (constants[10]+constants[12])/constants[11]
    constants[40] = (constants[17]+constants[19])/constants[18]
    constants[41] = (constants[25]+constants[27])/constants[26]
    constants[42] = (constants[30]+constants[32])/constants[31]
    constants[43] = (constants[25]+constants[27])/constants[26]
    constants[44] = (constants[30]+constants[32])/constants[31]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[9] = (constants[21]*states[11]+constants[8]*states[12]*(constants[2]-states[6])+constants[9]*states[12]*(constants[3]-states[4]))-(constants[20]*constants[1]*states[9]+constants[22]*states[9])
    rates[12] = (constants[21]*states[17]+constants[8]*states[9]*(constants[2]-states[6])+constants[9]*states[9]*(constants[3]-states[4]))-(constants[20]*constants[1]*states[12]+constants[22]*states[12])
    rates[11] = (constants[21]*states[9]+constants[8]*states[17]*(constants[2]-states[6])+constants[9]*states[17]*(constants[3]-states[4]))-(constants[20]*constants[1]*states[11]+constants[33]*states[11])
    rates[17] = (constants[21]*states[12]+constants[8]*states[17]*(constants[2]-states[6])+constants[9]*states[17]*(constants[3]-states[4]))-(constants[20]*constants[1]*states[17]+constants[33]*states[17])
    rates[2] = constants[11]*states[0]*constants[0]-(constants[10]*states[2]+constants[12]*states[2])
    rates[8] = constants[18]*states[3]*constants[0]-(constants[17]*states[8]+constants[19]*states[8])
    rates[14] = constants[11]*states[5]*constants[0]-(constants[10]*states[14]+constants[12]*states[14])
    rates[16] = constants[18]*states[10]*constants[0]-(constants[17]*states[16]+constants[19]*states[16])
    rates[18] = constants[26]*states[1]*constants[0]-(constants[25]*states[18]+constants[27]*states[18])
    rates[19] = constants[31]*states[7]*constants[0]-(constants[30]*states[19]+constants[32]*states[19])
    rates[20] = constants[26]*states[13]*constants[0]-(constants[25]*states[20]+constants[27]*states[20])
    rates[21] = constants[31]*states[15]*constants[0]-(constants[30]*states[21]+constants[32]*states[21])
    algebraic[0] = constants[12]*(constants[0]+states[2])
    rates[0] = (constants[5]*states[1]+constants[7]*states[4]*states[3]+constants[8]*states[5]*(constants[2]-states[6])+constants[9]*states[5]*(constants[3]-states[4]))-(constants[4]*constants[1]*states[0]+constants[6]*states[0]+algebraic[0]*(states[0]/(constants[37]+states[0])))
    algebraic[1] = constants[19]*(constants[0]+states[8])
    rates[3] = (constants[14]*states[7]+constants[16]*states[4]*states[9]+constants[8]*states[10]*(constants[2]-states[6])+constants[9]*states[10]*(constants[3]-states[4]))-(constants[13]*constants[1]*states[3]+constants[15]*states[3]+algebraic[1]*(states[3]/(constants[38]+states[3])))
    algebraic[2] = constants[12]*(constants[0]+states[14])
    rates[5] = (constants[5]*states[13]+constants[7]*states[4]*states[10]+constants[8]*states[0]*(constants[2]-states[6])+constants[9]*states[0]*(constants[3]-states[4]))-(constants[4]*constants[1]*states[5]+constants[6]*states[5]+algebraic[2]*(states[5]/(constants[39]+states[5])))
    algebraic[3] = constants[19]*(constants[0]+states[16])
    rates[10] = (constants[14]*states[15]+constants[16]*states[4]*states[12]+constants[8]*states[3]*(constants[2]-states[6])+constants[9]*states[3]*(constants[3]-states[4]))-(constants[13]*constants[1]*states[10]+constants[15]*states[10]+algebraic[3]*(states[10]/(constants[40]+states[10])))
    algebraic[4] = constants[27]*(constants[0]+states[18])
    rates[1] = (constants[5]*states[0]+constants[24]*states[4]*states[7]+constants[8]*states[13]*(constants[2]-states[6])+constants[9]*states[13]*(constants[3]-states[4]))-(constants[4]*constants[1]*states[1]+constants[23]*states[1]+algebraic[4]*(states[1]/(constants[41]+states[1])))
    algebraic[5] = constants[32]*(constants[0]+states[19])
    rates[7] = (constants[14]*states[3]+constants[29]*states[4]*states[11]+constants[8]*states[15]*(constants[2]-states[6])+constants[9]*states[15]*(constants[3]-states[4]))-(constants[13]*constants[1]*states[7]+constants[28]*states[7]+algebraic[5]*(states[7]/(constants[42]+states[7])))
    algebraic[6] = constants[27]*(constants[0]+states[20])
    rates[13] = (constants[5]*states[5]+constants[24]*states[4]*states[15]+constants[8]*states[1]*(constants[2]-states[6])+constants[9]*states[1]*(constants[3]-states[4]))-(constants[4]*constants[1]*states[13]+constants[23]*states[13]+algebraic[6]*(states[13]/(constants[43]+states[13])))
    algebraic[7] = constants[32]*(constants[0]+states[21])
    rates[15] = (constants[14]*states[10]+constants[29]*states[4]*states[17]+constants[8]*states[7]*(constants[2]-states[6])+constants[9]*states[7]*(constants[3]-states[4]))-(constants[13]*constants[1]*states[15]+constants[28]*states[15]+algebraic[7]*(states[15]/(constants[44]+states[15])))
    algebraic[8] = states[5]+states[13]+states[10]+states[15]+states[12]+states[17]
    rates[6] = constants[8]*algebraic[8]*(constants[2]-states[6])-constants[35]*constants[34]*states[6]
    rates[4] = constants[9]*algebraic[8]*(constants[3]-states[4])-constants[36]*states[4]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[12]*(constants[0]+states[2])
    algebraic[1] = constants[19]*(constants[0]+states[8])
    algebraic[2] = constants[12]*(constants[0]+states[14])
    algebraic[3] = constants[19]*(constants[0]+states[16])
    algebraic[4] = constants[27]*(constants[0]+states[18])
    algebraic[5] = constants[32]*(constants[0]+states[19])
    algebraic[6] = constants[27]*(constants[0]+states[20])
    algebraic[7] = constants[32]*(constants[0]+states[21])
    algebraic[8] = states[5]+states[13]+states[10]+states[15]+states[12]+states[17]
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