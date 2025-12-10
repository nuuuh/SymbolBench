# Size of variable arrays:
sizeAlgebraic = 21
sizeStates = 11
sizeConstants = 33
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_constants[0] = "k1 in component rate_constants (first_order_rate_constant)"
    legend_constants[1] = "k1_ in component rate_constants (first_order_rate_constant)"
    legend_constants[2] = "k1__ in component rate_constants (first_order_rate_constant)"
    legend_constants[3] = "kn1 in component rate_constants (first_order_rate_constant)"
    legend_constants[4] = "k2 in component rate_constants (first_order_rate_constant)"
    legend_constants[5] = "kn2 in component rate_constants (first_order_rate_constant)"
    legend_constants[6] = "k3 in component rate_constants (first_order_rate_constant)"
    legend_constants[7] = "k3_ in component rate_constants (first_order_rate_constant)"
    legend_constants[8] = "kn4 in component rate_constants (first_order_rate_constant)"
    legend_constants[9] = "k5 in component rate_constants (first_order_rate_constant)"
    legend_constants[10] = "kn6 in component rate_constants (first_order_rate_constant)"
    legend_constants[11] = "k8 in component rate_constants (first_order_rate_constant)"
    legend_constants[12] = "k9 in component rate_constants (first_order_rate_constant)"
    legend_constants[13] = "k10 in component rate_constants (first_order_rate_constant)"
    legend_constants[14] = "k17 in component rate_constants (first_order_rate_constant)"
    legend_constants[15] = "k18 in component rate_constants (first_order_rate_constant)"
    legend_constants[16] = "k19 in component rate_constants (first_order_rate_constant)"
    legend_constants[17] = "k20 in component rate_constants (first_order_rate_constant)"
    legend_constants[18] = "k21 in component rate_constants (first_order_rate_constant)"
    legend_constants[19] = "k22 in component rate_constants (first_order_rate_constant)"
    legend_constants[20] = "k24 in component rate_constants (first_order_rate_constant)"
    legend_constants[21] = "k25 in component rate_constants (first_order_rate_constant)"
    legend_constants[22] = "k25_ in component rate_constants (first_order_rate_constant)"
    legend_constants[23] = "k26 in component rate_constants (first_order_rate_constant)"
    legend_constants[24] = "k26_ in component rate_constants (first_order_rate_constant)"
    legend_constants[25] = "k28 in component rate_constants (first_order_rate_constant)"
    legend_constants[26] = "k29 in component rate_constants (first_order_rate_constant)"
    legend_states[0] = "a_CyclinE_Cdk2 in component a_CyclinE_Cdk2 (dimensionless)"
    legend_algebraic[0] = "V_2 in component V_2 (first_order_rate_constant)"
    legend_algebraic[13] = "V_10 in component V_10 (first_order_rate_constant)"
    legend_algebraic[2] = "V_n2 in component V_n2 (first_order_rate_constant)"
    legend_algebraic[12] = "V_9 in component V_9 (first_order_rate_constant)"
    legend_algebraic[14] = "V_21 in component V_21 (first_order_rate_constant)"
    legend_states[1] = "i_CyclinE_Cdk2 in component i_CyclinE_Cdk2 (dimensionless)"
    legend_algebraic[5] = "V_3 in component V_3 (first_order_rate_constant)"
    legend_algebraic[9] = "V_5 in component V_5 (first_order_rate_constant)"
    legend_states[2] = "pRB_E2F in component pRB_E2F (dimensionless)"
    legend_constants[32] = "V_n1 in component V_n1 (first_order_rate_constant)"
    legend_algebraic[1] = "V_1 in component V_1 (first_order_rate_constant)"
    legend_states[3] = "E2F in component E2F (dimensionless)"
    legend_constants[27] = "V_4 in component E2F (first_order_rate_constant)"
    legend_algebraic[6] = "V_18 in component V_18 (first_order_rate_constant)"
    legend_algebraic[3] = "V_n4 in component V_n4 (first_order_rate_constant)"
    legend_states[4] = "pRB in component pRB (dimensionless)"
    legend_constants[28] = "V_27 in component pRB (first_order_rate_constant)"
    legend_algebraic[4] = "V_26 in component V_26 (first_order_rate_constant)"
    legend_algebraic[10] = "V_29 in component V_29 (first_order_rate_constant)"
    legend_algebraic[7] = "V_28 in component V_28 (first_order_rate_constant)"
    legend_states[5] = "CycD_Cdk4 in component CycD_Cdk4 (dimensionless)"
    legend_constants[29] = "V_6 in component CycD_Cdk4 (first_order_rate_constant)"
    legend_algebraic[18] = "V_20 in component V_20 (first_order_rate_constant)"
    legend_algebraic[8] = "V_n6 in component V_n6 (first_order_rate_constant)"
    legend_algebraic[16] = "V_19 in component V_19 (first_order_rate_constant)"
    legend_algebraic[15] = "V_17 in component V_17 (first_order_rate_constant)"
    legend_states[6] = "p27 in component p27 (dimensionless)"
    legend_constants[30] = "V_7 in component p27 (first_order_rate_constant)"
    legend_algebraic[11] = "V_8 in component V_8 (first_order_rate_constant)"
    legend_algebraic[20] = "V_22 in component V_22 (first_order_rate_constant)"
    legend_states[7] = "CycE_Cdk2_p27 in component CycE_Cdk2_p27 (dimensionless)"
    legend_states[8] = "CycD_Cdk4_p27 in component CycD_Cdk4_p27 (dimensionless)"
    legend_states[9] = "p16 in component p16 (dimensionless)"
    legend_constants[31] = "V_23 in component p16 (first_order_rate_constant)"
    legend_algebraic[19] = "V_25 in component V_25 (first_order_rate_constant)"
    legend_algebraic[17] = "V_24 in component V_24 (first_order_rate_constant)"
    legend_states[10] = "pRB_P in component pRB_P (dimensionless)"
    legend_rates[0] = "d/dt a_CyclinE_Cdk2 in component a_CyclinE_Cdk2 (dimensionless)"
    legend_rates[1] = "d/dt i_CyclinE_Cdk2 in component i_CyclinE_Cdk2 (dimensionless)"
    legend_rates[2] = "d/dt pRB_E2F in component pRB_E2F (dimensionless)"
    legend_rates[3] = "d/dt E2F in component E2F (dimensionless)"
    legend_rates[4] = "d/dt pRB in component pRB (dimensionless)"
    legend_rates[5] = "d/dt CycD_Cdk4 in component CycD_Cdk4 (dimensionless)"
    legend_rates[6] = "d/dt p27 in component p27 (dimensionless)"
    legend_rates[7] = "d/dt CycE_Cdk2_p27 in component CycE_Cdk2_p27 (dimensionless)"
    legend_rates[8] = "d/dt CycD_Cdk4_p27 in component CycD_Cdk4_p27 (dimensionless)"
    legend_rates[9] = "d/dt p16 in component p16 (dimensionless)"
    legend_rates[10] = "d/dt pRB_P in component pRB_P (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 0.1
    constants[1] = 0.5
    constants[2] = 0.5
    constants[3] = 0.001
    constants[4] = 0.1
    constants[5] = 1
    constants[6] = 1.42
    constants[7] = 0
    constants[8] = 0.016
    constants[9] = 0.02
    constants[10] = 5
    constants[11] = 2
    constants[12] = 2
    constants[13] = 0.035
    constants[14] = 3.5
    constants[15] = 0.0001
    constants[16] = 0.05
    constants[17] = 0.01
    constants[18] = 0.1
    constants[19] = 0.001
    constants[20] = 0.1
    constants[21] = 0.01
    constants[22] = 0.02
    constants[23] = 0.01
    constants[24] = 0.1
    constants[25] = 0.01
    constants[26] = 0.001
    states[0] = 0
    states[1] = 0.01
    states[2] = 1.95
    states[3] = 0
    constants[27] = 0.000001
    states[4] = 0.05
    constants[28] = 0.01
    states[5] = 0
    constants[29] = 0.018
    states[6] = 5
    constants[30] = 0.0001
    states[7] = 1
    states[8] = 0
    states[9] = 5
    constants[31] = 0.2
    states[10] = 0.01
    constants[32] = constants[3]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = constants[1]*states[5]*states[2]+constants[2]*states[8]*states[2]+constants[0]*states[0]*states[2]
    rates[2] = constants[32]-algebraic[1]
    algebraic[6] = constants[15]*states[3]
    algebraic[3] = constants[8]*states[3]
    rates[3] = (algebraic[1]+constants[27]+algebraic[6])-(constants[32]+algebraic[3])
    algebraic[0] = constants[4]*states[0]*states[1]
    algebraic[2] = constants[5]*states[0]
    algebraic[5] = constants[6]*states[3]+constants[7]
    algebraic[9] = constants[9]*states[1]
    rates[1] = (algebraic[5]+algebraic[2])-(algebraic[0]+algebraic[9])
    algebraic[4] = constants[23]/((1.00000+constants[24]*states[9])*1.00000)
    algebraic[10] = constants[26]*states[10]
    algebraic[7] = constants[25]*states[4]
    rates[4] = (constants[28]+algebraic[4]+algebraic[10])-(constants[32]+algebraic[7])
    rates[10] = algebraic[1]-algebraic[10]
    algebraic[13] = constants[13]*states[7]
    algebraic[12] = constants[12]*states[6]*states[0]
    rates[7] = algebraic[12]-algebraic[13]
    algebraic[14] = constants[18]*(power(states[0], 2.00000))
    rates[0] = (algebraic[0]+algebraic[13])-(algebraic[2]+algebraic[12]+algebraic[14])
    algebraic[18] = constants[17]*states[8]
    algebraic[8] = constants[10]*states[5]
    algebraic[16] = constants[16]*states[6]*states[5]
    algebraic[15] = constants[14]*states[9]*states[5]
    rates[5] = (constants[29]+algebraic[18])-(algebraic[8]+algebraic[16]+algebraic[15])
    rates[8] = algebraic[16]-algebraic[18]
    algebraic[19] = constants[21]/((1.00000+constants[22]*states[4])*1.00000)
    algebraic[17] = constants[20]*states[9]
    rates[9] = (constants[31]+algebraic[19])-(algebraic[15]+algebraic[17])
    algebraic[11] = constants[11]*states[6]*states[0]
    algebraic[20] = constants[19]*states[6]
    rates[6] = (constants[30]+algebraic[13]+algebraic[18])-(algebraic[11]+algebraic[12]+algebraic[16]+algebraic[20])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = constants[1]*states[5]*states[2]+constants[2]*states[8]*states[2]+constants[0]*states[0]*states[2]
    algebraic[6] = constants[15]*states[3]
    algebraic[3] = constants[8]*states[3]
    algebraic[0] = constants[4]*states[0]*states[1]
    algebraic[2] = constants[5]*states[0]
    algebraic[5] = constants[6]*states[3]+constants[7]
    algebraic[9] = constants[9]*states[1]
    algebraic[4] = constants[23]/((1.00000+constants[24]*states[9])*1.00000)
    algebraic[10] = constants[26]*states[10]
    algebraic[7] = constants[25]*states[4]
    algebraic[13] = constants[13]*states[7]
    algebraic[12] = constants[12]*states[6]*states[0]
    algebraic[14] = constants[18]*(power(states[0], 2.00000))
    algebraic[18] = constants[17]*states[8]
    algebraic[8] = constants[10]*states[5]
    algebraic[16] = constants[16]*states[6]*states[5]
    algebraic[15] = constants[14]*states[9]*states[5]
    algebraic[19] = constants[21]/((1.00000+constants[22]*states[4])*1.00000)
    algebraic[17] = constants[20]*states[9]
    algebraic[11] = constants[11]*states[6]*states[0]
    algebraic[20] = constants[19]*states[6]
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