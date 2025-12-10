# Size of variable arrays:
sizeAlgebraic = 29
sizeStates = 15
sizeConstants = 46
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
    legend_constants[1] = "kn1 in component rate_constants (first_order_rate_constant)"
    legend_constants[2] = "k2 in component rate_constants (first_order_rate_constant)"
    legend_constants[3] = "kn2 in component rate_constants (first_order_rate_constant)"
    legend_constants[4] = "k2_ in component rate_constants (first_order_rate_constant)"
    legend_constants[5] = "k3_ in component rate_constants (first_order_rate_constant)"
    legend_constants[6] = "kn3_ in component rate_constants (first_order_rate_constant)"
    legend_constants[7] = "k4 in component rate_constants (first_order_rate_constant)"
    legend_constants[8] = "k5 in component rate_constants (first_order_rate_constant)"
    legend_constants[9] = "k6 in component rate_constants (first_order_rate_constant)"
    legend_constants[10] = "k7 in component rate_constants (first_order_rate_constant)"
    legend_constants[11] = "kn7 in component rate_constants (first_order_rate_constant)"
    legend_constants[12] = "k8 in component rate_constants (first_order_rate_constant)"
    legend_constants[13] = "kn8 in component rate_constants (first_order_rate_constant)"
    legend_constants[14] = "k8_ in component rate_constants (first_order_rate_constant)"
    legend_constants[15] = "k9 in component rate_constants (first_order_rate_constant)"
    legend_constants[16] = "kn9 in component rate_constants (first_order_rate_constant)"
    legend_constants[17] = "k9_ in component rate_constants (first_order_rate_constant)"
    legend_constants[18] = "k10 in component rate_constants (first_order_rate_constant)"
    legend_constants[19] = "k11 in component rate_constants (first_order_rate_constant)"
    legend_constants[20] = "k12 in component rate_constants (first_order_rate_constant)"
    legend_constants[21] = "k13 in component rate_constants (first_order_rate_constant)"
    legend_constants[22] = "k14 in component rate_constants (first_order_rate_constant)"
    legend_constants[23] = "k14_ in component rate_constants (first_order_rate_constant)"
    legend_constants[24] = "k15 in component rate_constants (first_order_rate_constant)"
    legend_constants[25] = "k16 in component rate_constants (first_order_rate_constant)"
    legend_constants[26] = "k17 in component rate_constants (first_order_rate_constant)"
    legend_constants[27] = "k18 in component rate_constants (first_order_rate_constant)"
    legend_constants[28] = "kn18 in component rate_constants (first_order_rate_constant)"
    legend_constants[29] = "k20 in component rate_constants (first_order_rate_constant)"
    legend_constants[30] = "k21 in component rate_constants (first_order_rate_constant)"
    legend_constants[31] = "k22 in component rate_constants (first_order_rate_constant)"
    legend_constants[32] = "k23 in component rate_constants (first_order_rate_constant)"
    legend_constants[33] = "kn23 in component rate_constants (first_order_rate_constant)"
    legend_constants[34] = "kPlk1 in component rate_constants (first_order_rate_constant)"
    legend_constants[35] = "kPlk1_ in component rate_constants (first_order_rate_constant)"
    legend_constants[36] = "kctak1 in component rate_constants (first_order_rate_constant)"
    legend_constants[37] = "kctak1_ in component rate_constants (first_order_rate_constant)"
    legend_constants[38] = "kex in component rate_constants (first_order_rate_constant)"
    legend_states[0] = "Chk1P in component Chk1P (dimensionless)"
    legend_algebraic[11] = "V_n1 in component V_n1 (first_order_rate_constant)"
    legend_algebraic[7] = "V_1 in component V_1 (first_order_rate_constant)"
    legend_states[1] = "Rad3 in component Rad3 (dimensionless)"
    legend_constants[41] = "V_4 in component V_4 (first_order_rate_constant)"
    legend_algebraic[1] = "V_5 in component V_5 (first_order_rate_constant)"
    legend_states[2] = "p53 in component p53 (dimensionless)"
    legend_constants[42] = "V_10 in component V_10 (first_order_rate_constant)"
    legend_algebraic[2] = "V_11 in component V_11 (first_order_rate_constant)"
    legend_states[3] = "preMPF in component preMPF (dimensionless)"
    legend_algebraic[12] = "V_14 in component V_14 (first_order_rate_constant)"
    legend_algebraic[8] = "V_n9 in component V_n9 (first_order_rate_constant)"
    legend_algebraic[3] = "V_9 in component V_9 (first_order_rate_constant)"
    legend_states[4] = "MPF in component MPF (dimensionless)"
    legend_algebraic[26] = "V_n23 in component V_n23 (first_order_rate_constant)"
    legend_algebraic[21] = "V_23 in component V_23 (first_order_rate_constant)"
    legend_algebraic[13] = "V_15 in component V_15 (first_order_rate_constant)"
    legend_states[5] = "p21 in component p21 (dimensionless)"
    legend_constants[43] = "V_21 in component V_21 (first_order_rate_constant)"
    legend_algebraic[16] = "V_20 in component V_20 (first_order_rate_constant)"
    legend_algebraic[18] = "V_22 in component V_22 (first_order_rate_constant)"
    legend_states[6] = "p21_MPF in component p21_MPF (dimensionless)"
    legend_states[7] = "iCdc25 in component iCdc25 (dimensionless)"
    legend_algebraic[27] = "V_n7 in component V_n7 (first_order_rate_constant)"
    legend_algebraic[19] = "V_n3_ in component V_n3_ (first_order_rate_constant)"
    legend_algebraic[22] = "V_7 in component V_7 (first_order_rate_constant)"
    legend_algebraic[14] = "V_2_ in component V_2_ (first_order_rate_constant)"
    legend_constants[39] = "V_in in component V_in (first_order_rate_constant)"
    legend_states[8] = "iCdc25Ps216 in component iCdc25Ps216 (dimensionless)"
    legend_algebraic[23] = "V_n18 in component V_n18 (first_order_rate_constant)"
    legend_algebraic[20] = "V_18 in component V_18 (first_order_rate_constant)"
    legend_algebraic[17] = "V_3_ in component V_3_ (first_order_rate_constant)"
    legend_states[9] = "iCdc25Ps216_protein1433 in component iCdc25Ps216_protein1433 (dimensionless)"
    legend_algebraic[24] = "V_ex in component V_ex (first_order_rate_constant)"
    legend_states[10] = "aCdc25 in component aCdc25 (dimensionless)"
    legend_algebraic[9] = "V_n2 in component V_n2 (first_order_rate_constant)"
    legend_algebraic[4] = "V_2 in component V_2 (first_order_rate_constant)"
    legend_states[11] = "aCdc25Ps216 in component aCdc25Ps216 (dimensionless)"
    legend_states[12] = "protein1433 in component protein1433 (dimensionless)"
    legend_algebraic[25] = "V_6 in component V_6 (first_order_rate_constant)"
    legend_constants[44] = "V_13 in component V_13 (first_order_rate_constant)"
    legend_algebraic[28] = "V_12 in component V_12 (first_order_rate_constant)"
    legend_states[13] = "Wee1 in component Wee1 (dimensionless)"
    legend_constants[45] = "V_16 in component V_16 (first_order_rate_constant)"
    legend_algebraic[10] = "V_n8 in component V_n8 (first_order_rate_constant)"
    legend_algebraic[5] = "V_8 in component V_8 (first_order_rate_constant)"
    legend_states[14] = "Wee1P in component Wee1P (dimensionless)"
    legend_algebraic[15] = "V_17 in component V_17 (first_order_rate_constant)"
    legend_algebraic[6] = "aCdc25_T in component aCdc25_T (dimensionless)"
    legend_algebraic[0] = "Chk1 in component V_1 (dimensionless)"
    legend_constants[40] = "Chk1_T in component V_1 (dimensionless)"
    legend_rates[0] = "d/dt Chk1P in component Chk1P (dimensionless)"
    legend_rates[1] = "d/dt Rad3 in component Rad3 (dimensionless)"
    legend_rates[2] = "d/dt p53 in component p53 (dimensionless)"
    legend_rates[3] = "d/dt preMPF in component preMPF (dimensionless)"
    legend_rates[4] = "d/dt MPF in component MPF (dimensionless)"
    legend_rates[5] = "d/dt p21 in component p21 (dimensionless)"
    legend_rates[6] = "d/dt p21_MPF in component p21_MPF (dimensionless)"
    legend_rates[7] = "d/dt iCdc25 in component iCdc25 (dimensionless)"
    legend_rates[8] = "d/dt iCdc25Ps216 in component iCdc25Ps216 (dimensionless)"
    legend_rates[9] = "d/dt iCdc25Ps216_protein1433 in component iCdc25Ps216_protein1433 (dimensionless)"
    legend_rates[10] = "d/dt aCdc25 in component aCdc25 (dimensionless)"
    legend_rates[11] = "d/dt aCdc25Ps216 in component aCdc25Ps216 (dimensionless)"
    legend_rates[12] = "d/dt protein1433 in component protein1433 (dimensionless)"
    legend_rates[13] = "d/dt Wee1 in component Wee1 (dimensionless)"
    legend_rates[14] = "d/dt Wee1P in component Wee1P (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1
    constants[1] = 10
    constants[2] = 0.1
    constants[3] = 0.01
    constants[4] = 0.1
    constants[5] = 100
    constants[6] = 0
    constants[7] = 0
    constants[8] = 1
    constants[9] = 0.01
    constants[10] = 1
    constants[11] = 0.01
    constants[12] = 0.1
    constants[13] = 0.01
    constants[14] = 0
    constants[15] = 1
    constants[16] = 1
    constants[17] = 0
    constants[18] = 0
    constants[19] = 1
    constants[20] = 0.1
    constants[21] = 1
    constants[22] = 0.0005
    constants[23] = 1
    constants[24] = 0.01
    constants[25] = 0.0002
    constants[26] = 0.1
    constants[27] = 1
    constants[28] = 0.01
    constants[29] = 0.1
    constants[30] = 0.01
    constants[31] = 0.1
    constants[32] = 0.1
    constants[33] = 1
    constants[34] = 0
    constants[35] = 0
    constants[36] = 0
    constants[37] = 0
    constants[38] = 1
    states[0] = 0.000001
    states[1] = 0
    states[2] = 0
    states[3] = 0.000001
    states[4] = 0.00000001
    states[5] = 0
    states[6] = 0
    states[7] = 0.000001
    constants[39] = 0.00001
    states[8] = 0.00002
    states[9] = 0.03
    states[10] = 0.000001
    states[11] = 0
    states[12] = 2
    states[13] = 0.001
    states[14] = 0
    constants[40] = 1
    constants[41] = constants[7]
    constants[42] = constants[18]
    constants[43] = constants[30]
    constants[44] = constants[21]
    constants[45] = constants[25]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = constants[8]*states[1]
    rates[1] = constants[41]-algebraic[1]
    algebraic[2] = constants[19]*states[2]
    rates[2] = constants[42]-algebraic[2]
    algebraic[10] = constants[13]*states[14]
    algebraic[5] = constants[12]*states[4]*states[13]+constants[14]*states[13]
    rates[13] = (constants[45]+algebraic[10])-algebraic[5]
    algebraic[11] = constants[1]*states[0]
    algebraic[0] = constants[40]-states[0]
    algebraic[7] = constants[0]*algebraic[0]*states[1]
    rates[0] = algebraic[7]-algebraic[11]
    algebraic[12] = (constants[22]/(1.00000+constants[23]*states[2]))*1.00000
    algebraic[8] = constants[16]*states[4]*states[13]
    algebraic[3] = constants[15]*(states[10]+states[11])*states[3]+constants[17]*states[3]
    rates[3] = (algebraic[12]+algebraic[8])-algebraic[3]
    algebraic[15] = constants[26]*states[14]
    rates[14] = algebraic[5]-(algebraic[10]+algebraic[15])
    algebraic[14] = constants[4]*states[0]*states[7]+constants[37]*states[7]
    algebraic[23] = constants[28]*states[11]
    algebraic[20] = constants[27]*states[4]*states[8]+constants[35]*states[8]
    algebraic[17] = constants[5]*states[8]*states[12]
    rates[8] = (algebraic[14]+algebraic[23])-(algebraic[20]+algebraic[17])
    algebraic[19] = constants[6]*states[9]
    algebraic[24] = constants[38]*states[9]
    rates[9] = algebraic[17]-(algebraic[19]+algebraic[24])
    algebraic[9] = constants[3]*states[11]
    algebraic[4] = constants[2]*states[0]*states[10]+constants[36]*states[10]
    rates[11] = (algebraic[4]+algebraic[20])-(algebraic[9]+algebraic[23])
    algebraic[26] = constants[33]*states[6]
    algebraic[21] = constants[32]*states[4]*states[5]
    algebraic[13] = constants[24]*(power(states[4], 2.00000))
    rates[4] = (algebraic[3]+algebraic[26])-(algebraic[8]+algebraic[21]+algebraic[13])
    algebraic[16] = constants[29]*states[2]
    algebraic[18] = constants[31]*states[5]
    rates[5] = (algebraic[16]+constants[43]+algebraic[26])-(algebraic[18]+algebraic[21])
    rates[6] = algebraic[21]-algebraic[26]
    algebraic[27] = constants[11]*states[10]
    algebraic[22] = constants[10]*states[4]*states[7]+constants[34]*states[7]
    rates[7] = (algebraic[27]+algebraic[19]+constants[39])-(algebraic[22]+algebraic[14])
    rates[10] = (algebraic[22]+algebraic[9])-(algebraic[27]+algebraic[4])
    algebraic[25] = constants[9]*states[2]
    algebraic[28] = constants[20]*states[12]
    rates[12] = (algebraic[25]+constants[44]+algebraic[19])-(algebraic[17]+algebraic[28])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = constants[8]*states[1]
    algebraic[2] = constants[19]*states[2]
    algebraic[10] = constants[13]*states[14]
    algebraic[5] = constants[12]*states[4]*states[13]+constants[14]*states[13]
    algebraic[11] = constants[1]*states[0]
    algebraic[0] = constants[40]-states[0]
    algebraic[7] = constants[0]*algebraic[0]*states[1]
    algebraic[12] = (constants[22]/(1.00000+constants[23]*states[2]))*1.00000
    algebraic[8] = constants[16]*states[4]*states[13]
    algebraic[3] = constants[15]*(states[10]+states[11])*states[3]+constants[17]*states[3]
    algebraic[15] = constants[26]*states[14]
    algebraic[14] = constants[4]*states[0]*states[7]+constants[37]*states[7]
    algebraic[23] = constants[28]*states[11]
    algebraic[20] = constants[27]*states[4]*states[8]+constants[35]*states[8]
    algebraic[17] = constants[5]*states[8]*states[12]
    algebraic[19] = constants[6]*states[9]
    algebraic[24] = constants[38]*states[9]
    algebraic[9] = constants[3]*states[11]
    algebraic[4] = constants[2]*states[0]*states[10]+constants[36]*states[10]
    algebraic[26] = constants[33]*states[6]
    algebraic[21] = constants[32]*states[4]*states[5]
    algebraic[13] = constants[24]*(power(states[4], 2.00000))
    algebraic[16] = constants[29]*states[2]
    algebraic[18] = constants[31]*states[5]
    algebraic[27] = constants[11]*states[10]
    algebraic[22] = constants[10]*states[4]*states[7]+constants[34]*states[7]
    algebraic[25] = constants[9]*states[2]
    algebraic[28] = constants[20]*states[12]
    algebraic[6] = states[10]+states[11]
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