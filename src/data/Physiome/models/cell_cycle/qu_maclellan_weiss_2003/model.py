# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 12
sizeConstants = 34
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
    legend_constants[1] = "k2 in component rate_constants (first_order_rate_constant)"
    legend_constants[2] = "k3 in component rate_constants (first_order_rate_constant)"
    legend_constants[3] = "k4 in component rate_constants (first_order_rate_constant)"
    legend_constants[4] = "k5 in component rate_constants (first_order_rate_constant)"
    legend_constants[5] = "k6 in component rate_constants (first_order_rate_constant)"
    legend_constants[6] = "k7 in component rate_constants (first_order_rate_constant)"
    legend_constants[7] = "k8 in component rate_constants (first_order_rate_constant)"
    legend_constants[8] = "k9 in component rate_constants (first_order_rate_constant)"
    legend_constants[9] = "k10 in component rate_constants (first_order_rate_constant)"
    legend_constants[10] = "k11 in component rate_constants (first_order_rate_constant)"
    legend_constants[11] = "k12 in component rate_constants (first_order_rate_constant)"
    legend_constants[12] = "k13 in component rate_constants (first_order_rate_constant)"
    legend_constants[13] = "k14 in component rate_constants (first_order_rate_constant)"
    legend_constants[14] = "k15 in component rate_constants (first_order_rate_constant)"
    legend_constants[15] = "k16 in component rate_constants (first_order_rate_constant)"
    legend_constants[16] = "k2u in component rate_constants (first_order_rate_constant)"
    legend_constants[17] = "k7u in component rate_constants (first_order_rate_constant)"
    legend_constants[18] = "k16u in component rate_constants (first_order_rate_constant)"
    legend_constants[19] = "az in component rate_constants (first_order_rate_constant)"
    legend_constants[20] = "bz in component rate_constants (first_order_rate_constant)"
    legend_constants[21] = "cz in component rate_constants (first_order_rate_constant)"
    legend_constants[22] = "aw in component rate_constants (first_order_rate_constant)"
    legend_constants[23] = "bw in component rate_constants (first_order_rate_constant)"
    legend_constants[24] = "cw in component rate_constants (first_order_rate_constant)"
    legend_constants[25] = "ai in component rate_constants (first_order_rate_constant)"
    legend_constants[26] = "bi in component rate_constants (first_order_rate_constant)"
    legend_constants[27] = "ci in component rate_constants (first_order_rate_constant)"
    legend_states[0] = "cyclin_CDK_active in component cyclin_CDK_active (dimensionless)"
    legend_algebraic[0] = "kplus_z in component rate_constants (first_order_rate_constant)"
    legend_algebraic[1] = "kplus_w in component rate_constants (first_order_rate_constant)"
    legend_algebraic[2] = "kplus_i in component rate_constants (first_order_rate_constant)"
    legend_constants[31] = "kminus_z in component rate_constants (first_order_rate_constant)"
    legend_constants[32] = "kminus_w in component rate_constants (first_order_rate_constant)"
    legend_constants[33] = "kminus_i in component rate_constants (first_order_rate_constant)"
    legend_states[1] = "cyclin in component cyclin (dimensionless)"
    legend_states[2] = "cyclin_CDK_inactive in component cyclin_CDK_inactive (dimensionless)"
    legend_algebraic[3] = "CDK in component CDK (dimensionless)"
    legend_states[3] = "APC in component APC (dimensionless)"
    legend_states[4] = "wee1 in component wee1 (dimensionless)"
    legend_states[5] = "Cdc25_2P in component Cdc25_2P (dimensionless)"
    legend_states[6] = "CKI in component CKI (dimensionless)"
    legend_states[7] = "cyclin_CDK_CKI in component cyclin_CDK_CKI (dimensionless)"
    legend_states[8] = "cyclin_CDK_CKI_P in component cyclin_CDK_CKI_P (dimensionless)"
    legend_states[9] = "Cdc25 in component Cdc25 (dimensionless)"
    legend_states[10] = "Cdc25_P in component Cdc25_P (dimensionless)"
    legend_states[11] = "wee1_P in component wee1_P (dimensionless)"
    legend_constants[28] = "tau in component APC (minute)"
    legend_algebraic[4] = "h_x in component APC (dimensionless)"
    legend_constants[29] = "a in component APC (dimensionless)"
    legend_constants[30] = "CDK_T in component CDK (dimensionless)"
    legend_algebraic[5] = "cyclin_T in component cyclin_T (dimensionless)"
    legend_rates[1] = "d/dt cyclin in component cyclin (dimensionless)"
    legend_rates[2] = "d/dt cyclin_CDK_inactive in component cyclin_CDK_inactive (dimensionless)"
    legend_rates[0] = "d/dt cyclin_CDK_active in component cyclin_CDK_active (dimensionless)"
    legend_rates[9] = "d/dt Cdc25 in component Cdc25 (dimensionless)"
    legend_rates[10] = "d/dt Cdc25_P in component Cdc25_P (dimensionless)"
    legend_rates[5] = "d/dt Cdc25_2P in component Cdc25_2P (dimensionless)"
    legend_rates[4] = "d/dt wee1 in component wee1 (dimensionless)"
    legend_rates[11] = "d/dt wee1_P in component wee1_P (dimensionless)"
    legend_rates[3] = "d/dt APC in component APC (dimensionless)"
    legend_rates[6] = "d/dt CKI in component CKI (dimensionless)"
    legend_rates[7] = "d/dt cyclin_CDK_CKI in component cyclin_CDK_CKI (dimensionless)"
    legend_rates[8] = "d/dt cyclin_CDK_CKI_P in component cyclin_CDK_CKI_P (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 300
    constants[1] = 5
    constants[2] = 30
    constants[3] = 30
    constants[4] = 0.1
    constants[5] = 1
    constants[6] = 10
    constants[7] = 100
    constants[8] = 1
    constants[9] = 10
    constants[10] = 1
    constants[11] = 0
    constants[12] = 1
    constants[13] = 1
    constants[14] = 1
    constants[15] = 2
    constants[16] = 50
    constants[17] = 0
    constants[18] = 25
    constants[19] = 10
    constants[20] = 0.1
    constants[21] = 1
    constants[22] = 10
    constants[23] = 0.1
    constants[24] = 1
    constants[25] = 10
    constants[26] = 0.1
    constants[27] = 1
    states[0] = 0.1
    states[1] = 0
    states[2] = 0.1
    states[3] = 0
    states[4] = 0
    states[5] = 0
    states[6] = 0
    states[7] = 0
    states[8] = 0
    states[9] = 0
    states[10] = 0
    states[11] = 0
    constants[28] = 25
    constants[29] = 4
    constants[30] = 200
    constants[31] = constants[19]
    constants[32] = constants[22]
    constants[33] = constants[25]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (((((constants[4]*states[2]+states[5]*states[2]*1.00000)-constants[5]*states[0])-states[4]*states[0]*1.00000)-(constants[6]+constants[17]*states[3])*states[0])-constants[13]*states[0]*states[6])+constants[14]*states[7]+(constants[15]+constants[18]*states[3])*states[8]
    rates[6] = ((constants[11]-constants[12]*states[6])-constants[11]*states[0]*states[6])+constants[14]*states[7]
    algebraic[3] = ((((constants[30]-states[0])-states[2])-states[7])-states[8])/constants[30]
    rates[1] = ((constants[0]+constants[3]*states[2])-constants[2]*states[1]*algebraic[3])-(constants[1]+constants[16]*states[3])*states[1]
    rates[2] = (((constants[2]*states[1]*algebraic[3]+constants[5]*states[0]+states[4]*states[0]*1.00000)-constants[3]*states[2])-constants[4]*states[2])-states[5]*states[2]*1.00000
    algebraic[0] = constants[20]+constants[21]*states[0]
    rates[9] = ((constants[7]+constants[31]*states[10])-algebraic[0]*states[9])-constants[8]*states[9]
    rates[10] = (algebraic[0]*(states[9]-states[10])+constants[31]*(states[5]-states[10]))-constants[8]*states[10]
    rates[5] = (algebraic[0]*states[10]-constants[31]*states[5])-constants[8]*states[5]
    algebraic[1] = constants[23]+constants[24]*states[0]
    rates[4] = ((constants[9]+constants[32]*states[11])-algebraic[1]*states[4])-constants[10]*states[4]
    rates[11] = (algebraic[1]*states[4]-constants[32]*states[11])-constants[10]*states[11]
    algebraic[4] = (power(states[0], 2.00000))/(power(constants[29], 2.00000)+power(states[0], 2.00000))
    rates[3] = (algebraic[4]-states[3])/constants[28]
    algebraic[2] = constants[26]+constants[27]*states[0]
    rates[7] = ((constants[13]*states[0]*states[6]-constants[14]*states[7])+constants[33]*states[8])-algebraic[2]*states[7]
    rates[8] = (algebraic[2]*states[7]-constants[33]*states[8])-(constants[15]+constants[18]*states[3])*states[8]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = ((((constants[30]-states[0])-states[2])-states[7])-states[8])/constants[30]
    algebraic[0] = constants[20]+constants[21]*states[0]
    algebraic[1] = constants[23]+constants[24]*states[0]
    algebraic[4] = (power(states[0], 2.00000))/(power(constants[29], 2.00000)+power(states[0], 2.00000))
    algebraic[2] = constants[26]+constants[27]*states[0]
    algebraic[5] = states[1]+states[0]+states[2]
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