# Size of variable arrays:
sizeAlgebraic = 34
sizeStates = 16
sizeConstants = 47
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_algebraic[0] = "v1 in component v1 (flux)"
    legend_constants[0] = "Vmax in component v1 (flux)"
    legend_constants[1] = "Km in component v1 (umol_per_g)"
    legend_states[0] = "GL in component GL (umol_per_g)"
    legend_algebraic[1] = "v2 in component v2 (flux)"
    legend_constants[2] = "Vmax in component v2 (flux)"
    legend_constants[3] = "Km in component v2 (umol_per_g)"
    legend_states[1] = "O2 in component O2 (umol_per_g)"
    legend_algebraic[15] = "v3 in component v3 (flux)"
    legend_constants[4] = "Vmax in component v3 (first_order_rate_constant)"
    legend_algebraic[14] = "k3 in component v3 (first_order_rate_constant)"
    legend_states[2] = "FA in component FA (umol_per_g)"
    legend_algebraic[13] = "PS in component PS (dimensionless)"
    legend_constants[5] = "PS0 in component PS0 (dimensionless)"
    legend_algebraic[30] = "v4 in component v4 (flux)"
    legend_constants[6] = "Vmax in component v4 (first_order_rate_constant)"
    legend_algebraic[29] = "k4 in component v4 (first_order_rate_constant)"
    legend_constants[7] = "epsilon in component v4 (dimensionless)"
    legend_algebraic[4] = "RS in component RS (dimensionless)"
    legend_constants[8] = "RS0 in component RS0 (dimensionless)"
    legend_algebraic[28] = "AF in component AF (dimensionless)"
    legend_constants[9] = "AF0 in component AF0 (dimensionless)"
    legend_algebraic[17] = "v5 in component v5 (flux)"
    legend_constants[10] = "Vmax in component v5 (first_order_rate_constant)"
    legend_algebraic[16] = "k5 in component v5 (first_order_rate_constant)"
    legend_constants[11] = "epsilon in component v5 (dimensionless)"
    legend_states[3] = "GP in component GP (umol_per_g)"
    legend_algebraic[19] = "v6 in component v6 (flux)"
    legend_constants[12] = "Vmax in component v6 (first_order_rate_constant)"
    legend_algebraic[18] = "k6 in component v6 (first_order_rate_constant)"
    legend_constants[13] = "epsilon in component v6 (dimensionless)"
    legend_algebraic[3] = "CS in component CS (dimensionless)"
    legend_constants[14] = "CS0 in component CS0 (dimensionless)"
    legend_algebraic[23] = "v7 in component v7 (flux)"
    legend_constants[15] = "Vmax in component v7 (first_order_rate_constant)"
    legend_algebraic[21] = "k7 in component v7 (first_order_rate_constant)"
    legend_constants[16] = "epsilon in component v7 (dimensionless)"
    legend_states[4] = "GY in component GY (umol_per_g)"
    legend_algebraic[2] = "v8 in component v8 (flux)"
    legend_constants[17] = "Vmax in component v8 (first_order_rate_constant)"
    legend_constants[45] = "k8 in component v8 (first_order_rate_constant)"
    legend_states[5] = "TG in component TG (umol_per_g)"
    legend_algebraic[6] = "v9 in component v9 (flux)"
    legend_constants[18] = "Vmax in component v9 (first_order_rate_constant)"
    legend_algebraic[5] = "k9 in component v9 (first_order_rate_constant)"
    legend_states[6] = "PY in component PY (umol_per_g)"
    legend_algebraic[32] = "v10 in component v10 (flux)"
    legend_constants[19] = "Vmax in component v10 (first_order_rate_constant)"
    legend_algebraic[31] = "k10 in component v10 (first_order_rate_constant)"
    legend_constants[20] = "epsilon in component v10 (dimensionless)"
    legend_algebraic[8] = "v11 in component v11 (flux)"
    legend_constants[21] = "Vmax in component v11 (first_order_rate_constant)"
    legend_algebraic[7] = "k11 in component v11 (first_order_rate_constant)"
    legend_states[7] = "LA in component LA (umol_per_g)"
    legend_algebraic[22] = "v12 in component v12 (flux)"
    legend_constants[22] = "Vmax in component v12 (first_order_rate_constant)"
    legend_algebraic[20] = "k12 in component v12 (first_order_rate_constant)"
    legend_constants[23] = "epsilon in component v12 (dimensionless)"
    legend_states[8] = "AC in component AC (umol_per_g)"
    legend_algebraic[25] = "v13 in component v13 (flux)"
    legend_constants[24] = "Vmax in component v13 (first_order_rate_constant)"
    legend_algebraic[24] = "k13 in component v13 (first_order_rate_constant)"
    legend_states[9] = "CR in component CR (umol_per_g)"
    legend_algebraic[27] = "v14 in component v14 (flux)"
    legend_constants[25] = "Vmax in component v14 (first_order_rate_constant)"
    legend_algebraic[26] = "k14 in component v14 (first_order_rate_constant)"
    legend_states[10] = "PC in component PC (umol_per_g)"
    legend_algebraic[10] = "v15 in component v15 (flux)"
    legend_constants[26] = "Vmax in component v15 (first_order_rate_constant)"
    legend_algebraic[9] = "k15 in component v15 (first_order_rate_constant)"
    legend_constants[27] = "epsilon in component v15 (dimensionless)"
    legend_states[11] = "CoA_pool in component CoA_pool (umol_per_g)"
    legend_states[12] = "FC in component FC (umol_per_g)"
    legend_constants[28] = "FC0 in component FC0 (umol_per_g)"
    legend_algebraic[12] = "v16 in component v16 (flux)"
    legend_constants[29] = "Vmax in component v16 (first_order_rate_constant)"
    legend_algebraic[11] = "k16 in component v16 (first_order_rate_constant)"
    legend_constants[30] = "epsilon in component v16 (dimensionless)"
    legend_constants[46] = "v17 in component v17 (flux)"
    legend_constants[31] = "Vmax in component v17 (first_order_rate_constant)"
    legend_constants[44] = "k17 in component v17 (first_order_rate_constant)"
    legend_constants[32] = "ATP in component ATP (umol_per_g)"
    legend_constants[33] = "aGL in component GL (umol_per_ml)"
    legend_constants[34] = "sigmaGL in component GL (g_per_ml)"
    legend_algebraic[33] = "F in component model_parameters (ml_per_g_min)"
    legend_constants[35] = "aFA in component FA (umol_per_ml)"
    legend_constants[36] = "sigmaFA in component FA (g_per_ml)"
    legend_constants[37] = "aLA in component LA (umol_per_ml)"
    legend_constants[38] = "sigmaLA in component LA (g_per_ml)"
    legend_constants[39] = "aO2 in component O2 (umol_per_ml)"
    legend_constants[40] = "sigmaO2 in component O2 (g_per_ml)"
    legend_states[13] = "CO2 in component CO2 (umol_per_g)"
    legend_constants[41] = "aCO2 in component CO2 (umol_per_ml)"
    legend_constants[42] = "sigmaCO2 in component CO2 (g_per_ml)"
    legend_states[14] = "NAD in component NAD (umol_per_g)"
    legend_states[15] = "ADP in component ADP (umol_per_g)"
    legend_constants[43] = "NADH in component NADH (umol_per_g)"
    legend_rates[0] = "d/dt GL in component GL (umol_per_g)"
    legend_rates[2] = "d/dt FA in component FA (umol_per_g)"
    legend_rates[3] = "d/dt GP in component GP (umol_per_g)"
    legend_rates[4] = "d/dt GY in component GY (umol_per_g)"
    legend_rates[5] = "d/dt TG in component TG (umol_per_g)"
    legend_rates[6] = "d/dt PY in component PY (umol_per_g)"
    legend_rates[7] = "d/dt LA in component LA (umol_per_g)"
    legend_rates[8] = "d/dt AC in component AC (umol_per_g)"
    legend_rates[12] = "d/dt FC in component FC (umol_per_g)"
    legend_rates[11] = "d/dt CoA_pool in component CoA_pool (umol_per_g)"
    legend_rates[1] = "d/dt O2 in component O2 (umol_per_g)"
    legend_rates[13] = "d/dt CO2 in component CO2 (umol_per_g)"
    legend_rates[14] = "d/dt NAD in component NAD (umol_per_g)"
    legend_rates[15] = "d/dt ADP in component ADP (umol_per_g)"
    legend_rates[10] = "d/dt PC in component PC (umol_per_g)"
    legend_rates[9] = "d/dt CR in component CR (umol_per_g)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 5.90
    constants[1] = 0.05
    states[0] = 0.998
    constants[2] = 67.6
    constants[3] = 0.01
    states[1] = 0.963
    constants[4] = 4.90
    states[2] = 0.021
    constants[5] = 0.2
    constants[6] = 21.3
    constants[7] = 0.6
    constants[8] = 0.111
    constants[9] = 0.523
    constants[10] = 2.82
    constants[11] = 0.254
    states[3] = 0.171
    constants[12] = 3.14
    constants[13] = 0.5
    constants[14] = 1.0
    constants[15] = 0.0162
    constants[16] = 0.5
    states[4] = 33.0
    constants[17] = 0.005
    states[5] = 3.96
    constants[18] = 1.8
    states[6] = 0.20
    constants[19] = 12.6
    constants[20] = 0.98
    constants[21] = 0.96
    states[7] = 1.98
    constants[22] = 695.7
    constants[23] = 0.75
    states[8] = 0.0046
    constants[24] = 0.455
    states[9] = 3.5
    constants[25] = 0.455
    states[10] = 8.80
    constants[26] = 626.1
    constants[27] = 0.669
    states[11] = 0.043
    states[12] = 0.0088
    constants[28] = 0.0088
    constants[29] = 67.0
    constants[30] = 0.775
    constants[31] = 7.76
    constants[32] = 4.5
    constants[33] = 4.0
    constants[34] = 3.76
    constants[35] = 0.5
    constants[36] = 13.2
    constants[37] = 1.8
    constants[38] = 0.51
    constants[39] = 6.53
    constants[40] = 1.0
    states[13] = 20.0
    constants[41] = 15.5
    constants[42] = 1.0
    states[14] = 1.81
    states[15] = 0.90
    constants[43] = 0.19
    constants[44] = constants[31]
    constants[45] = constants[17]
    constants[46] = constants[44]*constants[32]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[4] = constants[43]/states[14]
    algebraic[9] = constants[26]*(constants[27]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[27])*((power(states[12], -1.00000))/(power(constants[28], -1.00000)+power(states[12], -1.00000))))
    algebraic[10] = algebraic[9]*states[11]
    algebraic[11] = constants[29]*(constants[30]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[30])*((power(states[12], -1.00000))/(power(constants[28], -1.00000)+power(states[12], -1.00000))))
    algebraic[12] = algebraic[11]*states[8]
    rates[11] = algebraic[12]-algebraic[10]
    algebraic[13] = states[15]/constants[32]
    algebraic[14] = constants[4]*((power(algebraic[13], -1.00000))/(power(constants[5], -1.00000)+power(algebraic[13], -1.00000)))
    algebraic[15] = algebraic[14]*states[2]
    algebraic[2] = constants[45]*states[5]
    rates[5] = (1.00000/3.00000)*algebraic[15]-algebraic[2]
    algebraic[0] = constants[0]*(states[0]/(constants[1]+states[0]))
    algebraic[16] = constants[10]*(constants[11]*(algebraic[13]/(constants[5]+algebraic[13]))+(1.00000-constants[11])*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000))))
    algebraic[17] = algebraic[16]*states[3]
    algebraic[3] = states[9]/states[10]
    algebraic[18] = constants[12]*(constants[13]*((power(algebraic[13], -1.00000))/(power(constants[5], -1.00000)+power(algebraic[13], -1.00000)))+(1.00000-constants[13])*(power((power(algebraic[3], -1.00000))/(power(constants[14], -1.00000)+power(algebraic[3], -1.00000)), 2.00000)))
    algebraic[19] = algebraic[18]*states[3]
    algebraic[21] = constants[15]*(constants[16]*(algebraic[13]/(constants[5]+algebraic[13]))+(1.00000-constants[16])*(power(algebraic[3]/(constants[14]+algebraic[3]), 2.00000)))
    algebraic[23] = algebraic[21]*states[4]
    rates[3] = (algebraic[0]+algebraic[23])-(algebraic[17]+algebraic[19])
    rates[4] = algebraic[19]-algebraic[23]
    algebraic[24] = constants[24]*((power(algebraic[13], -1.00000))/(power(constants[5], -1.00000)+power(algebraic[13], -1.00000)))
    algebraic[25] = algebraic[24]*states[9]
    algebraic[26] = constants[25]*(algebraic[13]/(constants[5]+algebraic[13]))
    algebraic[27] = algebraic[26]*states[10]
    rates[10] = algebraic[25]-algebraic[27]
    rates[9] = algebraic[27]-algebraic[25]
    algebraic[1] = constants[2]*(states[1]/(constants[3]+states[1]))
    algebraic[28] = states[8]/states[12]
    algebraic[29] = constants[6]*(constants[7]*((power(algebraic[28], -1.00000))/(power(constants[9], -1.00000)+power(algebraic[28], -1.00000)))+(1.00000-constants[7])*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000))))
    algebraic[30] = algebraic[29]*states[2]
    algebraic[20] = constants[22]*(constants[23]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[23])*(algebraic[13]/(constants[5]+algebraic[13])))
    algebraic[22] = algebraic[20]*states[8]
    rates[15] = (algebraic[0]+algebraic[19]+2.00000*algebraic[30]+2.00000*algebraic[15]+algebraic[25]+constants[46])-(3.00000*algebraic[17]+algebraic[22]+6.00000*algebraic[1]+algebraic[27])
    algebraic[5] = constants[18]*(algebraic[4]/(constants[8]+algebraic[4]))
    algebraic[6] = algebraic[5]*states[6]
    algebraic[31] = constants[19]*(constants[20]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[20])*((power(algebraic[28], -1.00000))/(power(constants[9], -1.00000)+power(algebraic[28], -1.00000))))
    algebraic[32] = algebraic[31]*states[6]
    algebraic[7] = constants[21]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))
    algebraic[8] = algebraic[7]*states[7]
    rates[6] = (2.00000*algebraic[17]+algebraic[8])-(algebraic[6]+algebraic[32])
    rates[8] = (algebraic[32]+algebraic[30])-(algebraic[22]+algebraic[12])
    rates[12] = (algebraic[22]+algebraic[10]+algebraic[12])-(algebraic[32]+algebraic[30])
    rates[14] = (algebraic[6]+2.00000*algebraic[1])-(2.00000*algebraic[17]+algebraic[32]+algebraic[8]+(11.0000/3.00000)*algebraic[22]+(35.0000/3.00000)*algebraic[30])
    algebraic[33] = custom_piecewise([greater(voi , 0.00000) & less(voi , 5.00000), 1.00000 , True, 0.400000])
    rates[0] = algebraic[33]*(constants[33]-constants[34]*states[0])-algebraic[0]
    rates[2] = (3.00000*algebraic[2]+algebraic[33]*(constants[35]-constants[36]*states[2]))-(algebraic[15]+algebraic[30])
    rates[7] = (algebraic[6]+algebraic[33]*(constants[37]-constants[38]*states[7]))-algebraic[8]
    rates[1] = algebraic[33]*(constants[39]-constants[40]*states[1])-algebraic[1]
    rates[13] = algebraic[32]+2.00000*algebraic[22]+algebraic[33]*(constants[41]-constants[42]*states[13])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[4] = constants[43]/states[14]
    algebraic[9] = constants[26]*(constants[27]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[27])*((power(states[12], -1.00000))/(power(constants[28], -1.00000)+power(states[12], -1.00000))))
    algebraic[10] = algebraic[9]*states[11]
    algebraic[11] = constants[29]*(constants[30]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[30])*((power(states[12], -1.00000))/(power(constants[28], -1.00000)+power(states[12], -1.00000))))
    algebraic[12] = algebraic[11]*states[8]
    algebraic[13] = states[15]/constants[32]
    algebraic[14] = constants[4]*((power(algebraic[13], -1.00000))/(power(constants[5], -1.00000)+power(algebraic[13], -1.00000)))
    algebraic[15] = algebraic[14]*states[2]
    algebraic[2] = constants[45]*states[5]
    algebraic[0] = constants[0]*(states[0]/(constants[1]+states[0]))
    algebraic[16] = constants[10]*(constants[11]*(algebraic[13]/(constants[5]+algebraic[13]))+(1.00000-constants[11])*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000))))
    algebraic[17] = algebraic[16]*states[3]
    algebraic[3] = states[9]/states[10]
    algebraic[18] = constants[12]*(constants[13]*((power(algebraic[13], -1.00000))/(power(constants[5], -1.00000)+power(algebraic[13], -1.00000)))+(1.00000-constants[13])*(power((power(algebraic[3], -1.00000))/(power(constants[14], -1.00000)+power(algebraic[3], -1.00000)), 2.00000)))
    algebraic[19] = algebraic[18]*states[3]
    algebraic[21] = constants[15]*(constants[16]*(algebraic[13]/(constants[5]+algebraic[13]))+(1.00000-constants[16])*(power(algebraic[3]/(constants[14]+algebraic[3]), 2.00000)))
    algebraic[23] = algebraic[21]*states[4]
    algebraic[24] = constants[24]*((power(algebraic[13], -1.00000))/(power(constants[5], -1.00000)+power(algebraic[13], -1.00000)))
    algebraic[25] = algebraic[24]*states[9]
    algebraic[26] = constants[25]*(algebraic[13]/(constants[5]+algebraic[13]))
    algebraic[27] = algebraic[26]*states[10]
    algebraic[1] = constants[2]*(states[1]/(constants[3]+states[1]))
    algebraic[28] = states[8]/states[12]
    algebraic[29] = constants[6]*(constants[7]*((power(algebraic[28], -1.00000))/(power(constants[9], -1.00000)+power(algebraic[28], -1.00000)))+(1.00000-constants[7])*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000))))
    algebraic[30] = algebraic[29]*states[2]
    algebraic[20] = constants[22]*(constants[23]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[23])*(algebraic[13]/(constants[5]+algebraic[13])))
    algebraic[22] = algebraic[20]*states[8]
    algebraic[5] = constants[18]*(algebraic[4]/(constants[8]+algebraic[4]))
    algebraic[6] = algebraic[5]*states[6]
    algebraic[31] = constants[19]*(constants[20]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))+(1.00000-constants[20])*((power(algebraic[28], -1.00000))/(power(constants[9], -1.00000)+power(algebraic[28], -1.00000))))
    algebraic[32] = algebraic[31]*states[6]
    algebraic[7] = constants[21]*((power(algebraic[4], -1.00000))/(power(constants[8], -1.00000)+power(algebraic[4], -1.00000)))
    algebraic[8] = algebraic[7]*states[7]
    algebraic[33] = custom_piecewise([greater(voi , 0.00000) & less(voi , 5.00000), 1.00000 , True, 0.400000])
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

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