# Size of variable arrays:
sizeAlgebraic = 35
sizeStates = 29
sizeConstants = 45
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "BMK1 in component BMK1 (micromolar)"
    legend_algebraic[9] = "v10a in component v10 (flux)"
    legend_algebraic[10] = "v10b in component v10 (flux)"
    legend_states[1] = "MRNA in component MRNA (micromolar)"
    legend_algebraic[28] = "v18 in component v18 (flux)"
    legend_algebraic[30] = "v19 in component v19 (flux)"
    legend_states[2] = "NFATc in component NFATc (micromolar)"
    legend_algebraic[0] = "v1 in component v1 (flux)"
    legend_algebraic[31] = "v14 in component v14 (flux)"
    legend_algebraic[34] = "v15 in component v15 (flux)"
    legend_states[3] = "NFATn in component NFATn (micromolar)"
    legend_algebraic[19] = "v16 in component v16 (flux)"
    legend_algebraic[25] = "v17 in component v17 (flux)"
    legend_states[4] = "NFATpc in component NFATpc (micromolar)"
    legend_algebraic[1] = "v2 in component v2 (flux)"
    legend_algebraic[8] = "v9 in component v9 (flux)"
    legend_states[5] = "NFATpn in component NFATpn (micromolar)"
    legend_states[6] = "GSK3betac in component GSK3betac (micromolar)"
    legend_algebraic[2] = "v3 in component v3 (flux)"
    legend_algebraic[16] = "v12a in component v12 (flux)"
    legend_algebraic[18] = "v12b in component v12 (flux)"
    legend_algebraic[32] = "v15a in component v15 (flux)"
    legend_algebraic[33] = "v15b in component v15 (flux)"
    legend_states[7] = "GSK3betan in component GSK3betan (micromolar)"
    legend_algebraic[21] = "v17a in component v17 (flux)"
    legend_algebraic[23] = "v17b in component v17 (flux)"
    legend_states[8] = "CaNc in component CaNc (micromolar)"
    legend_algebraic[5] = "v6 in component v6 (flux)"
    legend_states[9] = "CaNc_star in component CaNc_star (micromolar)"
    legend_algebraic[3] = "v4 in component v4 (flux)"
    legend_algebraic[6] = "v7 in component v7 (flux)"
    legend_algebraic[12] = "v11a in component v11 (flux)"
    legend_algebraic[13] = "v11b in component v11 (flux)"
    legend_algebraic[22] = "v13a in component v13 (flux)"
    legend_algebraic[24] = "v13b in component v13 (flux)"
    legend_algebraic[26] = "v14a in component v14 (flux)"
    legend_algebraic[29] = "v14b in component v14 (flux)"
    legend_states[10] = "CaNn_star in component CaNn_star (micromolar)"
    legend_algebraic[15] = "v16a in component v16 (flux)"
    legend_algebraic[17] = "v16b in component v16 (flux)"
    legend_states[11] = "CaNn in component CaNn (micromolar)"
    legend_states[12] = "CaM in component CaM (micromolar)"
    legend_algebraic[4] = "v5 in component v5 (flux)"
    legend_states[13] = "CaMCa in component CaMCa (micromolar)"
    legend_states[14] = "MCIP in component MCIP (micromolar)"
    legend_algebraic[11] = "v10 in component v10 (flux)"
    legend_algebraic[14] = "v11 in component v11 (flux)"
    legend_states[15] = "MCIPp in component MCIPp (micromolar)"
    legend_algebraic[20] = "v12 in component v12 (flux)"
    legend_algebraic[27] = "v13 in component v13 (flux)"
    legend_states[16] = "MCIPpp in component MCIPpp (micromolar)"
    legend_algebraic[7] = "v8 in component v8 (flux)"
    legend_states[17] = "Comp1 in component Comp1 (micromolar)"
    legend_states[18] = "Comp2 in component Comp2 (micromolar)"
    legend_states[19] = "Comp3 in component Comp3 (micromolar)"
    legend_states[20] = "P1433 in component P1433 (micromolar)"
    legend_states[21] = "MCIP_BMK1 in component MCIP_BMK1 (micromolar)"
    legend_states[22] = "MCIPp_CaNc_star in component MCIPp_CaNc_star (micromolar)"
    legend_states[23] = "MCIPp_GSK3betac in component MCIPp_GSK3betac (micromolar)"
    legend_states[24] = "MCIPpp_CaNc_star in component MCIPpp_CaNc_star (micromolar)"
    legend_states[25] = "NFATpc_CaNc_star in component NFATpc_CaNc_star (micromolar)"
    legend_states[26] = "NFATc_GSK3betac in component NFATc_GSK3betac (micromolar)"
    legend_states[27] = "NFATpn_CaNn_star in component NFATpn_CaNn_star (micromolar)"
    legend_states[28] = "NFATn_GSK3betan in component NFATn_GSK3betan (micromolar)"
    legend_constants[0] = "k29 in component model_parameters (first_order_rate_constant)"
    legend_constants[1] = "k30 in component model_parameters (first_order_rate_constant)"
    legend_constants[2] = "k31 in component model_parameters (first_order_rate_constant)"
    legend_constants[3] = "k32 in component model_parameters (first_order_rate_constant)"
    legend_constants[4] = "k33 in component model_parameters (first_order_rate_constant)"
    legend_constants[5] = "k34 in component model_parameters (first_order_rate_constant)"
    legend_constants[6] = "Ca in component model_parameters (micromolar)"
    legend_constants[7] = "k1 in component model_parameters (fifth_order_rate_constant)"
    legend_constants[8] = "k2 in component model_parameters (first_order_rate_constant)"
    legend_constants[9] = "k3 in component model_parameters (second_order_rate_constant)"
    legend_constants[10] = "k4 in component model_parameters (first_order_rate_constant)"
    legend_constants[11] = "k5 in component model_parameters (second_order_rate_constant)"
    legend_constants[12] = "k6 in component model_parameters (first_order_rate_constant)"
    legend_constants[13] = "k19 in component model_parameters (second_order_rate_constant)"
    legend_constants[14] = "k20 in component model_parameters (first_order_rate_constant)"
    legend_constants[15] = "k27 in component model_parameters (second_order_rate_constant)"
    legend_constants[16] = "k28 in component model_parameters (first_order_rate_constant)"
    legend_constants[17] = "k7 in component model_parameters (second_order_rate_constant)"
    legend_constants[18] = "k8 in component model_parameters (first_order_rate_constant)"
    legend_constants[19] = "k9 in component model_parameters (first_order_rate_constant)"
    legend_constants[20] = "k10 in component model_parameters (second_order_rate_constant)"
    legend_constants[21] = "k11 in component model_parameters (first_order_rate_constant)"
    legend_constants[22] = "k12 in component model_parameters (first_order_rate_constant)"
    legend_constants[23] = "k13 in component model_parameters (second_order_rate_constant)"
    legend_constants[24] = "k14 in component model_parameters (first_order_rate_constant)"
    legend_constants[25] = "k15 in component model_parameters (first_order_rate_constant)"
    legend_constants[26] = "k16 in component model_parameters (second_order_rate_constant)"
    legend_constants[27] = "k17 in component model_parameters (first_order_rate_constant)"
    legend_constants[28] = "k18 in component model_parameters (first_order_rate_constant)"
    legend_constants[29] = "k21 in component model_parameters (second_order_rate_constant)"
    legend_constants[30] = "k22 in component model_parameters (first_order_rate_constant)"
    legend_constants[31] = "k23 in component model_parameters (first_order_rate_constant)"
    legend_constants[32] = "k24 in component model_parameters (second_order_rate_constant)"
    legend_constants[33] = "k25 in component model_parameters (first_order_rate_constant)"
    legend_constants[34] = "k26 in component model_parameters (first_order_rate_constant)"
    legend_constants[35] = "k35 in component model_parameters (second_order_rate_constant)"
    legend_constants[36] = "k36 in component model_parameters (first_order_rate_constant)"
    legend_constants[37] = "k37 in component model_parameters (first_order_rate_constant)"
    legend_constants[38] = "k38 in component model_parameters (second_order_rate_constant)"
    legend_constants[39] = "k39 in component model_parameters (first_order_rate_constant)"
    legend_constants[40] = "k40 in component model_parameters (first_order_rate_constant)"
    legend_constants[41] = "k41 in component model_parameters (first_order_rate_constant)"
    legend_constants[42] = "k42 in component model_parameters (first_order_rate_constant)"
    legend_constants[43] = "t_half in component model_parameters (minute)"
    legend_constants[44] = "k43 in component model_parameters (first_order_rate_constant)"
    legend_rates[0] = "d/dt BMK1 in component BMK1 (micromolar)"
    legend_rates[1] = "d/dt MRNA in component MRNA (micromolar)"
    legend_rates[2] = "d/dt NFATc in component NFATc (micromolar)"
    legend_rates[3] = "d/dt NFATn in component NFATn (micromolar)"
    legend_rates[4] = "d/dt NFATpc in component NFATpc (micromolar)"
    legend_rates[5] = "d/dt NFATpn in component NFATpn (micromolar)"
    legend_rates[6] = "d/dt GSK3betac in component GSK3betac (micromolar)"
    legend_rates[7] = "d/dt GSK3betan in component GSK3betan (micromolar)"
    legend_rates[8] = "d/dt CaNc in component CaNc (micromolar)"
    legend_rates[9] = "d/dt CaNc_star in component CaNc_star (micromolar)"
    legend_rates[10] = "d/dt CaNn_star in component CaNn_star (micromolar)"
    legend_rates[11] = "d/dt CaNn in component CaNn (micromolar)"
    legend_rates[12] = "d/dt CaM in component CaM (micromolar)"
    legend_rates[13] = "d/dt CaMCa in component CaMCa (micromolar)"
    legend_rates[14] = "d/dt MCIP in component MCIP (micromolar)"
    legend_rates[15] = "d/dt MCIPp in component MCIPp (micromolar)"
    legend_rates[16] = "d/dt MCIPpp in component MCIPpp (micromolar)"
    legend_rates[17] = "d/dt Comp1 in component Comp1 (micromolar)"
    legend_rates[18] = "d/dt Comp2 in component Comp2 (micromolar)"
    legend_rates[19] = "d/dt Comp3 in component Comp3 (micromolar)"
    legend_rates[20] = "d/dt P1433 in component P1433 (micromolar)"
    legend_rates[21] = "d/dt MCIP_BMK1 in component MCIP_BMK1 (micromolar)"
    legend_rates[22] = "d/dt MCIPp_CaNc_star in component MCIPp_CaNc_star (micromolar)"
    legend_rates[23] = "d/dt MCIPp_GSK3betac in component MCIPp_GSK3betac (micromolar)"
    legend_rates[24] = "d/dt MCIPpp_CaNc_star in component MCIPpp_CaNc_star (micromolar)"
    legend_rates[25] = "d/dt NFATpc_CaNc_star in component NFATpc_CaNc_star (micromolar)"
    legend_rates[26] = "d/dt NFATc_GSK3betac in component NFATc_GSK3betac (micromolar)"
    legend_rates[27] = "d/dt NFATpn_CaNn_star in component NFATpn_CaNn_star (micromolar)"
    legend_rates[28] = "d/dt NFATn_GSK3betan in component NFATn_GSK3betan (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.012
    states[1] = 3.33E-4
    states[2] = 2E-5
    states[3] = 4.99E-4
    states[4] = 4.94E-3
    states[5] = 8.01E-5
    states[6] = 0.17
    states[7] = 0.339
    states[8] = 0.91
    states[9] = 0.0275
    states[10] = 0.0568
    states[11] = 0.0057
    states[12] = 25.2
    states[13] = 7.88E-7
    states[14] = 2.15E-4
    states[15] = 7.76E-3
    states[16] = 0.0798
    states[17] = 5.21E-3
    states[18] = 0.283
    states[19] = 0.014
    states[20] = 0.708
    states[21] = 2.14E-5
    states[22] = 1.07E-4
    states[23] = 1.1E-3
    states[24] = 1.1E-3
    states[25] = 8.15E-5
    states[26] = 1.36E-6
    states[27] = 2.27E-6
    states[28] = 8.46E-5
    constants[0] = 0.4
    constants[1] = 0.1
    constants[2] = 0.1
    constants[3] = 0.05
    constants[4] = 0.114
    constants[5] = 0.0552
    constants[6] = 0.2
    constants[7] = 5
    constants[8] = 100
    constants[9] = 2760
    constants[10] = 0.072
    constants[11] = 50
    constants[12] = 0.0567
    constants[13] = 0.5
    constants[14] = 0.1
    constants[15] = 0.4
    constants[16] = 0.1
    constants[17] = 5
    constants[18] = 0.1
    constants[19] = 0.5
    constants[20] = 0.1
    constants[21] = 0.1
    constants[22] = 0.1
    constants[23] = 0.5
    constants[24] = 0.5
    constants[25] = 0.1
    constants[26] = 0.1
    constants[27] = 0.1
    constants[28] = 0.1
    constants[29] = 0.15
    constants[30] = 0.15
    constants[31] = 0.15
    constants[32] = 0.1
    constants[33] = 0.15
    constants[34] = 0.1
    constants[35] = 0.15
    constants[36] = 0.1
    constants[37] = 0.2
    constants[38] = 0.1
    constants[39] = 0.1
    constants[40] = 0.1
    constants[41] = 0.02
    constants[42] = 0.03
    constants[43] = 15
    constants[44] = 0.03
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[4] = constants[7]*states[12]*(power(constants[6], 4.00000))-constants[8]*states[13]
    rates[12] = -algebraic[4]
    algebraic[5] = constants[9]*states[13]*states[8]-constants[10]*states[9]
    rates[8] = -algebraic[5]
    rates[11] = -algebraic[5]
    rates[13] = algebraic[4]-algebraic[5]
    algebraic[6] = constants[11]*states[9]*states[14]-constants[12]*states[17]
    rates[17] = algebraic[6]
    algebraic[7] = constants[13]*states[20]*states[16]-constants[14]*states[18]
    rates[18] = algebraic[7]
    algebraic[8] = constants[15]*states[4]*states[20]-constants[16]*states[19]
    rates[19] = algebraic[8]
    rates[20] = -(algebraic[7]+algebraic[8])
    algebraic[9] = constants[17]*states[14]*states[0]-constants[18]*states[21]
    algebraic[10] = constants[19]*states[21]
    rates[0] = algebraic[10]-algebraic[9]
    rates[21] = algebraic[9]-algebraic[10]
    algebraic[12] = constants[20]*states[15]*states[9]-constants[21]*states[22]
    algebraic[13] = constants[22]*states[22]
    rates[22] = algebraic[12]-algebraic[13]
    algebraic[3] = constants[4]*states[9]-constants[5]*states[10]
    algebraic[15] = constants[35]*states[5]*states[10]-constants[36]*states[27]
    algebraic[17] = constants[37]*states[27]
    rates[10] = (algebraic[3]+algebraic[17])-algebraic[15]
    algebraic[16] = constants[23]*states[15]*states[6]-constants[24]*states[23]
    algebraic[18] = constants[25]*states[23]
    rates[23] = algebraic[16]-algebraic[18]
    rates[27] = algebraic[15]-algebraic[17]
    algebraic[2] = constants[2]*states[6]-constants[3]*states[7]
    algebraic[21] = constants[38]*states[3]*states[7]-constants[39]*states[28]
    algebraic[23] = constants[40]*states[28]
    rates[7] = (algebraic[2]+algebraic[23])-algebraic[21]
    algebraic[22] = constants[26]*states[16]*states[9]-constants[27]*states[24]
    algebraic[24] = constants[28]*states[24]
    rates[24] = algebraic[22]-algebraic[24]
    rates[28] = algebraic[21]-algebraic[23]
    algebraic[19] = algebraic[15]+algebraic[17]
    algebraic[25] = algebraic[21]+algebraic[23]
    algebraic[1] = constants[1]*states[5]
    rates[5] = algebraic[25]-(algebraic[1]+algebraic[19])
    algebraic[11] = algebraic[9]+algebraic[10]
    algebraic[14] = algebraic[12]+algebraic[13]
    algebraic[20] = algebraic[16]+algebraic[18]
    algebraic[27] = algebraic[22]+algebraic[24]
    rates[15] = (algebraic[11]+algebraic[27])-(algebraic[14]+algebraic[20])
    rates[16] = algebraic[20]-(algebraic[7]+algebraic[27])
    algebraic[28] = constants[41]*states[3]-constants[42]*states[1]
    algebraic[0] = constants[0]*states[2]
    rates[3] = (algebraic[0]+algebraic[19])-(algebraic[25]+algebraic[28])
    algebraic[26] = constants[29]*states[4]*states[9]-constants[30]*states[25]
    algebraic[29] = constants[31]*states[25]
    rates[9] = (algebraic[5]+algebraic[13]+algebraic[24]+algebraic[29])-(algebraic[3]+algebraic[6]+algebraic[12]+algebraic[22]+algebraic[26])
    rates[25] = algebraic[26]-algebraic[29]
    algebraic[30] = constants[44]*states[1]-(log(2.00000)/constants[43])*states[14]
    rates[1] = algebraic[30]-algebraic[28]
    rates[14] = algebraic[14]-(algebraic[6]+algebraic[11]+algebraic[30])
    algebraic[32] = constants[32]*states[2]*states[6]-constants[33]*states[26]
    algebraic[33] = constants[34]*states[26]
    rates[6] = (algebraic[18]+algebraic[33])-(algebraic[2]+algebraic[16]+algebraic[32])
    rates[26] = algebraic[32]-algebraic[33]
    algebraic[31] = algebraic[26]+algebraic[29]
    algebraic[34] = algebraic[32]+algebraic[33]
    rates[2] = algebraic[31]-(algebraic[0]+algebraic[34])
    rates[4] = (algebraic[1]+algebraic[34])-(algebraic[8]+algebraic[31])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[4] = constants[7]*states[12]*(power(constants[6], 4.00000))-constants[8]*states[13]
    algebraic[5] = constants[9]*states[13]*states[8]-constants[10]*states[9]
    algebraic[6] = constants[11]*states[9]*states[14]-constants[12]*states[17]
    algebraic[7] = constants[13]*states[20]*states[16]-constants[14]*states[18]
    algebraic[8] = constants[15]*states[4]*states[20]-constants[16]*states[19]
    algebraic[9] = constants[17]*states[14]*states[0]-constants[18]*states[21]
    algebraic[10] = constants[19]*states[21]
    algebraic[12] = constants[20]*states[15]*states[9]-constants[21]*states[22]
    algebraic[13] = constants[22]*states[22]
    algebraic[3] = constants[4]*states[9]-constants[5]*states[10]
    algebraic[15] = constants[35]*states[5]*states[10]-constants[36]*states[27]
    algebraic[17] = constants[37]*states[27]
    algebraic[16] = constants[23]*states[15]*states[6]-constants[24]*states[23]
    algebraic[18] = constants[25]*states[23]
    algebraic[2] = constants[2]*states[6]-constants[3]*states[7]
    algebraic[21] = constants[38]*states[3]*states[7]-constants[39]*states[28]
    algebraic[23] = constants[40]*states[28]
    algebraic[22] = constants[26]*states[16]*states[9]-constants[27]*states[24]
    algebraic[24] = constants[28]*states[24]
    algebraic[19] = algebraic[15]+algebraic[17]
    algebraic[25] = algebraic[21]+algebraic[23]
    algebraic[1] = constants[1]*states[5]
    algebraic[11] = algebraic[9]+algebraic[10]
    algebraic[14] = algebraic[12]+algebraic[13]
    algebraic[20] = algebraic[16]+algebraic[18]
    algebraic[27] = algebraic[22]+algebraic[24]
    algebraic[28] = constants[41]*states[3]-constants[42]*states[1]
    algebraic[0] = constants[0]*states[2]
    algebraic[26] = constants[29]*states[4]*states[9]-constants[30]*states[25]
    algebraic[29] = constants[31]*states[25]
    algebraic[30] = constants[44]*states[1]-(log(2.00000)/constants[43])*states[14]
    algebraic[32] = constants[32]*states[2]*states[6]-constants[33]*states[26]
    algebraic[33] = constants[34]*states[26]
    algebraic[31] = algebraic[26]+algebraic[29]
    algebraic[34] = algebraic[32]+algebraic[33]
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