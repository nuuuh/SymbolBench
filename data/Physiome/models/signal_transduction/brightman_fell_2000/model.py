# Size of variable arrays:
sizeAlgebraic = 29
sizeStates = 29
sizeConstants = 67
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "Rs in component Rs (dimensionless)"
    legend_algebraic[0] = "v1 in component v1 (first_order_rate_constant)"
    legend_algebraic[7] = "v3 in component v3 (first_order_rate_constant)"
    legend_states[1] = "RL in component RL (dimensionless)"
    legend_algebraic[4] = "v2 in component v2 (first_order_rate_constant)"
    legend_algebraic[10] = "v4 in component v4 (first_order_rate_constant)"
    legend_states[2] = "Ri in component Ri (dimensionless)"
    legend_algebraic[15] = "v6 in component v6 (first_order_rate_constant)"
    legend_states[3] = "L in component L (dimensionless)"
    legend_states[4] = "R2L2 in component R2L2 (dimensionless)"
    legend_algebraic[13] = "v5 in component v5 (first_order_rate_constant)"
    legend_algebraic[18] = "v7 in component v7 (first_order_rate_constant)"
    legend_states[5] = "R2_CPP in component R2_CPP (dimensionless)"
    legend_algebraic[22] = "v8 in component v8 (first_order_rate_constant)"
    legend_states[6] = "Li in component Li (dimensionless)"
    legend_states[7] = "R2i in component R2i (dimensionless)"
    legend_states[8] = "Shc in component Shc (dimensionless)"
    legend_algebraic[1] = "v9 in component v9 (first_order_rate_constant)"
    legend_algebraic[5] = "v10 in component v10 (first_order_rate_constant)"
    legend_states[9] = "ShcP in component ShcP (dimensionless)"
    legend_algebraic[8] = "v11 in component v11 (first_order_rate_constant)"
    legend_algebraic[16] = "v27 in component v27 (first_order_rate_constant)"
    legend_states[10] = "ShcGS in component ShcGS (dimensionless)"
    legend_algebraic[14] = "v13 in component v13 (first_order_rate_constant)"
    legend_algebraic[11] = "v12 in component v12 (first_order_rate_constant)"
    legend_states[11] = "GS in component GS (dimensionless)"
    legend_algebraic[19] = "v28 in component v28 (first_order_rate_constant)"
    legend_states[12] = "GSP in component GSP (dimensionless)"
    legend_states[13] = "RasGDP in component RasGDP (dimensionless)"
    legend_algebraic[20] = "v15 in component v15 (first_order_rate_constant)"
    legend_states[14] = "Ras_ShcGS in component Ras_ShcGS (dimensionless)"
    legend_states[15] = "RasGTP in component RasGTP (dimensionless)"
    legend_algebraic[23] = "v17 in component v17 (first_order_rate_constant)"
    legend_algebraic[17] = "v14 in component v14 (first_order_rate_constant)"
    legend_algebraic[21] = "v16 in component v16 (first_order_rate_constant)"
    legend_states[16] = "GAP in component GAP (dimensionless)"
    legend_states[17] = "Ras_GAP in component Ras_GAP (dimensionless)"
    legend_states[18] = "Raf in component Raf (dimensionless)"
    legend_algebraic[24] = "v18 in component v18 (first_order_rate_constant)"
    legend_states[19] = "Ras_Raf in component Ras_Raf (dimensionless)"
    legend_states[20] = "Rafa in component Rafa (dimensionless)"
    legend_algebraic[25] = "v19 in component v19 (first_order_rate_constant)"
    legend_algebraic[27] = "v21 in component v21 (first_order_rate_constant)"
    legend_states[21] = "MEK in component MEK (dimensionless)"
    legend_algebraic[26] = "v20 in component v20 (first_order_rate_constant)"
    legend_states[22] = "MEKP in component MEKP (dimensionless)"
    legend_algebraic[28] = "v22 in component v22 (first_order_rate_constant)"
    legend_states[23] = "MEKPP in component MEKPP (dimensionless)"
    legend_states[24] = "ERK in component ERK (dimensionless)"
    legend_algebraic[6] = "v24 in component v24 (first_order_rate_constant)"
    legend_algebraic[2] = "v23 in component v23 (first_order_rate_constant)"
    legend_states[25] = "ERKP in component ERKP (dimensionless)"
    legend_algebraic[12] = "v26 in component v26 (first_order_rate_constant)"
    legend_algebraic[9] = "v25 in component v25 (first_order_rate_constant)"
    legend_states[26] = "ERKPP in component ERKPP (dimensionless)"
    legend_states[27] = "t in component t (dimensionless)"
    legend_constants[65] = "v29 in component v29 (first_order_rate_constant)"
    legend_states[28] = "X in component X (dimensionless)"
    legend_algebraic[3] = "v30 in component v30 (first_order_rate_constant)"
    legend_constants[0] = "k1 in component v1 (first_order_rate_constant)"
    legend_constants[1] = "kn1 in component v1 (first_order_rate_constant)"
    legend_constants[2] = "DT in component v2 (dimensionless)"
    legend_constants[3] = "E in component v2 (dimensionless)"
    legend_constants[4] = "k2 in component v2 (first_order_rate_constant)"
    legend_constants[5] = "f in component v2 (dimensionless)"
    legend_constants[6] = "DT in component v3 (dimensionless)"
    legend_constants[7] = "E in component v3 (dimensionless)"
    legend_constants[8] = "kn3 in component v3 (first_order_rate_constant)"
    legend_constants[9] = "f in component v3 (dimensionless)"
    legend_constants[10] = "k3 in component v3 (first_order_rate_constant)"
    legend_constants[11] = "k2_4 in component v4 (first_order_rate_constant)"
    legend_constants[12] = "k4 in component v4 (first_order_rate_constant)"
    legend_constants[13] = "DT in component v5 (dimensionless)"
    legend_constants[14] = "E in component v5 (dimensionless)"
    legend_constants[15] = "k5 in component v5 (first_order_rate_constant)"
    legend_constants[16] = "f in component v5 (dimensionless)"
    legend_constants[17] = "DT in component v6 (dimensionless)"
    legend_constants[18] = "E in component v6 (dimensionless)"
    legend_constants[19] = "k6 in component v6 (first_order_rate_constant)"
    legend_constants[20] = "k7 in component v7 (first_order_rate_constant)"
    legend_constants[21] = "f in component v7 (dimensionless)"
    legend_constants[22] = "kn7 in component v7 (first_order_rate_constant)"
    legend_constants[23] = "DT in component v8 (dimensionless)"
    legend_constants[24] = "E in component v8 (dimensionless)"
    legend_constants[25] = "k8 in component v8 (first_order_rate_constant)"
    legend_constants[26] = "k9 in component v9 (first_order_rate_constant)"
    legend_constants[27] = "K_9 in component v9 (dimensionless)"
    legend_constants[28] = "V_10 in component v10 (first_order_rate_constant)"
    legend_constants[29] = "K_10 in component v10 (dimensionless)"
    legend_constants[30] = "k11 in component v11 (first_order_rate_constant)"
    legend_constants[31] = "kn11 in component v11 (first_order_rate_constant)"
    legend_constants[32] = "k12 in component v12 (first_order_rate_constant)"
    legend_constants[33] = "kn12 in component v12 (first_order_rate_constant)"
    legend_constants[34] = "k_13 in component v13 (first_order_rate_constant)"
    legend_constants[35] = "k14 in component v14 (first_order_rate_constant)"
    legend_constants[36] = "kn14 in component v14 (first_order_rate_constant)"
    legend_constants[37] = "k15 in component v15 (first_order_rate_constant)"
    legend_constants[38] = "k16 in component v16 (first_order_rate_constant)"
    legend_constants[39] = "kn16 in component v16 (first_order_rate_constant)"
    legend_constants[40] = "k17 in component v17 (first_order_rate_constant)"
    legend_constants[41] = "V_18 in component v18 (first_order_rate_constant)"
    legend_constants[42] = "K_18 in component v18 (dimensionless)"
    legend_constants[43] = "k19 in component v19 (first_order_rate_constant)"
    legend_constants[44] = "K_19 in component v19 (dimensionless)"
    legend_constants[45] = "V_20 in component v20 (first_order_rate_constant)"
    legend_constants[46] = "K_20 in component v20 (dimensionless)"
    legend_constants[47] = "k21 in component v21 (first_order_rate_constant)"
    legend_constants[48] = "K_21 in component v21 (dimensionless)"
    legend_constants[49] = "V_22 in component v22 (first_order_rate_constant)"
    legend_constants[50] = "K_22 in component v22 (dimensionless)"
    legend_constants[51] = "k23 in component v23 (first_order_rate_constant)"
    legend_constants[52] = "K_23 in component v23 (dimensionless)"
    legend_constants[53] = "V_24 in component v24 (first_order_rate_constant)"
    legend_constants[54] = "K_24 in component v24 (dimensionless)"
    legend_constants[55] = "k25 in component v25 (first_order_rate_constant)"
    legend_constants[56] = "K_25 in component v25 (dimensionless)"
    legend_constants[57] = "V_26 in component v26 (first_order_rate_constant)"
    legend_constants[58] = "K_26 in component v26 (dimensionless)"
    legend_constants[59] = "k27 in component v27 (first_order_rate_constant)"
    legend_constants[60] = "K_27 in component v27 (dimensionless)"
    legend_constants[61] = "V_28 in component v28 (first_order_rate_constant)"
    legend_constants[62] = "K_28 in component v28 (dimensionless)"
    legend_constants[63] = "v_1 in component v29 (first_order_rate_constant)"
    legend_constants[64] = "k_11 in component v30 (first_order_rate_constant)"
    legend_rates[0] = "d/dt Rs in component Rs (dimensionless)"
    legend_rates[1] = "d/dt RL in component RL (dimensionless)"
    legend_rates[2] = "d/dt Ri in component Ri (dimensionless)"
    legend_rates[3] = "d/dt L in component L (dimensionless)"
    legend_rates[4] = "d/dt R2L2 in component R2L2 (dimensionless)"
    legend_rates[5] = "d/dt R2_CPP in component R2_CPP (dimensionless)"
    legend_rates[6] = "d/dt Li in component Li (dimensionless)"
    legend_rates[7] = "d/dt R2i in component R2i (dimensionless)"
    legend_rates[8] = "d/dt Shc in component Shc (dimensionless)"
    legend_rates[9] = "d/dt ShcP in component ShcP (dimensionless)"
    legend_rates[10] = "d/dt ShcGS in component ShcGS (dimensionless)"
    legend_rates[11] = "d/dt GS in component GS (dimensionless)"
    legend_rates[12] = "d/dt GSP in component GSP (dimensionless)"
    legend_rates[13] = "d/dt RasGDP in component RasGDP (dimensionless)"
    legend_rates[14] = "d/dt Ras_ShcGS in component Ras_ShcGS (dimensionless)"
    legend_rates[15] = "d/dt RasGTP in component RasGTP (dimensionless)"
    legend_rates[16] = "d/dt GAP in component GAP (dimensionless)"
    legend_rates[17] = "d/dt Ras_GAP in component Ras_GAP (dimensionless)"
    legend_rates[18] = "d/dt Raf in component Raf (dimensionless)"
    legend_rates[19] = "d/dt Ras_Raf in component Ras_Raf (dimensionless)"
    legend_rates[20] = "d/dt Rafa in component Rafa (dimensionless)"
    legend_rates[21] = "d/dt MEK in component MEK (dimensionless)"
    legend_rates[22] = "d/dt MEKP in component MEKP (dimensionless)"
    legend_rates[23] = "d/dt MEKPP in component MEKPP (dimensionless)"
    legend_rates[24] = "d/dt ERK in component ERK (dimensionless)"
    legend_rates[25] = "d/dt ERKP in component ERKP (dimensionless)"
    legend_rates[26] = "d/dt ERKPP in component ERKPP (dimensionless)"
    legend_rates[27] = "d/dt t in component t (dimensionless)"
    legend_rates[28] = "d/dt X in component X (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 11100
    states[1] = 0
    states[2] = 3900
    states[3] = 0.0000001
    states[4] = 0
    states[5] = 0
    states[6] = 0
    states[7] = 0
    states[8] = 30000
    states[9] = 0
    states[10] = 0
    states[11] = 20000
    states[12] = 0
    states[13] = 19800
    states[14] = 0
    states[15] = 200
    states[16] = 15000
    states[17] = 0
    states[18] = 10000
    states[19] = 0
    states[20] = 0
    states[21] = 360000
    states[22] = 0
    states[23] = 0
    states[24] = 750000
    states[25] = 0
    states[26] = 0
    states[27] = 0
    states[28] = 0
    constants[0] = 384210000
    constants[1] = 0.73
    constants[2] = 6.5
    constants[3] = 0.12
    constants[4] = 0.7
    constants[5] = 0.2
    constants[6] = 6.5
    constants[7] = 0.12
    constants[8] = 0.7
    constants[9] = 0.2
    constants[10] = 0.0484
    constants[11] = 0.0000001
    constants[12] = 0.001383
    constants[13] = 6.5
    constants[14] = 0.12
    constants[15] = 0.35
    constants[16] = 0.2
    constants[17] = 6.5
    constants[18] = 0.12
    constants[19] = 0.35
    constants[20] = 1
    constants[21] = 0.2
    constants[22] = 0.000347
    constants[23] = 6.5
    constants[24] = 0.12
    constants[25] = 0.35
    constants[26] = 12
    constants[27] = 6000
    constants[28] = 300000
    constants[29] = 6000
    constants[30] = 0.002
    constants[31] = 3.81
    constants[32] = 0.0163
    constants[33] = 10
    constants[34] = 15
    constants[35] = 0.005
    constants[36] = 60
    constants[37] = 720
    constants[38] = 0.0012
    constants[39] = 3
    constants[40] = 27
    constants[41] = 97000
    constants[42] = 6000
    constants[43] = 50
    constants[44] = 9000
    constants[45] = 920000
    constants[46] = 600000
    constants[47] = 50
    constants[48] = 9000
    constants[49] = 920000
    constants[50] = 600000
    constants[51] = 8.3
    constants[52] = 90000
    constants[53] = 200000
    constants[54] = 600000
    constants[55] = 8.3
    constants[56] = 90000
    constants[57] = 400000
    constants[58] = 600000
    constants[59] = 1.6
    constants[60] = 600000
    constants[61] = 75
    constants[62] = 20000
    constants[63] = 1
    constants[64] = 0
    constants[65] = constants[63]
    constants[66] = constants[65]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[27] = constants[66]
    algebraic[0] = constants[0]*states[0]*states[3]-constants[1]*states[1]
    rates[3] = -algebraic[0]
    algebraic[3] = constants[64]*states[27]
    rates[28] = algebraic[3]
    algebraic[1] = (constants[26]*2.00000*(states[4]+states[7]+states[5])*states[8])/(constants[27]+states[8])
    algebraic[5] = (constants[28]*states[9])/(constants[29]+states[9])
    rates[8] = algebraic[5]-algebraic[1]
    algebraic[6] = (constants[53]*states[25])/(constants[54]+states[25])
    algebraic[2] = (constants[51]*states[24]*(states[22]+states[23]))/(constants[52]+states[24])
    rates[24] = algebraic[6]-algebraic[2]
    algebraic[7] = constants[10]*states[2]-constants[9]*constants[8]*(constants[7]+(1.00000-constants[7])*(1.00000-exp(-(power(states[27]/constants[6], 3.00000)))))*states[0]
    rates[0] = algebraic[7]-algebraic[0]
    algebraic[4] = constants[5]*constants[4]*(constants[3]+(1.00000-constants[3])*(1.00000-exp(-(power(states[27]/constants[2], 3.00000)))))*states[1]
    algebraic[10] = constants[12]*states[1]*states[1]-constants[11]*states[4]
    rates[1] = (algebraic[0]-algebraic[4])-algebraic[10]
    algebraic[12] = (constants[57]*states[26])/(constants[58]+states[26])
    algebraic[9] = (constants[55]*states[25]*(states[22]+states[23]))/(constants[56]+states[25])
    rates[25] = ((algebraic[2]+algebraic[12])-algebraic[6])-algebraic[9]
    rates[26] = algebraic[9]-algebraic[12]
    algebraic[14] = constants[34]*states[14]
    algebraic[11] = constants[32]*states[13]*states[10]-constants[33]*states[14]
    rates[14] = algebraic[11]-algebraic[14]
    algebraic[15] = constants[19]*(constants[18]+(1.00000-constants[18])*(1.00000-exp(-(power(states[27]/constants[17], 3.00000)))))*states[7]
    rates[2] = algebraic[4]+algebraic[7]+algebraic[15]
    rates[6] = algebraic[4]+algebraic[15]
    algebraic[8] = constants[30]*states[9]*states[11]-constants[31]*states[10]
    algebraic[16] = (states[26]*states[10]*constants[59])/(constants[60]+states[10])
    rates[9] = (algebraic[16]-algebraic[5])-algebraic[8]
    rates[10] = ((algebraic[14]+algebraic[8])-algebraic[16])-algebraic[11]
    algebraic[13] = constants[16]*constants[15]*(constants[14]+(1.00000-constants[14])*(1.00000-exp(-(power(states[27]/constants[13], 3.00000)))))*states[4]
    algebraic[18] = constants[20]*constants[21]*states[4]-constants[22]*states[5]
    rates[4] = (algebraic[10]-algebraic[13])-algebraic[18]
    algebraic[19] = (constants[61]*states[12])/(constants[62]+states[12])
    rates[11] = algebraic[19]-algebraic[8]
    rates[12] = algebraic[16]-algebraic[19]
    algebraic[20] = constants[37]*states[17]
    rates[13] = algebraic[20]-algebraic[11]
    algebraic[17] = constants[35]*states[15]*states[16]-constants[36]*states[17]
    rates[16] = algebraic[20]-algebraic[17]
    rates[17] = algebraic[17]-algebraic[20]
    algebraic[22] = constants[25]*(constants[24]+(1.00000-constants[24])*(1.00000-exp(-(power(states[27]/constants[23], 3.00000)))))*states[5]
    rates[5] = algebraic[18]-algebraic[22]
    rates[7] = (algebraic[13]+algebraic[22])-algebraic[15]
    algebraic[23] = constants[40]*states[19]
    algebraic[21] = constants[38]*states[15]*states[18]-constants[39]*states[19]
    rates[15] = ((algebraic[14]+algebraic[23])-algebraic[17])-algebraic[21]
    rates[19] = algebraic[21]-algebraic[23]
    algebraic[24] = (constants[41]*states[20])/(constants[42]+states[20])
    rates[18] = algebraic[24]-algebraic[21]
    algebraic[25] = (states[20]*states[21]*constants[43])/(constants[44]+states[21])
    algebraic[26] = (constants[45]*states[22])/(constants[46]+states[22])
    rates[21] = algebraic[26]-algebraic[25]
    algebraic[27] = (states[20]*states[22]*constants[47])/(constants[48]+states[22])
    rates[20] = ((algebraic[23]-algebraic[25])-algebraic[27])-algebraic[24]
    algebraic[28] = (constants[49]*states[23])/(constants[50]+states[23])
    rates[22] = ((algebraic[25]+algebraic[28])-algebraic[26])-algebraic[27]
    rates[23] = algebraic[27]-algebraic[28]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[0]*states[0]*states[3]-constants[1]*states[1]
    algebraic[3] = constants[64]*states[27]
    algebraic[1] = (constants[26]*2.00000*(states[4]+states[7]+states[5])*states[8])/(constants[27]+states[8])
    algebraic[5] = (constants[28]*states[9])/(constants[29]+states[9])
    algebraic[6] = (constants[53]*states[25])/(constants[54]+states[25])
    algebraic[2] = (constants[51]*states[24]*(states[22]+states[23]))/(constants[52]+states[24])
    algebraic[7] = constants[10]*states[2]-constants[9]*constants[8]*(constants[7]+(1.00000-constants[7])*(1.00000-exp(-(power(states[27]/constants[6], 3.00000)))))*states[0]
    algebraic[4] = constants[5]*constants[4]*(constants[3]+(1.00000-constants[3])*(1.00000-exp(-(power(states[27]/constants[2], 3.00000)))))*states[1]
    algebraic[10] = constants[12]*states[1]*states[1]-constants[11]*states[4]
    algebraic[12] = (constants[57]*states[26])/(constants[58]+states[26])
    algebraic[9] = (constants[55]*states[25]*(states[22]+states[23]))/(constants[56]+states[25])
    algebraic[14] = constants[34]*states[14]
    algebraic[11] = constants[32]*states[13]*states[10]-constants[33]*states[14]
    algebraic[15] = constants[19]*(constants[18]+(1.00000-constants[18])*(1.00000-exp(-(power(states[27]/constants[17], 3.00000)))))*states[7]
    algebraic[8] = constants[30]*states[9]*states[11]-constants[31]*states[10]
    algebraic[16] = (states[26]*states[10]*constants[59])/(constants[60]+states[10])
    algebraic[13] = constants[16]*constants[15]*(constants[14]+(1.00000-constants[14])*(1.00000-exp(-(power(states[27]/constants[13], 3.00000)))))*states[4]
    algebraic[18] = constants[20]*constants[21]*states[4]-constants[22]*states[5]
    algebraic[19] = (constants[61]*states[12])/(constants[62]+states[12])
    algebraic[20] = constants[37]*states[17]
    algebraic[17] = constants[35]*states[15]*states[16]-constants[36]*states[17]
    algebraic[22] = constants[25]*(constants[24]+(1.00000-constants[24])*(1.00000-exp(-(power(states[27]/constants[23], 3.00000)))))*states[5]
    algebraic[23] = constants[40]*states[19]
    algebraic[21] = constants[38]*states[15]*states[18]-constants[39]*states[19]
    algebraic[24] = (constants[41]*states[20])/(constants[42]+states[20])
    algebraic[25] = (states[20]*states[21]*constants[43])/(constants[44]+states[21])
    algebraic[26] = (constants[45]*states[22])/(constants[46]+states[22])
    algebraic[27] = (states[20]*states[22]*constants[47])/(constants[48]+states[22])
    algebraic[28] = (constants[49]*states[23])/(constants[50]+states[23])
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