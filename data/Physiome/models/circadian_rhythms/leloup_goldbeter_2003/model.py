# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 19
sizeConstants = 63
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_states[0] = "MP in component MP (nanomolar)"
    legend_constants[0] = "vsP in component MP (flux)"
    legend_constants[1] = "vmP in component MP (flux)"
    legend_constants[2] = "kdmp in component MP (first_order_rate_constant)"
    legend_constants[3] = "KAP in component MP (nanomolar)"
    legend_constants[4] = "KmP in component MP (nanomolar)"
    legend_constants[5] = "n in component model_parameters (dimensionless)"
    legend_states[1] = "BN in component BN (nanomolar)"
    legend_states[2] = "MC in component MC (nanomolar)"
    legend_constants[6] = "vsC in component MC (flux)"
    legend_constants[7] = "vmC in component MC (flux)"
    legend_constants[8] = "kdmc in component MC (first_order_rate_constant)"
    legend_constants[9] = "KAC in component MC (nanomolar)"
    legend_constants[10] = "KmC in component MC (nanomolar)"
    legend_states[3] = "MB in component MB (nanomolar)"
    legend_constants[11] = "vsB in component MB (flux)"
    legend_constants[12] = "vmB in component MB (flux)"
    legend_constants[13] = "kdmb in component MB (first_order_rate_constant)"
    legend_constants[14] = "KIB in component MB (nanomolar)"
    legend_constants[15] = "KmB in component MB (nanomolar)"
    legend_constants[16] = "m in component model_parameters (dimensionless)"
    legend_states[4] = "RN in component RN (nanomolar)"
    legend_states[5] = "MR in component MR (nanomolar)"
    legend_constants[17] = "vsR in component MR (flux)"
    legend_constants[18] = "vmR in component MR (flux)"
    legend_constants[19] = "kdmr in component MR (first_order_rate_constant)"
    legend_constants[20] = "KAR in component MR (nanomolar)"
    legend_constants[21] = "KmR in component MR (nanomolar)"
    legend_constants[22] = "h in component model_parameters (dimensionless)"
    legend_states[6] = "PC in component PC (nanomolar)"
    legend_constants[23] = "ksP in component model_parameters (first_order_rate_constant)"
    legend_constants[24] = "Kp in component model_parameters (nanomolar)"
    legend_constants[25] = "Kdp in component model_parameters (nanomolar)"
    legend_constants[26] = "k3 in component model_parameters (second_order_rate_constant)"
    legend_constants[27] = "k4 in component model_parameters (first_order_rate_constant)"
    legend_constants[28] = "kdn in component model_parameters (first_order_rate_constant)"
    legend_constants[29] = "V1P in component model_parameters (flux)"
    legend_constants[30] = "V2P in component model_parameters (flux)"
    legend_states[7] = "PCP in component PCP (nanomolar)"
    legend_states[8] = "PCC in component PCC (nanomolar)"
    legend_states[9] = "CC in component CC (nanomolar)"
    legend_constants[31] = "ksC in component model_parameters (first_order_rate_constant)"
    legend_constants[32] = "kdnc in component model_parameters (first_order_rate_constant)"
    legend_constants[33] = "V1C in component model_parameters (flux)"
    legend_constants[34] = "V2C in component model_parameters (flux)"
    legend_states[10] = "CCP in component CCP (nanomolar)"
    legend_states[11] = "RC in component RC (nanomolar)"
    legend_constants[35] = "ksR in component model_parameters (first_order_rate_constant)"
    legend_constants[36] = "Kd in component model_parameters (nanomolar)"
    legend_constants[37] = "k9 in component model_parameters (first_order_rate_constant)"
    legend_constants[38] = "k10 in component model_parameters (first_order_rate_constant)"
    legend_constants[39] = "vdRC in component model_parameters (flux)"
    legend_constants[40] = "vdPC in component model_parameters (flux)"
    legend_constants[41] = "vdCC in component model_parameters (flux)"
    legend_constants[42] = "k1 in component model_parameters (first_order_rate_constant)"
    legend_constants[43] = "k2 in component model_parameters (first_order_rate_constant)"
    legend_constants[44] = "V1PC in component model_parameters (flux)"
    legend_constants[45] = "V2PC in component model_parameters (flux)"
    legend_states[12] = "PCCP in component PCCP (nanomolar)"
    legend_states[13] = "PCN in component PCN (nanomolar)"
    legend_constants[46] = "k7 in component model_parameters (second_order_rate_constant)"
    legend_constants[47] = "k8 in component model_parameters (first_order_rate_constant)"
    legend_constants[48] = "V3PC in component model_parameters (flux)"
    legend_constants[49] = "V4PC in component model_parameters (flux)"
    legend_states[14] = "PCNP in component PCNP (nanomolar)"
    legend_states[15] = "IN in component IN (nanomolar)"
    legend_constants[50] = "vdRN in component model_parameters (flux)"
    legend_constants[51] = "vdPCC in component model_parameters (flux)"
    legend_constants[52] = "vdPCN in component model_parameters (flux)"
    legend_states[16] = "BC in component BC (nanomolar)"
    legend_constants[53] = "ksB in component model_parameters (first_order_rate_constant)"
    legend_constants[54] = "k5 in component model_parameters (first_order_rate_constant)"
    legend_constants[55] = "k6 in component model_parameters (first_order_rate_constant)"
    legend_constants[56] = "V1B in component model_parameters (flux)"
    legend_constants[57] = "V2B in component model_parameters (flux)"
    legend_states[17] = "BCP in component BCP (nanomolar)"
    legend_constants[58] = "vdBC in component model_parameters (flux)"
    legend_constants[59] = "V3B in component model_parameters (flux)"
    legend_constants[60] = "V4B in component model_parameters (flux)"
    legend_states[18] = "BNP in component BNP (nanomolar)"
    legend_constants[61] = "vdBN in component model_parameters (flux)"
    legend_constants[62] = "vdIN in component model_parameters (flux)"
    legend_rates[0] = "d/dt MP in component MP (nanomolar)"
    legend_rates[2] = "d/dt MC in component MC (nanomolar)"
    legend_rates[3] = "d/dt MB in component MB (nanomolar)"
    legend_rates[5] = "d/dt MR in component MR (nanomolar)"
    legend_rates[6] = "d/dt PC in component PC (nanomolar)"
    legend_rates[9] = "d/dt CC in component CC (nanomolar)"
    legend_rates[11] = "d/dt RC in component RC (nanomolar)"
    legend_rates[7] = "d/dt PCP in component PCP (nanomolar)"
    legend_rates[10] = "d/dt CCP in component CCP (nanomolar)"
    legend_rates[8] = "d/dt PCC in component PCC (nanomolar)"
    legend_rates[13] = "d/dt PCN in component PCN (nanomolar)"
    legend_rates[4] = "d/dt RN in component RN (nanomolar)"
    legend_rates[12] = "d/dt PCCP in component PCCP (nanomolar)"
    legend_rates[14] = "d/dt PCNP in component PCNP (nanomolar)"
    legend_rates[16] = "d/dt BC in component BC (nanomolar)"
    legend_rates[17] = "d/dt BCP in component BCP (nanomolar)"
    legend_rates[1] = "d/dt BN in component BN (nanomolar)"
    legend_rates[18] = "d/dt BNP in component BNP (nanomolar)"
    legend_rates[15] = "d/dt IN in component IN (nanomolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.1
    constants[0] = 2.4
    constants[1] = 2.2
    constants[2] = 0.02
    constants[3] = 0.6
    constants[4] = 0.3
    constants[5] = 2.0
    states[1] = 0.1
    states[2] = 1.2
    constants[6] = 2.2
    constants[7] = 2.0
    constants[8] = 0.02
    constants[9] = 0.6
    constants[10] = 0.4
    states[3] = 9
    constants[11] = 1.8
    constants[12] = 1.3
    constants[13] = 0.02
    constants[14] = 2.2
    constants[15] = 0.4
    constants[16] = 2.0
    states[4] = 0.1
    states[5] = 1.5
    constants[17] = 1.6
    constants[18] = 1.6
    constants[19] = 0.02
    constants[20] = 0.6
    constants[21] = 0.4
    constants[22] = 2.0
    states[6] = 0.1
    constants[23] = 1.2
    constants[24] = 1.006
    constants[25] = 0.1
    constants[26] = 0.8
    constants[27] = 0.4
    constants[28] = 0.02
    constants[29] = 9.6
    constants[30] = 0.6
    states[7] = 0.1
    states[8] = 0.1
    states[9] = 0.1
    constants[31] = 3.2
    constants[32] = 0.02
    constants[33] = 1.2
    constants[34] = 0.2
    states[10] = 0.1
    states[11] = 0.1
    constants[35] = 1.7
    constants[36] = 0.3
    constants[37] = 0.8
    constants[38] = 0.4
    constants[39] = 4.4
    constants[40] = 3.4
    constants[41] = 1.4
    constants[42] = 0.8
    constants[43] = 0.4
    constants[44] = 2.4
    constants[45] = 0.2
    states[12] = 0.1
    states[13] = 0.1
    constants[46] = 1.0
    constants[47] = 0.2
    constants[48] = 2.4
    constants[49] = 0.2
    states[14] = 0.1
    states[15] = 0.1
    constants[50] = 0.8
    constants[51] = 1.4
    constants[52] = 1.4
    states[16] = 0.1
    constants[53] = 0.32
    constants[54] = 0.8
    constants[55] = 0.4
    constants[56] = 1.4
    constants[57] = 0.2
    states[17] = 0.1
    constants[58] = 3.0
    constants[59] = 1.4
    constants[60] = 0.4
    states[18] = 0.1
    constants[61] = 3.0
    constants[62] = 1.6
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[0]*((power(states[1], constants[5]))/(power(constants[3], constants[5])+power(states[1], constants[5])))-(constants[1]*(states[0]/(constants[4]+states[0]))+constants[2]*states[0])
    rates[2] = constants[6]*((power(states[1], constants[5]))/(power(constants[9], constants[5])+power(states[1], constants[5])))-(constants[7]*(states[2]/(constants[10]+states[2]))+constants[8]*states[2])
    rates[3] = constants[11]*((power(constants[14], constants[16]))/(power(constants[14], constants[16])+power(states[4], constants[16])))-(constants[12]*(states[3]/(constants[15]+states[3]))+constants[13]*states[3])
    rates[5] = constants[17]*((power(states[1], constants[22]))/(power(constants[20], constants[22])+power(states[1], constants[22])))-(constants[18]*(states[5]/(constants[21]+states[5]))+constants[19]*states[5])
    rates[6] = (constants[23]*states[0]+constants[30]*(states[7]/(constants[25]+states[7]))+constants[27]*states[8])-(constants[29]*(states[6]/(constants[24]+states[6]))+constants[26]*states[6]*states[9]+constants[28]*states[6])
    rates[9] = (constants[31]*states[2]+constants[34]*(states[10]/(constants[25]+states[10]))+constants[27]*states[8])-(constants[33]*(states[9]/(constants[24]+states[9]))+constants[26]*states[6]*states[9]+constants[32]*states[9])
    rates[11] = (constants[35]*states[5]+constants[38]*states[4])-(constants[37]*states[11]+constants[39]*(states[11]/(constants[36]+states[11]))+constants[28]*states[11])
    rates[7] = constants[29]*(states[6]/(constants[24]+states[6]))-(constants[30]*(states[7]/(constants[25]+states[7]))+constants[40]*(states[7]/(constants[36]+states[7]))+constants[28]*states[7])
    rates[10] = constants[33]*(states[9]/(constants[24]+states[9]))-(constants[34]*(states[10]/(constants[25]+states[10]))+constants[41]*(states[10]/(constants[36]+states[10]))+constants[28]*states[10])
    rates[8] = (constants[45]*(states[12]/(constants[25]+states[12]))+constants[26]*states[6]*states[9]+constants[43]*states[13])-(constants[44]*(states[8]/(constants[24]+states[8]))+constants[27]*states[8]+constants[42]*states[8]+constants[28]*states[8])
    rates[13] = (constants[49]*(states[14]/(constants[25]+states[14]))+constants[42]*states[8]+constants[47]*states[15])-(constants[48]*(states[13]/(constants[24]+states[13]))+constants[43]*states[13]+constants[46]*states[1]*states[13]+constants[28]*states[13])
    rates[4] = constants[37]*states[11]-(constants[38]*states[4]+constants[50]*(states[4]/(constants[36]+states[4]))+constants[28]*states[4])
    rates[12] = constants[44]*(states[8]/(constants[24]+states[8]))-(constants[45]*(states[12]/(constants[25]+states[12]))+constants[51]*(states[12]/(constants[36]+states[12]))+constants[28]*states[12])
    rates[14] = constants[48]*(states[13]/(constants[24]+states[13]))-(constants[49]*(states[14]/(constants[25]+states[14]))+constants[52]*(states[14]/(constants[36]+states[14]))+constants[28]*states[14])
    rates[16] = (constants[57]*(states[17]/(constants[25]+states[17]))+constants[55]*states[1]+constants[53]*states[3])-(constants[56]*(states[16]/(constants[24]+states[16]))+constants[54]*states[16]+constants[28]*states[16])
    rates[17] = constants[56]*(states[16]/(constants[24]+states[16]))-(constants[57]*(states[17]/(constants[25]+states[17]))+constants[58]*(states[17]/(constants[36]+states[17]))+constants[28]*states[17])
    rates[1] = (constants[60]*(states[18]/(constants[25]+states[18]))+constants[54]*states[16]+constants[47]*states[15])-(constants[59]*(states[1]/(constants[24]+states[1]))+constants[55]*states[1]+constants[46]*states[1]*states[13]+constants[28]*states[1])
    rates[18] = constants[59]*(states[1]/(constants[24]+states[1]))-(constants[60]*(states[18]/(constants[25]+states[18]))+constants[61]*(states[18]/(constants[36]+states[18]))+constants[28]*states[18])
    rates[15] = constants[46]*states[1]*states[13]-(constants[47]*states[15]+constants[62]*(states[15]/(constants[36]+states[15]))+constants[28]*states[15])
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