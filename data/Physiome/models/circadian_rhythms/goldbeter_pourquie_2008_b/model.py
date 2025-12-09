# Size of variable arrays:
sizeAlgebraic = 17
sizeStates = 16
sizeConstants = 72
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "N in component N (nanomolar)"
    legend_constants[0] = "KdN in component N (nanomolar)"
    legend_constants[1] = "vsN in component N (flux)"
    legend_constants[2] = "vdN in component N (flux)"
    legend_constants[3] = "kc in component model_parameters (first_order_rate_constant)"
    legend_constants[4] = "KIF in component model_parameters (nanomolar)"
    legend_states[1] = "F in component F (nanomolar)"
    legend_constants[5] = "epsilon in component model_parameters (dimensionless)"
    legend_constants[6] = "j in component model_parameters (dimensionless)"
    legend_states[2] = "Na in component Na (nanomolar)"
    legend_algebraic[0] = "Vtr in component Na (flux)"
    legend_constants[7] = "KdNa in component Na (nanomolar)"
    legend_constants[8] = "VdNa in component Na (flux)"
    legend_constants[9] = "kt1 in component Na (first_order_rate_constant)"
    legend_constants[10] = "kt2 in component Na (first_order_rate_constant)"
    legend_states[3] = "Nan in component Nan (nanomolar)"
    legend_constants[11] = "KdNan in component Nan (nanomolar)"
    legend_constants[12] = "VdNan in component Nan (flux)"
    legend_states[4] = "MF in component MF (nanomolar)"
    legend_constants[13] = "KdMF in component MF (nanomolar)"
    legend_constants[14] = "vsF in component MF (flux)"
    legend_constants[15] = "vmF in component MF (flux)"
    legend_constants[16] = "KA in component MF (nanomolar)"
    legend_constants[17] = "p in component model_parameters (dimensionless)"
    legend_constants[18] = "KdF in component F (nanomolar)"
    legend_constants[19] = "vdF in component F (flux)"
    legend_constants[20] = "ksF in component F (first_order_rate_constant)"
    legend_states[5] = "K in component K (nanomolar)"
    legend_algebraic[6] = "V1 in component Wnt_parameters (flux)"
    legend_constants[21] = "theta in component model_parameters (dimensionless)"
    legend_states[6] = "B in component B (nanomolar)"
    legend_constants[22] = "kd1 in component B (first_order_rate_constant)"
    legend_constants[23] = "vsB in component B (flux)"
    legend_algebraic[12] = "VK in component Wnt_parameters (flux)"
    legend_algebraic[16] = "VP in component Wnt_parameters (flux)"
    legend_algebraic[7] = "V2 in component Wnt_parameters (flux)"
    legend_constants[24] = "Kt in component Wnt_parameters (nanomolar)"
    legend_algebraic[1] = "AK in component Wnt_parameters (nanomolar)"
    legend_states[7] = "Bp in component Bp (nanomolar)"
    legend_constants[25] = "kd2 in component Bp (first_order_rate_constant)"
    legend_states[8] = "BN in component BN (nanomolar)"
    legend_states[9] = "MAx in component MAx (nanomolar)"
    legend_constants[26] = "v0 in component MAx (flux)"
    legend_constants[27] = "vMB in component MAx (flux)"
    legend_constants[28] = "vmd in component MAx (flux)"
    legend_constants[29] = "KaB in component MAx (nanomolar)"
    legend_constants[30] = "Kmd in component MAx (nanomolar)"
    legend_constants[31] = "n in component MAx (dimensionless)"
    legend_states[10] = "A in component A (nanomolar)"
    legend_constants[32] = "ksAx in component A (first_order_rate_constant)"
    legend_constants[33] = "vdAx in component A (flux)"
    legend_constants[34] = "KdAx in component A (nanomolar)"
    legend_constants[35] = "d1 in component Wnt_parameters (first_order_rate_constant)"
    legend_constants[36] = "a1 in component Wnt_parameters (second_order_rate_constant)"
    legend_constants[37] = "K1 in component Wnt_parameters (nanomolar)"
    legend_constants[38] = "K2 in component Wnt_parameters (nanomolar)"
    legend_constants[39] = "D in component Wnt_parameters (nanomolar)"
    legend_constants[40] = "KID in component Wnt_parameters (nanomolar)"
    legend_constants[41] = "kt3 in component Wnt_parameters (first_order_rate_constant)"
    legend_constants[42] = "kt4 in component Wnt_parameters (first_order_rate_constant)"
    legend_constants[43] = "VMK in component Wnt_parameters (flux)"
    legend_constants[44] = "VMP in component Wnt_parameters (flux)"
    legend_states[11] = "Rasa in component Rasa (nanomolar)"
    legend_algebraic[8] = "VaRas in component FGF_parameters (flux)"
    legend_algebraic[13] = "VdRas in component FGF_parameters (flux)"
    legend_constants[45] = "eta in component model_parameters (dimensionless)"
    legend_states[12] = "ERKa in component ERKa (nanomolar)"
    legend_algebraic[9] = "VaErk in component FGF_parameters (flux)"
    legend_algebraic[14] = "VdErk in component FGF_parameters (flux)"
    legend_states[13] = "Xa in component Xa (nanomolar)"
    legend_algebraic[10] = "VaX in component FGF_parameters (flux)"
    legend_algebraic[15] = "VdX in component FGF_parameters (flux)"
    legend_states[14] = "MDusp in component MDusp (nanomolar)"
    legend_algebraic[5] = "VsMDusp in component FGF_parameters (flux)"
    legend_algebraic[11] = "VdMDusp in component FGF_parameters (flux)"
    legend_states[15] = "Dusp in component Dusp (nanomolar)"
    legend_constants[46] = "ksDusp in component Dusp (first_order_rate_constant)"
    legend_constants[47] = "vdDusp in component Dusp (flux)"
    legend_constants[48] = "KdDusp in component Dusp (nanomolar)"
    legend_algebraic[2] = "Rasi in component FGF_parameters (nanomolar)"
    legend_algebraic[3] = "ERKi in component FGF_parameters (nanomolar)"
    legend_algebraic[4] = "Xi in component FGF_parameters (nanomolar)"
    legend_constants[49] = "Rast in component FGF_parameters (nanomolar)"
    legend_constants[50] = "ERKt in component FGF_parameters (nanomolar)"
    legend_constants[51] = "Xt in component FGF_parameters (nanomolar)"
    legend_constants[52] = "kcDusp in component FGF_parameters (first_order_rate_constant)"
    legend_constants[53] = "VMaRas in component FGF_parameters (flux)"
    legend_constants[54] = "VMdRas in component FGF_parameters (flux)"
    legend_constants[55] = "VMaErk in component FGF_parameters (flux)"
    legend_constants[56] = "VMaX in component FGF_parameters (flux)"
    legend_constants[57] = "VMdX in component FGF_parameters (flux)"
    legend_constants[58] = "VMsMDusp in component FGF_parameters (flux)"
    legend_constants[59] = "VMdMDusp in component FGF_parameters (flux)"
    legend_constants[60] = "Fgf in component FGF_parameters (nanomolar)"
    legend_constants[61] = "KaFgf in component FGF_parameters (nanomolar)"
    legend_constants[62] = "KaRas in component FGF_parameters (nanomolar)"
    legend_constants[63] = "KdRas in component FGF_parameters (nanomolar)"
    legend_constants[64] = "KdErk in component FGF_parameters (nanomolar)"
    legend_constants[65] = "KaErk in component FGF_parameters (nanomolar)"
    legend_constants[66] = "KaX in component FGF_parameters (nanomolar)"
    legend_constants[67] = "KdX in component FGF_parameters (nanomolar)"
    legend_constants[68] = "KaMDusp in component FGF_parameters (nanomolar)"
    legend_constants[69] = "KdMDusp in component FGF_parameters (nanomolar)"
    legend_constants[70] = "q in component FGF_parameters (dimensionless)"
    legend_constants[71] = "r in component FGF_parameters (dimensionless)"
    legend_rates[0] = "d/dt N in component N (nanomolar)"
    legend_rates[2] = "d/dt Na in component Na (nanomolar)"
    legend_rates[3] = "d/dt Nan in component Nan (nanomolar)"
    legend_rates[4] = "d/dt MF in component MF (nanomolar)"
    legend_rates[1] = "d/dt F in component F (nanomolar)"
    legend_rates[5] = "d/dt K in component K (nanomolar)"
    legend_rates[6] = "d/dt B in component B (nanomolar)"
    legend_rates[7] = "d/dt Bp in component Bp (nanomolar)"
    legend_rates[8] = "d/dt BN in component BN (nanomolar)"
    legend_rates[9] = "d/dt MAx in component MAx (nanomolar)"
    legend_rates[10] = "d/dt A in component A (nanomolar)"
    legend_rates[11] = "d/dt Rasa in component Rasa (nanomolar)"
    legend_rates[12] = "d/dt ERKa in component ERKa (nanomolar)"
    legend_rates[13] = "d/dt Xa in component Xa (nanomolar)"
    legend_rates[14] = "d/dt MDusp in component MDusp (nanomolar)"
    legend_rates[15] = "d/dt Dusp in component Dusp (nanomolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.5
    constants[0] = 1.4
    constants[1] = 0.23
    constants[2] = 2.82
    constants[3] = 3.45
    constants[4] = 0.5
    states[1] = 0.001
    constants[5] = 0.3
    constants[6] = 2.0
    states[2] = 0.2
    constants[7] = 0.001
    constants[8] = 0.01
    constants[9] = 0.1
    constants[10] = 0.1
    states[3] = 0.0
    constants[11] = 0.001
    constants[12] = 0.1
    states[4] = 0.1
    constants[13] = 0.768
    constants[14] = 3.0
    constants[15] = 1.92
    constants[16] = 0.05
    constants[17] = 2.0
    constants[18] = 0.37
    constants[19] = 0.39
    constants[20] = 0.3
    states[5] = 3.0
    constants[21] = 1.5
    states[6] = 0.1
    constants[22] = 0.0
    constants[23] = 0.087
    constants[24] = 3.0
    states[7] = 0.1
    constants[25] = 7.062
    states[8] = 0.001
    states[9] = 0.1
    constants[26] = 0.06
    constants[27] = 1.64
    constants[28] = 0.8
    constants[29] = 0.7
    constants[30] = 0.48
    constants[31] = 2.0
    states[10] = 0.1
    constants[32] = 0.02
    constants[33] = 0.6
    constants[34] = 0.63
    constants[35] = 0.1
    constants[36] = 1.8
    constants[37] = 0.28
    constants[38] = 0.03
    constants[39] = 2.0
    constants[40] = 0.5
    constants[41] = 0.7
    constants[42] = 1.5
    constants[43] = 5.08
    constants[44] = 1.0
    states[11] = 0.5
    constants[45] = 0.3
    states[12] = 0.2
    states[13] = 0.1
    states[14] = 0.1
    states[15] = 0.1
    constants[46] = 0.5
    constants[47] = 2.0
    constants[48] = 0.5
    constants[49] = 2.0
    constants[50] = 2.0
    constants[51] = 2.0
    constants[52] = 1.35
    constants[53] = 4.968
    constants[54] = 0.41
    constants[55] = 3.30
    constants[56] = 1.6
    constants[57] = 0.5
    constants[58] = 0.9
    constants[59] = 0.5
    constants[60] = 1.0
    constants[61] = 0.5
    constants[62] = 0.103
    constants[63] = 0.1
    constants[64] = 0.05
    constants[65] = 0.05
    constants[66] = 0.05
    constants[67] = 0.05
    constants[68] = 0.5
    constants[69] = 0.5
    constants[70] = 2.0
    constants[71] = 2.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[5]*(constants[1]-(constants[2]*(states[0]/(constants[0]+states[0]))+constants[3]*states[0]*((power(constants[4], constants[6]))/(power(constants[4], constants[6])+power(states[1], constants[6])))))
    rates[4] = constants[5]*(constants[14]*((power(states[3], constants[17]))/(power(constants[16], constants[17])+power(states[3], constants[17])))-constants[15]*(states[4]/(constants[13]+states[4])))
    rates[1] = constants[5]*(constants[20]*states[4]-constants[19]*(states[1]/(constants[18]+states[1])))
    rates[9] = constants[21]*((constants[26]+constants[27]*((power(states[8], constants[31]))/(power(constants[29], constants[31])+power(states[8], constants[31]))))-constants[28]*(states[9]/(constants[30]+states[9])))
    rates[15] = constants[45]*(constants[46]*states[14]-constants[47]*(states[15]/(constants[48]+states[15])))
    algebraic[0] = constants[9]*states[2]-constants[10]*states[3]
    rates[2] = constants[5]*(constants[3]*states[0]*((power(constants[4], constants[6]))/(power(constants[4], constants[6])+power(states[1], constants[6])))-(constants[8]*(states[2]/(constants[7]+states[2]))+algebraic[0]))
    rates[3] = constants[5]*(algebraic[0]-constants[12]*(states[3]/(constants[11]+states[3])))
    algebraic[1] = constants[24]-states[5]
    algebraic[6] = constants[35]*algebraic[1]-constants[36]*states[10]*states[5]
    rates[5] = constants[21]*algebraic[6]
    algebraic[7] = constants[42]*states[8]-constants[41]*states[6]
    rates[8] = -(constants[21]*algebraic[7])
    rates[10] = constants[21]*((constants[32]*states[9]+algebraic[6])-constants[33]*(states[10]/(constants[34]+states[10])))
    algebraic[5] = constants[58]*((power(states[13], constants[70]))/(power(constants[68], constants[70])+power(states[13], constants[70])))
    algebraic[11] = constants[59]*(states[14]/(constants[69]+states[14]))
    rates[14] = constants[45]*(algebraic[5]-algebraic[11])
    algebraic[2] = constants[49]-states[11]
    algebraic[8] = constants[53]*((power(constants[60], constants[71]))/(power(constants[61], constants[71])+power(constants[60], constants[71])))*(algebraic[2]/(constants[62]+algebraic[2]))
    algebraic[13] = constants[54]*(states[11]/(constants[63]+states[11]))
    rates[11] = constants[45]*(algebraic[8]-algebraic[13])
    algebraic[3] = constants[50]-states[12]
    algebraic[9] = constants[55]*(states[11]/constants[49])*(algebraic[3]/(constants[65]+algebraic[3]))
    algebraic[14] = constants[52]*states[15]*(states[12]/(constants[64]+states[12]))
    rates[12] = constants[45]*(algebraic[9]-algebraic[14])
    algebraic[4] = constants[51]-states[13]
    algebraic[10] = constants[56]*(states[12]/constants[50])*(algebraic[4]/(constants[66]+algebraic[4]))
    algebraic[15] = constants[57]*(states[13]/(constants[67]+states[13]))
    rates[13] = constants[45]*(algebraic[10]-algebraic[15])
    algebraic[12] = constants[43]*(constants[40]/(constants[40]+constants[39]))*(states[6]/(constants[37]+states[6]))
    algebraic[16] = constants[44]*(states[7]/(constants[38]+states[7]))
    rates[6] = constants[21]*((constants[23]+algebraic[16]+algebraic[7])-(algebraic[12]*(algebraic[1]/constants[24])+constants[22]*states[6]))
    rates[7] = constants[21]*(algebraic[12]*(algebraic[1]/constants[24])-(algebraic[16]+constants[25]*states[7]))
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[9]*states[2]-constants[10]*states[3]
    algebraic[1] = constants[24]-states[5]
    algebraic[6] = constants[35]*algebraic[1]-constants[36]*states[10]*states[5]
    algebraic[7] = constants[42]*states[8]-constants[41]*states[6]
    algebraic[5] = constants[58]*((power(states[13], constants[70]))/(power(constants[68], constants[70])+power(states[13], constants[70])))
    algebraic[11] = constants[59]*(states[14]/(constants[69]+states[14]))
    algebraic[2] = constants[49]-states[11]
    algebraic[8] = constants[53]*((power(constants[60], constants[71]))/(power(constants[61], constants[71])+power(constants[60], constants[71])))*(algebraic[2]/(constants[62]+algebraic[2]))
    algebraic[13] = constants[54]*(states[11]/(constants[63]+states[11]))
    algebraic[3] = constants[50]-states[12]
    algebraic[9] = constants[55]*(states[11]/constants[49])*(algebraic[3]/(constants[65]+algebraic[3]))
    algebraic[14] = constants[52]*states[15]*(states[12]/(constants[64]+states[12]))
    algebraic[4] = constants[51]-states[13]
    algebraic[10] = constants[56]*(states[12]/constants[50])*(algebraic[4]/(constants[66]+algebraic[4]))
    algebraic[15] = constants[57]*(states[13]/(constants[67]+states[13]))
    algebraic[12] = constants[43]*(constants[40]/(constants[40]+constants[39]))*(states[6]/(constants[37]+states[6]))
    algebraic[16] = constants[44]*(states[7]/(constants[38]+states[7]))
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