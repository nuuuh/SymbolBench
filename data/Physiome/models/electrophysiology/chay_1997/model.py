# Size of variable arrays:
sizeAlgebraic = 15
sizeStates = 6
sizeConstants = 33
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "R in component membrane (millijoule_per_mole_kelvin)"
    legend_constants[1] = "T in component membrane (kelvin)"
    legend_constants[2] = "F in component membrane (coulomb_per_mole)"
    legend_constants[3] = "Cm in component membrane (microF_per_cm2)"
    legend_algebraic[11] = "i_K_dr in component delayed_rectifier_K_channel_current (nanoA_per_cm2)"
    legend_algebraic[12] = "i_K_Ca in component Ca_sensitive_K_current (nanoA_per_cm2)"
    legend_algebraic[13] = "i_K_ATP in component ATP_sensitive_K_current (nanoA_per_cm2)"
    legend_algebraic[4] = "i_fast in component fast_current (nanoA_per_cm2)"
    legend_algebraic[9] = "i_Ca in component calcium_current (nanoA_per_cm2)"
    legend_algebraic[10] = "i_NS in component cationic_nonselective_inward_current (nanoA_per_cm2)"
    legend_algebraic[14] = "i_NaL in component Na_leak_current (nanoA_per_cm2)"
    legend_constants[4] = "g_fast in component fast_current (microS_per_cm2)"
    legend_constants[5] = "V_fast in component fast_current (millivolt)"
    legend_algebraic[0] = "m_infinity in component fast_current_m_gate (dimensionless)"
    legend_states[1] = "h in component fast_current_h_gate (dimensionless)"
    legend_constants[6] = "Vm in component fast_current_m_gate (millivolt)"
    legend_constants[7] = "Sm in component fast_current_m_gate (millivolt)"
    legend_constants[8] = "lamda_h in component fast_current_h_gate (per_second)"
    legend_algebraic[5] = "tau_h in component fast_current_h_gate (second)"
    legend_algebraic[1] = "h_infinity in component fast_current_h_gate (dimensionless)"
    legend_constants[9] = "Vh in component fast_current_h_gate (millivolt)"
    legend_constants[10] = "Sh in component fast_current_h_gate (millivolt)"
    legend_constants[11] = "K_Ca in component calcium_current (micromolar)"
    legend_constants[12] = "P_Ca in component calcium_current (nanoA_per_micromolar_per_cm2)"
    legend_constants[13] = "Ca_o in component calcium_current (micromolar)"
    legend_states[2] = "Ca_i in component cytosolic_calcium (micromolar)"
    legend_algebraic[8] = "f_infinity in component calcium_current_f_gate (dimensionless)"
    legend_states[3] = "d in component calcium_current_d_gate (dimensionless)"
    legend_constants[14] = "lamda_d in component calcium_current_d_gate (per_second)"
    legend_algebraic[6] = "tau_d in component calcium_current_d_gate (second)"
    legend_algebraic[2] = "d_infinity in component calcium_current_d_gate (dimensionless)"
    legend_constants[15] = "Vd in component calcium_current_d_gate (millivolt)"
    legend_constants[16] = "Sd in component calcium_current_d_gate (millivolt)"
    legend_constants[17] = "g_NS in component cationic_nonselective_inward_current (microS_per_cm2)"
    legend_constants[18] = "K_NS in component cationic_nonselective_inward_current (micromolar)"
    legend_constants[19] = "VNS in component cationic_nonselective_inward_current (millivolt)"
    legend_states[4] = "Ca_lum in component cytosolic_calcium (micromolar)"
    legend_constants[20] = "V_K in component delayed_rectifier_K_channel_current (millivolt)"
    legend_constants[21] = "g_K_dr in component delayed_rectifier_K_channel_current (microS_per_cm2)"
    legend_states[5] = "n in component delayed_rectifier_K_channel_current_n_gate (dimensionless)"
    legend_constants[22] = "lamda_n in component delayed_rectifier_K_channel_current_n_gate (per_second)"
    legend_constants[23] = "Vn in component delayed_rectifier_K_channel_current_n_gate (millivolt)"
    legend_constants[24] = "Sn in component delayed_rectifier_K_channel_current_n_gate (millivolt)"
    legend_algebraic[3] = "n_infinity in component delayed_rectifier_K_channel_current_n_gate (dimensionless)"
    legend_algebraic[7] = "tau_n in component delayed_rectifier_K_channel_current_n_gate (second)"
    legend_constants[25] = "g_K_Ca in component Ca_sensitive_K_current (microS_per_cm2)"
    legend_constants[26] = "g_K_ATP in component ATP_sensitive_K_current (microS_per_cm2)"
    legend_constants[27] = "g_NaL in component Na_leak_current (microS_per_cm2)"
    legend_constants[28] = "V_Na in component Na_leak_current (millivolt)"
    legend_constants[29] = "k_rel in component cytosolic_calcium (per_second)"
    legend_constants[30] = "k_Ca in component cytosolic_calcium (per_second)"
    legend_constants[31] = "k_pump in component cytosolic_calcium (per_second)"
    legend_constants[32] = "omega in component cytosolic_calcium (micromolar_cm2_per_nanoA_per_second)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt h in component fast_current_h_gate (dimensionless)"
    legend_rates[3] = "d/dt d in component calcium_current_d_gate (dimensionless)"
    legend_rates[5] = "d/dt n in component delayed_rectifier_K_channel_current_n_gate (dimensionless)"
    legend_rates[2] = "d/dt Ca_i in component cytosolic_calcium (micromolar)"
    legend_rates[4] = "d/dt Ca_lum in component cytosolic_calcium (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -38.34146
    constants[0] = 8314
    constants[1] = 310
    constants[2] = 96485
    constants[3] = 1
    constants[4] = 600
    constants[5] = 80
    states[1] = 0.214723
    constants[6] = -25
    constants[7] = 9
    constants[8] = 12.5
    constants[9] = -48
    constants[10] = -7
    constants[11] = 1
    constants[12] = 2
    constants[13] = 2500
    states[2] = 0.6959466
    states[3] = 0.0031711238
    constants[14] = 2.5
    constants[15] = -10
    constants[16] = 5
    constants[17] = 5
    constants[18] = 50
    constants[19] = -20
    states[4] = 102.686
    constants[20] = -75
    constants[21] = 600
    states[5] = 0.1836403
    constants[22] = 12.5
    constants[23] = -18
    constants[24] = 14
    constants[25] = 5
    constants[26] = 2
    constants[27] = 0.3
    constants[28] = 80
    constants[29] = 0.2
    constants[30] = 7
    constants[31] = 30
    constants[32] = 0.2
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[4] = -constants[29]*(states[4]-states[2])+constants[31]*states[2]
    algebraic[5] = 1.00000/(constants[8]*(exp((constants[9]-states[0])/(2.00000*constants[10]))+exp((states[0]-constants[9])/(2.00000*constants[10]))))
    algebraic[1] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    rates[1] = (algebraic[1]-states[1])/algebraic[5]
    algebraic[6] = 1.00000/(constants[14]*(exp((constants[15]-states[0])/(2.00000*constants[16]))+exp((states[0]-constants[15])/(2.00000*constants[16]))))
    algebraic[2] = 1.00000/(1.00000+exp((constants[15]-states[0])/constants[16]))
    rates[3] = (algebraic[2]-states[3])/algebraic[6]
    algebraic[3] = 1.00000/(1.00000+exp((constants[23]-states[0])/constants[24]))
    algebraic[7] = 1.00000/(constants[22]*(exp((constants[23]-states[0])/(2.00000*constants[24]))+exp((states[0]-constants[23])/(2.00000*constants[24]))))
    rates[5] = (algebraic[3]-states[5])/algebraic[7]
    algebraic[8] = constants[11]/(constants[11]+states[2])
    algebraic[9] = (((constants[12]*states[3]*algebraic[8]*2.00000*constants[2]*states[0])/(constants[0]*constants[1]))*(constants[13]-states[2]*exp((2.00000*constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((2.00000*constants[2]*states[0])/(constants[0]*constants[1])))
    rates[2] = constants[29]*(states[4]-states[2])-(constants[32]*algebraic[9]+constants[30]*states[2]+constants[31]*states[2])
    algebraic[11] = constants[21]*(power(states[5], 4.00000))*(states[0]-constants[20])
    algebraic[12] = ((constants[25]*(power(states[2], 3.00000)))/(power(constants[11], 3.00000)+power(states[2], 3.00000)))*(states[0]-constants[20])
    algebraic[13] = constants[26]*(states[0]-constants[20])
    algebraic[0] = 1.00000/(1.00000+exp((constants[6]-states[0])/constants[7]))
    algebraic[4] = constants[4]*(power(algebraic[0], 3.00000))*states[1]*(states[0]-constants[5])
    algebraic[10] = ((constants[17]*(power(constants[18], 2.00000)))/(power(constants[18], 2.00000)+power(states[4], 2.00000)))*((states[0]-constants[19])/(1.00000-exp(0.100000*(constants[19]-states[0])))-10.0000)
    algebraic[14] = constants[27]*(states[0]-constants[28])
    rates[0] = -(algebraic[11]+algebraic[12]+algebraic[13]+algebraic[4]+algebraic[9]+algebraic[10]+algebraic[14])/constants[3]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[5] = 1.00000/(constants[8]*(exp((constants[9]-states[0])/(2.00000*constants[10]))+exp((states[0]-constants[9])/(2.00000*constants[10]))))
    algebraic[1] = 1.00000/(1.00000+exp((constants[9]-states[0])/constants[10]))
    algebraic[6] = 1.00000/(constants[14]*(exp((constants[15]-states[0])/(2.00000*constants[16]))+exp((states[0]-constants[15])/(2.00000*constants[16]))))
    algebraic[2] = 1.00000/(1.00000+exp((constants[15]-states[0])/constants[16]))
    algebraic[3] = 1.00000/(1.00000+exp((constants[23]-states[0])/constants[24]))
    algebraic[7] = 1.00000/(constants[22]*(exp((constants[23]-states[0])/(2.00000*constants[24]))+exp((states[0]-constants[23])/(2.00000*constants[24]))))
    algebraic[8] = constants[11]/(constants[11]+states[2])
    algebraic[9] = (((constants[12]*states[3]*algebraic[8]*2.00000*constants[2]*states[0])/(constants[0]*constants[1]))*(constants[13]-states[2]*exp((2.00000*constants[2]*states[0])/(constants[0]*constants[1]))))/(1.00000-exp((2.00000*constants[2]*states[0])/(constants[0]*constants[1])))
    algebraic[11] = constants[21]*(power(states[5], 4.00000))*(states[0]-constants[20])
    algebraic[12] = ((constants[25]*(power(states[2], 3.00000)))/(power(constants[11], 3.00000)+power(states[2], 3.00000)))*(states[0]-constants[20])
    algebraic[13] = constants[26]*(states[0]-constants[20])
    algebraic[0] = 1.00000/(1.00000+exp((constants[6]-states[0])/constants[7]))
    algebraic[4] = constants[4]*(power(algebraic[0], 3.00000))*states[1]*(states[0]-constants[5])
    algebraic[10] = ((constants[17]*(power(constants[18], 2.00000)))/(power(constants[18], 2.00000)+power(states[4], 2.00000)))*((states[0]-constants[19])/(1.00000-exp(0.100000*(constants[19]-states[0])))-10.0000)
    algebraic[14] = constants[27]*(states[0]-constants[28])
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