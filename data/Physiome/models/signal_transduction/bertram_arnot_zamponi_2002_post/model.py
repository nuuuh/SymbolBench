# Size of variable arrays:
sizeAlgebraic = 21
sizeStates = 11
sizeConstants = 19
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (ms)"
    legend_constants[0] = "R_gas_const in component parameters (millijoule_per_mole_kelvin)"
    legend_constants[1] = "Temp in component parameters (kelvin)"
    legend_constants[2] = "F in component parameters (coulomb_per_mole)"
    legend_states[0] = "R in component transmitter_release (dimensionless)"
    legend_constants[3] = "kr_plus in component transmitter_release (per_uM_per_ms)"
    legend_constants[4] = "kr_minus in component transmitter_release (per_ms)"
    legend_algebraic[19] = "Ca in component calcium_concentration (uM)"
    legend_algebraic[0] = "T in component transmitter_release (uM)"
    legend_constants[5] = "T_bar in component transmitter_release (uM)"
    legend_constants[6] = "Ca_ex in component calcium_concentration (uM)"
    legend_algebraic[8] = "Ca_open in component calcium_concentration (uM)"
    legend_constants[7] = "Dc in component calcium_concentration (micrometre2_per_second)"
    legend_constants[8] = "r in component calcium_concentration (nanometre)"
    legend_algebraic[5] = "sigma in component calcium_concentration (uM_per_ms)"
    legend_algebraic[1] = "i_V in component calcium_concentration (uA)"
    legend_constants[9] = "g_Ca in component calcium_concentration (pS)"
    legend_constants[10] = "P in component calcium_concentration (mV_per_uM)"
    legend_states[1] = "V_post in component membrane_post (mV)"
    legend_algebraic[16] = "O in component O (dimensionless)"
    legend_algebraic[10] = "alpha in component rate_constants (per_ms)"
    legend_algebraic[12] = "alpha_ in component rate_constants (per_ms)"
    legend_algebraic[14] = "beta in component rate_constants (per_ms)"
    legend_algebraic[17] = "beta_ in component rate_constants (per_ms)"
    legend_algebraic[20] = "kG_plus in component rate_constants (per_ms)"
    legend_states[2] = "b in component rate_constants (dimensionless)"
    legend_constants[11] = "kb_plus in component rate_constants (per_uM_per_ms)"
    legend_constants[12] = "kb_minus in component rate_constants (per_ms)"
    legend_constants[13] = "kG_minus in component rate_constants (per_ms)"
    legend_constants[14] = "kG2_minus in component rate_constants (per_ms)"
    legend_constants[15] = "kG3_minus in component rate_constants (per_ms)"
    legend_states[3] = "C1 in component C1 (dimensionless)"
    legend_states[4] = "C2 in component C2 (dimensionless)"
    legend_states[5] = "C_G1 in component C_G1 (dimensionless)"
    legend_states[6] = "C3 in component C3 (dimensionless)"
    legend_states[7] = "C_G2 in component C_G2 (dimensionless)"
    legend_states[8] = "C4 in component C4 (dimensionless)"
    legend_states[9] = "C_G3 in component C_G3 (dimensionless)"
    legend_algebraic[2] = "C_G in component O (dimensionless)"
    legend_constants[16] = "Cm in component membrane_post (uF_per_cm2)"
    legend_algebraic[3] = "i_syn in component synaptic_current (uA_per_cm2)"
    legend_algebraic[15] = "i_Na_post in component sodium_current_post (uA_per_cm2)"
    legend_algebraic[18] = "i_K_post in component potassium_current_post (uA_per_cm2)"
    legend_algebraic[6] = "i_leak_post in component leak_current_post (uA_per_cm2)"
    legend_constants[17] = "g_syn in component synaptic_current (mS_per_cm2)"
    legend_constants[18] = "V_syn in component synaptic_current (mV)"
    legend_algebraic[13] = "x_infinity in component sodium_current_post (dimensionless)"
    legend_algebraic[9] = "alpha_x in component sodium_current_post (dimensionless)"
    legend_algebraic[11] = "beta_x in component sodium_current_post (dimensionless)"
    legend_states[10] = "n_post in component potassium_current_n_gate_post (dimensionless)"
    legend_algebraic[4] = "alpha_n in component potassium_current_n_gate_post (per_ms)"
    legend_algebraic[7] = "beta_n in component potassium_current_n_gate_post (per_ms)"
    legend_rates[0] = "d/dt R in component transmitter_release (dimensionless)"
    legend_rates[2] = "d/dt b in component rate_constants (dimensionless)"
    legend_rates[3] = "d/dt C1 in component C1 (dimensionless)"
    legend_rates[4] = "d/dt C2 in component C2 (dimensionless)"
    legend_rates[6] = "d/dt C3 in component C3 (dimensionless)"
    legend_rates[8] = "d/dt C4 in component C4 (dimensionless)"
    legend_rates[5] = "d/dt C_G1 in component C_G1 (dimensionless)"
    legend_rates[7] = "d/dt C_G2 in component C_G2 (dimensionless)"
    legend_rates[9] = "d/dt C_G3 in component C_G3 (dimensionless)"
    legend_rates[1] = "d/dt V_post in component membrane_post (mV)"
    legend_rates[10] = "d/dt n_post in component potassium_current_n_gate_post (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 8314.41
    constants[1] = 310
    constants[2] = 96485
    states[0] = 0
    constants[3] = 0.15
    constants[4] = 2.5
    constants[5] = 4000.0
    constants[6] = 2000.0
    constants[7] = 220
    constants[8] = 10
    constants[9] = 1.2
    constants[10] = 0.006
    states[1] = -65
    states[2] = 0
    constants[11] = 2000.0
    constants[12] = 1.0
    constants[13] = 0.00025
    constants[14] = 0.016
    constants[15] = 1.024
    states[3] = 1
    states[4] = 0
    states[5] = 0
    states[6] = 0
    states[7] = 0
    states[8] = 0
    states[9] = 0
    constants[16] = 1.0
    constants[17] = 0.2
    constants[18] = 0
    states[10] = 0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = constants[5]*states[0]
    rates[2] = constants[11]*algebraic[0]*(1.00000-states[2])-constants[12]*states[2]
    algebraic[4] = (0.0200000*(states[1]+55.0000))/(1.00000-exp(-(states[1]+55.0000)/10.0000))
    algebraic[7] = 0.250000*exp(-(states[1]+65.0000)/80.0000)
    rates[10] = algebraic[4]*(1.00000-states[10])-algebraic[7]*states[10]
    algebraic[16] = ((((((1.00000-states[3])-states[4])-states[6])-states[8])-states[5])-states[7])-states[9]
    algebraic[10] = 0.450000*exp(states[1]/22.0000)
    algebraic[14] = 0.0150000*exp(-states[1]/14.0000)
    rates[8] = (2.00000*algebraic[10]*states[6]+4.00000*algebraic[14]*algebraic[16])-states[8]*(3.00000*algebraic[14]+algebraic[10])
    algebraic[3] = constants[17]*states[2]*(states[1]-constants[18])
    algebraic[9] = (0.200000*(states[1]+40.0000))/(1.00000-1.00000*exp(-(states[1]+40.0000)/10.0000))
    algebraic[11] = 8.00000*exp(1.00000/-(states[1]+65.0000/18.0000))
    algebraic[13] = algebraic[9]/(algebraic[9]+algebraic[11])
    algebraic[15] = 120.000*(power(algebraic[13], 3.00000))*(1.00000-states[10])*(states[1]-120.000)
    algebraic[18] = 36.0000*(power(states[10], 4.00000))*(states[1]+77.0000)
    algebraic[6] = 0.300000*(states[1]+54.0000)
    rates[1] = -(algebraic[15]+algebraic[18]+algebraic[6]+algebraic[3])/constants[16]
    algebraic[1] = (((constants[9]*constants[10]*2.00000*constants[2]*states[1])/(constants[0]*constants[1]))*constants[6])/(1.00000-exp((2.00000*constants[2]*states[1])/(constants[0]*constants[1])))
    algebraic[5] = -5.18200*algebraic[1]
    algebraic[8] = algebraic[5]/(2.00000*constants[7]*constants[8]* pi)
    algebraic[19] = algebraic[16]*algebraic[8]+0.100000
    rates[0] = constants[3]*algebraic[19]*(1.00000-states[0])-constants[4]*states[0]
    algebraic[20] = (3.00000*states[2])/(680.000+320.000*states[2])
    rates[3] = (algebraic[14]*states[4]+constants[13]*states[5])-states[3]*(4.00000*algebraic[10]+algebraic[20])
    rates[4] = (4.00000*algebraic[10]*states[3]+2.00000*algebraic[14]*states[6]+constants[14]*states[7])-states[4]*(algebraic[14]+3.00000*algebraic[10]+algebraic[20])
    rates[6] = (3.00000*algebraic[10]*states[4]+3.00000*algebraic[14]*states[8]+constants[15]*states[9])-states[6]*(2.00000*algebraic[14]+2.00000*algebraic[10]+algebraic[20])
    algebraic[12] = algebraic[10]/8.00000
    algebraic[17] = algebraic[14]*8.00000
    rates[5] = (algebraic[17]*states[7]+algebraic[20]*states[3])-states[5]*(4.00000*algebraic[12]+constants[13])
    rates[7] = (4.00000*algebraic[12]*states[5]+2.00000*algebraic[17]*states[9]+algebraic[20]*states[4])-states[7]*(algebraic[17]+3.00000*algebraic[12]+constants[14])
    rates[9] = (3.00000*algebraic[12]*states[7]+algebraic[20]*states[6])-states[9]*(2.00000*algebraic[17]+constants[15])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[5]*states[0]
    algebraic[4] = (0.0200000*(states[1]+55.0000))/(1.00000-exp(-(states[1]+55.0000)/10.0000))
    algebraic[7] = 0.250000*exp(-(states[1]+65.0000)/80.0000)
    algebraic[16] = ((((((1.00000-states[3])-states[4])-states[6])-states[8])-states[5])-states[7])-states[9]
    algebraic[10] = 0.450000*exp(states[1]/22.0000)
    algebraic[14] = 0.0150000*exp(-states[1]/14.0000)
    algebraic[3] = constants[17]*states[2]*(states[1]-constants[18])
    algebraic[9] = (0.200000*(states[1]+40.0000))/(1.00000-1.00000*exp(-(states[1]+40.0000)/10.0000))
    algebraic[11] = 8.00000*exp(1.00000/-(states[1]+65.0000/18.0000))
    algebraic[13] = algebraic[9]/(algebraic[9]+algebraic[11])
    algebraic[15] = 120.000*(power(algebraic[13], 3.00000))*(1.00000-states[10])*(states[1]-120.000)
    algebraic[18] = 36.0000*(power(states[10], 4.00000))*(states[1]+77.0000)
    algebraic[6] = 0.300000*(states[1]+54.0000)
    algebraic[1] = (((constants[9]*constants[10]*2.00000*constants[2]*states[1])/(constants[0]*constants[1]))*constants[6])/(1.00000-exp((2.00000*constants[2]*states[1])/(constants[0]*constants[1])))
    algebraic[5] = -5.18200*algebraic[1]
    algebraic[8] = algebraic[5]/(2.00000*constants[7]*constants[8]* pi)
    algebraic[19] = algebraic[16]*algebraic[8]+0.100000
    algebraic[20] = (3.00000*states[2])/(680.000+320.000*states[2])
    algebraic[12] = algebraic[10]/8.00000
    algebraic[17] = algebraic[14]*8.00000
    algebraic[2] = states[5]+states[7]+states[9]
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