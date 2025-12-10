# Size of variable arrays:
sizeAlgebraic = 21
sizeStates = 11
sizeConstants = 20
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_constants[0] = "Cm in component membrane (picoF)"
    legend_constants[1] = "i_app in component membrane (picoA)"
    legend_algebraic[8] = "i_Na in component sodium_current (picoA)"
    legend_algebraic[14] = "i_Ca_L in component L_type_calcium_current (picoA)"
    legend_algebraic[15] = "i_Ca_T in component T_type_calcium_current (picoA)"
    legend_algebraic[16] = "i_K_dr in component delayed_rectifier_K_channel_current (picoA)"
    legend_algebraic[17] = "i_M in component M_like_K_current (picoA)"
    legend_algebraic[19] = "i_ir in component inward_rectifier_K_current (picoA)"
    legend_algebraic[20] = "i_d in component inward_leak_current (picoA)"
    legend_constants[2] = "g_Na in component sodium_current (nanoS)"
    legend_constants[3] = "V_Na in component sodium_current (millivolt)"
    legend_algebraic[0] = "O in component sodium_current (dimensionless)"
    legend_states[1] = "A in component A (dimensionless)"
    legend_constants[4] = "k1 in component reaction_constants (first_order_rate_constant)"
    legend_constants[5] = "k1_ in component reaction_constants (first_order_rate_constant)"
    legend_algebraic[1] = "alpha in component reaction_constants (first_order_rate_constant)"
    legend_algebraic[9] = "beta in component reaction_constants (first_order_rate_constant)"
    legend_states[2] = "D in component D (dimensionless)"
    legend_states[3] = "A_ in component A_ (dimensionless)"
    legend_constants[19] = "a in component reaction_constants (dimensionless)"
    legend_states[4] = "D_ in component D_ (dimensionless)"
    legend_constants[6] = "k2 in component reaction_constants (first_order_rate_constant)"
    legend_constants[7] = "k2_ in component reaction_constants (first_order_rate_constant)"
    legend_constants[8] = "V_Ca in component L_type_calcium_current (millivolt)"
    legend_constants[9] = "g_Ca_L in component L_type_calcium_current (nanoS)"
    legend_states[5] = "m in component L_type_calcium_current_m_gate (dimensionless)"
    legend_algebraic[2] = "m_infinity in component L_type_calcium_current_m_gate (dimensionless)"
    legend_algebraic[10] = "tau_m in component L_type_calcium_current_m_gate (millisecond)"
    legend_constants[10] = "Vh in component L_type_calcium_current_m_gate (millivolt)"
    legend_constants[11] = "g_Ca_T in component T_type_calcium_current (nanoS)"
    legend_states[6] = "m in component T_type_calcium_current_m_gate (dimensionless)"
    legend_states[7] = "h in component T_type_calcium_current_h_gate (dimensionless)"
    legend_algebraic[3] = "m_infinity in component T_type_calcium_current_m_gate (dimensionless)"
    legend_algebraic[11] = "tau_m in component T_type_calcium_current_m_gate (millisecond)"
    legend_algebraic[4] = "h_infinity in component T_type_calcium_current_h_gate (dimensionless)"
    legend_constants[12] = "tau_h in component T_type_calcium_current_h_gate (millisecond)"
    legend_constants[13] = "V_K in component delayed_rectifier_K_channel_current (millivolt)"
    legend_constants[14] = "g_K_dr in component delayed_rectifier_K_channel_current (nanoS)"
    legend_states[8] = "n in component delayed_rectifier_K_channel_current_n_gate (dimensionless)"
    legend_states[9] = "h in component delayed_rectifier_K_channel_current_h_gate (dimensionless)"
    legend_algebraic[5] = "n_infinity in component delayed_rectifier_K_channel_current_n_gate (dimensionless)"
    legend_algebraic[12] = "tau_n in component delayed_rectifier_K_channel_current_n_gate (millisecond)"
    legend_algebraic[6] = "h_infinity in component delayed_rectifier_K_channel_current_h_gate (dimensionless)"
    legend_constants[15] = "tau_h in component delayed_rectifier_K_channel_current_h_gate (millisecond)"
    legend_constants[16] = "g_M in component M_like_K_current (nanoS)"
    legend_states[10] = "n in component M_like_K_current_n_gate (dimensionless)"
    legend_algebraic[7] = "n_infinity in component M_like_K_current_n_gate (dimensionless)"
    legend_algebraic[13] = "tau_n in component M_like_K_current_n_gate (millisecond)"
    legend_constants[17] = "g_ir in component inward_rectifier_K_current (nanoS)"
    legend_algebraic[18] = "n_infinity in component inward_rectifier_K_current_n_gate (dimensionless)"
    legend_constants[18] = "g_d in component inward_leak_current (nanoS)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt A in component A (dimensionless)"
    legend_rates[3] = "d/dt A_ in component A_ (dimensionless)"
    legend_rates[2] = "d/dt D in component D (dimensionless)"
    legend_rates[4] = "d/dt D_ in component D_ (dimensionless)"
    legend_rates[5] = "d/dt m in component L_type_calcium_current_m_gate (dimensionless)"
    legend_rates[6] = "d/dt m in component T_type_calcium_current_m_gate (dimensionless)"
    legend_rates[7] = "d/dt h in component T_type_calcium_current_h_gate (dimensionless)"
    legend_rates[8] = "d/dt n in component delayed_rectifier_K_channel_current_n_gate (dimensionless)"
    legend_rates[9] = "d/dt h in component delayed_rectifier_K_channel_current_h_gate (dimensionless)"
    legend_rates[10] = "d/dt n in component M_like_K_current_n_gate (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -60
    constants[0] = 14
    constants[1] = 15
    constants[2] = 60
    constants[3] = 60
    states[1] = 1
    constants[4] = 0.3
    constants[5] = 0.03
    states[2] = 0
    states[3] = 0
    states[4] = 0
    constants[6] = 0.001
    constants[7] = 0.01
    constants[8] = 100
    constants[9] = 1.3
    states[5] = 0
    constants[10] = 40
    constants[11] = 0.94
    states[6] = 0
    states[7] = 0
    constants[12] = 22
    constants[13] = -80
    constants[14] = 20
    states[8] = 0
    states[9] = 0
    constants[15] = 1000
    constants[16] = 4
    states[10] = 0
    constants[17] = 1.71
    constants[18] = 0.044
    constants[19] = power((constants[4]*constants[7])/(constants[5]*constants[6]), 1.0/2)
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[4] = 1.00000/(1.00000+exp((states[0]+86.4000)/4.70000))
    rates[7] = (algebraic[4]-states[7])/constants[12]
    algebraic[6] = 0.700000/(1.00000+exp(-(states[0]+35.0000)/10.0000))+0.300000
    rates[9] = (algebraic[6]-states[9])/constants[15]
    algebraic[1] = 10.0000/(1.00000+exp(-(states[0]+6.00000)/10.0000))
    algebraic[9] = 10.0000/(1.00000+exp((states[0]+54.4000)/4.60000))
    rates[1] = (algebraic[1]*states[2]+constants[5]*states[3])-(algebraic[9]*states[1]+constants[4]*states[1])
    rates[3] = (algebraic[1]*constants[19]*states[4]+constants[4]*states[1])-((algebraic[9]/constants[19])*states[3]+constants[5]*states[3])
    rates[2] = (algebraic[9]*states[1]+constants[7]*states[4])-(algebraic[1]*states[2]+constants[6]*states[2])
    rates[4] = ((algebraic[9]/constants[19])*states[3]+constants[6]*states[2])-(algebraic[1]*constants[19]*states[4]+constants[7]*states[4])
    algebraic[2] = 1.00000/(1.00000+exp(-(states[0]-constants[10])/12.0000))
    algebraic[10] = 5.00000/(exp((states[0]+15.0000)/25.0000)+exp(-(states[0]+15.0000)/25.0000))
    rates[5] = (algebraic[2]-states[5])/algebraic[10]
    algebraic[3] = 1.00000/(1.00000+exp(-(states[0]-56.1000)/10.0000))
    algebraic[11] = 7.00000/(exp((states[0]+50.0000)/9.00000)+exp(-(states[0]+50.0000)/9.00000))+0.800000
    rates[6] = (algebraic[3]-states[6])/algebraic[11]
    algebraic[5] = 1.00000/(1.00000+exp(-(states[0]+25.0000)/15.0000))
    algebraic[12] = 15.0000/(exp((states[0]+30.0000)/15.0000)+exp(-(states[0]+30.0000)/15.0000))+1.00000
    rates[8] = (algebraic[5]-states[8])/algebraic[12]
    algebraic[7] = 1.00000/(1.00000+exp(-(states[0]+37.0000)/4.00000))
    algebraic[13] = 80.0000/(exp((states[0]+30.0000)/15.0000)+exp(-(states[0]+30.0000)/15.0000))
    rates[10] = (algebraic[7]-states[10])/algebraic[13]
    algebraic[0] = power(states[1], 3.00000)
    algebraic[8] = constants[2]*algebraic[0]*(states[0]-constants[3])
    algebraic[14] = constants[9]*(power(states[5], 2.00000))*(states[0]-constants[8])
    algebraic[15] = constants[11]*(power(states[6], 2.00000))*states[7]*(states[0]-constants[8])
    algebraic[16] = constants[14]*(power(states[8], 4.00000))*states[9]*(states[0]-constants[13])
    algebraic[17] = constants[16]*states[10]*(states[0]-constants[13])
    algebraic[18] = (0.800000*1.00000)/(1.00000+exp((states[0]+80.0000)/12.0000))+0.200000
    algebraic[19] = constants[17]*algebraic[18]*(states[0]-constants[13])
    algebraic[20] = constants[18]*(states[0]-constants[8])
    rates[0] = (constants[1]-(algebraic[8]+algebraic[14]+algebraic[15]+algebraic[16]+algebraic[17]+algebraic[19]+algebraic[20]))/constants[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[4] = 1.00000/(1.00000+exp((states[0]+86.4000)/4.70000))
    algebraic[6] = 0.700000/(1.00000+exp(-(states[0]+35.0000)/10.0000))+0.300000
    algebraic[1] = 10.0000/(1.00000+exp(-(states[0]+6.00000)/10.0000))
    algebraic[9] = 10.0000/(1.00000+exp((states[0]+54.4000)/4.60000))
    algebraic[2] = 1.00000/(1.00000+exp(-(states[0]-constants[10])/12.0000))
    algebraic[10] = 5.00000/(exp((states[0]+15.0000)/25.0000)+exp(-(states[0]+15.0000)/25.0000))
    algebraic[3] = 1.00000/(1.00000+exp(-(states[0]-56.1000)/10.0000))
    algebraic[11] = 7.00000/(exp((states[0]+50.0000)/9.00000)+exp(-(states[0]+50.0000)/9.00000))+0.800000
    algebraic[5] = 1.00000/(1.00000+exp(-(states[0]+25.0000)/15.0000))
    algebraic[12] = 15.0000/(exp((states[0]+30.0000)/15.0000)+exp(-(states[0]+30.0000)/15.0000))+1.00000
    algebraic[7] = 1.00000/(1.00000+exp(-(states[0]+37.0000)/4.00000))
    algebraic[13] = 80.0000/(exp((states[0]+30.0000)/15.0000)+exp(-(states[0]+30.0000)/15.0000))
    algebraic[0] = power(states[1], 3.00000)
    algebraic[8] = constants[2]*algebraic[0]*(states[0]-constants[3])
    algebraic[14] = constants[9]*(power(states[5], 2.00000))*(states[0]-constants[8])
    algebraic[15] = constants[11]*(power(states[6], 2.00000))*states[7]*(states[0]-constants[8])
    algebraic[16] = constants[14]*(power(states[8], 4.00000))*states[9]*(states[0]-constants[13])
    algebraic[17] = constants[16]*states[10]*(states[0]-constants[13])
    algebraic[18] = (0.800000*1.00000)/(1.00000+exp((states[0]+80.0000)/12.0000))+0.200000
    algebraic[19] = constants[17]*algebraic[18]*(states[0]-constants[13])
    algebraic[20] = constants[18]*(states[0]-constants[8])
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