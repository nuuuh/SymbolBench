# Size of variable arrays:
sizeAlgebraic = 9
sizeStates = 5
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
    legend_constants[0] = "tau_c in component nucleotides (second)"
    legend_constants[1] = "eta in component nucleotides (dimensionless)"
    legend_constants[2] = "v in component nucleotides (dimensionless)"
    legend_constants[3] = "k in component nucleotides (dimensionless)"
    legend_algebraic[0] = "phi in component nucleotides (dimensionless)"
    legend_states[0] = "ADP in component nucleotides (dimensionless)"
    legend_states[1] = "ATP in component nucleotides (dimensionless)"
    legend_constants[4] = "C_m in component membrane (femtofarad)"
    legend_algebraic[3] = "I_Ca in component Ca_current (femtoampere)"
    legend_algebraic[4] = "I_K in component K_current (femtoampere)"
    legend_algebraic[7] = "I_KCa in component Ca_activated_K_current (femtoampere)"
    legend_algebraic[8] = "I_KATP in component ATP_sensitive_K_current (femtoampere)"
    legend_states[2] = "V in component membrane (millivolt)"
    legend_constants[5] = "g_Ca_ in component Ca_current (picosiemens)"
    legend_constants[6] = "V_Ca in component Ca_current (millivolt)"
    legend_constants[7] = "v_m in component Ca_current (millivolt)"
    legend_constants[8] = "s_m in component Ca_current (millivolt)"
    legend_algebraic[1] = "m_infinity in component Ca_current (dimensionless)"
    legend_constants[9] = "g_K_ in component K_current (picosiemens)"
    legend_constants[10] = "V_K in component K_current (millivolt)"
    legend_states[3] = "n in component K_channel_activation (dimensionless)"
    legend_constants[11] = "g_KCa_ in component Ca_activated_K_current (picosiemens)"
    legend_constants[12] = "k_D in component Ca_activated_K_current (micromolar)"
    legend_states[4] = "c in component cytosolic_Ca (micromolar)"
    legend_algebraic[6] = "omega in component Ca_activated_K_current (dimensionless)"
    legend_constants[13] = "g_KATP_ in component ATP_sensitive_K_current (picosiemens)"
    legend_constants[14] = "tau_n in component K_channel_activation (millisecond)"
    legend_constants[15] = "v_n in component K_channel_activation (millivolt)"
    legend_constants[16] = "s_n in component K_channel_activation (millivolt)"
    legend_algebraic[2] = "n_infinity in component K_channel_activation (dimensionless)"
    legend_algebraic[5] = "J_mem in component Ca_influx (micromolar_per_ms)"
    legend_constants[17] = "f in component Ca_influx (dimensionless)"
    legend_constants[18] = "alpha in component Ca_influx (micromolar_per_fA_ms)"
    legend_constants[19] = "k_c in component Ca_influx (per_millisecond)"
    legend_rates[1] = "d/dt ATP in component nucleotides (dimensionless)"
    legend_rates[0] = "d/dt ADP in component nucleotides (dimensionless)"
    legend_rates[2] = "d/dt V in component membrane (millivolt)"
    legend_rates[3] = "d/dt n in component K_channel_activation (dimensionless)"
    legend_rates[4] = "d/dt c in component cytosolic_Ca (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1200
    constants[1] = 185
    constants[2] = 10
    constants[3] = 20
    states[0] = 0.085817
    states[1] = 2.1047
    constants[4] = 5300
    states[2] = -67.018
    constants[5] = 1200
    constants[6] = 25
    constants[7] = -20
    constants[8] = 12
    constants[9] = 3000
    constants[10] = -75
    states[3] = 0.00011
    constants[11] = 300
    constants[12] = 0.3
    states[4] = 0.15666
    constants[13] = 350
    constants[14] = 16
    constants[15] = -16
    constants[16] = 5.6
    constants[17] = 0.001
    constants[18] = 0.00000225
    constants[19] = 0.1
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = states[1]*(power(1.00000+constants[3]*states[0], 2.00000))
    rates[1] = (constants[2]-algebraic[0])/(1000.00*constants[0])
    rates[0] = (algebraic[0]-constants[1]*states[0])/(1000.00*constants[0])
    algebraic[2] = 1.00000/(1.00000+exp((constants[15]-states[2])/constants[16]))
    rates[3] = (algebraic[2]-states[3])/constants[14]
    algebraic[1] = 1.00000/(1.00000+exp((constants[7]-states[2])/constants[8]))
    algebraic[3] = constants[5]*algebraic[1]*(states[2]-constants[6])
    algebraic[5] = -constants[17]*(constants[18]*algebraic[3]+constants[19]*states[4])
    rates[4] = algebraic[5]
    algebraic[4] = constants[9]*states[3]*(states[2]-constants[10])
    algebraic[6] = 1.00000/(1.00000+constants[12]/states[4])
    algebraic[7] = constants[11]*algebraic[6]*(states[2]-constants[10])
    algebraic[8] = ((states[2]-constants[10])*constants[13])/states[1]
    rates[2] = -(algebraic[3]+algebraic[4]+algebraic[7]+algebraic[8])/constants[4]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[1]*(power(1.00000+constants[3]*states[0], 2.00000))
    algebraic[2] = 1.00000/(1.00000+exp((constants[15]-states[2])/constants[16]))
    algebraic[1] = 1.00000/(1.00000+exp((constants[7]-states[2])/constants[8]))
    algebraic[3] = constants[5]*algebraic[1]*(states[2]-constants[6])
    algebraic[5] = -constants[17]*(constants[18]*algebraic[3]+constants[19]*states[4])
    algebraic[4] = constants[9]*states[3]*(states[2]-constants[10])
    algebraic[6] = 1.00000/(1.00000+constants[12]/states[4])
    algebraic[7] = constants[11]*algebraic[6]*(states[2]-constants[10])
    algebraic[8] = ((states[2]-constants[10])*constants[13])/states[1]
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