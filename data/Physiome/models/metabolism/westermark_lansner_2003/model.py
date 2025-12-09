# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 3
sizeConstants = 28
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[27] = "v_GK in component v_GK (flux)"
    legend_constants[0] = "V_GK_min in component v_GK (enzyme_activity)"
    legend_constants[24] = "V_GK in component v_GK (enzyme_activity)"
    legend_constants[1] = "Sgk in component v_GK (millimolar)"
    legend_constants[2] = "h_GK in component v_GK (dimensionless)"
    legend_constants[26] = "Glc in component Glc (millimolar)"
    legend_constants[3] = "min_to_sec in component model_parameters (dimensionless)"
    legend_constants[4] = "dw_per_ml in component model_parameters (dimensionless)"
    legend_algebraic[1] = "v_PFK in component v_PFK (flux)"
    legend_constants[5] = "V_PFK_min in component v_PFK (enzyme_activity)"
    legend_constants[22] = "V_PFK in component v_PFK (enzyme_activity)"
    legend_constants[6] = "Spfk in component v_PFK (millimolar)"
    legend_constants[7] = "Sfba in component v_PFK (millimolar)"
    legend_constants[8] = "Xpfk in component v_PFK (millimolar)"
    legend_constants[9] = "hx in component v_PFK (dimensionless)"
    legend_constants[10] = "alpha in component v_PFK (dimensionless)"
    legend_constants[11] = "h_PFK in component v_PFK (dimensionless)"
    legend_constants[12] = "h_act in component v_PFK (dimensionless)"
    legend_states[0] = "FBP in component FBP (millimolar)"
    legend_algebraic[0] = "F6P in component F6P (millimolar)"
    legend_algebraic[5] = "v_FBA in component v_FBA (flux)"
    legend_constants[13] = "V_FBA_min in component v_FBA (enzyme_activity)"
    legend_constants[23] = "V_FBA in component v_FBA (enzyme_activity)"
    legend_constants[14] = "Qfba in component v_FBA (millimolar)"
    legend_constants[15] = "Sfba in component v_FBA (millimolar)"
    legend_constants[16] = "Pfba in component v_FBA (millimolar)"
    legend_constants[17] = "Keq_FBA in component v_FBA (millimolar)"
    legend_algebraic[2] = "G3P in component G3P (millimolar)"
    legend_algebraic[4] = "DHAP in component DHAP (millimolar)"
    legend_algebraic[3] = "v_GAPDH in component v_GAPDH (flux)"
    legend_constants[18] = "V_GAPDH_min in component v_GAPDH (enzyme_activity)"
    legend_constants[25] = "V_GAPDH in component v_GAPDH (enzyme_activity)"
    legend_constants[19] = "Sgapdh in component v_GAPDH (millimolar)"
    legend_states[1] = "G6P_F6P in component G6P_F6P (millimolar)"
    legend_constants[20] = "Keq_GPI in component F6P (dimensionless)"
    legend_states[2] = "DHAP_G3P in component DHAP_G3P (millimolar)"
    legend_constants[21] = "Keq_TPI in component G3P (dimensionless)"
    legend_rates[1] = "d/dt G6P_F6P in component G6P_F6P (millimolar)"
    legend_rates[0] = "d/dt FBP in component FBP (millimolar)"
    legend_rates[2] = "d/dt DHAP_G3P in component DHAP_G3P (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 10.0
    constants[1] = 8.0
    constants[2] = 1.7
    constants[3] = 60.0
    constants[4] = 0.3333
    constants[5] = 100.0
    constants[6] = 4.0
    constants[7] = 0.005
    constants[8] = 0.01
    constants[9] = 2.5
    constants[10] = 5.0
    constants[11] = 2.5
    constants[12] = 1.0
    states[0] = 0.00063612
    constants[13] = 25.0
    constants[14] = 0.275
    constants[15] = 0.005
    constants[16] = 0.5
    constants[17] = 0.1
    constants[18] = 250.0
    constants[19] = 0.005
    states[1] = 3.71728
    constants[20] = 0.3
    states[2] = 0.00262966
    constants[21] = 0.045455
    constants[22] = (constants[5]*constants[4])/constants[3]
    constants[23] = (constants[13]*constants[4])/constants[3]
    constants[24] = (constants[0]*constants[4])/constants[3]
    constants[25] = (constants[18]*constants[4])/constants[3]
    constants[26] = 10.0000
    constants[27] = (constants[24]*(power(constants[26]/constants[1], constants[2])))/(1.00000+power(constants[26]/constants[1], constants[2]))
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = (states[1]*constants[20])/(1.00000+constants[20])
    algebraic[1] = (constants[22]*(power(algebraic[0]/constants[6], constants[11]-(constants[11]-constants[12])*((states[0]/constants[7])/(1.00000+states[0]/constants[7])))))/(power(algebraic[0]/constants[6], constants[11]-(constants[11]-constants[12])*((states[0]/constants[7])/(1.00000+states[0]/constants[7])))+(1.00000+power(states[0]/constants[8], constants[9]))/(1.00000+(power(constants[10], constants[11]-(constants[11]-constants[12])*((states[0]/constants[7])/(1.00000+states[0]/constants[7]))))*(power(states[0]/constants[8], constants[9]))))
    rates[1] = constants[27]-algebraic[1]
    algebraic[2] = (states[2]*constants[21])/(1.00000+constants[21])
    algebraic[4] = states[2]-algebraic[2]
    algebraic[5] = (constants[23]*(states[0]/constants[15]-(algebraic[2]*algebraic[4])/(constants[16]*constants[14]*constants[17])))/(1.00000+states[0]/constants[15]+algebraic[4]/constants[14]+(algebraic[2]*algebraic[4])/(constants[16]*constants[14]))
    rates[0] = algebraic[1]-algebraic[5]
    algebraic[3] = (constants[25]*algebraic[2])/(constants[19]+algebraic[2])
    rates[2] = 2.00000*algebraic[5]-algebraic[3]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = (states[1]*constants[20])/(1.00000+constants[20])
    algebraic[1] = (constants[22]*(power(algebraic[0]/constants[6], constants[11]-(constants[11]-constants[12])*((states[0]/constants[7])/(1.00000+states[0]/constants[7])))))/(power(algebraic[0]/constants[6], constants[11]-(constants[11]-constants[12])*((states[0]/constants[7])/(1.00000+states[0]/constants[7])))+(1.00000+power(states[0]/constants[8], constants[9]))/(1.00000+(power(constants[10], constants[11]-(constants[11]-constants[12])*((states[0]/constants[7])/(1.00000+states[0]/constants[7]))))*(power(states[0]/constants[8], constants[9]))))
    algebraic[2] = (states[2]*constants[21])/(1.00000+constants[21])
    algebraic[4] = states[2]-algebraic[2]
    algebraic[5] = (constants[23]*(states[0]/constants[15]-(algebraic[2]*algebraic[4])/(constants[16]*constants[14]*constants[17])))/(1.00000+states[0]/constants[15]+algebraic[4]/constants[14]+(algebraic[2]*algebraic[4])/(constants[16]*constants[14]))
    algebraic[3] = (constants[25]*algebraic[2])/(constants[19]+algebraic[2])
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