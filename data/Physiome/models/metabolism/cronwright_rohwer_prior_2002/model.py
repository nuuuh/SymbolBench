# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 1
sizeConstants = 19
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_constants[12] = "F16BP in component F16BP (millimolar)"
    legend_states[0] = "G3P in component G3P (millimolar)"
    legend_algebraic[0] = "V_Gpd_p in component V_Gpd_p (flux)"
    legend_algebraic[1] = "V_Gpp_p in component V_Gpp_p (flux)"
    legend_constants[13] = "DHAP in component DHAP (millimolar)"
    legend_constants[14] = "ATP in component ATP (millimolar)"
    legend_constants[15] = "ADP in component ADP (millimolar)"
    legend_constants[16] = "NADH in component NADH (millimolar)"
    legend_constants[17] = "NAD in component NAD (millimolar)"
    legend_constants[18] = "Pi_ in component Pi (millimolar)"
    legend_constants[0] = "K_F16BP in component V_Gpd_p (millimolar)"
    legend_constants[1] = "K_ATP in component V_Gpd_p (millimolar)"
    legend_constants[2] = "K_ADP in component V_Gpd_p (millimolar)"
    legend_constants[3] = "K_NAD in component V_Gpd_p (millimolar)"
    legend_constants[4] = "K_NADH in component V_Gpd_p (millimolar)"
    legend_constants[5] = "K_G3P in component V_Gpd_p (millimolar)"
    legend_constants[6] = "K_DHAP in component V_Gpd_p (millimolar)"
    legend_constants[7] = "K_eq in component V_Gpd_p (dimensionless)"
    legend_constants[8] = "Vf in component V_Gpd_p (flux)"
    legend_constants[9] = "K_G3P in component V_Gpp_p (millimolar)"
    legend_constants[10] = "K_Pi in component V_Gpp_p (millimolar)"
    legend_constants[11] = "V in component V_Gpp_p (flux)"
    legend_rates[0] = "d/dt G3P in component G3P (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 24
    constants[0] = 4.8
    constants[1] = 0.73
    constants[2] = 2
    constants[3] = 0.93
    constants[4] = 0.023
    constants[5] = 1.2
    constants[6] = 0.54
    constants[7] = 1e4
    constants[8] = 36
    constants[9] = 3.5
    constants[10] = 1
    constants[11] = 18
    constants[12] = 0.00000
    constants[13] = 0.590000
    constants[14] = 2.37000
    constants[15] = 2.17000
    constants[16] = 1.87000
    constants[17] = 1.45000
    constants[18] = 2.17000
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = ((constants[8]/(constants[4]*constants[6]))*(constants[16]*constants[13]-(constants[17]*states[0])/constants[7]))/((1.00000+constants[12]/constants[0]+constants[14]/constants[1]+constants[15]/constants[2])*(1.00000+constants[16]/constants[4]+constants[17]/constants[3])*(1.00000+constants[13]/constants[6]+states[0]/constants[5]))
    algebraic[1] = ((constants[11]*states[0])/constants[9])/((1.00000+states[0]/constants[9])*(1.00000+constants[18]/constants[10]))
    rates[0] = -algebraic[1]+algebraic[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = ((constants[8]/(constants[4]*constants[6]))*(constants[16]*constants[13]-(constants[17]*states[0])/constants[7]))/((1.00000+constants[12]/constants[0]+constants[14]/constants[1]+constants[15]/constants[2])*(1.00000+constants[16]/constants[4]+constants[17]/constants[3])*(1.00000+constants[13]/constants[6]+states[0]/constants[5]))
    algebraic[1] = ((constants[11]*states[0])/constants[9])/((1.00000+states[0]/constants[9])*(1.00000+constants[18]/constants[10]))
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