# Size of variable arrays:
sizeAlgebraic = 13
sizeStates = 5
sizeConstants = 15
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
    legend_algebraic[0] = "Vs in component membrane (millivolt)"
    legend_constants[0] = "V_I in component membrane (millivolt)"
    legend_constants[1] = "V_K in component membrane (millivolt)"
    legend_constants[2] = "V_L in component membrane (millivolt)"
    legend_constants[3] = "V_H_Na in component membrane (millivolt)"
    legend_constants[4] = "V_H_K in component membrane (millivolt)"
    legend_constants[13] = "g_I in component membrane (milliS_per_microF)"
    legend_constants[5] = "g_K in component membrane (milliS_per_microF)"
    legend_constants[6] = "g_L in component membrane (milliS_per_microF)"
    legend_constants[14] = "g_T in component membrane (milliS_per_microF)"
    legend_constants[7] = "g_P in component membrane (milliS_per_microF)"
    legend_constants[8] = "Kp in component membrane (millimolar)"
    legend_states[1] = "c in component calcium_concentration (millimolar)"
    legend_algebraic[8] = "sI in component sI_gate (dimensionless)"
    legend_states[2] = "yI in component yI_gate (dimensionless)"
    legend_states[3] = "xT in component xT_gate (dimensionless)"
    legend_states[4] = "xK in component xK_gate (dimensionless)"
    legend_algebraic[1] = "alpha_m in component sI_gate (per_millisecond)"
    legend_algebraic[5] = "beta_m in component sI_gate (per_millisecond)"
    legend_algebraic[9] = "ZI in component yI_gate (dimensionless)"
    legend_algebraic[2] = "alpha_h in component yI_gate (per_millisecond)"
    legend_algebraic[6] = "beta_h in component yI_gate (per_millisecond)"
    legend_algebraic[11] = "tau_yI in component yI_gate (millisecond)"
    legend_algebraic[3] = "sT in component xT_gate (dimensionless)"
    legend_constants[9] = "tau_xT in component xT_gate (millisecond)"
    legend_constants[10] = "V_Ca in component calcium_concentration (millivolt)"
    legend_constants[11] = "rho in component calcium_concentration (per_millisecond)"
    legend_constants[12] = "K_c in component calcium_concentration (millimolar_per_millivolt)"
    legend_algebraic[4] = "alpha_n in component xK_gate (per_millisecond)"
    legend_algebraic[7] = "beta_n in component xK_gate (per_millisecond)"
    legend_algebraic[12] = "tau_xK in component xK_gate (millisecond)"
    legend_algebraic[10] = "sK in component xK_gate (dimensionless)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[2] = "d/dt yI in component yI_gate (dimensionless)"
    legend_rates[3] = "d/dt xT in component xT_gate (dimensionless)"
    legend_rates[1] = "d/dt c in component calcium_concentration (millimolar)"
    legend_rates[4] = "d/dt xK in component xK_gate (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -54
    constants[0] = 30.0
    constants[1] = -75.0
    constants[2] = -40.0
    constants[3] = 115.0
    constants[4] = -12.0
    constants[5] = 0.3
    constants[6] = 0.003
    constants[7] = 0.03
    constants[8] = 0.5
    states[1] = 0.1
    states[2] = 0.1
    states[3] = 0.1
    states[4] = 0.1
    constants[9] = 235.0
    constants[10] = 140.0
    constants[11] = 0.0003
    constants[12] = 0.0085
    constants[13] = 1.00000*((constants[3]-constants[4])/(constants[0]-constants[1]))
    constants[14] = 1.00000*((constants[3]*constants[1]-constants[0]*constants[4])/(constants[0]-constants[1]))
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = constants[11]*(constants[12]*states[3]*(constants[10]-states[0])-states[1])
    algebraic[0] = 1.00000*constants[13]*states[0]+1.00000*constants[14]
    algebraic[3] = 1.00000/(exp(0.150000*(-50.0000-algebraic[0]))+1.00000)
    rates[3] = (algebraic[3]-states[3])/constants[9]
    algebraic[1] = (0.100000*(50.0000-algebraic[0]))/-exp((50.0000-algebraic[0])/10.0000)
    algebraic[5] = 4.00000*exp((25.0000-algebraic[0])/18.0000)
    algebraic[8] = algebraic[1]/(algebraic[1]+algebraic[5])
    rates[0] = (constants[13]*(power(algebraic[8], 3.00000))*states[2]+constants[14]*states[3])*(constants[0]-states[0])+(constants[5]*(power(states[4], 4.00000))+constants[7]*states[1]*(power(constants[8]+states[1], -1.00000)))*(constants[1]-states[0])+constants[6]*(constants[2]-states[0])
    algebraic[2] = 0.0700000*exp((25.0000-algebraic[0])/20.0000)
    algebraic[6] = 1.00000/(exp((55.0000-algebraic[0])/10.0000)+1.00000)
    algebraic[9] = algebraic[2]/(algebraic[2]+algebraic[6])
    algebraic[11] = 12.5000/(algebraic[2]+algebraic[6])
    rates[2] = (algebraic[9]-states[2])/algebraic[11]
    algebraic[4] = (0.0100000*(55.0000-algebraic[0]))/(exp((55.0000-algebraic[0])/10.0000)-1.00000)
    algebraic[7] = 0.125000*exp((45.0000-algebraic[0])/80.0000)
    algebraic[12] = 12.5000/(algebraic[4]+algebraic[7])
    algebraic[10] = algebraic[4]/(algebraic[4]+algebraic[7])
    rates[4] = (algebraic[10]-states[4])/algebraic[12]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = 1.00000*constants[13]*states[0]+1.00000*constants[14]
    algebraic[3] = 1.00000/(exp(0.150000*(-50.0000-algebraic[0]))+1.00000)
    algebraic[1] = (0.100000*(50.0000-algebraic[0]))/-exp((50.0000-algebraic[0])/10.0000)
    algebraic[5] = 4.00000*exp((25.0000-algebraic[0])/18.0000)
    algebraic[8] = algebraic[1]/(algebraic[1]+algebraic[5])
    algebraic[2] = 0.0700000*exp((25.0000-algebraic[0])/20.0000)
    algebraic[6] = 1.00000/(exp((55.0000-algebraic[0])/10.0000)+1.00000)
    algebraic[9] = algebraic[2]/(algebraic[2]+algebraic[6])
    algebraic[11] = 12.5000/(algebraic[2]+algebraic[6])
    algebraic[4] = (0.0100000*(55.0000-algebraic[0]))/(exp((55.0000-algebraic[0])/10.0000)-1.00000)
    algebraic[7] = 0.125000*exp((45.0000-algebraic[0])/80.0000)
    algebraic[12] = 12.5000/(algebraic[4]+algebraic[7])
    algebraic[10] = algebraic[4]/(algebraic[4]+algebraic[7])
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