# Size of variable arrays:
sizeAlgebraic = 4
sizeStates = 3
sizeConstants = 15
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_algebraic[0] = "Qv in component v (per_second)"
    legend_states[0] = "Vv in component v (mV)"
    legend_constants[0] = "tau_v in component v (second)"
    legend_constants[1] = "v_vm in component v (mV_second)"
    legend_algebraic[1] = "Qm in component m (per_second)"
    legend_constants[2] = "Qmax in component model_parameters (per_second)"
    legend_algebraic[3] = "D in component D (mV)"
    legend_constants[3] = "theta in component model_parameters (mV)"
    legend_constants[4] = "sigma in component model_parameters (mV)"
    legend_states[1] = "Vm in component m (mV)"
    legend_constants[5] = "tau_m in component m (second)"
    legend_constants[6] = "v_mv in component m (mV_second)"
    legend_constants[7] = "A in component m (mV)"
    legend_states[2] = "H in component H (nM)"
    legend_constants[8] = "chi in component H (hour)"
    legend_constants[9] = "mu in component H (nM_second)"
    legend_algebraic[2] = "C in component D (dimensionless)"
    legend_constants[10] = "c0 in component D (dimensionless)"
    legend_constants[14] = "omega in component D (per_hour)"
    legend_constants[11] = "alpha in component D (hour)"
    legend_constants[12] = "v_vc in component D (mV)"
    legend_constants[13] = "v_vh in component D (mV_per_nM)"
    legend_rates[0] = "d/dt Vv in component v (mV)"
    legend_rates[1] = "d/dt Vm in component m (mV)"
    legend_rates[2] = "d/dt H in component H (nM)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -6.0
    constants[0] = 10.0
    constants[1] = -1.9
    constants[2] = 100.0
    constants[3] = 10.0
    constants[4] = 3.0
    states[1] = 0.0
    constants[5] = 10.0
    constants[6] = -1.9
    constants[7] = 1.0
    states[2] = 15.0
    constants[8] = 11.0
    constants[9] = 3.6
    constants[10] = 1.0
    constants[11] = 16.0
    constants[12] = -3.2
    constants[13] = 0.19
    constants[14] = (2.00000* pi)/24.0000
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = constants[2]/(1.00000+exp(-(states[0]-constants[3])/constants[4]))
    rates[1] = ((constants[7]+constants[6]*algebraic[0])-states[1])/(constants[5]/3600.00)
    algebraic[1] = constants[2]/(1.00000+exp(-(states[1]-constants[3])/constants[4]))
    rates[2] = (constants[9]*algebraic[1]-states[2])/constants[8]
    algebraic[2] = constants[10]+cos(constants[14]*(voi-constants[11]))
    algebraic[3] = constants[12]*algebraic[2]+constants[13]*states[2]
    rates[0] = ((constants[1]*algebraic[1]+algebraic[3])-states[0])/(constants[0]/3600.00)
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[2]/(1.00000+exp(-(states[0]-constants[3])/constants[4]))
    algebraic[1] = constants[2]/(1.00000+exp(-(states[1]-constants[3])/constants[4]))
    algebraic[2] = constants[10]+cos(constants[14]*(voi-constants[11]))
    algebraic[3] = constants[12]*algebraic[2]+constants[13]*states[2]
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