# Size of variable arrays:
sizeAlgebraic = 9
sizeStates = 3
sizeConstants = 8
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component model (minute)"
    legend_states[0] = "NFATP_cyt in component model (molecule)"
    legend_states[1] = "NFAT_cyt in component model (molecule)"
    legend_states[2] = "NFAT_nuc in component model (molecule)"
    legend_algebraic[0] = "NFAT_tot in component model (molecule)"
    legend_constants[0] = "k1_unstim in component model (per_minute)"
    legend_constants[1] = "k1_stim in component model (per_minute)"
    legend_algebraic[2] = "k1 in component model (per_minute)"
    legend_constants[2] = "k2 in component model (per_minute)"
    legend_constants[3] = "k3 in component model (per_minute)"
    legend_constants[4] = "k4 in component model (per_minute)"
    legend_constants[5] = "stim_wavelength in component model (minute)"
    legend_constants[6] = "stim_duration in component model (minute)"
    legend_algebraic[1] = "stim_on in component model (dimensionless)"
    legend_constants[7] = "time_before_stim in component model (minute)"
    legend_algebraic[6] = "Jdephosphorylation in component model (molecules_per_minute)"
    legend_algebraic[7] = "Jtranslocate in component model (molecules_per_minute)"
    legend_algebraic[8] = "Jexport in component model (molecules_per_minute)"
    legend_algebraic[3] = "percentage_NFAT_cyt in component model (dimensionless)"
    legend_algebraic[4] = "percentage_NFATP_cyt in component model (dimensionless)"
    legend_algebraic[5] = "percentage_NFAT_nuc in component model (dimensionless)"
    legend_rates[0] = "d/dt NFATP_cyt in component model (molecule)"
    legend_rates[1] = "d/dt NFAT_cyt in component model (molecule)"
    legend_rates[2] = "d/dt NFAT_nuc in component model (molecule)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 5000
    states[1] = 0
    states[2] = 0
    constants[0] = 0
    constants[1] = 0.359
    constants[2] = 0.147
    constants[3] = 0.06
    constants[4] = 0.035
    constants[5] = 3
    constants[6] = 0.5
    constants[7] = 1
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = custom_piecewise([greater_equal(voi , constants[7]) & less_equal( voi-constants[7] % constants[5] , constants[6]), 1.00000 , True, 0.00000])
    algebraic[2] = custom_piecewise([equal(algebraic[1] , 1.00000), constants[1] , True, constants[0]])
    algebraic[6] = algebraic[2]*states[0]-constants[2]*states[1]
    algebraic[7] = constants[3]*states[1]
    rates[1] = algebraic[6]-algebraic[7]
    algebraic[8] = constants[4]*states[2]
    rates[0] = algebraic[8]-algebraic[6]
    rates[2] = algebraic[7]-algebraic[8]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = custom_piecewise([greater_equal(voi , constants[7]) & less_equal( voi-constants[7] % constants[5] , constants[6]), 1.00000 , True, 0.00000])
    algebraic[2] = custom_piecewise([equal(algebraic[1] , 1.00000), constants[1] , True, constants[0]])
    algebraic[6] = algebraic[2]*states[0]-constants[2]*states[1]
    algebraic[7] = constants[3]*states[1]
    algebraic[8] = constants[4]*states[2]
    algebraic[0] = states[0]+states[1]+states[2]
    algebraic[3] = (states[1]*100.000)/algebraic[0]
    algebraic[4] = (states[0]*100.000)/algebraic[0]
    algebraic[5] = (states[2]*100.000)/algebraic[0]
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

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