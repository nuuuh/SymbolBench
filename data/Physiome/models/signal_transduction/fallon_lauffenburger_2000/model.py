# Size of variable arrays:
sizeAlgebraic = 0
sizeStates = 8
sizeConstants = 13
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "Rs in component Rs (number_per_cell)"
    legend_states[1] = "L in component L (picomolar)"
    legend_states[2] = "Cs in component Cs (number_per_cell)"
    legend_constants[0] = "Vs in component model_parameters (number_per_cell_minute)"
    legend_constants[10] = "kf in component model_parameters (second_order_rate_constant)"
    legend_constants[1] = "kr in component model_parameters (first_order_rate_constant)"
    legend_constants[2] = "kt in component model_parameters (first_order_rate_constant)"
    legend_constants[3] = "ksyn in component model_parameters (first_order_rate_constant)"
    legend_constants[4] = "ke in component model_parameters (first_order_rate_constant)"
    legend_states[3] = "Ri in component Ri (number_per_cell)"
    legend_states[4] = "Li in component Li (picomolar)"
    legend_states[5] = "Ci in component Ci (number_per_cell)"
    legend_constants[12] = "kfe in component model_parameters (second_order_rate_constant)"
    legend_constants[11] = "kre in component model_parameters (first_order_rate_constant)"
    legend_constants[5] = "kh in component model_parameters (first_order_rate_constant)"
    legend_constants[6] = "kx in component model_parameters (first_order_rate_constant)"
    legend_constants[7] = "Ve in component model_parameters (L_per_cell)"
    legend_constants[8] = "NA in component model_parameters (number_per_picomole)"
    legend_states[6] = "Ld in component Ld (number_per_cell)"
    legend_states[7] = "Y in component Y (cell_per_L)"
    legend_constants[9] = "IL2 in component model_parameters (dimensionless)"
    legend_rates[0] = "d/dt Rs in component Rs (number_per_cell)"
    legend_rates[2] = "d/dt Cs in component Cs (number_per_cell)"
    legend_rates[3] = "d/dt Ri in component Ri (number_per_cell)"
    legend_rates[5] = "d/dt Ci in component Ci (number_per_cell)"
    legend_rates[4] = "d/dt Li in component Li (picomolar)"
    legend_rates[6] = "d/dt Ld in component Ld (number_per_cell)"
    legend_rates[1] = "d/dt L in component L (picomolar)"
    legend_rates[7] = "d/dt Y in component Y (cell_per_L)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 1500
    states[1] = 10
    states[2] = 1
    constants[0] = 11
    constants[1] = 0.0138
    constants[2] = 0.007
    constants[3] = 0.0011
    constants[4] = 0.04
    states[3] = 300
    states[4] = 0.01
    states[5] = 1
    constants[5] = 0.035
    constants[6] = 0.15
    constants[7] = 1e-14
    constants[8] = 6.022e11
    states[6] = 1
    states[7] = 2.5e8
    constants[9] = 1
    constants[10] = custom_piecewise([equal(constants[9] , 1.00000), constants[1]/11.1000 , True, constants[1]/8.20000])
    constants[11] = custom_piecewise([equal(constants[9] , 1.00000), constants[1]*8.00000 , True, constants[1]*5.00000])
    constants[12] = custom_piecewise([equal(constants[9] , 1.00000), constants[11]/1000.00 , True, constants[11]/3000.00])
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = ((constants[1]+constants[3])*states[2]+constants[0])-(constants[10]*states[1]*states[0]+constants[2]*states[0])
    rates[2] = constants[10]*states[1]*states[0]-(constants[1]+constants[4])*states[2]
    rates[3] = (constants[11]*states[5]+constants[2]*states[0])-(constants[12]*states[4]*states[3]+constants[5]*states[3])
    rates[5] = (constants[12]*states[4]*states[3]+constants[4]*states[2])-(constants[11]+constants[5])*states[5]
    rates[4] = (constants[11]*states[5]-constants[12]*states[4]*states[3])/(constants[7]*constants[8])-constants[6]*states[4]
    rates[6] = constants[5]*states[5]
    rates[1] = (((constants[1]*states[2]+constants[6]*states[4]*constants[7]*constants[8])-constants[10]*states[1]*states[0])*states[7])/constants[8]
    rates[7] = custom_piecewise([greater((600.000*states[2])/(250.000+states[2])-200.000 , 0.00000), ((600.000*states[2])/(250.000+states[2])-200.000)*1000.00 , True, 0.00000])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
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