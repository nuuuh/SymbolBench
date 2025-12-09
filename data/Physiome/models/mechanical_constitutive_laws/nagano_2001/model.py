# Size of variable arrays:
sizeAlgebraic = 8
sizeStates = 1
sizeConstants = 12
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_algebraic[0] = "F_isom in component contraction (newton)"
    legend_constants[11] = "c in component contraction (newton)"
    legend_states[0] = "L_ce in component contraction (metre)"
    legend_constants[0] = "L_ce_opt in component contraction (metre)"
    legend_algebraic[1] = "L in component contraction (metre)"
    legend_constants[1] = "width in component contraction (metre)"
    legend_constants[2] = "Factor in component contraction (per_second)"
    legend_constants[3] = "A_REL in component contraction (newton)"
    legend_constants[4] = "B_REL in component contraction (dimensionless)"
    legend_algebraic[7] = "v_ce in component contraction (metre_per_second)"
    legend_algebraic[5] = "F in component contraction (newton)"
    legend_constants[5] = "F_max in component contraction (newton)"
    legend_constants[6] = "q in component contraction (dimensionless)"
    legend_algebraic[4] = "c1 in component contraction (per_second)"
    legend_algebraic[2] = "c2 in component contraction (newton)"
    legend_algebraic[6] = "c3 in component contraction (per_newton_second)"
    legend_constants[7] = "slope in component contraction (newton)"
    legend_constants[8] = "F_asympt in component contraction (dimensionless)"
    legend_algebraic[3] = "L_see in component contraction (metre)"
    legend_constants[9] = "L_slack in component contraction (metre)"
    legend_constants[10] = "alpha in component contraction (newton_per_metre)"
    legend_rates[0] = "d/dt L_ce in component contraction (metre)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.05
    constants[0] = 0.055
    constants[1] = 0.888
    constants[2] = 1
    constants[3] = 0.41
    constants[4] = 5.2
    constants[5] = 3277.4
    constants[6] = 1
    constants[7] = 2
    constants[8] = 1.5
    constants[9] = 0.42
    constants[10] = 1449.027
    constants[11] = -1.00000/(power(constants[1], 2.00000))
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = (constants[11]*(power(states[0]/constants[0], 2.00000))-(2.00000*constants[11]*states[0])/constants[0])+constants[11]+1.00000
    algebraic[1] = custom_piecewise([less_equal(voi , 1.00000), 1.00000 , greater(voi , 1.00000) & less(voi , 5.00000), 0.920000 , True, 0.900000])
    algebraic[3] = algebraic[1]-states[0]
    algebraic[5] = constants[10]*(algebraic[3]-constants[9])
    algebraic[7] = -constants[2]*states[0]*(((algebraic[0]+constants[3])*constants[4])/(1.00000*(algebraic[5]/(constants[5]*constants[6]))+constants[3])-constants[4])
    rates[0] = algebraic[7]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = (constants[11]*(power(states[0]/constants[0], 2.00000))-(2.00000*constants[11]*states[0])/constants[0])+constants[11]+1.00000
    algebraic[1] = custom_piecewise([less_equal(voi , 1.00000), 1.00000 , greater(voi , 1.00000) & less(voi , 5.00000), 0.920000 , True, 0.900000])
    algebraic[3] = algebraic[1]-states[0]
    algebraic[5] = constants[10]*(algebraic[3]-constants[9])
    algebraic[7] = -constants[2]*states[0]*(((algebraic[0]+constants[3])*constants[4])/(1.00000*(algebraic[5]/(constants[5]*constants[6]))+constants[3])-constants[4])
    algebraic[2] = algebraic[0]*constants[8]
    algebraic[4] = (constants[2]*constants[4]*(power(algebraic[0]+algebraic[2], 2.00000)))/((algebraic[0]+constants[3])*constants[7])
    algebraic[6] = algebraic[4]/(algebraic[0]+algebraic[2])
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