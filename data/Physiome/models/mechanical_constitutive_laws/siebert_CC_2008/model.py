# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 1
sizeConstants = 18
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "L_1 in component contraction (mm)"
    legend_constants[1] = "L_2 in component contraction (mm)"
    legend_constants[2] = "L_3 in component contraction (mm)"
    legend_constants[3] = "L_4 in component contraction (mm)"
    legend_constants[4] = "f_c in component contraction (newton)"
    legend_constants[5] = "v_max in component contraction (mm_per_second)"
    legend_constants[6] = "curv in component contraction (dimensionless)"
    legend_constants[7] = "k_1 in component contraction (newton)"
    legend_constants[8] = "k_2 in component contraction (per_mm)"
    legend_constants[9] = "F_1 in component contraction (newton)"
    legend_constants[10] = "d_LSEC1 in component contraction (mm)"
    legend_constants[11] = "k_sh in component contraction (dimensionless)"
    legend_constants[12] = "L_m in component contraction (mm)"
    legend_constants[13] = "F_im in component contraction (newton)"
    legend_constants[14] = "tau in component contraction (second)"
    legend_algebraic[4] = "v_cc in component contraction (mm_per_second)"
    legend_algebraic[5] = "f_v in component contraction (dimensionless)"
    legend_algebraic[0] = "f_L in component contraction (newton)"
    legend_algebraic[3] = "f_sec in component contraction (newton)"
    legend_algebraic[2] = "delta_L_sec in component contraction (mm)"
    legend_constants[15] = "delta_L_sec1 in component contraction (mm)"
    legend_constants[16] = "k in component contraction (newton_per_mm)"
    legend_states[0] = "L_cc in component contraction (mm)"
    legend_algebraic[1] = "L_mtc in component contraction (mm)"
    legend_constants[17] = "A in component contraction (dimensionless)"
    legend_rates[0] = "d/dt L_cc in component contraction (mm)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = -24
    constants[1] = -14
    constants[2] = 4
    constants[3] = 14
    constants[4] = 0.51
    constants[5] = -152
    constants[6] = 9.5
    constants[7] = 0.0064
    constants[8] = 0.414
    constants[9] = 8.9
    constants[10] = 2.9
    constants[11] = 2.5
    constants[12] = 2.6
    constants[13] = 20.6
    constants[14] = 0.065
    constants[15] = 4.1
    constants[16] = 8.2
    states[0] = 0.2
    constants[17] = 1
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = custom_piecewise([greater_equal(states[0] , constants[0]) & less_equal(states[0] , constants[1]), (constants[4]/(constants[1]-constants[0]))*(states[0]-constants[0]) , greater(states[0] , constants[1]) & less_equal(states[0] , 0.00000), ((1.00000-constants[4])/-constants[1])*(states[0]-constants[1]) , greater(states[0] , 0.00000) & less_equal(states[0] , constants[2]), 1.00000 , greater(states[0] , constants[2]) & less_equal(states[0] , constants[3]), (-1.00000/(constants[3]-constants[2]))*(states[0]-constants[2]) , True, float('nan')])
    algebraic[1] = custom_piecewise([less_equal(voi , 1.00000), 0.290000 , greater(voi , 1.00000) & less(voi , 5.00000), 0.220000 , True, 0.190000])
    algebraic[2] = (algebraic[1]-states[0])-constants[12]
    algebraic[3] = custom_piecewise([greater(algebraic[2] , 0.00000) & less(algebraic[2] , constants[15]), (constants[9]/(exp(constants[11])-1.00000))*(exp((constants[11]*algebraic[2])/constants[15])-1.00000) , less_equal(algebraic[2] , constants[15]), constants[9]+constants[16]*(algebraic[2]-constants[15]) , True, float('nan')])
    rootfind_0(voi, constants, rates, states, algebraic)
    rates[0] = algebraic[4]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(states[0] , constants[0]) & less_equal(states[0] , constants[1]), (constants[4]/(constants[1]-constants[0]))*(states[0]-constants[0]) , greater(states[0] , constants[1]) & less_equal(states[0] , 0.00000), ((1.00000-constants[4])/-constants[1])*(states[0]-constants[1]) , greater(states[0] , 0.00000) & less_equal(states[0] , constants[2]), 1.00000 , greater(states[0] , constants[2]) & less_equal(states[0] , constants[3]), (-1.00000/(constants[3]-constants[2]))*(states[0]-constants[2]) , True, float('nan')])
    algebraic[1] = custom_piecewise([less_equal(voi , 1.00000), 0.290000 , greater(voi , 1.00000) & less(voi , 5.00000), 0.220000 , True, 0.190000])
    algebraic[2] = (algebraic[1]-states[0])-constants[12]
    algebraic[3] = custom_piecewise([greater(algebraic[2] , 0.00000) & less(algebraic[2] , constants[15]), (constants[9]/(exp(constants[11])-1.00000))*(exp((constants[11]*algebraic[2])/constants[15])-1.00000) , less_equal(algebraic[2] , constants[15]), constants[9]+constants[16]*(algebraic[2]-constants[15]) , True, float('nan')])
    return algebraic

initialGuess0 = None
def rootfind_0(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess0
    if initialGuess0 is None: initialGuess0 = ones(2)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_0, initialGuess0, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess0 = soln
        algebraic[4] = soln[0]
        algebraic[5] = soln[1]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess0 = soln
            algebraic[4][i] = soln[0]
            algebraic[5][i] = soln[1]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 2)
    algebraic[4] = algebraicCandidate[0]
    algebraic[5] = algebraicCandidate[1]
    resid[0] = (algebraic[4]-(1.00000/algebraic[5])*(algebraic[3]/(constants[17]*algebraic[0]*constants[13])))
    resid[1] = (algebraic[5]-(constants[5]-algebraic[4])/(constants[5]+algebraic[4]*constants[6]))
    return resid

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