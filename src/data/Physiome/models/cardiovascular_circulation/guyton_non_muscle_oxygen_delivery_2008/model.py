# Size of variable arrays:
sizeAlgebraic = 9
sizeStates = 1
sizeConstants = 6
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_constants[0] = "BFN in component non_muscle_O2_delivery (L_per_minute)"
    legend_constants[1] = "OVA in component non_muscle_O2_delivery (mL_per_L)"
    legend_constants[2] = "HM in component non_muscle_O2_delivery (dimensionless)"
    legend_constants[3] = "AOM in component non_muscle_O2_delivery (dimensionless)"
    legend_constants[5] = "O2ARTN in component NM_O2_blood_supply (mL_per_minute)"
    legend_algebraic[4] = "DOB in component delivery_of_O2_to_NM_tissues (mL_per_minute)"
    legend_algebraic[5] = "POV in component NM_venous_O2_content (mmHg)"
    legend_algebraic[6] = "OSV in component NM_venous_O2_content (dimensionless)"
    legend_algebraic[1] = "POT in component pressure_of_O2_in_NM_tissue_cells (mmHg)"
    legend_algebraic[3] = "MO2 in component O2_consumption_by_NM_tissue (mL_per_minute)"
    legend_constants[4] = "O2M in component parameter_values (mL_per_minute)"
    legend_algebraic[2] = "P1O in component O2_consumption_by_NM_tissue (mmHg)"
    legend_algebraic[0] = "QO2 in component volume_of_O2_in_NM_tissue (mL)"
    legend_algebraic[8] = "DO2N in component volume_of_O2_in_NM_tissue (mL_per_minute)"
    legend_algebraic[7] = "DO2N1 in component volume_of_O2_in_NM_tissue (mL_per_minute)"
    legend_states[0] = "QO2T in component volume_of_O2_in_NM_tissue (mL)"
    legend_rates[0] = "d/dt QO2T in component volume_of_O2_in_NM_tissue (mL)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 2.79521
    constants[1] = 204.497
    constants[2] = 40.0381
    constants[3] = 1.00002
    constants[4] = 164
    states[0] = 72.2362
    constants[5] = constants[1]*constants[0]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = custom_piecewise([less(states[0] , 0.00000), 0.00000 , True, states[0]])
    algebraic[1] = algebraic[0]*0.486110
    rootfind_0(voi, constants, rates, states, algebraic)
    algebraic[2] = custom_piecewise([greater(algebraic[1] , 35.0000), 35.0000 , True, algebraic[1]])
    algebraic[3] = constants[3]*constants[4]*(1.00000-(power(35.0001-algebraic[2], 3.00000))/42875.0)
    algebraic[7] = algebraic[4]-algebraic[3]
    algebraic[8] = custom_piecewise([less(algebraic[0] , 6.00000) & less(algebraic[7] , 0.00000), algebraic[7]*0.100000 , True, algebraic[7]])
    rates[0] = algebraic[8]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less(states[0] , 0.00000), 0.00000 , True, states[0]])
    algebraic[1] = algebraic[0]*0.486110
    algebraic[2] = custom_piecewise([greater(algebraic[1] , 35.0000), 35.0000 , True, algebraic[1]])
    algebraic[3] = constants[3]*constants[4]*(1.00000-(power(35.0001-algebraic[2], 3.00000))/42875.0)
    algebraic[7] = algebraic[4]-algebraic[3]
    algebraic[8] = custom_piecewise([less(algebraic[0] , 6.00000) & less(algebraic[7] , 0.00000), algebraic[7]*0.100000 , True, algebraic[7]])
    return algebraic

initialGuess0 = None
def rootfind_0(voi, constants, rates, states, algebraic):
    """Calculate values of algebraic variables for DAE"""
    from scipy.optimize import fsolve
    global initialGuess0
    if initialGuess0 is None: initialGuess0 = ones(3)*0.1
    if not iterable(voi):
        soln = fsolve(residualSN_0, initialGuess0, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess0 = soln
        algebraic[4] = soln[0]
        algebraic[5] = soln[1]
        algebraic[6] = soln[2]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess0 = soln
            algebraic[4][i] = soln[0]
            algebraic[5][i] = soln[1]
            algebraic[6][i] = soln[2]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 3)
    algebraic[4] = algebraicCandidate[0]
    algebraic[5] = algebraicCandidate[1]
    algebraic[6] = algebraicCandidate[2]
    resid[0] = (algebraic[6]-(constants[5]-algebraic[4])/(constants[2]*5.25000*constants[0]))
    resid[1] = (algebraic[5]-algebraic[6]*57.1400)
    resid[2] = (algebraic[4]-(algebraic[5]-algebraic[1])*12.8570*constants[0])
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