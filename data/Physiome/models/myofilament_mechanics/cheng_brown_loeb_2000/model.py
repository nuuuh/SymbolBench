# Size of variable arrays:
sizeAlgebraic = 13
sizeStates = 7
sizeConstants = 41
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_constants[40] = "F_SE in component F_SE (newton)"
    legend_constants[0] = "cT in component F_SE (dimensionless)"
    legend_constants[1] = "kT in component F_SE (dimensionless)"
    legend_constants[2] = "LT_r in component F_SE (dimensionless)"
    legend_constants[3] = "LT in component user_defined_constants (dimensionless)"
    legend_constants[4] = "F_max in component user_defined_constants (newton)"
    legend_algebraic[1] = "F_PE1 in component F_PE1 (dimensionless)"
    legend_constants[5] = "c1 in component F_PE1 (dimensionless)"
    legend_constants[6] = "k1 in component F_PE1 (dimensionless)"
    legend_constants[7] = "L_r1 in component F_PE1 (dimensionless)"
    legend_constants[8] = "eta in component F_PE1 (millisecond)"
    legend_states[0] = "L in component L (dimensionless)"
    legend_constants[9] = "L_max in component user_defined_constants (dimensionless)"
    legend_states[1] = "V in component V (first_order_rate_constant)"
    legend_algebraic[4] = "F_PE2 in component F_PE2 (dimensionless)"
    legend_constants[10] = "c2 in component F_PE2 (dimensionless)"
    legend_constants[11] = "k2 in component F_PE2 (dimensionless)"
    legend_constants[12] = "L_r2 in component F_PE2 (dimensionless)"
    legend_algebraic[5] = "FL in component FL (dimensionless)"
    legend_constants[13] = "beta in component FL (dimensionless)"
    legend_constants[14] = "omega in component FL (dimensionless)"
    legend_constants[15] = "rho in component FL (dimensionless)"
    legend_algebraic[6] = "FV in component FV (dimensionless)"
    legend_constants[16] = "av0 in component FV (dimensionless)"
    legend_constants[17] = "av1 in component FV (dimensionless)"
    legend_constants[18] = "av2 in component FV (dimensionless)"
    legend_constants[19] = "cv0 in component FV (dimensionless)"
    legend_constants[20] = "cv1 in component FV (dimensionless)"
    legend_constants[21] = "bv in component FV (first_order_rate_constant)"
    legend_constants[22] = "V_max in component FV (first_order_rate_constant)"
    legend_algebraic[7] = "Af in component Af (dimensionless)"
    legend_constants[23] = "af in component Af (dimensionless)"
    legend_constants[24] = "nf0 in component Af (dimensionless)"
    legend_constants[25] = "nf1 in component Af (dimensionless)"
    legend_constants[26] = "nf in component Af (dimensionless)"
    legend_states[2] = "Y in component Y (dimensionless)"
    legend_states[3] = "S in component S (dimensionless)"
    legend_states[4] = "f_eff in component rise_and_fall_time (dimensionless)"
    legend_states[5] = "L_eff in component L_eff (dimensionless)"
    legend_algebraic[8] = "F0 in component F0 (dimensionless)"
    legend_algebraic[11] = "F_CE in component F_CE (newton)"
    legend_algebraic[12] = "F_total in component F_total (newton)"
    legend_constants[27] = "T_L in component L_eff (millisecond)"
    legend_constants[28] = "T_s in component S (millisecond)"
    legend_constants[29] = "as1 in component S (dimensionless)"
    legend_constants[30] = "as2 in component S (dimensionless)"
    legend_algebraic[0] = "as_ in component S (dimensionless)"
    legend_constants[31] = "c_Y in component Y (dimensionless)"
    legend_constants[32] = "V_Y in component Y (first_order_rate_constant)"
    legend_constants[33] = "T_Y in component Y (millisecond)"
    legend_states[6] = "f_int in component rise_and_fall_time (dimensionless)"
    legend_algebraic[9] = "df_eff_dt in component rise_and_fall_time (first_order_rate_constant)"
    legend_algebraic[10] = "T_f in component rise_and_fall_time (millisecond)"
    legend_constants[34] = "T_f1 in component rise_and_fall_time (millisecond)"
    legend_constants[35] = "T_f2 in component rise_and_fall_time (millisecond)"
    legend_constants[36] = "T_f3 in component rise_and_fall_time (millisecond)"
    legend_constants[37] = "T_f4 in component rise_and_fall_time (millisecond)"
    legend_constants[38] = "f_env in component user_defined_constants (dimensionless)"
    legend_constants[39] = "mass in component V (kilogram)"
    legend_algebraic[2] = "V0 in component V0 (first_order_rate_constant)"
    legend_algebraic[3] = "L0 in component L0 (dimensionless)"
    legend_rates[5] = "d/dt L_eff in component L_eff (dimensionless)"
    legend_rates[3] = "d/dt S in component S (dimensionless)"
    legend_rates[2] = "d/dt Y in component Y (dimensionless)"
    legend_rates[6] = "d/dt f_int in component rise_and_fall_time (dimensionless)"
    legend_rates[4] = "d/dt f_eff in component rise_and_fall_time (dimensionless)"
    legend_rates[1] = "d/dt V in component V (first_order_rate_constant)"
    legend_rates[0] = "d/dt L in component L (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 27.8
    constants[1] = 0.0047
    constants[2] = 0.964
    constants[3] = 0.02
    constants[4] = 23
    constants[5] = 23
    constants[6] = 0.046
    constants[7] = 1.17
    constants[8] = 0.001
    states[0] = 0.15
    constants[9] = 0.13
    states[1] = 0.09314
    constants[10] = 23
    constants[11] = 0.046
    constants[12] = 1.17
    constants[13] = 1.55
    constants[14] = 0.75
    constants[15] = 2.12
    constants[16] = -1.53
    constants[17] = 0
    constants[18] = 0
    constants[19] = -5.7
    constants[20] = 9.18
    constants[21] = 0.69
    constants[22] = -9.15
    constants[23] = 0.56
    constants[24] = 2.1
    constants[25] = 3.3
    constants[26] = 1
    states[2] = 1
    states[3] = 1
    states[4] = 0
    states[5] = 0.1497
    constants[27] = 0.088
    constants[28] = 43
    constants[29] = 1.76
    constants[30] = 0.96
    constants[31] = 0.35
    constants[32] = 0.1
    constants[33] = 200
    states[6] = 0
    constants[34] = 0.35
    constants[35] = 0.1
    constants[36] = 200
    constants[37] = 200
    constants[38] = 1
    constants[39] = 0.005
    constants[40] = constants[0]*constants[4]*constants[1]*log(exp((constants[3]-constants[2])/constants[1])+1.00000)
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[2] = (1.00000-(constants[31]*(1.00000-exp(-fabs(states[1])/constants[32]))+states[2]))/constants[33]
    rates[0] = states[1]
    algebraic[0] = custom_piecewise([less(states[4] , 0.100000), constants[29] , True, constants[30]])
    rates[3] = (algebraic[0]-states[3])/constants[28]
    algebraic[7] = 1.00000-exp(-(power((states[2]*states[3]*states[4])/(constants[23]*constants[26]), constants[26])))
    rates[5] = (power(states[0]-states[5], 3.00000))/(constants[27]*(1.00000-algebraic[7]))
    rootfind_0(voi, constants, rates, states, algebraic)
    rates[6] = (constants[38]-states[6])/algebraic[10]
    rates[4] = algebraic[9]
    algebraic[1] = constants[5]*constants[6]*log(exp((states[0]/constants[9]-constants[7])/constants[6])+1.00000)+constants[8]*states[1]
    algebraic[4] = constants[10]*(exp(constants[11]*(states[0]-constants[12]))-1.00000)
    algebraic[5] = exp(-(power(fabs((power(states[0], constants[13])-1.00000)/constants[14]), constants[15])))
    algebraic[6] = custom_piecewise([less_equal(states[1] , 0.00000), (constants[22]-states[1])/(constants[22]+(constants[19]+constants[20]*states[0])*states[1]) , True, (constants[21]-(constants[16]+constants[17]*states[0]+constants[18]*(power(states[0], 2.00000)))*states[1])/(constants[21]+states[1])])
    algebraic[8] = algebraic[7]*(algebraic[5]+algebraic[6]+algebraic[4])+algebraic[1]
    algebraic[11] = algebraic[8]*constants[4]
    algebraic[12] = constants[40]-algebraic[11]
    rates[1] = algebraic[12]/(1.00000*constants[39])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less(states[4] , 0.100000), constants[29] , True, constants[30]])
    algebraic[7] = 1.00000-exp(-(power((states[2]*states[3]*states[4])/(constants[23]*constants[26]), constants[26])))
    algebraic[1] = constants[5]*constants[6]*log(exp((states[0]/constants[9]-constants[7])/constants[6])+1.00000)+constants[8]*states[1]
    algebraic[4] = constants[10]*(exp(constants[11]*(states[0]-constants[12]))-1.00000)
    algebraic[5] = exp(-(power(fabs((power(states[0], constants[13])-1.00000)/constants[14]), constants[15])))
    algebraic[6] = custom_piecewise([less_equal(states[1] , 0.00000), (constants[22]-states[1])/(constants[22]+(constants[19]+constants[20]*states[0])*states[1]) , True, (constants[21]-(constants[16]+constants[17]*states[0]+constants[18]*(power(states[0], 2.00000)))*states[1])/(constants[21]+states[1])])
    algebraic[8] = algebraic[7]*(algebraic[5]+algebraic[6]+algebraic[4])+algebraic[1]
    algebraic[11] = algebraic[8]*constants[4]
    algebraic[12] = constants[40]-algebraic[11]
    algebraic[2] = states[1]/constants[9]
    algebraic[3] = states[0]/constants[9]
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
        algebraic[9] = soln[0]
        algebraic[10] = soln[1]
    else:
        for (i,t) in enumerate(voi):
            soln = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates[:i], states[:,i]), xtol=1E-6)
            initialGuess0 = soln
            algebraic[9][i] = soln[0]
            algebraic[10][i] = soln[1]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    resid = array([0.0] * 2)
    algebraic[9] = algebraicCandidate[0]
    algebraic[10] = algebraicCandidate[1]
    resid[0] = (algebraic[9]-(states[6]-states[4])/algebraic[10])
    resid[1] = (algebraic[10]-(custom_piecewise([greater_equal(algebraic[9] , 0.00000), constants[34]*(power(states[0], 2.00000))+constants[35]*constants[38] , True, (constants[36]+constants[37]*algebraic[7])/states[0]])))
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