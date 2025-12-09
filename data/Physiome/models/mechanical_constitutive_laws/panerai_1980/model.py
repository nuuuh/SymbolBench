# Size of variable arrays:
sizeAlgebraic = 1
sizeStates = 2
sizeConstants = 34
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "t in component environment (second)"
    legend_constants[0] = "x in component environment (nanometer)"
    legend_states[0] = "n in component Crossbridges_attached (dimensionless)"
    legend_states[1] = "A_c in component Actin_free (dimensionless)"
    legend_constants[23] = "f in component f (per_second)"
    legend_constants[24] = "g in component g (per_second)"
    legend_constants[1] = "h in component Crossbridges_attached (nanometer)"
    legend_constants[2] = "f_1 in component f (per_second)"
    legend_constants[3] = "g_1 in component g (per_second)"
    legend_constants[4] = "g_2 in component g (per_second)"
    legend_algebraic[0] = "Ca_f in component Ca_sarcoplasm (molar)"
    legend_constants[5] = "t_d in component Ca_sarcoplasm (second)"
    legend_constants[6] = "a_1 in component Ca_sarcoplasm (per_second_squared)"
    legend_constants[7] = "b_1 in component Ca_sarcoplasm (per_second_squared)"
    legend_constants[8] = "Ca_0 in component Ca_sarcoplasm (molar)"
    legend_constants[9] = "c_1 in component Actin_free (per_second)"
    legend_constants[26] = "c_2 in component Actin_free (per_second)"
    legend_constants[10] = "c_2_0 in component Actin_free (per_second)"
    legend_constants[11] = "k_i in component Actin_free (dimensionless)"
    legend_constants[25] = "s_h in component s_h (muscle_length)"
    legend_constants[12] = "q in component Actin_free (dimensionless)"
    legend_constants[13] = "AT_0 in component Actin_free (dimensionless)"
    legend_constants[32] = "F_SE in component Series_Elastic_Element (force)"
    legend_constants[14] = "alpha_s in component Series_Elastic_Element (force)"
    legend_constants[15] = "beta_s in component Series_Elastic_Element (muscle_length)"
    legend_constants[31] = "x_s in component SE_constants (muscle_length)"
    legend_constants[16] = "x_so in component Series_Elastic_Element (muscle_length)"
    legend_constants[30] = "X_M_0 in component X_0 (muscle_length)"
    legend_constants[17] = "L_max in component Series_Elastic_Element (muscle_length)"
    legend_constants[28] = "F_PE in component Parallel_Elastic_Element (force)"
    legend_constants[18] = "alpha_p in component Parallel_Elastic_Element (force)"
    legend_constants[19] = "beta_p in component Parallel_Elastic_Element (muscle_length)"
    legend_constants[27] = "x_p in component PE_constants (muscle_length)"
    legend_constants[20] = "x_po in component Parallel_Elastic_Element (muscle_length)"
    legend_constants[33] = "F_CE in component Contractile_Element (force)"
    legend_constants[21] = "F_PL in component s_h (force)"
    legend_constants[29] = "X_S_0 in component X_0 (muscle_length)"
    legend_constants[22] = "F_PL in component X_0 (force)"
    legend_rates[0] = "d/dt n in component Crossbridges_attached (dimensionless)"
    legend_rates[1] = "d/dt A_c in component Actin_free (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 10
    states[0] = 0
    states[1] = 1
    constants[1] = 12
    constants[2] = 70
    constants[3] = 40
    constants[4] = 240
    constants[5] = 0.3
    constants[6] = 200
    constants[7] = 5
    constants[8] = 0.45e-6
    constants[9] = 200e12
    constants[10] = 20
    constants[11] = 30.9
    constants[12] = 1.45
    constants[13] = 2
    constants[14] = 0.1027
    constants[15] = 20
    constants[16] = 0.0387
    constants[17] = 1
    constants[18] = 0.00224
    constants[19] = 20
    constants[20] = 0.221
    constants[21] = 3
    constants[22] = 3
    constants[23] = custom_piecewise([less(constants[0] , 0.00000), 0.00000 , greater_equal(constants[0] , 0.00000) & less(constants[0] , constants[1]), (constants[2]*constants[0])/constants[1] , True, 0.00000])
    constants[24] = custom_piecewise([less(constants[0] , 0.00000), constants[4] , greater_equal(constants[0] , 0.00000) & less(constants[0] , constants[1]), (constants[3]*constants[0])/constants[1] , True, (constants[3]*constants[0])/constants[1]])
    constants[25] = constants[20]-(1.00000*1.00000*log(1.00000+constants[21]/constants[18], 10))/constants[19]
    constants[26] = constants[10]*exp(constants[11]*(power(constants[25]/1.00000, constants[12])))
    constants[27] = constants[20]-constants[25]
    constants[28] = constants[18]*(exp((constants[19]*constants[27])/(1.00000*1.00000))-1.00000)
    constants[29] = (1.00000*1.00000*log(1.00000+constants[22]/constants[14], 10))/constants[15]
    constants[30] = ((constants[29]+constants[17])-constants[25])-constants[16]
    constants[31] = (constants[16]+constants[25]+constants[30])-constants[17]
    constants[32] = constants[14]*(exp((constants[15]*constants[31])/(1.00000*1.00000))-1.00000)
    constants[33] = constants[32]-constants[28]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[23]*(states[1]-states[0])-constants[24]*states[0]
    algebraic[0] = constants[8]*fabs(1.00000-exp(-constants[6]*(power(voi, 2.00000))))*exp(-constants[7]*(power(voi-constants[5], 2.00000)))
    rates[1] = constants[9]*(power(algebraic[0]/1.00000, 2.00000))*(constants[13]-states[1])-constants[26]*states[1]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[8]*fabs(1.00000-exp(-constants[6]*(power(voi, 2.00000))))*exp(-constants[7]*(power(voi-constants[5], 2.00000)))
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