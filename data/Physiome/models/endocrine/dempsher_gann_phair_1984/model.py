# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 13
sizeConstants = 30
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "cAMP in component cAMP (micromolar)"
    legend_constants[0] = "Ko in component cAMP (dimensionless)"
    legend_constants[1] = "Ka in component cAMP (per_micromolar)"
    legend_constants[2] = "Kb in component cAMP (per_micromolar)"
    legend_constants[3] = "Kdsm in component cAMP (micromolar)"
    legend_constants[4] = "Vmsm in component cAMP (flux)"
    legend_constants[5] = "Vdsm in component cAMP (flux)"
    legend_algebraic[0] = "ACTH in component cAMP (micromolar)"
    legend_algebraic[1] = "IS in component IS (dimensionless)"
    legend_constants[6] = "Crpt in component IS (dimensionless)"
    legend_constants[7] = "K in component IS (dimensionless)"
    legend_constants[8] = "Kd in component IS (micromolar)"
    legend_constants[9] = "n in component IS (dimensionless)"
    legend_states[1] = "V in component V (flux)"
    legend_constants[10] = "P in component V (micromolar_per_minute2)"
    legend_constants[11] = "Q in component V (first_order_rate_constant)"
    legend_states[2] = "W in component W (flux)"
    legend_constants[12] = "T in component W (micromolar2_per_minute2)"
    legend_constants[13] = "U in component W (first_order_rate_constant)"
    legend_states[3] = "CHOC in component CHOC (micromolar)"
    legend_constants[14] = "Lmtr in component model_parameters (first_order_rate_constant)"
    legend_states[4] = "Kmtr in component Kmtr (first_order_rate_constant)"
    legend_states[5] = "CHOM in component CHOM (micromolar)"
    legend_constants[15] = "Kbac in component model_parameters (first_order_rate_constant)"
    legend_states[6] = "Kfor in component Kfor (first_order_rate_constant)"
    legend_constants[16] = "Kcb in component model_parameters (first_order_rate_constant)"
    legend_constants[17] = "Kcf in component model_parameters (first_order_rate_constant)"
    legend_states[7] = "CHON in component CHON (micromolar)"
    legend_states[8] = "CHOL in component CHOL (micromolar)"
    legend_constants[18] = "C in component Kmtr (per_minute2)"
    legend_constants[19] = "D in component Kmtr (first_order_rate_constant)"
    legend_constants[20] = "R in component Kfor (per_minute2)"
    legend_constants[21] = "S in component Kfor (first_order_rate_constant)"
    legend_constants[22] = "Vm in component model_parameters (flux)"
    legend_constants[23] = "Km in component model_parameters (micromolar)"
    legend_states[9] = "PREG in component PREG (micromolar)"
    legend_constants[24] = "Vmptr in component model_parameters (flux)"
    legend_constants[25] = "Kmptr in component model_parameters (micromolar)"
    legend_states[10] = "PRO in component PRO (micromolar)"
    legend_constants[26] = "HA in component PRO (dimensionless)"
    legend_constants[27] = "AH in component model_parameters (first_order_rate_constant)"
    legend_states[11] = "HYPR in component HYPR (micromolar)"
    legend_constants[28] = "HY in component model_parameters (first_order_rate_constant)"
    legend_states[12] = "CORT in component CORT (micromolar)"
    legend_constants[29] = "LH in component CORT (first_order_rate_constant)"
    legend_rates[0] = "d/dt cAMP in component cAMP (micromolar)"
    legend_rates[1] = "d/dt V in component V (flux)"
    legend_rates[2] = "d/dt W in component W (flux)"
    legend_rates[3] = "d/dt CHOC in component CHOC (micromolar)"
    legend_rates[5] = "d/dt CHOM in component CHOM (micromolar)"
    legend_rates[8] = "d/dt CHOL in component CHOL (micromolar)"
    legend_rates[4] = "d/dt Kmtr in component Kmtr (first_order_rate_constant)"
    legend_rates[6] = "d/dt Kfor in component Kfor (first_order_rate_constant)"
    legend_rates[7] = "d/dt CHON in component CHON (micromolar)"
    legend_rates[9] = "d/dt PREG in component PREG (micromolar)"
    legend_rates[10] = "d/dt PRO in component PRO (micromolar)"
    legend_rates[11] = "d/dt HYPR in component HYPR (micromolar)"
    legend_rates[12] = "d/dt CORT in component CORT (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.95
    constants[0] = 0.013
    constants[1] = 10
    constants[2] = 1000000.0
    constants[3] = 10.0
    constants[4] = 6.0
    constants[5] = 10.0
    constants[6] = 3.0
    constants[7] = 80.0
    constants[8] = 2.11
    constants[9] = 4.0
    states[1] = 11.3
    constants[10] = 0.052
    constants[11] = 0.042
    states[2] = 10.0
    constants[12] = 8.0
    constants[13] = 0.0015
    states[3] = 532.0
    constants[14] = 1.65
    states[4] = 0.446
    states[5] = 139.0
    constants[15] = 10.0
    states[6] = 0.370
    constants[16] = 0.01
    constants[17] = 0.00033
    states[7] = 3.03
    states[8] = 3000.0
    constants[18] = 6.25
    constants[19] = 125.0
    constants[20] = 3.0
    constants[21] = 76.0
    constants[22] = 1890.0
    constants[23] = 270.0
    states[9] = 6.56
    constants[24] = 500.0
    constants[25] = 150.0
    states[10] = 0.64
    constants[26] = 0.5
    constants[27] = 16.4
    states[11] = 0.64
    constants[28] = 16.4
    states[12] = 5.2
    constants[29] = 0.724
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[2] = constants[12]*(power(states[3], -1.00000))-constants[13]*states[2]
    rates[3] = (states[1]+states[2]+constants[14]*states[5])-states[4]*states[3]
    rates[5] = (states[4]*states[3]+constants[15]*states[7]+constants[17]*states[8])-(constants[14]*states[5]+constants[16]*states[5]+states[6]*states[5])
    rates[8] = constants[16]*states[5]-constants[17]*states[8]
    rates[7] = states[6]*states[5]-(constants[15]*states[7]+(constants[22]*states[7])/(constants[23]+states[7]))
    rates[9] = (constants[22]*states[7])/(constants[23]+states[7])-(constants[24]*states[9])/(constants[25]+states[9])
    rates[10] = constants[26]*((constants[24]*states[9])/(constants[25]+states[9]))-constants[27]*states[10]
    rates[11] = constants[27]*states[10]-constants[28]*states[11]
    rates[12] = constants[28]*states[11]-constants[29]*states[12]
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 10.0000), 1.20000e-05 , greater_equal(voi , 10.0000) & less(voi , 20.0000), 1.60000e-05 , greater_equal(voi , 20.0000) & less(voi , 35.0000), 1.20000e-05 , greater_equal(voi , 35.0000) & less(voi , 45.0000), 1.60000e-05 , True, 1.20000e-05])
    rates[0] = (constants[4]*constants[0]*(1.00000+constants[2]*algebraic[0]))/((1.00000+constants[1]*algebraic[0])+constants[0]*(1.00000+constants[2]*algebraic[0]))-(constants[5]*states[0])/(constants[3]+states[0])
    algebraic[1] = (constants[7]*constants[6]*(power(states[0], constants[9])))/(power(constants[8], constants[9])+power(states[0], constants[9]))
    rates[1] = constants[10]*algebraic[1]-constants[11]*states[1]
    rates[4] = constants[18]*algebraic[1]-constants[19]*states[4]
    rates[6] = constants[20]*algebraic[1]-constants[21]*states[6]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 10.0000), 1.20000e-05 , greater_equal(voi , 10.0000) & less(voi , 20.0000), 1.60000e-05 , greater_equal(voi , 20.0000) & less(voi , 35.0000), 1.20000e-05 , greater_equal(voi , 35.0000) & less(voi , 45.0000), 1.60000e-05 , True, 1.20000e-05])
    algebraic[1] = (constants[7]*constants[6]*(power(states[0], constants[9])))/(power(constants[8], constants[9])+power(states[0], constants[9]))
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