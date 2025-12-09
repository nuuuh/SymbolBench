# Size of variable arrays:
sizeAlgebraic = 8
sizeStates = 7
sizeConstants = 35
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "Y in component Y (dimensionless)"
    legend_constants[0] = "k1 in component Y (first_order_rate_constant)"
    legend_algebraic[6] = "k2 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[1] = "k3 in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[0] = "C in component C (dimensionless)"
    legend_states[1] = "M in component M (dimensionless)"
    legend_algebraic[7] = "kwee in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[5] = "kcdc25 in component kinetic_parameters (first_order_rate_constant)"
    legend_states[2] = "preMPF in component preMPF (dimensionless)"
    legend_states[3] = "Cdc25P in component Cdc25P (dimensionless)"
    legend_constants[2] = "k25 in component Cdc25P (first_order_rate_constant)"
    legend_constants[3] = "Km25 in component Cdc25P (dimensionless)"
    legend_constants[4] = "k25r in component Cdc25P (first_order_rate_constant)"
    legend_constants[33] = "k25ro in component Cdc25P (first_order_rate_constant)"
    legend_constants[5] = "Km25r in component Cdc25P (dimensionless)"
    legend_constants[32] = "p in component p (dimensionless)"
    legend_algebraic[2] = "Cdc25 in component Cdc25 (dimensionless)"
    legend_states[4] = "Wee1 in component Wee1 (dimensionless)"
    legend_constants[6] = "kw in component Wee1 (first_order_rate_constant)"
    legend_constants[7] = "Kmw in component Wee1 (dimensionless)"
    legend_constants[8] = "kwr in component Wee1 (first_order_rate_constant)"
    legend_constants[34] = "kwro in component Wee1 (first_order_rate_constant)"
    legend_constants[9] = "Kmwr in component Wee1 (dimensionless)"
    legend_algebraic[3] = "Wee1P in component Wee1P (dimensionless)"
    legend_states[5] = "IEP in component IEP (dimensionless)"
    legend_constants[10] = "kie in component IEP (first_order_rate_constant)"
    legend_constants[11] = "Kmie in component IEP (dimensionless)"
    legend_constants[12] = "kier in component IEP (first_order_rate_constant)"
    legend_constants[13] = "Kmier in component IEP (dimensionless)"
    legend_algebraic[1] = "IE in component IE (dimensionless)"
    legend_states[6] = "APC_ in component APC_ (dimensionless)"
    legend_constants[14] = "kap in component APC_ (first_order_rate_constant)"
    legend_constants[15] = "Kmap in component APC_ (dimensionless)"
    legend_constants[16] = "kapr in component APC_ (first_order_rate_constant)"
    legend_constants[17] = "Kmapr in component APC_ (dimensionless)"
    legend_algebraic[4] = "APC in component APC (dimensionless)"
    legend_constants[18] = "Cdc2tot in component C (dimensionless)"
    legend_constants[19] = "Cdc25tot in component Cdc25 (dimensionless)"
    legend_constants[20] = "Wee1tot in component Wee1P (dimensionless)"
    legend_constants[21] = "IEtot in component IE (dimensionless)"
    legend_constants[22] = "APCtot in component APC (dimensionless)"
    legend_constants[23] = "R in component p (dimensionless)"
    legend_constants[24] = "s in component p (dimensionless)"
    legend_constants[25] = "q in component p (dimensionless)"
    legend_constants[26] = "V2 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[27] = "V2_ in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[28] = "V25 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[29] = "V25_ in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[30] = "Vwee in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[31] = "Vwee_ in component kinetic_parameters (first_order_rate_constant)"
    legend_rates[0] = "d/dt Y in component Y (dimensionless)"
    legend_rates[1] = "d/dt M in component M (dimensionless)"
    legend_rates[2] = "d/dt preMPF in component preMPF (dimensionless)"
    legend_rates[3] = "d/dt Cdc25P in component Cdc25P (dimensionless)"
    legend_rates[4] = "d/dt Wee1 in component Wee1 (dimensionless)"
    legend_rates[5] = "d/dt IEP in component IEP (dimensionless)"
    legend_rates[6] = "d/dt APC_ in component APC_ (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.01577552
    constants[0] = 0.024
    constants[1] = 1.15
    states[1] = 0.0142077
    states[2] = 0.02541723
    states[3] = 0.4844362
    constants[2] = 18
    constants[3] = 0.1
    constants[4] = 0.8
    constants[5] = 1
    states[4] = 0.5155638
    constants[6] = 18
    constants[7] = 0.1
    constants[8] = 0.8
    constants[9] = 1
    states[5] = 0.002287817
    constants[10] = 4.5
    constants[11] = 0.01
    constants[12] = 0.34
    constants[13] = 0.01
    states[6] = 0.5051103
    constants[14] = 0.3
    constants[15] = 0.01
    constants[16] = 0.3
    constants[17] = 1
    constants[18] = 1
    constants[19] = 1
    constants[20] = 1
    constants[21] = 1
    constants[22] = 1
    constants[23] = 1
    constants[24] = 0.021
    constants[25] = 1
    constants[26] = 0.01
    constants[27] = 0.6
    constants[28] = 0.04
    constants[29] = 0.4
    constants[30] = 0.025
    constants[31] = 2.5
    constants[32] = constants[24]*constants[23]+constants[25]
    constants[33] = constants[4]/constants[32]
    constants[34] = constants[8]/constants[32]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = constants[21]-states[5]
    rates[5] = (constants[10]*states[1]*algebraic[1])/(constants[11]+algebraic[1])-(constants[12]*states[5])/(constants[13]+states[5])
    algebraic[2] = constants[19]-states[3]
    rates[3] = (constants[2]*states[1]*algebraic[2])/(constants[3]+algebraic[2])-(constants[4]*states[3])/(constants[5]+states[3])
    algebraic[3] = constants[20]-states[4]
    rates[4] = (-constants[6]*states[1]*states[4])/(constants[7]+states[4])+(constants[8]*algebraic[3])/(constants[9]+algebraic[3])
    algebraic[4] = constants[22]-states[6]
    rates[6] = (constants[14]*states[5]*algebraic[4])/(constants[15]+algebraic[4])-(constants[16]*states[6])/(constants[17]+states[6])
    algebraic[6] = constants[26]*algebraic[4]+constants[27]*states[6]
    algebraic[0] = constants[18]-(states[1]+states[2])
    rates[0] = constants[0]-(algebraic[6]*states[0]+constants[1]*states[0]*algebraic[0])
    algebraic[7] = constants[30]*algebraic[3]+constants[31]*states[4]
    algebraic[5] = constants[28]*algebraic[2]+constants[29]*states[3]
    rates[1] = (algebraic[5]*states[2]+constants[1]*states[0]*algebraic[0])-(algebraic[6]*states[1]+algebraic[7]*states[1])
    rates[2] = algebraic[7]*states[1]-(algebraic[6]*states[2]+algebraic[5]*states[2])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = constants[21]-states[5]
    algebraic[2] = constants[19]-states[3]
    algebraic[3] = constants[20]-states[4]
    algebraic[4] = constants[22]-states[6]
    algebraic[6] = constants[26]*algebraic[4]+constants[27]*states[6]
    algebraic[0] = constants[18]-(states[1]+states[2])
    algebraic[7] = constants[30]*algebraic[3]+constants[31]*states[4]
    algebraic[5] = constants[28]*algebraic[2]+constants[29]*states[3]
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