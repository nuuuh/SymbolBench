# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 13
sizeConstants = 23
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "Cdk2_CycE in component Cdk2_CycE (dimensionless)"
    legend_constants[0] = "kwee in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[1] = "kon in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[2] = "k25A in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[3] = "koff in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[4] = "kassoc in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[5] = "kdissoc in component kinetic_parameters (first_order_rate_constant)"
    legend_algebraic[4] = "phi in component kinetic_parameters (dimensionless)"
    legend_states[1] = "Wee1_a in component Wee1_a (dimensionless)"
    legend_states[2] = "PCdk2_CycE in component PCdk2_CycE (dimensionless)"
    legend_states[3] = "Xic in component Xic (dimensionless)"
    legend_states[4] = "Cdk2_CycErem in component Cdk2_CycErem (dimensionless)"
    legend_states[5] = "Xic_Cdk2_CycE in component Xic_Cdk2_CycE (dimensionless)"
    legend_states[6] = "PCdk2_CycErem in component PCdk2_CycErem (dimensionless)"
    legend_states[7] = "Xic_PCdk2_CycE in component Xic_PCdk2_CycE (dimensionless)"
    legend_constants[6] = "Jwact in component Wee1_a (dimensionless)"
    legend_constants[7] = "Jwinact in component Wee1_a (dimensionless)"
    legend_constants[8] = "Wee1_total in component Wee1_a (dimensionless)"
    legend_constants[9] = "kwact in component Wee1_a (first_order_rate_constant)"
    legend_constants[10] = "kwinact in component Wee1_a (first_order_rate_constant)"
    legend_states[8] = "Kin_a in component Kin_a (dimensionless)"
    legend_constants[11] = "Jiact in component Kin_a (dimensionless)"
    legend_constants[12] = "Jiinact in component Kin_a (dimensionless)"
    legend_constants[13] = "kiact in component Kin_a (first_order_rate_constant)"
    legend_constants[14] = "kiinact in component Kin_a (first_order_rate_constant)"
    legend_constants[15] = "kedeg in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[16] = "kxdeg in component kinetic_parameters (first_order_rate_constant)"
    legend_states[9] = "Deg_a in component Deg_a (dimensionless)"
    legend_states[10] = "Xic_Cdk2_CycErem in component Xic_Cdk2_CycErem (dimensionless)"
    legend_states[11] = "Xic_PCdk2_CycErem in component Xic_PCdk2_CycErem (dimensionless)"
    legend_algebraic[5] = "Heav in component Deg_a (dimensionless)"
    legend_constants[17] = "kdact in component Deg_a (first_order_rate_constant)"
    legend_constants[18] = "theta in component Deg_a (dimensionless)"
    legend_algebraic[1] = "x in component Deg_a (dimensionless)"
    legend_constants[19] = "Deg in component kinetic_parameters (dimensionless)"
    legend_states[12] = "Xicrem in component Xicrem (dimensionless)"
    legend_algebraic[2] = "Cyc_total in component Cyc_total (dimensionless)"
    legend_algebraic[3] = "Xic_total in component Xic_total (dimensionless)"
    legend_constants[20] = "epsilon in component kinetic_parameters (dimensionless)"
    legend_algebraic[0] = "pool in component kinetic_parameters (dimensionless)"
    legend_constants[21] = "n in component kinetic_parameters (dimensionless)"
    legend_constants[22] = "L in component kinetic_parameters (dimensionless)"
    legend_rates[0] = "d/dt Cdk2_CycE in component Cdk2_CycE (dimensionless)"
    legend_rates[2] = "d/dt PCdk2_CycE in component PCdk2_CycE (dimensionless)"
    legend_rates[1] = "d/dt Wee1_a in component Wee1_a (dimensionless)"
    legend_rates[8] = "d/dt Kin_a in component Kin_a (dimensionless)"
    legend_rates[4] = "d/dt Cdk2_CycErem in component Cdk2_CycErem (dimensionless)"
    legend_rates[6] = "d/dt PCdk2_CycErem in component PCdk2_CycErem (dimensionless)"
    legend_rates[9] = "d/dt Deg_a in component Deg_a (dimensionless)"
    legend_rates[3] = "d/dt Xic in component Xic (dimensionless)"
    legend_rates[5] = "d/dt Xic_Cdk2_CycE in component Xic_Cdk2_CycE (dimensionless)"
    legend_rates[7] = "d/dt Xic_PCdk2_CycE in component Xic_PCdk2_CycE (dimensionless)"
    legend_rates[10] = "d/dt Xic_Cdk2_CycErem in component Xic_Cdk2_CycErem (dimensionless)"
    legend_rates[11] = "d/dt Xic_PCdk2_CycErem in component Xic_PCdk2_CycErem (dimensionless)"
    legend_rates[12] = "d/dt Xicrem in component Xicrem (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.06
    constants[0] = 1.5
    constants[1] = 0.02
    constants[2] = 0.1
    constants[3] = 0.0001
    constants[4] = 0.1
    constants[5] = 0.001
    states[1] = 1.02
    states[2] = 0.94
    states[3] = 3
    states[4] = 0
    states[5] = 0
    states[6] = 0
    states[7] = 0
    constants[6] = 0.01
    constants[7] = 0.01
    constants[8] = 8
    constants[9] = 0.75
    constants[10] = 1.5
    states[8] = 0.6
    constants[11] = 0.01
    constants[12] = 0.01
    constants[13] = 0.15
    constants[14] = 0.6
    constants[15] = 0.017
    constants[16] = 0.01
    states[9] = 0
    states[10] = 0
    states[11] = 0
    constants[17] = 0.023
    constants[18] = 0.3
    constants[19] = 0
    states[12] = 0
    constants[20] = 0.001
    constants[21] = 4
    constants[22] = 0.4
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[1] = (constants[9]*(constants[8]-states[1]))/((constants[6]+constants[8])-states[1])-(constants[10]*states[8]*states[1])/(constants[7]+states[1])
    rates[8] = (constants[13]*(1.00000-states[8]))/((constants[11]+1.00000)-states[8])-(constants[14]*states[0]*states[8])/(constants[12]+states[8])
    rates[3] = constants[5]*(states[5]+states[7]+states[10]+states[11])-constants[4]*states[3]*(states[0]+states[2]+states[4]+states[6])
    rates[12] = (constants[15]*states[9]*states[10]+constants[15]*states[9]*states[11])-constants[16]*states[12]
    algebraic[0] = states[4]+states[6]+states[10]+states[11]
    algebraic[4] = (constants[20]+power(algebraic[0], constants[21]))/(power(constants[22], constants[21])+power(algebraic[0], constants[21]))
    rates[0] = (constants[2]*states[2]+constants[3]*states[4]+constants[5]*states[5])-(constants[0]*states[1]*states[0]+constants[1]*algebraic[4]*states[0]+constants[4]*states[3]*states[0])
    rates[2] = (constants[0]*states[1]*states[0]+constants[3]*states[6]+constants[5]*states[7])-(constants[2]*states[2]+constants[1]*algebraic[4]*states[2]+constants[4]*states[3]*states[2])
    rates[4] = (constants[1]*algebraic[4]*states[0]+constants[16]*states[10]+constants[5]*states[10])-(constants[3]*states[4]+constants[15]*states[9]*states[4]+constants[4]*states[3]*states[4])
    rates[6] = (constants[1]*algebraic[4]*states[2]+constants[16]*states[11]+constants[5]*states[11])-(constants[3]*states[6]+constants[15]*states[9]*states[6]+constants[4]*states[3]*states[6])
    algebraic[1] = states[4]-constants[18]
    algebraic[5] = custom_piecewise([less(algebraic[1] , 0.00000), 0.00000 , True, 1.00000])
    rates[9] = constants[17]*algebraic[5]
    rates[5] = (constants[4]*states[3]*states[0]+constants[3]*states[10]+constants[2]*states[7])-(constants[5]*states[5]+constants[1]*algebraic[4]*states[5]+constants[0]*states[1]*states[5])
    rates[7] = (constants[4]*states[3]*states[2]+constants[3]*states[11]+constants[0]*states[1]*states[5])-(constants[5]*states[7]+constants[1]*algebraic[4]*states[7]+constants[2]*states[7])
    rates[10] = (constants[4]*states[3]*states[4]+constants[1]*algebraic[4]*states[5])-(constants[5]*states[10]+constants[3]*states[10]+constants[15]*states[9]*states[10]+constants[16]*states[10])
    rates[11] = (constants[4]*states[3]*states[6]+constants[1]*algebraic[4]*states[7])-(constants[5]*states[11]+constants[3]*states[11]+constants[15]*constants[19]*states[11]+constants[16]*states[11])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[4]+states[6]+states[10]+states[11]
    algebraic[4] = (constants[20]+power(algebraic[0], constants[21]))/(power(constants[22], constants[21])+power(algebraic[0], constants[21]))
    algebraic[1] = states[4]-constants[18]
    algebraic[5] = custom_piecewise([less(algebraic[1] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[2] = states[7]+states[5]+states[11]+states[10]+states[6]+states[4]+states[0]+states[2]
    algebraic[3] = (states[3]+states[7]+states[5]+states[11]+states[10]+states[12])/3.00000
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