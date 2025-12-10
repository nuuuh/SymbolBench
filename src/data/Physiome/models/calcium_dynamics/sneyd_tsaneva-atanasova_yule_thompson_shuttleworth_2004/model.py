# Size of variable arrays:
sizeAlgebraic = 10
sizeStates = 8
sizeConstants = 29
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "c in component c (micromolar)"
    legend_constants[0] = "delta in component c (dimensionless)"
    legend_algebraic[2] = "J_serca in component J_serca (flux)"
    legend_constants[28] = "J_in in component J_in (flux)"
    legend_algebraic[4] = "J_pm in component J_pm (flux)"
    legend_algebraic[0] = "J_IPR in component J_IPR (flux)"
    legend_states[1] = "ce in component ce (micromolar)"
    legend_constants[1] = "gamma in component ce (dimensionless)"
    legend_constants[2] = "p in component p (micromolar)"
    legend_constants[3] = "kf in component J_IPR (first_order_rate_constant)"
    legend_constants[4] = "g1 in component J_IPR (first_order_rate_constant)"
    legend_states[2] = "A in component A (dimensionless)"
    legend_states[3] = "O in component O (dimensionless)"
    legend_states[4] = "R in component R (dimensionless)"
    legend_algebraic[1] = "phi_1 in component IPR_parameters (first_order_rate_constant)"
    legend_algebraic[3] = "phi_2 in component IPR_parameters (second_order_rate_constant)"
    legend_algebraic[5] = "phi_2_ in component IPR_parameters (first_order_rate_constant)"
    legend_constants[5] = "k_1_ in component IPR_parameters (first_order_rate_constant)"
    legend_constants[6] = "l_2_ in component IPR_parameters (first_order_rate_constant)"
    legend_states[5] = "I_1 in component I_1 (dimensionless)"
    legend_states[6] = "S in component S (dimensionless)"
    legend_algebraic[6] = "phi_3 in component IPR_parameters (second_order_rate_constant)"
    legend_algebraic[7] = "phi_4 in component IPR_parameters (first_order_rate_constant)"
    legend_algebraic[8] = "phi_4_ in component IPR_parameters (first_order_rate_constant)"
    legend_constants[7] = "k_3_ in component IPR_parameters (first_order_rate_constant)"
    legend_states[7] = "I_2 in component I_2 (dimensionless)"
    legend_algebraic[9] = "phi_5 in component IPR_parameters (first_order_rate_constant)"
    legend_constants[8] = "k_1 in component IPR_parameters (second_order_rate_constant)"
    legend_constants[9] = "k_2 in component IPR_parameters (second_order_rate_constant)"
    legend_constants[10] = "k_2_ in component IPR_parameters (first_order_rate_constant)"
    legend_constants[11] = "k_3 in component IPR_parameters (second_order_rate_constant)"
    legend_constants[12] = "k_4 in component IPR_parameters (second_order_rate_constant)"
    legend_constants[13] = "k_4_ in component IPR_parameters (first_order_rate_constant)"
    legend_constants[14] = "l_2 in component IPR_parameters (first_order_rate_constant)"
    legend_constants[15] = "l_4 in component IPR_parameters (second_order_rate_constant)"
    legend_constants[16] = "l_4_ in component IPR_parameters (second_order_rate_constant)"
    legend_constants[17] = "l_6 in component IPR_parameters (first_order_rate_constant)"
    legend_constants[18] = "l_6_ in component IPR_parameters (first_order_rate_constant)"
    legend_constants[19] = "L_1 in component IPR_parameters (micromolar)"
    legend_constants[20] = "L_3 in component IPR_parameters (micromolar)"
    legend_constants[21] = "L_5 in component IPR_parameters (micromolar)"
    legend_constants[22] = "Vs in component J_serca (micromolar2_per_second)"
    legend_constants[23] = "Ks in component J_serca (micromolar)"
    legend_constants[24] = "Vp in component J_pm (flux)"
    legend_constants[25] = "Kp in component J_pm (micromolar)"
    legend_constants[26] = "alpha1 in component J_in (flux)"
    legend_constants[27] = "alpha2 in component J_in (first_order_rate_constant)"
    legend_rates[0] = "d/dt c in component c (micromolar)"
    legend_rates[1] = "d/dt ce in component ce (micromolar)"
    legend_rates[4] = "d/dt R in component R (dimensionless)"
    legend_rates[3] = "d/dt O in component O (dimensionless)"
    legend_rates[5] = "d/dt I_1 in component I_1 (dimensionless)"
    legend_rates[7] = "d/dt I_2 in component I_2 (dimensionless)"
    legend_rates[6] = "d/dt S in component S (dimensionless)"
    legend_rates[2] = "d/dt A in component A (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.01
    constants[0] = 0.1
    states[1] = 0.01
    constants[1] = 5.4
    constants[2] = 10.0
    constants[3] = 0.96
    constants[4] = 0.002
    states[2] = 0.16
    states[3] = 0.01
    states[4] = 0.16
    constants[5] = 0.04
    constants[6] = 0.8
    states[5] = 0.16
    states[6] = 0.16
    constants[7] = 29.8
    states[7] = 0.16
    constants[8] = 0.64
    constants[9] = 37.4
    constants[10] = 1.4
    constants[11] = 0.11
    constants[12] = 4.0
    constants[13] = 0.54
    constants[14] = 1.7
    constants[15] = 1.7
    constants[16] = 2.5
    constants[17] = 4707.0
    constants[18] = 11.4
    constants[19] = 0.12
    constants[20] = 0.025
    constants[21] = 54.7
    constants[22] = 120.0
    constants[23] = 0.18
    constants[24] = 28.0
    constants[25] = 0.42
    constants[26] = 0.03
    constants[27] = 0.2
    constants[28] = constants[26]+constants[27]*constants[2]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = ((constants[8]*constants[19]+constants[14])*states[0])/(constants[19]+states[0]*(1.00000+constants[19]/constants[20]))
    rates[5] = algebraic[1]*states[4]-(constants[5]+constants[6])*states[5]
    algebraic[2] = ((constants[22]*states[0])/(constants[23]+states[0]))*(1.00000/states[1])
    algebraic[0] = (constants[3]*(power(0.100000*states[3]+0.900000*states[2], 4.00000))+constants[4])*(states[1]-states[0])
    rates[1] = constants[1]*(algebraic[2]-algebraic[0])
    algebraic[4] = (constants[24]*(power(states[0], 2.00000)))/(power(constants[25], 2.00000)+power(states[0], 2.00000))
    rates[0] = (algebraic[0]-algebraic[2])+constants[0]*(constants[28]-algebraic[4])
    algebraic[3] = (constants[9]*constants[20]+constants[15]*states[0])/(constants[20]+states[0]*(1.00000+constants[20]/constants[19]))
    algebraic[5] = (constants[10]+constants[16]*states[0])/(1.00000+states[0]/constants[21])
    rates[4] = (algebraic[5]*states[3]-(algebraic[3]*constants[2]*states[4]+algebraic[1]*states[4]))+(constants[6]+constants[5])*states[5]
    algebraic[6] = (constants[11]*constants[21])/(states[0]+constants[21])
    rates[6] = 1.00000*algebraic[6]*states[3]-constants[7]*states[6]
    algebraic[7] = ((constants[12]*constants[21]+constants[17])*states[0])/(states[0]+constants[21])
    algebraic[8] = (constants[19]*(constants[13]+constants[18]))/(states[0]+constants[19])
    rates[3] = (algebraic[3]*constants[2]*states[4]-(algebraic[5]+algebraic[7]+1.00000*algebraic[6])*states[3])+algebraic[8]*states[2]+constants[7]*states[6]
    algebraic[9] = ((constants[8]*constants[19]+constants[14])*states[0])/(states[0]+constants[19])
    rates[7] = algebraic[9]*states[2]-(constants[5]+constants[6])*states[7]
    rates[2] = (algebraic[7]*states[3]-(algebraic[8]*states[2]+algebraic[9]*states[2]))+(constants[5]+constants[6])*states[7]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = ((constants[8]*constants[19]+constants[14])*states[0])/(constants[19]+states[0]*(1.00000+constants[19]/constants[20]))
    algebraic[2] = ((constants[22]*states[0])/(constants[23]+states[0]))*(1.00000/states[1])
    algebraic[0] = (constants[3]*(power(0.100000*states[3]+0.900000*states[2], 4.00000))+constants[4])*(states[1]-states[0])
    algebraic[4] = (constants[24]*(power(states[0], 2.00000)))/(power(constants[25], 2.00000)+power(states[0], 2.00000))
    algebraic[3] = (constants[9]*constants[20]+constants[15]*states[0])/(constants[20]+states[0]*(1.00000+constants[20]/constants[19]))
    algebraic[5] = (constants[10]+constants[16]*states[0])/(1.00000+states[0]/constants[21])
    algebraic[6] = (constants[11]*constants[21])/(states[0]+constants[21])
    algebraic[7] = ((constants[12]*constants[21]+constants[17])*states[0])/(states[0]+constants[21])
    algebraic[8] = (constants[19]*(constants[13]+constants[18]))/(states[0]+constants[19])
    algebraic[9] = ((constants[8]*constants[19]+constants[14])*states[0])/(states[0]+constants[19])
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