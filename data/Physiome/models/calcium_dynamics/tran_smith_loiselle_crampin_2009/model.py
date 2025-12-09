# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 1
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
    legend_constants[0] = "k_p1 in component SERCA (second_order_rate_constant)"
    legend_constants[1] = "k_p2 in component SERCA (first_order_rate_constant)"
    legend_constants[2] = "k_p3 in component SERCA (first_order_rate_constant)"
    legend_constants[3] = "k_m1 in component SERCA (first_order_rate_constant)"
    legend_constants[4] = "k_m2 in component SERCA (second_order_rate_constant)"
    legend_constants[5] = "k_m3 in component SERCA (second_order_rate_constant)"
    legend_constants[6] = "kdcai in component SERCA (millimolar)"
    legend_constants[7] = "kdcasr in component SERCA (millimolar)"
    legend_constants[8] = "kdh1 in component SERCA (millimolar)"
    legend_constants[9] = "kdhi in component SERCA (millimolar_squared)"
    legend_constants[10] = "kdhsr in component SERCA (millimolar_squared)"
    legend_constants[11] = "kdh in component SERCA (millimolar)"
    legend_constants[12] = "n in component SERCA (dimensionless)"
    legend_constants[13] = "Ca_i in component SERCA (millimolar)"
    legend_states[0] = "Ca_sr in component SERCA (millimolar)"
    legend_constants[14] = "H_i in component SERCA (millimolar)"
    legend_constants[15] = "ATP in component SERCA (millimolar)"
    legend_constants[16] = "ADP in component SERCA (millimolar)"
    legend_constants[17] = "P_i in component SERCA (millimolar)"
    legend_constants[18] = "T_Cai in component SERCA (dimensionless)"
    legend_algebraic[0] = "T_Casr in component SERCA (dimensionless)"
    legend_constants[19] = "T_H1 in component SERCA (dimensionless)"
    legend_constants[20] = "T_Hi in component SERCA (dimensionless)"
    legend_constants[21] = "T_Hsr in component SERCA (dimensionless)"
    legend_constants[22] = "T_H in component SERCA (dimensionless)"
    legend_constants[23] = "a_p1 in component SERCA (first_order_rate_constant)"
    legend_constants[24] = "a_p2 in component SERCA (first_order_rate_constant)"
    legend_algebraic[1] = "a_p3 in component SERCA (first_order_rate_constant)"
    legend_constants[25] = "a_m1 in component SERCA (first_order_rate_constant)"
    legend_algebraic[2] = "a_m2 in component SERCA (first_order_rate_constant)"
    legend_constants[26] = "a_m3 in component SERCA (first_order_rate_constant)"
    legend_algebraic[3] = "s1 in component SERCA (per_second_squared)"
    legend_algebraic[4] = "s2 in component SERCA (per_second_squared)"
    legend_constants[27] = "s3 in component SERCA (per_second_squared)"
    legend_algebraic[5] = "v_cycle in component SERCA (first_order_rate_constant)"
    legend_rates[0] = "d/dt Ca_sr in component SERCA (millimolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 25900
    constants[1] = 2540
    constants[2] = 20.5
    constants[3] = 16
    constants[4] = 67200
    constants[5] = 149
    constants[6] = 0.9
    constants[7] = 2.24
    constants[8] = 1.09e-5
    constants[9] = 3.54e-3
    constants[10] = 1.05e-8
    constants[11] = 7.24e-5
    constants[12] = 2
    constants[13] = 150e-6
    states[0] = 0
    constants[14] = 1e-4
    constants[15] = 5
    constants[16] = 36.3e-3
    constants[17] = 1
    constants[18] = constants[13]/constants[6]
    constants[28] = 1.00000
    constants[19] = constants[14]/constants[8]
    constants[20] = (power(constants[14], constants[12]))/constants[9]
    constants[21] = (power(constants[14], constants[12]))/constants[10]
    constants[22] = constants[14]/constants[11]
    constants[23] = constants[0]*constants[15]
    constants[24] = (constants[1]*(power(constants[18], 2.00000)))/(power(constants[18], 2.00000)+(power(constants[18], 2.00000))*constants[20]+constants[20]*(1.00000+constants[19]))
    constants[25] = (constants[3]*constants[20])/(power(constants[18], 2.00000)+(power(constants[18], 2.00000))*constants[20]+constants[20]*(1.00000+constants[19]))
    constants[26] = constants[5]*constants[17]
    constants[27] = constants[23]*constants[24]+constants[26]*constants[25]+constants[26]*constants[24]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[28]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[0]/constants[7]
    algebraic[1] = (constants[2]*constants[21])/((power(algebraic[0], 2.00000))*constants[22]+constants[22]+constants[21]*(1.00000+constants[22]))
    algebraic[2] = (constants[4]*constants[16]*(power(algebraic[0], 2.00000))*constants[22])/((power(algebraic[0], 2.00000))*constants[22]+constants[22]+constants[21]*(1.00000+constants[22]))
    algebraic[3] = constants[24]*algebraic[1]+constants[25]*algebraic[1]+constants[25]*algebraic[2]
    algebraic[4] = constants[23]*algebraic[1]+algebraic[2]*constants[23]+algebraic[2]*constants[26]
    algebraic[5] = (constants[23]*constants[24]*algebraic[1]-constants[25]*algebraic[2]*constants[26])/(algebraic[3]+algebraic[4]+constants[27])
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