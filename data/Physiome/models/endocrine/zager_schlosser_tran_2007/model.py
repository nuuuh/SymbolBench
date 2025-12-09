# Size of variable arrays:
sizeAlgebraic = 12
sizeStates = 7
sizeConstants = 26
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (hour)"
    legend_states[0] = "Cgen_pp in component Cgen_pp (micromolar)"
    legend_constants[22] = "Qpp in component Cgen_pp (litre_hour)"
    legend_constants[0] = "Ppp in component Cgen_pp (dimensionless)"
    legend_algebraic[0] = "H_Cgen_pp in component Cgen_pp (dimensionless)"
    legend_constants[1] = "Vpp in component Cgen_pp (litre)"
    legend_algebraic[8] = "H_Cgen_B in component Cgen_B (dimensionless)"
    legend_constants[19] = "Qt in component model_parameters (litre_hour)"
    legend_algebraic[6] = "Cgen_B in component Cgen_B (micromolar)"
    legend_states[1] = "Cgen_rp in component Cgen_rp (micromolar)"
    legend_constants[2] = "Prp in component Cgen_rp (dimensionless)"
    legend_constants[23] = "Qrp in component Cgen_rp (litre_hour)"
    legend_algebraic[1] = "H_Cgen_rp in component Cgen_rp (dimensionless)"
    legend_constants[3] = "Vrp in component Cgen_rp (litre)"
    legend_constants[4] = "kgen_urine in component Cgen_rp (first_order_rate_constant)"
    legend_states[2] = "Cgen_l in component Cgen_l (micromolar)"
    legend_constants[24] = "Ql in component Cgen_l (litre_hour)"
    legend_algebraic[2] = "H_Cgen_l in component Cgen_l (dimensionless)"
    legend_constants[5] = "Vl in component Cgen_l (litre)"
    legend_constants[6] = "Pl in component Cgen_l (dimensionless)"
    legend_constants[7] = "Vmax in component model_parameters (micromolar_hour)"
    legend_constants[8] = "Km in component model_parameters (micromolar)"
    legend_algebraic[11] = "vinf in component model_parameters (micromole_hour)"
    legend_algebraic[3] = "H_Cgen_GI in component Cgen_GI (dimensionless)"
    legend_constants[25] = "QGI in component Cgen_GI (litre_hour)"
    legend_states[3] = "Cgen_GI in component Cgen_GI (micromolar)"
    legend_constants[9] = "VGI in component Cgen_GI (litre)"
    legend_states[4] = "Ccon_l in component Ccon_l (micromolar)"
    legend_algebraic[4] = "H_Ccon_l in component Ccon_l (dimensionless)"
    legend_constants[10] = "kbile in component model_parameters (first_order_rate_constant)"
    legend_algebraic[7] = "H_Ccon_ROB in component Ccon_ROB (dimensionless)"
    legend_states[5] = "Ccon_ROB in component Ccon_ROB (micromolar)"
    legend_constants[11] = "VROB in component Ccon_ROB (litre)"
    legend_constants[12] = "kcon_urine in component Ccon_ROB (first_order_rate_constant)"
    legend_states[6] = "Acon_b in component Acon_b (micromole)"
    legend_algebraic[5] = "d_t_Cgen_l in component Acon_b (hour)"
    legend_constants[13] = "tbasal in component Acon_b (hour)"
    legend_constants[14] = "ct in component Acon_b (hour)"
    legend_constants[15] = "cgen in component Acon_b (litre_hour_micromole)"
    legend_constants[16] = "cv in component model_parameters (micromole_hour)"
    legend_algebraic[9] = "F in component model_parameters (dimensionless)"
    legend_constants[20] = "G in component model_parameters (dimensionless)"
    legend_algebraic[10] = "H in component model_parameters (dimensionless)"
    legend_constants[21] = "I in component model_parameters (dimensionless)"
    legend_constants[17] = "epsilon in component model_parameters (hour)"
    legend_constants[18] = "m_rat in component model_parameters (dimensionless)"
    legend_rates[0] = "d/dt Cgen_pp in component Cgen_pp (micromolar)"
    legend_rates[1] = "d/dt Cgen_rp in component Cgen_rp (micromolar)"
    legend_rates[2] = "d/dt Cgen_l in component Cgen_l (micromolar)"
    legend_rates[3] = "d/dt Cgen_GI in component Cgen_GI (micromolar)"
    legend_rates[4] = "d/dt Ccon_l in component Ccon_l (micromolar)"
    legend_rates[5] = "d/dt Ccon_ROB in component Ccon_ROB (micromolar)"
    legend_rates[6] = "d/dt Acon_b in component Acon_b (micromole)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0
    constants[0] = 0.59
    constants[1] = 0.188
    states[1] = 0
    constants[2] = 1.94
    constants[3] = 0.012
    constants[4] = 9.5
    states[2] = 0
    constants[5] = 0.0092
    constants[6] = 3.61
    constants[7] = 25.25
    constants[8] = 0.6231
    states[3] = 0
    constants[9] = 0.011
    states[4] = 0
    constants[10] = 111.72
    states[5] = 0
    constants[11] = 0.2408
    constants[12] = 0.0
    states[6] = 0
    constants[13] = 0.0865
    constants[14] = 1.0107
    constants[15] = 0.06561
    constants[16] = 1.0
    constants[17] = 0.1
    constants[18] = 0.25
    constants[19] = 14.1000*(power(constants[18], 0.750000))
    constants[20] = 1.00000
    constants[21] = 0.00000
    constants[22] = 0.528000*constants[19]
    constants[23] = 0.289000*constants[19]
    constants[24] = 0.0300000*constants[19]
    constants[25] = 0.153000*constants[19]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[4] = custom_piecewise([less(states[4] , 0.00000), 0.00000 , True, 1.00000])
    rates[6] = constants[10]*algebraic[4]*states[4]*constants[5]
    algebraic[2] = custom_piecewise([less(states[2] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[7] = custom_piecewise([less(states[5] , 0.00000), 0.00000 , True, 1.00000])
    rates[4] = (constants[7]*algebraic[2]*states[2])/(constants[8]+algebraic[2]*states[2])-(constants[10]*algebraic[4]*states[4]+((constants[24]+constants[25])*(algebraic[4]*states[4]-algebraic[7]*states[5]))/constants[5])
    rates[5] = ((constants[24]+constants[25])*(algebraic[4]*states[4]-algebraic[7]*states[5]))/constants[11]-constants[12]*algebraic[7]*states[5]
    algebraic[0] = custom_piecewise([less(states[0] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[1] = custom_piecewise([less(states[1] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[6] = ((constants[23]*algebraic[1]*states[1])/constants[2]+(constants[22]*algebraic[0]*states[0])/constants[0]+((constants[24]+constants[25])*algebraic[2]*states[2])/constants[6])/constants[19]
    algebraic[8] = custom_piecewise([less(algebraic[6] , 0.00000), 0.00000 , True, 1.00000])
    rates[0] = (constants[22]*algebraic[8]*algebraic[6]-(constants[22]*algebraic[0]*states[0])/constants[0])/constants[1]
    rates[1] = (constants[23]*algebraic[8]*algebraic[6]-(constants[23]*algebraic[1]*states[1])/constants[2])/constants[3]-constants[4]*algebraic[1]*states[1]
    algebraic[3] = custom_piecewise([less(states[3] , 0.00000), 0.00000 , True, 1.00000])
    rates[3] = (constants[25]*algebraic[8]*algebraic[6]-(constants[25]*algebraic[3]*states[3])/constants[2])/constants[9]
    algebraic[9] = (-2.00000/(power(constants[17], 3.00000)))*(power(voi, 3.00000))+(3.00000/(power(constants[17], 2.00000)))*(power(voi, 2.00000))
    algebraic[10] = (3.00000*(power(1.00000-voi, 2.00000)))/(power(constants[17], 2.00000))-(2.00000*(power(1.00000-voi, 3.00000)))/(power(constants[17], 3.00000))
    algebraic[11] = custom_piecewise([less(voi , constants[17]) & greater_equal(voi , 0.00000), constants[16]*algebraic[9] , less_equal(voi , 1.00000-constants[17]) & greater_equal(voi , constants[17]), constants[16]*constants[20] , less_equal(voi , 1.00000) & greater(voi , 1.00000-constants[17]), constants[16]*algebraic[10] , less_equal(voi , 1.00000+5.00000/60.0000) & greater(voi , 1.00000), constants[16]*constants[21] , True, float('nan')])
    rates[2] = ((constants[24]*algebraic[8]*algebraic[6]-((constants[24]+constants[25])*algebraic[2]*states[2])/constants[6])+algebraic[11]+(constants[25]*algebraic[3]*states[3])/constants[2])/constants[5]-(constants[7]*algebraic[2]*states[2])/(constants[8]+algebraic[2]*states[2])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[4] = custom_piecewise([less(states[4] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[2] = custom_piecewise([less(states[2] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[7] = custom_piecewise([less(states[5] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[0] = custom_piecewise([less(states[0] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[1] = custom_piecewise([less(states[1] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[6] = ((constants[23]*algebraic[1]*states[1])/constants[2]+(constants[22]*algebraic[0]*states[0])/constants[0]+((constants[24]+constants[25])*algebraic[2]*states[2])/constants[6])/constants[19]
    algebraic[8] = custom_piecewise([less(algebraic[6] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[3] = custom_piecewise([less(states[3] , 0.00000), 0.00000 , True, 1.00000])
    algebraic[9] = (-2.00000/(power(constants[17], 3.00000)))*(power(voi, 3.00000))+(3.00000/(power(constants[17], 2.00000)))*(power(voi, 2.00000))
    algebraic[10] = (3.00000*(power(1.00000-voi, 2.00000)))/(power(constants[17], 2.00000))-(2.00000*(power(1.00000-voi, 3.00000)))/(power(constants[17], 3.00000))
    algebraic[11] = custom_piecewise([less(voi , constants[17]) & greater_equal(voi , 0.00000), constants[16]*algebraic[9] , less_equal(voi , 1.00000-constants[17]) & greater_equal(voi , constants[17]), constants[16]*constants[20] , less_equal(voi , 1.00000) & greater(voi , 1.00000-constants[17]), constants[16]*algebraic[10] , less_equal(voi , 1.00000+5.00000/60.0000) & greater(voi , 1.00000), constants[16]*constants[21] , True, float('nan')])
    algebraic[5] = constants[13]+constants[15]*algebraic[2]*states[2]*exp(-(voi/constants[14]))
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