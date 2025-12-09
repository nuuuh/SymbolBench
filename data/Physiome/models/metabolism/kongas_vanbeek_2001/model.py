# Size of variable arrays:
sizeAlgebraic = 11
sizeStates = 10
sizeConstants = 31
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "ATP in component ATP (micromolar)"
    legend_constants[0] = "V_hyd in component V_hyd (flux)"
    legend_algebraic[1] = "V_MMCK in component V_MMCK (flux)"
    legend_algebraic[8] = "J_ATP in component J_ATP (flux)"
    legend_constants[1] = "V_cyt in component fractional_volumes (dimensionless)"
    legend_states[1] = "ADP in component ADP (micromolar)"
    legend_algebraic[9] = "J_ADP in component J_ADP (flux)"
    legend_states[2] = "PCr in component PCr (micromolar)"
    legend_algebraic[5] = "J_PCr in component J_PCr (flux)"
    legend_states[3] = "Cr in component Cr (micromolar)"
    legend_algebraic[6] = "J_Cr in component J_Cr (flux)"
    legend_states[4] = "Pi in component Pi (micromolar)"
    legend_algebraic[10] = "J_Pi in component J_Pi (flux)"
    legend_states[5] = "ATP_i in component ATP_i (micromolar)"
    legend_algebraic[3] = "V_MiCK in component V_MiCK (flux)"
    legend_algebraic[7] = "V_syn in component V_syn (flux)"
    legend_constants[2] = "V_ims in component fractional_volumes (dimensionless)"
    legend_states[6] = "ADP_i in component ADP_i (micromolar)"
    legend_states[7] = "PCr_i in component PCr_i (micromolar)"
    legend_states[8] = "Cr_i in component Cr_i (micromolar)"
    legend_states[9] = "Pi_i in component Pi_i (micromolar)"
    legend_algebraic[0] = "Den_MMCK in component V_MMCK (dimensionless)"
    legend_constants[3] = "Kia in component V_MMCK (micromolar)"
    legend_constants[4] = "Kb in component V_MMCK (micromolar)"
    legend_constants[5] = "Kib in component V_MMCK (micromolar)"
    legend_constants[27] = "KIb in component V_MMCK (micromolar)"
    legend_constants[28] = "Kc in component V_MMCK (micromolar)"
    legend_constants[6] = "Kic in component V_MMCK (micromolar)"
    legend_constants[7] = "Kd in component V_MMCK (micromolar)"
    legend_constants[8] = "Kid in component V_MMCK (micromolar)"
    legend_constants[9] = "Vf in component V_MMCK (flux)"
    legend_constants[10] = "Vb in component V_MMCK (flux)"
    legend_algebraic[2] = "Den_MiCK in component V_MiCK (dimensionless)"
    legend_constants[11] = "Kia in component V_MiCK (micromolar)"
    legend_constants[12] = "Kb in component V_MiCK (micromolar)"
    legend_constants[13] = "Kib in component V_MiCK (micromolar)"
    legend_constants[29] = "KIb in component V_MiCK (micromolar)"
    legend_constants[30] = "Kc in component V_MiCK (micromolar)"
    legend_constants[14] = "Kic in component V_MiCK (micromolar)"
    legend_constants[15] = "Kd in component V_MiCK (micromolar)"
    legend_constants[16] = "Kid in component V_MiCK (micromolar)"
    legend_constants[17] = "Vf in component V_MiCK (flux)"
    legend_constants[18] = "Vb in component V_MiCK (flux)"
    legend_algebraic[4] = "Den_syn in component V_syn (dimensionless)"
    legend_constants[19] = "KPi in component V_syn (micromolar)"
    legend_constants[20] = "KADP in component V_syn (micromolar)"
    legend_constants[21] = "V_syn_max in component V_syn (flux)"
    legend_constants[22] = "R_ATP in component J_ATP (first_order_rate_constant)"
    legend_constants[23] = "R_ADP in component J_ADP (first_order_rate_constant)"
    legend_constants[24] = "R_PCr in component J_PCr (first_order_rate_constant)"
    legend_constants[25] = "R_Cr in component J_Cr (first_order_rate_constant)"
    legend_constants[26] = "R_Pi in component J_Pi (first_order_rate_constant)"
    legend_rates[0] = "d/dt ATP in component ATP (micromolar)"
    legend_rates[1] = "d/dt ADP in component ADP (micromolar)"
    legend_rates[2] = "d/dt PCr in component PCr (micromolar)"
    legend_rates[3] = "d/dt Cr in component Cr (micromolar)"
    legend_rates[4] = "d/dt Pi in component Pi (micromolar)"
    legend_rates[5] = "d/dt ATP_i in component ATP_i (micromolar)"
    legend_rates[6] = "d/dt ADP_i in component ADP_i (micromolar)"
    legend_rates[7] = "d/dt PCr_i in component PCr_i (micromolar)"
    legend_rates[8] = "d/dt Cr_i in component Cr_i (micromolar)"
    legend_rates[9] = "d/dt Pi_i in component Pi_i (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 9644.425
    constants[0] = 4.6E3
    constants[1] = 0.75
    states[1] = 60.0
    states[2] = 12500.0
    states[3] = 13500.0
    states[4] = 8000.0
    states[5] = 9644.425
    constants[2] = 0.0625
    states[6] = 2.5
    states[7] = 12500.0
    states[8] = 13500.0
    states[9] = 8000.0
    constants[3] = 9.0E2
    constants[4] = 1.55E4
    constants[5] = 3.49E4
    constants[6] = 2.224E2
    constants[7] = 1.67E3
    constants[8] = 4.73E3
    constants[9] = 6.966E3
    constants[10] = 2.925E4
    constants[11] = 7.5E2
    constants[12] = 5.2E3
    constants[13] = 2.88E4
    constants[14] = 2.048E2
    constants[15] = 5.0E2
    constants[16] = 1.6E3
    constants[17] = 2.658E3
    constants[18] = 1.116E4
    constants[19] = 20.0
    constants[20] = 8.0E2
    constants[21] = 4.6E3
    constants[22] = 8.16
    constants[23] = 8.16
    constants[24] = 14.6
    constants[25] = 14.6
    constants[26] = 18.4
    constants[27] = constants[5]
    constants[28] = (constants[6]*constants[7])/constants[8]
    constants[29] = constants[13]
    constants[30] = (constants[14]*constants[15])/constants[16]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = 1.00000+states[3]/constants[5]+states[2]/constants[8]+states[0]*(1.00000/constants[3]+states[3]/(constants[3]*constants[4]))+states[1]*(1.00000/constants[6]+states[2]/(constants[8]*constants[28])+states[3]/(constants[6]*constants[27]))
    algebraic[1] = (constants[9]*((states[0]*states[3])/(constants[3]*constants[4]))-constants[10]*((states[1]*states[2])/(constants[6]*constants[7])))/algebraic[0]
    algebraic[5] = constants[24]*(states[7]-states[2])
    rates[2] = (algebraic[5]+algebraic[1])/constants[1]
    algebraic[6] = constants[25]*(states[8]-states[3])
    rates[3] = (algebraic[6]-algebraic[1])/constants[1]
    algebraic[2] = 1.00000+states[8]/constants[13]+states[7]/constants[16]+states[5]*(1.00000/constants[11]+states[8]/(constants[11]*constants[12]))+states[6]*(1.00000/constants[14]+states[7]/(constants[16]*constants[30])+states[8]/(constants[14]*constants[29]))
    algebraic[3] = (constants[17]*((states[5]*states[8])/(constants[11]*constants[12]))-constants[18]*((states[6]*states[7])/(constants[14]*constants[15])))/algebraic[2]
    rates[7] = (algebraic[3]-algebraic[5])/constants[2]
    rates[8] = -(algebraic[3]+algebraic[6])/constants[2]
    algebraic[8] = constants[22]*(states[5]-states[0])
    rates[0] = (algebraic[8]-(constants[0]+algebraic[1]))/constants[1]
    algebraic[9] = constants[23]*(states[6]-states[1])
    rates[1] = (algebraic[9]+constants[0]+algebraic[1])/constants[1]
    algebraic[10] = constants[26]*(states[9]-states[4])
    rates[4] = (algebraic[10]+constants[0])/constants[1]
    algebraic[4] = 1.00000+states[6]/constants[20]+states[9]/constants[19]+(states[6]*states[9])/(constants[20]*constants[19])
    algebraic[7] = constants[21]*((states[6]*states[9])/(constants[19]*constants[20]*algebraic[4]))
    rates[5] = -(algebraic[8]+algebraic[7]+algebraic[3])/constants[2]
    rates[6] = ((algebraic[7]+algebraic[3])-algebraic[9])/constants[2]
    rates[9] = (algebraic[7]-algebraic[10])/constants[2]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = 1.00000+states[3]/constants[5]+states[2]/constants[8]+states[0]*(1.00000/constants[3]+states[3]/(constants[3]*constants[4]))+states[1]*(1.00000/constants[6]+states[2]/(constants[8]*constants[28])+states[3]/(constants[6]*constants[27]))
    algebraic[1] = (constants[9]*((states[0]*states[3])/(constants[3]*constants[4]))-constants[10]*((states[1]*states[2])/(constants[6]*constants[7])))/algebraic[0]
    algebraic[5] = constants[24]*(states[7]-states[2])
    algebraic[6] = constants[25]*(states[8]-states[3])
    algebraic[2] = 1.00000+states[8]/constants[13]+states[7]/constants[16]+states[5]*(1.00000/constants[11]+states[8]/(constants[11]*constants[12]))+states[6]*(1.00000/constants[14]+states[7]/(constants[16]*constants[30])+states[8]/(constants[14]*constants[29]))
    algebraic[3] = (constants[17]*((states[5]*states[8])/(constants[11]*constants[12]))-constants[18]*((states[6]*states[7])/(constants[14]*constants[15])))/algebraic[2]
    algebraic[8] = constants[22]*(states[5]-states[0])
    algebraic[9] = constants[23]*(states[6]-states[1])
    algebraic[10] = constants[26]*(states[9]-states[4])
    algebraic[4] = 1.00000+states[6]/constants[20]+states[9]/constants[19]+(states[6]*states[9])/(constants[20]*constants[19])
    algebraic[7] = constants[21]*((states[6]*states[9])/(constants[19]*constants[20]*algebraic[4]))
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