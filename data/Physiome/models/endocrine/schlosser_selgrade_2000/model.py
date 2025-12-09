# Size of variable arrays:
sizeAlgebraic = 12
sizeStates = 4
sizeConstants = 20
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (day)"
    legend_algebraic[6] = "rel_LH_E2_P4_RP_LH in component RP_LH (microg_day)"
    legend_states[0] = "RP_LH in component RP_LH (microg)"
    legend_algebraic[10] = "syn_LH_E2_P4 in component RP_LH (microg_day)"
    legend_constants[0] = "V0_LH in component RP_LH (microg_day)"
    legend_constants[1] = "V1_LH in component RP_LH (microg_day)"
    legend_constants[2] = "h in component RP_LH (dimensionless)"
    legend_constants[3] = "Km_LH in component RP_LH (ng_L)"
    legend_constants[4] = "Ki_LHP in component RP_LH (nmol_L)"
    legend_constants[5] = "kLH_rel in component RP_LH (first_order_rate_constant)"
    legend_constants[6] = "CLH_P in component RP_LH (L_nmol)"
    legend_constants[7] = "CLH_E in component RP_LH (L_ng)"
    legend_algebraic[3] = "E2 in component E2 (ng_L)"
    legend_algebraic[4] = "E2_dE in component E2_dE (ng_L)"
    legend_algebraic[5] = "P4 in component P4 (nmol_L)"
    legend_algebraic[8] = "P4_dP in component P4_dP (nmol_L)"
    legend_states[1] = "LH in component LH (microg_l)"
    legend_constants[8] = "v_dis in component LH (litre)"
    legend_algebraic[0] = "clear_LH in component LH (microg_l_day)"
    legend_constants[9] = "kLH_cl in component LH (first_order_rate_constant)"
    legend_algebraic[7] = "rel_FSH_E2_P4_RP_FSH in component RP_FSH (microg_day)"
    legend_states[2] = "RP_FSH in component RP_FSH (microg)"
    legend_algebraic[11] = "syn_FSH_Ih in component RP_FSH (microg_day)"
    legend_constants[10] = "V_FSH in component RP_FSH (microg_day)"
    legend_constants[11] = "Ki_FSH_Ih in component RP_FSH (U_L)"
    legend_constants[12] = "kFSH_rel in component RP_FSH (first_order_rate_constant)"
    legend_constants[13] = "CFSH_P in component RP_FSH (L_nmol)"
    legend_constants[14] = "CFSH_E in component RP_FSH (L_ng2)"
    legend_algebraic[9] = "Ih_dIh in component Ih_dIh (U_L)"
    legend_states[3] = "FSH in component FSH (microg_l)"
    legend_constants[15] = "v_dis in component FSH (litre)"
    legend_algebraic[2] = "clear_FSH in component FSH (microg_l_day)"
    legend_constants[16] = "kFSH_cl in component FSH (first_order_rate_constant)"
    legend_constants[17] = "dE in component E2_dE (day)"
    legend_constants[18] = "dP in component P4_dP (day)"
    legend_algebraic[1] = "Ih in component Ih (U_L)"
    legend_constants[19] = "dIh in component Ih_dIh (day)"
    legend_rates[0] = "d/dt RP_LH in component RP_LH (microg)"
    legend_rates[1] = "d/dt LH in component LH (microg_l)"
    legend_rates[2] = "d/dt RP_FSH in component RP_FSH (microg)"
    legend_rates[3] = "d/dt FSH in component FSH (microg_l)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 467.0
    constants[0] = 1400.0
    constants[1] = 95900.0
    constants[2] = 8.0
    constants[3] = 360.0
    constants[4] = 26.0
    constants[5] = 3.0
    constants[6] = 0.024
    constants[7] = 0.008
    states[1] = 40.0
    constants[8] = 2.5
    constants[9] = 14.0
    states[2] = 0.0
    constants[10] = 4400.0
    constants[11] = 1176.5
    constants[12] = 45.0
    constants[13] = 3.0
    constants[14] = 0.005
    states[3] = 150.0
    constants[15] = 2.5
    constants[16] = 8.21
    constants[17] = 0.42
    constants[18] = 2.9
    constants[19] = 2.0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[3] = (300.000-(240.000*(power(voi+1.00000, 2.00000)))/(3.00000+power(voi+1.00000, 2.00000)))+90.0000*exp(-((power(voi-8.00000, 2.00000))/10.0000))
    algebraic[5] = 52.0000*exp(-((power(voi-7.00000, 2.00000))/18.0000))
    algebraic[6] = (constants[5]*(1.00000+constants[6]*algebraic[5])*states[0])/(1.00000+constants[7]*algebraic[3])
    algebraic[0] = constants[9]*states[1]
    rates[1] = algebraic[6]/constants[8]-algebraic[0]
    algebraic[7] = (constants[12]*(1.00000+constants[13]*algebraic[5])*states[2])/(1.00000+constants[14]*(power(algebraic[3], 2.00000)))
    algebraic[2] = constants[16]*states[3]
    rates[3] = algebraic[7]/constants[15]-algebraic[2]
    algebraic[4] = (300.000-(240.000*(power((voi+1.00000)-constants[17], 2.00000)))/(3.00000+power((voi+1.00000)-constants[17], 2.00000)))+90.0000*exp(-((power(voi-(constants[17]+8.00000), 2.00000))/10.0000))
    algebraic[8] = 52.0000*exp(-((power(voi-(constants[18]+7.00000), 2.00000))/18.0000))
    algebraic[10] = (constants[0]+(constants[1]*(power(algebraic[4], constants[2])))/(power(constants[3], constants[2])+power(algebraic[4], constants[2])))/(1.00000+algebraic[8]/constants[4])
    rates[0] = algebraic[10]-algebraic[6]
    algebraic[9] = 300.000+1330.00*exp(-((power(voi-(7.00000+constants[19]), 2.00000))/19.0000))
    algebraic[11] = constants[10]/(1.00000+algebraic[9]/constants[11])
    rates[2] = algebraic[11]-algebraic[7]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = (300.000-(240.000*(power(voi+1.00000, 2.00000)))/(3.00000+power(voi+1.00000, 2.00000)))+90.0000*exp(-((power(voi-8.00000, 2.00000))/10.0000))
    algebraic[5] = 52.0000*exp(-((power(voi-7.00000, 2.00000))/18.0000))
    algebraic[6] = (constants[5]*(1.00000+constants[6]*algebraic[5])*states[0])/(1.00000+constants[7]*algebraic[3])
    algebraic[0] = constants[9]*states[1]
    algebraic[7] = (constants[12]*(1.00000+constants[13]*algebraic[5])*states[2])/(1.00000+constants[14]*(power(algebraic[3], 2.00000)))
    algebraic[2] = constants[16]*states[3]
    algebraic[4] = (300.000-(240.000*(power((voi+1.00000)-constants[17], 2.00000)))/(3.00000+power((voi+1.00000)-constants[17], 2.00000)))+90.0000*exp(-((power(voi-(constants[17]+8.00000), 2.00000))/10.0000))
    algebraic[8] = 52.0000*exp(-((power(voi-(constants[18]+7.00000), 2.00000))/18.0000))
    algebraic[10] = (constants[0]+(constants[1]*(power(algebraic[4], constants[2])))/(power(constants[3], constants[2])+power(algebraic[4], constants[2])))/(1.00000+algebraic[8]/constants[4])
    algebraic[9] = 300.000+1330.00*exp(-((power(voi-(7.00000+constants[19]), 2.00000))/19.0000))
    algebraic[11] = constants[10]/(1.00000+algebraic[9]/constants[11])
    algebraic[1] = 300.000+1330.00*exp(-((power(voi-7.00000, 2.00000))/19.0000))
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