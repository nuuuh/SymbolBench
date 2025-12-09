# Size of variable arrays:
sizeAlgebraic = 9
sizeStates = 7
sizeConstants = 37
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "R_T in component Parameters (dimensionless)"
    legend_constants[1] = "K_1 in component Parameters (micromolar)"
    legend_constants[2] = "K_2 in component Parameters (micromolar)"
    legend_constants[3] = "k_r in component Parameters (per_second)"
    legend_constants[4] = "k_p in component Parameters (per_second)"
    legend_constants[5] = "k_e in component Parameters (per_second)"
    legend_constants[6] = "xi in component Parameters (dimensionless)"
    legend_constants[7] = "G_T in component Parameters (dimensionless)"
    legend_constants[8] = "k_deg in component Parameters (per_second)"
    legend_constants[9] = "k_a in component Parameters (per_second)"
    legend_constants[10] = "k_d in component Parameters (per_second)"
    legend_constants[11] = "PIP_2_T in component Parameters (dimensionless)"
    legend_constants[12] = "r_r in component Parameters (per_second)"
    legend_constants[13] = "delta in component Parameters (dimensionless)"
    legend_constants[14] = "K_c in component Parameters (micromolar)"
    legend_constants[15] = "alpha in component Parameters (per_second)"
    legend_constants[16] = "N_a in component Parameters (per_micromole)"
    legend_constants[17] = "v in component Parameters (litre)"
    legend_constants[18] = "epsilon_r in component Parameters (dimensionless)"
    legend_constants[19] = "d_1 in component Parameters (micromolar)"
    legend_constants[20] = "d_2 in component Parameters (micromolar)"
    legend_constants[21] = "d_3 in component Parameters (micromolar)"
    legend_constants[22] = "d_5 in component Parameters (micromolar)"
    legend_constants[23] = "a_2 in component Parameters (per_micromolar_per_second)"
    legend_constants[24] = "B_e in component Parameters (micromolar)"
    legend_constants[25] = "K_e in component Parameters (micromolar)"
    legend_constants[26] = "B_ER in component Parameters (micromolar)"
    legend_constants[27] = "K_ER in component Parameters (micromolar)"
    legend_constants[28] = "B_x in component Parameters (micromolar)"
    legend_constants[29] = "K_x in component Parameters (micromolar)"
    legend_constants[30] = "k_3 in component Parameters (micromolar)"
    legend_constants[31] = "eta_1 in component Parameters (per_second)"
    legend_constants[32] = "eta_2 in component Parameters (per_second)"
    legend_constants[33] = "eta_3 in component Parameters (flux)"
    legend_constants[34] = "C_T in component Parameters (micromolar)"
    legend_constants[35] = "L in component ligand (micromolar)"
    legend_states[0] = "RS in component RS (dimensionless)"
    legend_states[1] = "RS_p in component RS_p (dimensionless)"
    legend_states[2] = "G in component G_GTP (dimensionless)"
    legend_algebraic[0] = "rho_r in component rho_r (dimensionless)"
    legend_states[3] = "IP_3 in component IP_3 (micromolar)"
    legend_states[4] = "PIP_2 in component PIP_2 (dimensionless)"
    legend_algebraic[1] = "r_h in component r_h (per_second)"
    legend_states[5] = "C in component C (micromolar)"
    legend_algebraic[8] = "C_ER in component C_ER (micromolar)"
    legend_algebraic[4] = "beta in component beta (dimensionless)"
    legend_states[6] = "h in component h (dimensionless)"
    legend_algebraic[2] = "m_infinit in component m_infinit (dimensionless)"
    legend_algebraic[7] = "h_infinit in component h_infinit (dimensionless)"
    legend_algebraic[5] = "tau_h in component tau_h (second)"
    legend_algebraic[3] = "zeta in component zeta (micromolar)"
    legend_algebraic[6] = "gamma in component gamma (dimensionless)"
    legend_constants[36] = "RS_E in component RS_E (dimensionless)"
    legend_rates[0] = "d/dt RS in component RS (dimensionless)"
    legend_rates[1] = "d/dt RS_p in component RS_p (dimensionless)"
    legend_rates[2] = "d/dt G in component G_GTP (dimensionless)"
    legend_rates[3] = "d/dt IP_3 in component IP_3 (micromolar)"
    legend_rates[4] = "d/dt PIP_2 in component PIP_2 (dimensionless)"
    legend_rates[5] = "d/dt C in component C (micromolar)"
    legend_rates[6] = "d/dt h in component h (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 2e4
    constants[1] = 5
    constants[2] = 100
    constants[3] = 1.75e-4
    constants[4] = 0.03
    constants[5] = 6e-3
    constants[6] = 0.85
    constants[7] = 1e5
    constants[8] = 1.25
    constants[9] = 0.017
    constants[10] = 0.15
    constants[11] = 5e7
    constants[12] = 10
    constants[13] = 1.238e-3
    constants[14] = 0.4
    constants[15] = 2.781e-5
    constants[16] = 6.02252e17
    constants[17] = 5e-13
    constants[18] = 0.185
    constants[19] = 0.13
    constants[20] = 1.05
    constants[21] = 0.943
    constants[22] = 0.0823
    constants[23] = 0.2
    constants[24] = 150
    constants[25] = 10
    constants[26] = 120000
    constants[27] = 1200
    constants[28] = 50
    constants[29] = 0.2
    constants[30] = 0.4
    constants[31] = 575
    constants[32] = 5.2
    constants[33] = 45
    constants[34] = 67
    constants[35] = 1000
    states[0] = 1.7e4
    states[1] = 0
    states[2] = 0
    states[3] = 0.01
    states[4] = 49997000
    states[5] = 0.0961
    states[6] = 0.6155
    constants[36] = ((constants[3]*(1.00000+((constants[4]/constants[5])*(constants[2]+constants[35]))/(constants[1]+constants[35])))/(constants[3]+(constants[4]*constants[35])/(constants[1]+constants[35])+(((constants[3]*constants[4])/constants[5])*(constants[2]+constants[35]))/(constants[1]+constants[35])))*constants[6]*constants[0]+(1.00000-constants[6])*constants[0]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = (constants[3]*constants[0]-(constants[3]+(constants[4]*constants[35])/(constants[1]+constants[35]))*states[0])-constants[3]*states[1]
    rates[1] = constants[35]*((constants[4]*states[0])/(constants[1]+constants[35])-(constants[5]*states[1])/(constants[2]+constants[35]))
    algebraic[0] = (constants[35]*states[0])/(constants[6]*constants[0]*(constants[1]+constants[35]))
    rates[2] = constants[9]*(constants[13]+algebraic[0])*(constants[7]-states[2])-constants[10]*states[2]
    algebraic[1] = ((constants[15]*states[5])/(constants[14]+states[5]))*states[2]
    rates[3] = (algebraic[1]*states[4])/(constants[16]*constants[17])-constants[8]*states[3]
    rates[4] = (-(algebraic[1]+constants[12])*states[4]-constants[12]*constants[16]*constants[17]*states[3])+constants[12]*constants[11]
    algebraic[3] = (constants[20]*(states[3]+constants[19]))/(states[3]+constants[21])
    algebraic[7] = algebraic[3]/(algebraic[3]+states[5])
    algebraic[5] = power(constants[23]*(algebraic[3]+states[5]), -1.00000)
    rates[6] = (algebraic[7]-states[6])/algebraic[5]
    algebraic[6] = power(1.00000+constants[24]/(constants[25]+states[5])+constants[28]/(constants[29]+states[5]), -1.00000)
    algebraic[8] = (constants[27]/(constants[26]*constants[18]))*(constants[34]-states[5]/algebraic[6])
    algebraic[4] = power(1.00000+(constants[25]*constants[24])/(power(constants[25]+states[5], 2.00000))+(constants[29]*constants[28])/(power(constants[29]+states[5], 2.00000)), -1.00000)
    algebraic[2] = ((states[3]/(constants[19]+states[3]))*states[5])/(constants[22]+states[5])
    rates[5] = algebraic[4]*(constants[18]*(constants[31]*(power(algebraic[2], 3.00000))*(power(states[6], 3.00000))+constants[32])*(algebraic[8]-states[5])-(constants[33]*(power(states[5], 2.00000)))/(power(constants[30], 2.00000)+power(states[5], 2.00000)))
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = (constants[35]*states[0])/(constants[6]*constants[0]*(constants[1]+constants[35]))
    algebraic[1] = ((constants[15]*states[5])/(constants[14]+states[5]))*states[2]
    algebraic[3] = (constants[20]*(states[3]+constants[19]))/(states[3]+constants[21])
    algebraic[7] = algebraic[3]/(algebraic[3]+states[5])
    algebraic[5] = power(constants[23]*(algebraic[3]+states[5]), -1.00000)
    algebraic[6] = power(1.00000+constants[24]/(constants[25]+states[5])+constants[28]/(constants[29]+states[5]), -1.00000)
    algebraic[8] = (constants[27]/(constants[26]*constants[18]))*(constants[34]-states[5]/algebraic[6])
    algebraic[4] = power(1.00000+(constants[25]*constants[24])/(power(constants[25]+states[5], 2.00000))+(constants[29]*constants[28])/(power(constants[29]+states[5], 2.00000)), -1.00000)
    algebraic[2] = ((states[3]/(constants[19]+states[3]))*states[5])/(constants[22]+states[5])
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