# Size of variable arrays:
sizeAlgebraic = 15
sizeStates = 7
sizeConstants = 40
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "cap in component parameters (nanofarad)"
    legend_constants[1] = "cc1lin in component parameters (per_second)"
    legend_constants[2] = "cc_2 in component parameters (per_second)"
    legend_constants[3] = "ck1lin in component parameters (per_second)"
    legend_constants[4] = "ck_2 in component parameters (per_second)"
    legend_constants[5] = "clmax in component parameters (nanosiemens)"
    legend_constants[6] = "cnmax in component parameters (nanosiemens)"
    legend_constants[7] = "cx1lin in component parameters (per_second)"
    legend_constants[8] = "cx2 in component parameters (per_second)"
    legend_constants[9] = "ef in component parameters (per_second)"
    legend_constants[10] = "gl in component parameters (nanosiemens)"
    legend_constants[11] = "hmc_1 in component parameters (uM)"
    legend_constants[12] = "hmc_2 in component parameters (uM)"
    legend_constants[13] = "inf in component parameters (uM_per_picocoulomb)"
    legend_constants[14] = "inhmax in component parameters (dimensionless)"
    legend_constants[15] = "k_1 in component parameters (per_uM_per_second)"
    legend_constants[16] = "k_2 in component parameters (per_second)"
    legend_constants[17] = "kI in component parameters (uM)"
    legend_constants[18] = "kinh in component parameters (uM)"
    legend_constants[19] = "kinhcng in component parameters (uM)"
    legend_constants[20] = "n_1 in component parameters (dimensionless)"
    legend_constants[21] = "n_2 in component parameters (dimensionless)"
    legend_constants[22] = "nI in component parameters (dimensionless)"
    legend_constants[23] = "ninh in component parameters (dimensionless)"
    legend_constants[24] = "ninhcng in component parameters (dimensionless)"
    legend_constants[25] = "pd in component parameters (per_second)"
    legend_constants[26] = "r_1 in component parameters (per_second)"
    legend_constants[27] = "r_2 in component parameters (per_second)"
    legend_constants[28] = "smax in component parameters (uM_per_second)"
    legend_constants[29] = "V_Cl in component parameters (millivolt)"
    legend_constants[30] = "V_cng in component parameters (millivolt)"
    legend_constants[31] = "V_l in component parameters (millivolt)"
    legend_constants[39] = "F_vol in component parameters (picocoulomb_per_uM)"
    legend_constants[32] = "F in component parameters (coulombs_per_mole)"
    legend_constants[33] = "C_vol in component parameters (liter)"
    legend_algebraic[9] = "O_stim in component O_stim (uM)"
    legend_constants[34] = "od in component O_stim (uM)"
    legend_constants[35] = "t_0 in component O_stim (second)"
    legend_constants[36] = "t_1 in component O_stim (second)"
    legend_algebraic[0] = "H_0 in component O_stim (dimensionless)"
    legend_algebraic[5] = "H_1 in component O_stim (dimensionless)"
    legend_states[0] = "bLR in component bLR (dimensionless)"
    legend_constants[37] = "R_tot in component bLR (dimensionless)"
    legend_states[1] = "aG in component aG (dimensionless)"
    legend_constants[38] = "G_tot in component aG (dimensionless)"
    legend_algebraic[1] = "k_G in component k_G (per_second)"
    legend_algebraic[6] = "r_G in component r_G (per_second)"
    legend_states[2] = "cAMP in component cAMP (uM)"
    legend_algebraic[2] = "synth in component synth (uM_per_second)"
    legend_algebraic[7] = "degrad in component degrad (uM_per_second)"
    legend_states[3] = "aCaMK in component aCaMK (uM)"
    legend_states[4] = "Ca in component Ca (uM)"
    legend_algebraic[10] = "I_CNG in component I_CNG (nanoampere)"
    legend_algebraic[12] = "J_NCX in component J_NCX (uM_per_second)"
    legend_algebraic[3] = "cc_1 in component cc_1 (uM_per_second)"
    legend_states[5] = "CaCaM in component CaCaM (uM)"
    legend_algebraic[4] = "ck_1 in component ck_1 (uM_per_second)"
    legend_states[6] = "V in component V (millivolt)"
    legend_algebraic[11] = "I_ClCa in component I_ClCa (nanoampere)"
    legend_algebraic[13] = "I_NCX in component I_NCX (nanoampere)"
    legend_algebraic[14] = "I_other in component I_other (nanoampere)"
    legend_algebraic[8] = "inhcng in component inhcng (dimensionless)"
    legend_rates[0] = "d/dt bLR in component bLR (dimensionless)"
    legend_rates[1] = "d/dt aG in component aG (dimensionless)"
    legend_rates[2] = "d/dt cAMP in component cAMP (uM)"
    legend_rates[4] = "d/dt Ca in component Ca (uM)"
    legend_rates[5] = "d/dt CaCaM in component CaCaM (uM)"
    legend_rates[3] = "d/dt aCaMK in component aCaMK (uM)"
    legend_rates[6] = "d/dt V in component V (millivolt)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 0.004
    constants[1] = 0.88
    constants[2] = 26
    constants[3] = 13
    constants[4] = 0.9
    constants[5] = 1
    constants[6] = 1
    constants[7] = 1
    constants[8] = 13
    constants[9] = 2
    constants[10] = 6
    constants[11] = 2
    constants[12] = 3
    constants[13] = 1.9
    constants[14] = 5
    constants[15] = 0.06
    constants[16] = 20
    constants[17] = 0.7
    constants[18] = 2
    constants[19] = 1
    constants[20] = 2
    constants[21] = 2
    constants[22] = 2
    constants[23] = 1.5
    constants[24] = 1.3
    constants[25] = 20
    constants[26] = 10
    constants[27] = 5
    constants[28] = 71
    constants[29] = -50
    constants[30] = 0
    constants[31] = -70
    constants[32] = 9.649e4
    constants[33] = 1e-13
    constants[34] = 20
    constants[35] = 0.5
    constants[36] = 1.5
    states[0] = 0
    constants[37] = 1
    states[1] = 0
    constants[38] = 1
    states[2] = 1.35648992164649e-88
    states[3] = 6.60756525051462e-8
    states[4] = 5.09073088043779e-12
    states[5] = 1.86113118246926e-13
    states[6] = -70
    constants[39] = (1.00000e+12/1.00000)*(1.00000/1000.00)*constants[32]*constants[33]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[3] = constants[1]*states[4]
    rates[5] = algebraic[3]-constants[2]*states[5]
    algebraic[4] = constants[3]*states[5]
    rates[3] = algebraic[4]-constants[4]*states[3]
    algebraic[1] = constants[16]*states[0]
    algebraic[6] = constants[27]*states[1]
    rates[1] = algebraic[1]*(constants[38]-states[1])-algebraic[6]
    algebraic[2] = (states[1]*constants[28])/(1.00000+power(states[3]/constants[18], constants[23]))
    algebraic[7] = constants[25]*states[2]
    rates[2] = algebraic[2]-algebraic[7]
    algebraic[0] = custom_piecewise([less(voi , constants[35]), 0.00000 , True, 1.00000])
    algebraic[5] = custom_piecewise([less(voi , constants[36]), 0.00000 , True, 1.00000])
    algebraic[9] = constants[34]*(algebraic[0]-algebraic[5])
    rates[0] = constants[15]*algebraic[9]*(constants[37]-states[0])-constants[26]*states[0]
    algebraic[8] = 1.00000+((constants[14]-1.00000)*(power(states[5], constants[24])))/(power(states[5], constants[24])+power(constants[19], constants[24]))
    algebraic[10] = ((constants[6]*(power(states[2], constants[20])))/(power(states[2], constants[20])+power(algebraic[8]*constants[11], constants[20])))*(1.00000/1000.00)*(constants[30]-states[6])
    algebraic[12] = constants[9]*states[4]
    rates[4] = ((1000.00/1.00000)*constants[13]*algebraic[10]-algebraic[12])-(algebraic[3]-constants[2]*states[5])
    algebraic[11] = ((constants[5]*(power(states[4], constants[21])))/(power(states[4], constants[21])+power(constants[12], constants[21])))*(1.00000/1000.00)*(constants[29]-states[6])
    algebraic[13] = (1.00000/1000.00)*constants[39]*algebraic[12]
    algebraic[14] = constants[10]*(1.00000/1000.00)*(constants[31]-states[6])
    rates[6] = (1000.00/1.00000)*(1.00000/constants[0])*(algebraic[10]+algebraic[11]+algebraic[13]+algebraic[14])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = constants[1]*states[4]
    algebraic[4] = constants[3]*states[5]
    algebraic[1] = constants[16]*states[0]
    algebraic[6] = constants[27]*states[1]
    algebraic[2] = (states[1]*constants[28])/(1.00000+power(states[3]/constants[18], constants[23]))
    algebraic[7] = constants[25]*states[2]
    algebraic[0] = custom_piecewise([less(voi , constants[35]), 0.00000 , True, 1.00000])
    algebraic[5] = custom_piecewise([less(voi , constants[36]), 0.00000 , True, 1.00000])
    algebraic[9] = constants[34]*(algebraic[0]-algebraic[5])
    algebraic[8] = 1.00000+((constants[14]-1.00000)*(power(states[5], constants[24])))/(power(states[5], constants[24])+power(constants[19], constants[24]))
    algebraic[10] = ((constants[6]*(power(states[2], constants[20])))/(power(states[2], constants[20])+power(algebraic[8]*constants[11], constants[20])))*(1.00000/1000.00)*(constants[30]-states[6])
    algebraic[12] = constants[9]*states[4]
    algebraic[11] = ((constants[5]*(power(states[4], constants[21])))/(power(states[4], constants[21])+power(constants[12], constants[21])))*(1.00000/1000.00)*(constants[29]-states[6])
    algebraic[13] = (1.00000/1000.00)*constants[39]*algebraic[12]
    algebraic[14] = constants[10]*(1.00000/1000.00)*(constants[31]-states[6])
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