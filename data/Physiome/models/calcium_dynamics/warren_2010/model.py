# Size of variable arrays:
sizeAlgebraic = 26
sizeStates = 8
sizeConstants = 46
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "volume_i0 in component volume_i (volume)"
    legend_constants[43] = "volume_i in component volume_i (volume)"
    legend_constants[1] = "gamma in component volume_er (dimensionless)"
    legend_constants[44] = "volume_er in component volume_er (volume)"
    legend_states[0] = "volumeCa_i in component volumeCa_i (microMVolume)"
    legend_algebraic[13] = "Jipt in component Jipt (microMvolumepersec)"
    legend_algebraic[15] = "Jryr in component Jryr (microMvolumepersec)"
    legend_algebraic[17] = "Jer in component Jer (microMvolumepersec)"
    legend_algebraic[18] = "Jserca in component Jserca (microMvolumepersec)"
    legend_algebraic[24] = "Jin in component Jin (microMvolumepersec)"
    legend_algebraic[25] = "Jpm in component Jpm (microMvolumepersec)"
    legend_algebraic[0] = "Ca_i in component Ca_i (microM)"
    legend_states[1] = "volumeCa_er in component volumeCa_er (microMVolume)"
    legend_algebraic[3] = "Ca_er in component Ca_er (microM)"
    legend_states[2] = "volume_iIP3 in component volume_iIP3 (microMVolume)"
    legend_algebraic[12] = "J_ip3P in component J_ip3P (microMpersec)"
    legend_algebraic[16] = "J_ip3D in component J_ip3D (microMpersec)"
    legend_algebraic[1] = "IP3 in component IP3 (microM)"
    legend_states[3] = "PIP2 in component PIP2 (dimensionless)"
    legend_constants[2] = "PIP2_Total in component PIP2 (microM)"
    legend_constants[45] = "J_ip3R in component J_ip3R (per_sec)"
    legend_algebraic[4] = "IPX in component IPX (dimensionless)"
    legend_constants[3] = "PIP2_Total in component IPX (microM)"
    legend_constants[4] = "J_ip3R0 in component J_ip3R (per_sec)"
    legend_states[4] = "ATP_e in component receptor (microM)"
    legend_algebraic[9] = "V_IP3 in component receptor (microMpersec)"
    legend_constants[5] = "MechanicalStimulation in component receptor (dimensionless)"
    legend_algebraic[6] = "V_mech in component receptor (dimensionless)"
    legend_constants[6] = "V_mech0 in component receptor (dimensionless)"
    legend_constants[7] = "C_1 in component receptor (microMpersec)"
    legend_constants[8] = "C_2 in component receptor (microM)"
    legend_constants[9] = "C_3 in component receptor (microMpersec)"
    legend_constants[10] = "C_4 in component receptor (microM)"
    legend_constants[11] = "V_ATP in component receptor (microMpersec)"
    legend_constants[12] = "K_ATP in component receptor (microM)"
    legend_constants[13] = "Beta_1 in component J_ip3D (microMpersec)"
    legend_constants[14] = "Beta_2 in component J_ip3D (microMpersec)"
    legend_algebraic[14] = "F_Ca in component J_ip3D (dimensionless)"
    legend_constants[15] = "K_rc in component J_ip3D (microM)"
    legend_constants[16] = "K_IP3 in component J_ip3D (microM)"
    legend_algebraic[10] = "Pipt in component Jipt (dimensionless)"
    legend_states[5] = "O in component O (dimensionless)"
    legend_constants[17] = "Jipt0 in component Jipt (per_sec)"
    legend_constants[18] = "Ka in component Jryr (microM)"
    legend_constants[19] = "Kb in component Jryr (microM)"
    legend_constants[20] = "Kc in component Jryr (dimensionless)"
    legend_constants[21] = "Jryr0 in component Jryr (per_sec)"
    legend_algebraic[5] = "Pryr in component Jryr (dimensionless)"
    legend_algebraic[2] = "W_inf in component Jryr (dimensionless)"
    legend_constants[22] = "Jer0 in component Jer (per_sec)"
    legend_constants[23] = "Vserca in component Jserca (microM2persec)"
    legend_constants[24] = "Kserca in component Jserca (microM)"
    legend_constants[25] = "Vm in component Jin (milliV)"
    legend_algebraic[19] = "V_Ca in component Jin (milliV)"
    legend_constants[26] = "R in component Jin (millijoule_per_mole_kelvin)"
    legend_constants[27] = "T in component Jin (kelvin)"
    legend_constants[28] = "F in component Jin (coulomb_per_mole)"
    legend_constants[29] = "z_Ca in component Jin (dimensionless)"
    legend_constants[30] = "Ca_ext in component Jin (microM)"
    legend_constants[31] = "P_pm_1 in component Jin (conductancepervolume)"
    legend_constants[32] = "P_pm_2 in component Jin (conductancepervolume)"
    legend_algebraic[20] = "Iin_1 in component Jin (current)"
    legend_algebraic[21] = "Iin_2 in component Jin (current)"
    legend_algebraic[22] = "Jin_1 in component Jin (microMvolumepersec)"
    legend_algebraic[23] = "Jin_2 in component Jin (microMvolumepersec)"
    legend_constants[33] = "Vpm in component Jpm (microMpersec)"
    legend_constants[34] = "Kpm in component Jpm (microM)"
    legend_algebraic[7] = "S in component O (dimensionless)"
    legend_states[6] = "I1 in component I1 (dimensionless)"
    legend_states[7] = "I2 in component I2 (dimensionless)"
    legend_algebraic[11] = "k1 in component constants (per_secmicroM)"
    legend_constants[35] = "k1_ in component constants (first_order_rate_constant)"
    legend_constants[36] = "k2 in component constants (first_order_rate_constant)"
    legend_constants[37] = "k3 in component constants (first_order_rate_constant)"
    legend_algebraic[8] = "k4 in component constants (first_order_rate_constant)"
    legend_constants[38] = "k5 in component constants (first_order_rate_constant)"
    legend_constants[39] = "alpha_1 in component constants (per_secmicroM)"
    legend_constants[40] = "beta_1 in component constants (microM)"
    legend_constants[41] = "alpha_4 in component constants (first_order_rate_constant)"
    legend_constants[42] = "beta_4 in component constants (microM)"
    legend_rates[0] = "d/dt volumeCa_i in component volumeCa_i (microMVolume)"
    legend_rates[1] = "d/dt volumeCa_er in component volumeCa_er (microMVolume)"
    legend_rates[2] = "d/dt volume_iIP3 in component volume_iIP3 (microMVolume)"
    legend_rates[3] = "d/dt PIP2 in component PIP2 (dimensionless)"
    legend_rates[4] = "d/dt ATP_e in component receptor (microM)"
    legend_rates[5] = "d/dt O in component O (dimensionless)"
    legend_rates[6] = "d/dt I1 in component I1 (dimensionless)"
    legend_rates[7] = "d/dt I2 in component I2 (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1450
    constants[1] = 4.9
    states[0] = 75.8335890358239
    states[1] = 1199.73747802144
    states[2] = -3.28287434380887e-93
    states[3] = 1
    constants[2] = 5
    constants[3] = 25
    constants[4] = 10
    states[4] = 0
    constants[5] = 0
    constants[6] = 1
    constants[7] = 2.6
    constants[8] = 0.02
    constants[9] = 2.3
    constants[10] = 10
    constants[11] = 0.001
    constants[12] = 0.0001
    constants[13] = 0
    constants[14] = 36
    constants[15] = 0.3
    constants[16] = 0.01
    states[5] = -2.42013500345263e-89
    constants[17] = 28
    constants[18] = 0.3722419
    constants[19] = 0.636005
    constants[20] = 0.0571428
    constants[21] = 0
    constants[22] = 0.0035
    constants[23] = 0.09
    constants[24] = 0.04
    constants[25] = -39
    constants[26] = 8314
    constants[27] = 310
    constants[28] = 96485
    constants[29] = 2
    constants[30] = 1300
    constants[31] = 0.6
    constants[32] = 28
    constants[33] = 0.072
    constants[34] = 0.6
    states[6] = 2.24343309680416e-81
    states[7] = 1.01276768045278e-40
    constants[35] = 0.88
    constants[36] = 0.5
    constants[37] = 0.5
    constants[38] = 0.02
    constants[39] = 40
    constants[40] = 0.8
    constants[41] = 0.06
    constants[42] = 0.01
    constants[43] = constants[0]
    constants[44] = constants[43]/constants[1]
    constants[45] = constants[4]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[4] = (-constants[11]*states[4])/(constants[12]+states[4])
    algebraic[1] = states[2]/constants[43]
    algebraic[8] = (constants[41]*algebraic[1])/(constants[42]+algebraic[1])
    rates[6] = constants[36]*states[5]-(constants[37]+algebraic[8])*states[6]
    rates[7] = algebraic[8]*states[6]-constants[38]*states[7]
    algebraic[7] = 1.00000-(states[5]+states[6]+states[7])
    algebraic[0] = states[0]/constants[43]
    algebraic[11] = (constants[39]*(power(algebraic[0], 3.00000)))/(power(constants[40], 3.00000)+power(algebraic[0], 3.00000))
    rates[5] = algebraic[11]*algebraic[1]*algebraic[7]-(constants[35]*states[5]+constants[36]*states[5])
    algebraic[9] = (constants[7]*states[4])/(constants[8]+states[4])+(constants[9]*(power(states[4], 2.00000)))/(power(constants[10], 2.00000)+power(states[4], 2.00000))
    algebraic[12] = algebraic[9]*states[3]
    algebraic[4] = (1.00000-algebraic[1]/constants[3])-states[3]
    rates[3] = constants[45]*algebraic[4]-algebraic[12]/constants[2]
    algebraic[14] = algebraic[0]/(constants[15]+algebraic[0])
    algebraic[16] = ((constants[13]+constants[14]*(power(algebraic[14], 2.00000)))*algebraic[1])/(constants[16]+algebraic[1])
    rates[2] = (algebraic[12]-algebraic[16])*constants[43]
    algebraic[3] = states[1]/constants[44]
    algebraic[10] = power(states[5], 4.00000)
    algebraic[13] = constants[17]*constants[0]*algebraic[10]*(algebraic[3]-algebraic[0])
    algebraic[15] = constants[21]*constants[0]*(algebraic[3]-algebraic[0])
    algebraic[17] = constants[22]*constants[0]*(algebraic[3]-algebraic[0])
    algebraic[18] = (((constants[0]*constants[23]*(power(algebraic[0], 2.00000)))/(power(constants[24], 2.00000)+power(algebraic[0], 2.00000)))*1.00000)/algebraic[3]
    rates[1] = -((algebraic[13]+algebraic[15]+algebraic[17])-algebraic[18])
    algebraic[19] = ((constants[26]*constants[27])/(constants[29]*constants[28]))*log(constants[30]/algebraic[0])
    algebraic[20] = constants[0]*constants[31]*(algebraic[19]-constants[25])
    algebraic[22] = algebraic[20]/(constants[29]*constants[28])
    algebraic[6] = custom_piecewise([greater(voi , 10.0000) & less(voi , 25.0000) & greater(constants[5] , 0.00000), constants[6] , True, 0.00000])
    algebraic[21] = constants[0]*constants[32]*algebraic[6]*(algebraic[19]-constants[25])
    algebraic[23] = algebraic[21]/(constants[29]*constants[28])
    algebraic[24] = algebraic[22]+algebraic[23]
    algebraic[25] = (constants[33]*constants[0]*(power(algebraic[0], 2.00000)))/(power(constants[34], 2.00000)+power(algebraic[0], 2.00000))
    rates[0] = ((algebraic[13]+algebraic[15]+algebraic[17])-algebraic[18])+(algebraic[24]-algebraic[25])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = states[2]/constants[43]
    algebraic[8] = (constants[41]*algebraic[1])/(constants[42]+algebraic[1])
    algebraic[7] = 1.00000-(states[5]+states[6]+states[7])
    algebraic[0] = states[0]/constants[43]
    algebraic[11] = (constants[39]*(power(algebraic[0], 3.00000)))/(power(constants[40], 3.00000)+power(algebraic[0], 3.00000))
    algebraic[9] = (constants[7]*states[4])/(constants[8]+states[4])+(constants[9]*(power(states[4], 2.00000)))/(power(constants[10], 2.00000)+power(states[4], 2.00000))
    algebraic[12] = algebraic[9]*states[3]
    algebraic[4] = (1.00000-algebraic[1]/constants[3])-states[3]
    algebraic[14] = algebraic[0]/(constants[15]+algebraic[0])
    algebraic[16] = ((constants[13]+constants[14]*(power(algebraic[14], 2.00000)))*algebraic[1])/(constants[16]+algebraic[1])
    algebraic[3] = states[1]/constants[44]
    algebraic[10] = power(states[5], 4.00000)
    algebraic[13] = constants[17]*constants[0]*algebraic[10]*(algebraic[3]-algebraic[0])
    algebraic[15] = constants[21]*constants[0]*(algebraic[3]-algebraic[0])
    algebraic[17] = constants[22]*constants[0]*(algebraic[3]-algebraic[0])
    algebraic[18] = (((constants[0]*constants[23]*(power(algebraic[0], 2.00000)))/(power(constants[24], 2.00000)+power(algebraic[0], 2.00000)))*1.00000)/algebraic[3]
    algebraic[19] = ((constants[26]*constants[27])/(constants[29]*constants[28]))*log(constants[30]/algebraic[0])
    algebraic[20] = constants[0]*constants[31]*(algebraic[19]-constants[25])
    algebraic[22] = algebraic[20]/(constants[29]*constants[28])
    algebraic[6] = custom_piecewise([greater(voi , 10.0000) & less(voi , 25.0000) & greater(constants[5] , 0.00000), constants[6] , True, 0.00000])
    algebraic[21] = constants[0]*constants[32]*algebraic[6]*(algebraic[19]-constants[25])
    algebraic[23] = algebraic[21]/(constants[29]*constants[28])
    algebraic[24] = algebraic[22]+algebraic[23]
    algebraic[25] = (constants[33]*constants[0]*(power(algebraic[0], 2.00000)))/(power(constants[34], 2.00000)+power(algebraic[0], 2.00000))
    algebraic[2] = (1.00000+power(constants[18]/algebraic[0], 4.00000)+power(algebraic[0]/constants[19], 3.00000))/(1.00000+1.00000/constants[20]+power(algebraic[0]/constants[19], 3.00000)+power(constants[18]/algebraic[0], 4.00000))
    algebraic[5] = (algebraic[2]*(1.00000+power(constants[18]/algebraic[0], 4.00000)+power(algebraic[0]/constants[19], 3.00000)))/(1.00000+1.00000/constants[20]+power(algebraic[0]/constants[19], 3.00000)+power(constants[18]/algebraic[0], 4.00000))
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