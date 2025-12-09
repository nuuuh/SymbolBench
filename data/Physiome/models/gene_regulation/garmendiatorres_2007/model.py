# Size of variable arrays:
sizeAlgebraic = 5
sizeStates = 11
sizeConstants = 55
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_constants[0] = "V_max1 in component V1 (micromolar_per_minute)"
    legend_constants[1] = "GEFt in component V1 (micromolar)"
    legend_constants[43] = "V_1 in component V1 (per_minute)"
    legend_constants[2] = "Str in component V2 (dimensionless)"
    legend_constants[3] = "V_max2 in component V2 (micromolar_per_minute)"
    legend_constants[44] = "V_2 in component V2 (per_minute)"
    legend_constants[4] = "k_c3 in component V3 (per_minute)"
    legend_constants[5] = "PKAt in component V3 (micromolar)"
    legend_constants[6] = "GAPt in component V3 (micromolar)"
    legend_constants[45] = "V_3 in component V3 (per_minute)"
    legend_constants[7] = "V_max4 in component V4 (micromolar_per_minute)"
    legend_constants[46] = "V_4 in component V4 (per_minute)"
    legend_constants[8] = "k_gef in component V5 (per_minute)"
    legend_constants[9] = "RASt in component V5 (micromolar)"
    legend_constants[47] = "V_5 in component V5 (per_minute)"
    legend_constants[10] = "k_gap in component V6 (per_minute)"
    legend_constants[48] = "V_6 in component V6 (per_minute)"
    legend_constants[11] = "k_c7 in component V7 (per_minute)"
    legend_constants[12] = "PDEt in component V7 (micromolar)"
    legend_constants[49] = "V_7 in component V7 (per_minute)"
    legend_constants[13] = "V_max8 in component V8 (micromolar_per_minute)"
    legend_constants[50] = "V_8 in component V8 (per_minute)"
    legend_constants[14] = "a in component VPKAact (per_micromolar_squared_minute)"
    legend_constants[15] = "r in component VPKAact (per_minute)"
    legend_states[0] = "R2C2 in component holoenzyme_R_C (dimensionless)"
    legend_states[1] = "cAMP in component cyclic_AMP (micromolar)"
    legend_algebraic[3] = "C in component C_subunit (dimensionless)"
    legend_algebraic[0] = "R2cAMP2 in component holoenzyme_R_cAMP (dimensionless)"
    legend_algebraic[4] = "V_PKAact in component VPKAact (per_minute)"
    legend_constants[16] = "K_1 in component active_GEF (dimensionless)"
    legend_constants[17] = "K_2 in component active_GEF (dimensionless)"
    legend_states[2] = "GEFa in component active_GEF (dimensionless)"
    legend_constants[18] = "K_3 in component active_GAP (dimensionless)"
    legend_constants[19] = "K_4 in component active_GAP (dimensionless)"
    legend_states[3] = "GAPa in component active_GAP (dimensionless)"
    legend_constants[20] = "K_5 in component RAS_to_GTP (dimensionless)"
    legend_constants[21] = "K_6 in component RAS_to_GTP (dimensionless)"
    legend_states[4] = "RGTP in component RAS_to_GTP (dimensionless)"
    legend_constants[22] = "k_a in component adenylate_cyclase (per_micromolar_minute)"
    legend_constants[23] = "k_i in component adenylate_cyclase (per_minute)"
    legend_states[5] = "CYCLa in component adenylate_cyclase (dimensionless)"
    legend_constants[24] = "K_7 in component active_PDE (dimensionless)"
    legend_constants[25] = "K_8 in component active_PDE (dimensionless)"
    legend_states[6] = "PDEa in component active_PDE (dimensionless)"
    legend_constants[26] = "k_s in component cyclic_AMP (per_minute)"
    legend_constants[27] = "k_d in component cyclic_AMP (per_minute)"
    legend_constants[28] = "CYCLt in component cyclic_AMP (micromolar)"
    legend_constants[29] = "K_md in component cyclic_AMP (micromolar)"
    legend_constants[30] = "k_c9 in component V9 (per_minute)"
    legend_constants[31] = "MSNt in component V9 (micromolar)"
    legend_constants[51] = "V_9 in component V9 (per_minute)"
    legend_constants[32] = "V_max10 in component V10 (micromolar_per_minute)"
    legend_constants[52] = "V_10 in component V10 (per_minute)"
    legend_constants[33] = "k_c11 in component V11 (per_minute)"
    legend_constants[53] = "V_11 in component V11 (per_minute)"
    legend_constants[34] = "V_max12 in component V12 (micromolar_per_minute)"
    legend_constants[54] = "V_12 in component V12 (per_minute)"
    legend_constants[35] = "k_t1 in component cytosol (per_minute)"
    legend_constants[36] = "k_t2 in component cytosol (per_minute)"
    legend_constants[37] = "K_11 in component cytosol (dimensionless)"
    legend_constants[38] = "K_12 in component cytosol (dimensionless)"
    legend_states[7] = "MN in component nucleus (dimensionless)"
    legend_states[8] = "MCP in component cytosol_phos (dimensionless)"
    legend_states[9] = "MC in component cytosol (dimensionless)"
    legend_constants[39] = "K_9 in component nucleus (dimensionless)"
    legend_constants[40] = "K_10 in component nucleus (dimensionless)"
    legend_states[10] = "MNP in component nucleus_phos (dimensionless)"
    legend_constants[41] = "k_t3 in component nucleus_phos (per_minute)"
    legend_constants[42] = "k_t4 in component nucleus_phos (per_minute)"
    legend_algebraic[1] = "M_cyto in component Mcyto (dimensionless)"
    legend_algebraic[2] = "M_nucl in component Mnucl (dimensionless)"
    legend_rates[2] = "d/dt GEFa in component active_GEF (dimensionless)"
    legend_rates[3] = "d/dt GAPa in component active_GAP (dimensionless)"
    legend_rates[4] = "d/dt RGTP in component RAS_to_GTP (dimensionless)"
    legend_rates[5] = "d/dt CYCLa in component adenylate_cyclase (dimensionless)"
    legend_rates[6] = "d/dt PDEa in component active_PDE (dimensionless)"
    legend_rates[1] = "d/dt cAMP in component cyclic_AMP (micromolar)"
    legend_rates[0] = "d/dt R2C2 in component holoenzyme_R_C (dimensionless)"
    legend_rates[9] = "d/dt MC in component cytosol (dimensionless)"
    legend_rates[7] = "d/dt MN in component nucleus (dimensionless)"
    legend_rates[10] = "d/dt MNP in component nucleus_phos (dimensionless)"
    legend_rates[8] = "d/dt MCP in component cytosol_phos (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 1
    constants[1] = 4
    constants[2] = 1
    constants[3] = 1
    constants[4] = 3.5
    constants[5] = 0.3
    constants[6] = 1.5
    constants[7] = 1.3
    constants[8] = 240
    constants[9] = 250
    constants[10] = 600
    constants[11] = 3.333
    constants[12] = 0.5
    constants[13] = 1.5
    constants[14] = 1
    constants[15] = 1
    states[0] = 0.5
    states[1] = 1
    constants[16] = 0.05
    constants[17] = 0.05
    states[2] = 0.36
    constants[18] = 0.01
    constants[19] = 0.01
    states[3] = 0.5
    constants[20] = 0.001
    constants[21] = 0.001
    states[4] = 0.1
    constants[22] = 0.01
    constants[23] = 1
    states[5] = 0.1
    constants[24] = 0.01
    constants[25] = 0.01
    states[6] = 0.5
    constants[26] = 4
    constants[27] = 100
    constants[28] = 0.7
    constants[29] = 20
    constants[30] = 3.333
    constants[31] = 1
    constants[32] = 0.6
    constants[33] = 3.333
    constants[34] = 2
    constants[35] = 10
    constants[36] = 0.001
    constants[37] = 0.05
    constants[38] = 0.05
    states[7] = 0.25
    states[8] = 0.25
    states[9] = 0.25
    constants[39] = 0.05
    constants[40] = 0.05
    states[10] = 0.25
    constants[41] = 0.001
    constants[42] = 10
    constants[43] = constants[0]/constants[1]
    constants[44] = (constants[2]*constants[3])/constants[1]
    constants[45] = (constants[4]*constants[5])/constants[6]
    constants[46] = constants[7]/constants[6]
    constants[47] = (constants[8]*constants[1])/constants[9]
    constants[48] = (constants[10]*constants[6])/constants[9]
    constants[49] = (constants[11]*constants[5])/constants[12]
    constants[50] = constants[13]/constants[12]
    constants[51] = (constants[30]*constants[5])/constants[31]
    constants[52] = (constants[2]*constants[32])/constants[31]
    constants[53] = (constants[33]*constants[5])/constants[31]
    constants[54] = (constants[2]*constants[34])/constants[31]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[2] = (constants[43]*(1.00000-states[2]))/(constants[16]+(1.00000-states[2]))-(constants[44]*states[2])/(constants[17]+states[2])
    rates[4] = (constants[47]*states[2]*(1.00000-states[4]))/(constants[20]+(1.00000-states[4]))-(constants[48]*states[3]*states[4])/(constants[21]+states[4])
    rates[5] = constants[22]*states[4]*constants[9]*(1.00000-states[5])-constants[23]*states[5]
    algebraic[3] = 2.00000*(1.00000-states[0])
    rates[3] = (constants[45]*algebraic[3]*(1.00000-states[3]))/(constants[18]+(1.00000-states[3]))-(constants[46]*states[3])/(constants[19]+states[3])
    rates[6] = (constants[49]*algebraic[3]*(1.00000-states[6]))/(constants[24]+(1.00000-states[6]))-(constants[50]*states[6])/(constants[25]+states[6])
    algebraic[0] = 1.00000-states[0]
    rates[0] = -constants[14]*states[0]*(power(states[1], 2.00000))+constants[15]*(power(algebraic[3], 2.00000))*algebraic[0]*(power(constants[5], 2.00000))*1.00000
    rates[9] = ((-constants[35]*states[9]+constants[36]*states[7])-(constants[53]*algebraic[3]*states[9])/(constants[37]+states[9]))+(constants[54]*states[8])/(constants[38]+states[8])
    rates[7] = ((constants[35]*states[9]-constants[36]*states[7])-(constants[51]*algebraic[3]*states[7])/(constants[39]+states[7]))+(constants[52]*states[10])/(constants[40]+states[10])
    rates[10] = (((constants[51]*algebraic[3]*states[7])/(constants[39]+states[7])-(constants[52]*states[10])/(constants[40]+states[10]))+constants[41]*states[8])-constants[42]*states[10]
    rates[8] = (-constants[41]*states[8]+constants[42]*states[10]+(constants[53]*algebraic[3]*states[9])/(constants[37]+states[9]))-(constants[54]*states[8])/(constants[38]+states[8])
    algebraic[4] = constants[14]*states[0]*(power(states[1], 2.00000))-constants[15]*algebraic[3]*algebraic[0]*(power(constants[5], 2.00000))*1.00000
    rates[1] = (constants[26]*states[5]*constants[28]-(constants[27]*constants[12]*states[6]*states[1])/(constants[29]+states[1]))-2.00000*algebraic[4]*constants[5]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[3] = 2.00000*(1.00000-states[0])
    algebraic[0] = 1.00000-states[0]
    algebraic[4] = constants[14]*states[0]*(power(states[1], 2.00000))-constants[15]*algebraic[3]*algebraic[0]*(power(constants[5], 2.00000))*1.00000
    algebraic[1] = states[9]+states[8]
    algebraic[2] = states[7]+states[10]
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