# Size of variable arrays:
sizeAlgebraic = 22
sizeStates = 10
sizeConstants = 43
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_constants[0] = "R_mt in component heart_parameters (kPa_second_per_mL)"
    legend_constants[1] = "R_av in component heart_parameters (kPa_second_per_mL)"
    legend_constants[2] = "R_tc in component heart_parameters (kPa_second_per_mL)"
    legend_constants[3] = "R_pv in component heart_parameters (kPa_second_per_mL)"
    legend_constants[4] = "R_pul in component heart_parameters (kPa_second_per_mL)"
    legend_constants[5] = "R_sys in component heart_parameters (kPa_second_per_mL)"
    legend_constants[6] = "L_tc in component heart_parameters (kPa_second2_per_mL)"
    legend_constants[7] = "L_pv in component heart_parameters (kPa_second2_per_mL)"
    legend_constants[8] = "L_mt in component heart_parameters (kPa_second2_per_mL)"
    legend_constants[9] = "L_av in component heart_parameters (kPa_second2_per_mL)"
    legend_constants[10] = "V_tot in component heart_parameters (mL)"
    legend_constants[11] = "P_th in component heart_parameters (kPa)"
    legend_algebraic[1] = "e_t in component driver_function (dimensionless)"
    legend_constants[12] = "A in component driver_function (dimensionless)"
    legend_constants[13] = "B in component driver_function (per_second2)"
    legend_constants[14] = "C in component driver_function (second)"
    legend_algebraic[0] = "tau in component driver_function (second)"
    legend_constants[15] = "period in component driver_function (second)"
    legend_algebraic[2] = "V_pcd in component pericardium (mL)"
    legend_algebraic[3] = "P_pcd in component pericardium (kPa)"
    legend_algebraic[4] = "P_peri in component pericardium (kPa)"
    legend_states[0] = "V_lv in component left_ventricle (mL)"
    legend_states[1] = "V_rv in component right_ventricle (mL)"
    legend_constants[16] = "P_0_pcd in component pericardium (kPa)"
    legend_constants[17] = "V_0_pcd in component pericardium (mL)"
    legend_constants[18] = "lambda_pcd in component pericardium (per_mL)"
    legend_algebraic[6] = "V_lvf in component left_ventricle (mL)"
    legend_algebraic[9] = "P_lvf in component left_ventricle (kPa)"
    legend_algebraic[10] = "P_lv in component left_ventricle (kPa)"
    legend_algebraic[5] = "V_spt in component septum (mL)"
    legend_algebraic[7] = "P_es_lvf in component lvf_calculator (kPa)"
    legend_algebraic[8] = "P_ed_lvf in component lvf_calculator (kPa)"
    legend_algebraic[18] = "P_pu in component pulmonary_vein (kPa)"
    legend_algebraic[17] = "P_ao in component aorta (kPa)"
    legend_constants[19] = "E_es_lvf in component lvf_calculator (kPa_per_mL)"
    legend_constants[20] = "lambda_lvf in component lvf_calculator (per_mL)"
    legend_constants[21] = "P_0_lvf in component lvf_calculator (kPa)"
    legend_states[2] = "Q_mt in component flow (mL_per_second)"
    legend_states[3] = "Q_av in component flow (mL_per_second)"
    legend_constants[22] = "V_d_lvf in component lvf_calculator (mL)"
    legend_constants[23] = "V_0_lvf in component lvf_calculator (mL)"
    legend_algebraic[11] = "V_rvf in component right_ventricle (mL)"
    legend_algebraic[14] = "P_rvf in component right_ventricle (kPa)"
    legend_algebraic[15] = "P_rv in component right_ventricle (kPa)"
    legend_algebraic[12] = "P_es_rvf in component rvf_calculator (kPa)"
    legend_algebraic[13] = "P_ed_rvf in component rvf_calculator (kPa)"
    legend_algebraic[16] = "P_pa in component pulmonary_artery (kPa)"
    legend_algebraic[19] = "P_vc in component vena_cava (kPa)"
    legend_constants[24] = "E_es_rvf in component rvf_calculator (kPa_per_mL)"
    legend_constants[25] = "lambda_rvf in component rvf_calculator (per_mL)"
    legend_constants[26] = "P_0_rvf in component rvf_calculator (kPa)"
    legend_states[4] = "Q_tc in component flow (mL_per_second)"
    legend_states[5] = "Q_pv in component flow (mL_per_second)"
    legend_constants[27] = "V_d_rvf in component rvf_calculator (mL)"
    legend_constants[28] = "V_0_rvf in component rvf_calculator (mL)"
    legend_constants[29] = "E_es_spt in component septum (kPa_per_mL)"
    legend_constants[30] = "V_d_spt in component septum (mL)"
    legend_constants[31] = "P_0_spt in component septum (kPa)"
    legend_constants[32] = "lambda_spt in component septum (per_mL)"
    legend_constants[33] = "V_0_spt in component septum (mL)"
    legend_constants[34] = "one in component septum (dimensionless)"
    legend_constants[35] = "E_es_pa in component pulmonary_artery (kPa_per_mL)"
    legend_states[6] = "V_pa in component pulmonary_artery (mL)"
    legend_constants[36] = "V_d_pa in component pulmonary_artery (mL)"
    legend_algebraic[20] = "Q_pul in component flow (mL_per_second)"
    legend_constants[37] = "E_es_pu in component pulmonary_vein (kPa_per_mL)"
    legend_states[7] = "V_pu in component pulmonary_vein (mL)"
    legend_constants[38] = "V_d_pu in component pulmonary_vein (mL)"
    legend_constants[39] = "E_es_ao in component aorta (kPa_per_mL)"
    legend_states[8] = "V_ao in component aorta (mL)"
    legend_constants[40] = "V_d_ao in component aorta (mL)"
    legend_algebraic[21] = "Q_sys in component flow (mL_per_second)"
    legend_constants[41] = "E_es_vc in component vena_cava (kPa_per_mL)"
    legend_states[9] = "V_vc in component vena_cava (mL)"
    legend_constants[42] = "V_d_vc in component vena_cava (mL)"
    legend_rates[0] = "d/dt V_lv in component left_ventricle (mL)"
    legend_rates[1] = "d/dt V_rv in component right_ventricle (mL)"
    legend_rates[6] = "d/dt V_pa in component pulmonary_artery (mL)"
    legend_rates[7] = "d/dt V_pu in component pulmonary_vein (mL)"
    legend_rates[8] = "d/dt V_ao in component aorta (mL)"
    legend_rates[9] = "d/dt V_vc in component vena_cava (mL)"
    legend_rates[2] = "d/dt Q_mt in component flow (mL_per_second)"
    legend_rates[3] = "d/dt Q_av in component flow (mL_per_second)"
    legend_rates[4] = "d/dt Q_tc in component flow (mL_per_second)"
    legend_rates[5] = "d/dt Q_pv in component flow (mL_per_second)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 0.0158
    constants[1] = 0.0180
    constants[2] = 0.0237
    constants[3] = 0.0055
    constants[4] = 0.1552
    constants[5] = 1.0889
    constants[6] = 8.0093e-5
    constants[7] = 1.4868e-4
    constants[8] = 7.6968e-5
    constants[9] = 1.2189e-4
    constants[10] = 5.5
    constants[11] = -4
    constants[12] = 1
    constants[13] = 80
    constants[14] = 0.375
    constants[15] = 0.75
    states[0] = 94.6812
    states[1] = 90.7302
    constants[16] = 0.5003
    constants[17] = 200
    constants[18] = 0.03
    constants[19] = 2.8798
    constants[20] = 0.033
    constants[21] = 0.1203
    states[2] = 245.5813
    states[3] = 0
    constants[22] = 0
    constants[23] = 0
    constants[24] = 0.585
    constants[25] = 0.023
    constants[26] = 0.2157
    states[4] = 190.0661
    states[5] = 0
    constants[27] = 0
    constants[28] = 0
    constants[29] = 48.754
    constants[30] = 2
    constants[31] = 1.1101
    constants[32] = 0.435
    constants[33] = 2
    constants[34] = 1
    constants[35] = 0.369
    states[6] = 43.0123
    constants[36] = 0
    constants[37] = 0.0073
    states[7] = 808.4579
    constants[38] = 0
    constants[39] = 0.6913
    states[8] = 133.3381
    constants[40] = 0
    constants[41] = 0.0059
    states[9] = 329.7803
    constants[42] = 0
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = custom_piecewise([less(states[2] , 0.00000) & less(states[3] , 0.00000), 0.00000 , less(states[2] , 0.00000), -states[3] , less(states[3] , 0.00000), states[2] , True, states[2]-states[3]])
    rates[1] = custom_piecewise([less(states[4] , 0.00000) & less(states[5] , 0.00000), 0.00000 , less(states[4] , 0.00000), -states[5] , less(states[5] , 0.00000), states[4] , True, states[4]-states[5]])
    algebraic[2] = states[0]+states[1]
    algebraic[3] = constants[16]*(exp(constants[18]*(algebraic[2]-constants[17]))-1.00000)
    algebraic[4] = algebraic[3]+constants[11]
    algebraic[0] = custom_piecewise([less_equal(voi , constants[15]), voi , less_equal(voi , constants[15]*2.00000), voi-constants[15] , less_equal(voi , constants[15]*3.00000), voi-constants[15]*2.00000 , less_equal(voi , constants[15]*4.00000), voi-constants[15]*3.00000 , less_equal(voi , constants[15]*5.00000), voi-constants[15]*4.00000 , less_equal(voi , constants[15]*6.00000), voi-constants[15]*5.00000 , less_equal(voi , constants[15]*7.00000), voi-constants[15]*6.00000 , less_equal(voi , constants[15]*8.00000), voi-constants[15]*7.00000 , less_equal(voi , constants[15]*9.00000), voi-constants[15]*8.00000 , less_equal(voi , constants[15]*10.0000), voi-constants[15]*9.00000 , less_equal(voi , constants[15]*11.0000), voi-constants[15]*10.0000 , less_equal(voi , constants[15]*12.0000), voi-constants[15]*11.0000 , less_equal(voi , constants[15]*13.0000), voi-constants[15]*12.0000 , True, float('nan')])
    algebraic[1] = constants[12]*exp(-constants[13]*(power(algebraic[0]-constants[14], 2.00000)))
    rootfind_0(voi, constants, rates, states, algebraic)
    algebraic[6] = states[0]-algebraic[5]
    algebraic[7] = constants[19]*(algebraic[6]-constants[22])
    algebraic[8] = constants[21]*(exp(constants[20]*(algebraic[6]-constants[23]))-1.00000)
    algebraic[9] = algebraic[1]*algebraic[7]+(1.00000-algebraic[1])*algebraic[8]
    algebraic[10] = algebraic[9]+algebraic[4]
    algebraic[17] = constants[39]*(states[8]-constants[40])
    rates[3] = custom_piecewise([less(algebraic[10]-algebraic[17] , 0.00000) & less(states[3] , 0.00000), 0.00000 , True, ((algebraic[10]-algebraic[17])-states[3]*constants[1])/constants[9]])
    algebraic[11] = states[1]+algebraic[5]
    algebraic[12] = constants[24]*(algebraic[11]-constants[27])
    algebraic[13] = constants[26]*(exp(constants[25]*(algebraic[11]-constants[28]))-1.00000)
    algebraic[14] = algebraic[1]*algebraic[12]+(1.00000-algebraic[1])*algebraic[13]
    algebraic[15] = algebraic[14]+algebraic[4]
    algebraic[16] = constants[35]*(states[6]-constants[36])+constants[11]
    rates[5] = custom_piecewise([less(algebraic[15]-algebraic[16] , 0.00000) & less(states[5] , 0.00000), 0.00000 , True, ((algebraic[15]-algebraic[16])-states[5]*constants[3])/constants[7]])
    algebraic[18] = constants[37]*(states[7]-constants[38])+constants[11]
    rates[2] = custom_piecewise([less(algebraic[18]-algebraic[10] , 0.00000) & less(states[2] , 0.00000), 0.00000 , True, ((algebraic[18]-algebraic[10])-states[2]*constants[0])/constants[8]])
    algebraic[19] = constants[41]*(states[9]-constants[42])
    rates[4] = custom_piecewise([less(algebraic[19]-algebraic[15] , 0.00000) & less(states[4] , 0.00000), 0.00000 , True, ((algebraic[19]-algebraic[15])-states[4]*constants[2])/constants[6]])
    algebraic[20] = (algebraic[16]-algebraic[18])/constants[4]
    rates[6] = custom_piecewise([less(states[5] , 0.00000), -algebraic[20] , True, states[5]-algebraic[20]])
    rates[7] = custom_piecewise([less(states[2] , 0.00000), algebraic[20] , True, algebraic[20]-states[2]])
    algebraic[21] = (algebraic[17]-algebraic[19])/constants[5]
    rates[8] = custom_piecewise([less(states[3] , 0.00000), -algebraic[21] , True, states[3]-algebraic[21]])
    rates[9] = custom_piecewise([less(states[4] , 0.00000), algebraic[21] , True, algebraic[21]-states[4]])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[2] = states[0]+states[1]
    algebraic[3] = constants[16]*(exp(constants[18]*(algebraic[2]-constants[17]))-1.00000)
    algebraic[4] = algebraic[3]+constants[11]
    algebraic[0] = custom_piecewise([less_equal(voi , constants[15]), voi , less_equal(voi , constants[15]*2.00000), voi-constants[15] , less_equal(voi , constants[15]*3.00000), voi-constants[15]*2.00000 , less_equal(voi , constants[15]*4.00000), voi-constants[15]*3.00000 , less_equal(voi , constants[15]*5.00000), voi-constants[15]*4.00000 , less_equal(voi , constants[15]*6.00000), voi-constants[15]*5.00000 , less_equal(voi , constants[15]*7.00000), voi-constants[15]*6.00000 , less_equal(voi , constants[15]*8.00000), voi-constants[15]*7.00000 , less_equal(voi , constants[15]*9.00000), voi-constants[15]*8.00000 , less_equal(voi , constants[15]*10.0000), voi-constants[15]*9.00000 , less_equal(voi , constants[15]*11.0000), voi-constants[15]*10.0000 , less_equal(voi , constants[15]*12.0000), voi-constants[15]*11.0000 , less_equal(voi , constants[15]*13.0000), voi-constants[15]*12.0000 , True, float('nan')])
    algebraic[1] = constants[12]*exp(-constants[13]*(power(algebraic[0]-constants[14], 2.00000)))
    algebraic[6] = states[0]-algebraic[5]
    algebraic[7] = constants[19]*(algebraic[6]-constants[22])
    algebraic[8] = constants[21]*(exp(constants[20]*(algebraic[6]-constants[23]))-1.00000)
    algebraic[9] = algebraic[1]*algebraic[7]+(1.00000-algebraic[1])*algebraic[8]
    algebraic[10] = algebraic[9]+algebraic[4]
    algebraic[17] = constants[39]*(states[8]-constants[40])
    algebraic[11] = states[1]+algebraic[5]
    algebraic[12] = constants[24]*(algebraic[11]-constants[27])
    algebraic[13] = constants[26]*(exp(constants[25]*(algebraic[11]-constants[28]))-1.00000)
    algebraic[14] = algebraic[1]*algebraic[12]+(1.00000-algebraic[1])*algebraic[13]
    algebraic[15] = algebraic[14]+algebraic[4]
    algebraic[16] = constants[35]*(states[6]-constants[36])+constants[11]
    algebraic[18] = constants[37]*(states[7]-constants[38])+constants[11]
    algebraic[19] = constants[41]*(states[9]-constants[42])
    algebraic[20] = (algebraic[16]-algebraic[18])/constants[4]
    algebraic[21] = (algebraic[17]-algebraic[19])/constants[5]
    return algebraic

initialGuess0 = None
def rootfind_0(voi, constants, states, algebraic):
    """Calculate value of algebraic variable for DAE"""
    from scipy.optimize import fsolve
    global initialGuess0
    if initialGuess0 is None: initialGuess0 = 0.1
    if not iterable(voi):
        algebraic[5] = fsolve(residualSN_0, initialGuess0, args=(algebraic, voi, constants, rates, states), xtol=1E-6)
        initialGuess0 = algebraic[5]
    else:
        for (i,t) in enumerate(voi):
            algebraic[5][i] = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], voi[i], constants, rates, states[:,i]), xtol=1E-6)
            initialGuess0 = algebraic[5][i]

def residualSN_0(algebraicCandidate, algebraic, voi, constants, rates, states):
    algebraic[5] = algebraicCandidate
    return (0.00000) - ((((algebraic[1]*constants[29]*(algebraic[5]-constants[30])+(constants[34]-algebraic[1])*constants[31]*(exp(constants[32]*(algebraic[5]-constants[33]))-constants[34]))-algebraic[1]*constants[19]*(states[0]-algebraic[5]))-(1.00000-algebraic[1])*constants[21]*(exp(constants[20]*(states[0]-algebraic[5]))-1.00000))+algebraic[1]*constants[24]*(states[1]+algebraic[5])+(1.00000-algebraic[1])*constants[26]*(exp(constants[25]*(states[1]+algebraic[5]))-1.00000))

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