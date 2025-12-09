# Size of variable arrays:
sizeAlgebraic = 24
sizeStates = 7
sizeConstants = 51
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "t in component environment (second)"
    legend_algebraic[11] = "Pi in component TempCDa (UnitP)"
    legend_states[0] = "Pi in component TempRLC (UnitP)"
    legend_algebraic[20] = "Qo in component TempRC (UnitQ)"
    legend_algebraic[18] = "Qvi in component TempDPump (UnitQ)"
    legend_algebraic[6] = "Qo in component TempCDv (UnitQ)"
    legend_algebraic[7] = "Pi in component TJoint (UnitP)"
    legend_algebraic[17] = "Qvo in component TempDPump (UnitQ)"
    legend_algebraic[0] = "Pi in component TJoint (UnitP)"
    legend_algebraic[4] = "Pi in component TempCDv (UnitP)"
    legend_algebraic[21] = "Qo in component TJoint (UnitQ)"
    legend_constants[0] = "CVao in component ParaHeartP (UnitCV)"
    legend_algebraic[3] = "E in component EVentricle (UnitE)"
    legend_states[1] = "V in component TempCDv (UnitV)"
    legend_constants[1] = "PlvIni in component ParaHeartP (UnitP)"
    legend_constants[2] = "VlvIni in component ParaHeartP (UnitV)"
    legend_algebraic[5] = "Tao in component TempCDv (dimensionless)"
    legend_constants[3] = "Vlv0 in component ParaHeartP (UnitV)"
    legend_algebraic[13] = "Qo in component TempCDa (UnitQ)"
    legend_constants[4] = "CVmi in component ParaHeartP (UnitCV)"
    legend_algebraic[10] = "E in component EAtrium (UnitE)"
    legend_states[2] = "V in component TempCDa (UnitV)"
    legend_constants[5] = "PlaIni in component ParaHeartP (UnitP)"
    legend_constants[6] = "VlaIni in component ParaHeartP (UnitV)"
    legend_algebraic[12] = "Tao in component TempCDa (dimensionless)"
    legend_constants[7] = "Vla0 in component ParaHeartP (UnitV)"
    legend_constants[8] = "ElvMax in component ParaHeartP (UnitE)"
    legend_constants[9] = "ElvMin in component ParaHeartP (UnitE)"
    legend_constants[10] = "T in component ParaHeartP (second)"
    legend_constants[11] = "Ts1 in component ParaHeartP (dimensionless)"
    legend_constants[12] = "Ts2 in component ParaHeartP (dimensionless)"
    legend_algebraic[1] = "mt in component EVentricle (second)"
    legend_algebraic[2] = "et in component EVentricle (dimensionless)"
    legend_constants[13] = "ElaMax in component ParaHeartP (UnitE)"
    legend_constants[14] = "ElaMin in component ParaHeartP (UnitE)"
    legend_constants[15] = "Tpwb in component ParaHeartP (dimensionless)"
    legend_constants[16] = "Tpww in component ParaHeartP (dimensionless)"
    legend_algebraic[8] = "mt in component EAtrium (second)"
    legend_algebraic[9] = "et in component EAtrium (dimensionless)"
    legend_states[3] = "Qo in component TempRLC (UnitQ)"
    legend_constants[17] = "Rsas in component ParaSys (UnitR)"
    legend_constants[18] = "Csas in component ParaSys (UnitC)"
    legend_constants[19] = "Lsas in component ParaSys (UnitL)"
    legend_constants[20] = "P0sas in component ParaSys (UnitP)"
    legend_constants[21] = "Q0sas in component ParaSys (UnitQ)"
    legend_states[4] = "Pi in component TempRLC (UnitP)"
    legend_algebraic[22] = "Qo in component TJoint (UnitQ)"
    legend_algebraic[23] = "Pi in component TempR (UnitP)"
    legend_states[5] = "Qo in component TempRLC (UnitQ)"
    legend_constants[22] = "Rsat in component ParaSys (UnitR)"
    legend_constants[23] = "Csat in component ParaSys (UnitC)"
    legend_constants[24] = "Lsat in component ParaSys (UnitL)"
    legend_constants[25] = "P0sat in component ParaSys (UnitP)"
    legend_constants[26] = "Q0sat in component ParaSys (UnitQ)"
    legend_algebraic[19] = "Pi in component TempR (UnitP)"
    legend_algebraic[14] = "Qo in component TempR (UnitQ)"
    legend_constants[27] = "Rsar in component ParaSys (UnitR)"
    legend_states[6] = "Pi in component TempRC (UnitP)"
    legend_algebraic[16] = "Qo in component TempR (UnitQ)"
    legend_constants[28] = "Rscp in component ParaSys (UnitR)"
    legend_constants[29] = "Rsvn in component ParaSys (UnitR)"
    legend_constants[30] = "Csvn in component ParaSys (UnitC)"
    legend_constants[31] = "P0svn in component ParaSys (UnitP)"
    legend_constants[32] = "Qai in component ParaDPump (UnitQ)"
    legend_constants[33] = "Qao in component ParaDPump (UnitQ)"
    legend_algebraic[15] = "mt in component TempDPump (second)"
    legend_constants[34] = "Ts2 in component ParaHeartP (dimensionless)"
    legend_constants[35] = "T in component ParaHeartP (second)"
    legend_constants[36] = "ElaMax in component ParaHeartP (UnitE)"
    legend_constants[37] = "ElaMin in component ParaHeartP (UnitE)"
    legend_constants[38] = "PlaIni in component ParaHeartP (UnitP)"
    legend_constants[39] = "VlaIni in component ParaHeartP (UnitV)"
    legend_constants[40] = "ElvMax in component ParaHeartP (UnitE)"
    legend_constants[41] = "ElvMin in component ParaHeartP (UnitE)"
    legend_constants[42] = "PlvIni in component ParaHeartP (UnitP)"
    legend_constants[43] = "VlvIni in component ParaHeartP (UnitV)"
    legend_constants[44] = "Tpwb in component ParaHeartP (dimensionless)"
    legend_constants[45] = "Tpww in component ParaHeartP (dimensionless)"
    legend_constants[46] = "Ts1 in component ParaHeartP (dimensionless)"
    legend_constants[47] = "CVao in component ParaHeartP (UnitCV)"
    legend_constants[48] = "CVmi in component ParaHeartP (UnitCV)"
    legend_constants[49] = "Vlv0 in component ParaHeartP (UnitV)"
    legend_constants[50] = "Vla0 in component ParaHeartP (UnitV)"
    legend_rates[1] = "d/dt V in component TempCDv (UnitV)"
    legend_rates[2] = "d/dt V in component TempCDa (UnitV)"
    legend_rates[0] = "d/dt Pi in component TempRLC (UnitP)"
    legend_rates[3] = "d/dt Qo in component TempRLC (UnitQ)"
    legend_rates[4] = "d/dt Pi in component TempRLC (UnitP)"
    legend_rates[5] = "d/dt Qo in component TempRLC (UnitQ)"
    legend_rates[6] = "d/dt Pi in component TempRC (UnitP)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 350.
    constants[1] = 1.0
    constants[2] = 5.0
    constants[3] = 500
    constants[4] = 400.
    constants[5] = 1.0
    constants[6] = 4.0
    constants[7] = 20
    constants[8] = 0.5
    constants[9] = 0.1
    constants[10] = 1.0
    constants[11] = 0.3
    constants[12] = 0.45
    constants[13] = 0.25
    constants[14] = 0.15
    constants[15] = 0.92
    constants[16] = 0.09
    constants[17] = 0.003
    constants[18] = 0.08
    constants[19] = 0.000062
    constants[20] = 100.
    constants[21] = 0.
    constants[22] = 0.05
    constants[23] = 1.6
    constants[24] = 0.0017
    constants[25] = 100.
    constants[26] = 0.
    constants[27] = 0.5
    constants[28] = 0.52
    constants[29] = 0.075
    constants[30] = 20.5
    constants[31] = 0.
    constants[32] = -300
    constants[33] = 200
    constants[34] = 0.45
    constants[35] = 1.0
    constants[36] = 0.25
    constants[37] = 0.15
    constants[38] = 1.0
    constants[39] = 4.0
    constants[40] = 0.5
    constants[41] = 0.1
    constants[42] = 1.0
    constants[43] = 5.0
    constants[44] = 0.92
    constants[45] = 0.09
    constants[46] = 0.3
    constants[47] = 350.
    constants[48] = 400.
    constants[49] = 500
    constants[50] = 20
    states[0] = constants[20]
    states[1] = constants[3]
    states[2] = constants[7]
    states[3] = constants[21]
    states[4] = constants[25]
    states[5] = constants[26]
    states[6] = constants[31]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = states[4]
    rates[3] = ((states[0]-algebraic[0])-constants[17]*states[3])/constants[19]
    algebraic[1] = voi-constants[10]*floor(voi/constants[10])
    algebraic[2] = custom_piecewise([greater_equal(algebraic[1] , 0.00000) & less_equal(algebraic[1] , constants[11]*constants[10]), 1.00000-cos((3.14159*algebraic[1])/(constants[11]*constants[10])) , greater(algebraic[1] , constants[11]*constants[10]) & less_equal(algebraic[1] , constants[12]*constants[10]), 1.00000+cos((3.14159*(algebraic[1]-constants[11]*constants[10]))/((constants[12]-constants[11])*constants[10])) , greater(algebraic[1] , constants[12]*constants[10]) & less(algebraic[1] , constants[10]), 0.00000 , True, float('nan')])
    algebraic[3] = constants[9]+(algebraic[2]*(constants[8]-constants[9]))/2.00000
    algebraic[4] = constants[1]+algebraic[3]*(states[1]-constants[2])
    algebraic[5] = custom_piecewise([greater_equal(algebraic[4] , states[0]), 1.00000 , less(algebraic[4] , states[0]), 0.00000 , True, float('nan')])
    algebraic[6] = custom_piecewise([greater_equal(algebraic[4] , states[0]), constants[0]*algebraic[5]*(power(fabs(algebraic[4]-states[0]), 0.500000)) , less(algebraic[4] , states[0]), -1.00000*constants[0]*algebraic[5]*(power(fabs(states[0]-algebraic[4]), 0.500000)) , True, float('nan')])
    rates[0] = (algebraic[6]-states[3])/constants[18]
    algebraic[15] = voi-constants[35]*floor(voi/constants[35])
    algebraic[18] = custom_piecewise([greater_equal(algebraic[15] , 0.00000) & less_equal(algebraic[15] , constants[34]*constants[35]), 0.500000*constants[32]*(1.00000-cos((2.00000*3.14159*algebraic[15])/(constants[34]*constants[35]))) , greater(algebraic[15] , constants[34]*constants[35]) & less(algebraic[15] , constants[35]), 0.00000 , True, float('nan')])
    algebraic[8] = voi-constants[10]*floor(voi/constants[10])
    algebraic[9] = custom_piecewise([greater_equal(algebraic[8] , 0.00000) & less_equal(algebraic[8] , ((constants[15]+constants[16])-1.00000)*constants[10]), 1.00000-cos((2.00000*3.14159*(algebraic[8]-(constants[15]-1.00000)*constants[10]))/(constants[16]*constants[10])) , greater(algebraic[8] , ((constants[15]+constants[16])-1.00000)*constants[10]) & less_equal(algebraic[8] , constants[15]*constants[10]), 0.00000 , greater(algebraic[8] , constants[15]*constants[10]) & less_equal(algebraic[8] , constants[10]), 1.00000-cos((2.00000*3.14159*(algebraic[8]-constants[15]*constants[10]))/(constants[16]*constants[10])) , True, float('nan')])
    algebraic[10] = constants[14]+(algebraic[9]*(constants[13]-constants[14]))/2.00000
    algebraic[11] = constants[5]+algebraic[10]*(states[2]-constants[6])
    algebraic[7] = algebraic[4]
    algebraic[12] = custom_piecewise([greater_equal(algebraic[11] , algebraic[7]), 1.00000 , less(algebraic[11] , algebraic[7]), 0.00000 , True, float('nan')])
    algebraic[13] = custom_piecewise([greater_equal(algebraic[11] , algebraic[7]), constants[4]*algebraic[12]*(power(fabs(algebraic[11]-algebraic[7]), 0.500000)) , less(algebraic[11] , algebraic[7]), -1.00000*constants[4]*algebraic[12]*(power(fabs(algebraic[7]-algebraic[11]), 0.500000)) , True, float('nan')])
    algebraic[21] = algebraic[13]+algebraic[18]
    rates[1] = algebraic[21]-algebraic[6]
    algebraic[20] = (states[6]-algebraic[11])/constants[29]
    rates[2] = algebraic[20]-algebraic[13]
    algebraic[17] = custom_piecewise([greater_equal(algebraic[15] , 0.00000) & less_equal(algebraic[15] , constants[34]*constants[35]), 0.00000 , greater(algebraic[15] , constants[34]*constants[35]) & less(algebraic[15] , constants[35]), 0.500000*constants[33]*(1.00000-cos((2.00000*3.14159*(algebraic[15]-constants[34]*constants[35]))/(constants[35]-constants[34]*constants[35]))) , True, float('nan')])
    algebraic[22] = states[3]+algebraic[17]
    rates[4] = (algebraic[22]-states[5])/constants[23]
    algebraic[14] = states[5]
    algebraic[19] = states[6]+constants[28]*algebraic[14]
    algebraic[23] = algebraic[19]+constants[27]*states[5]
    rates[5] = ((states[4]-algebraic[23])-constants[22]*states[5])/constants[24]
    algebraic[16] = algebraic[14]
    rates[6] = (algebraic[16]-algebraic[20])/constants[30]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[4]
    algebraic[1] = voi-constants[10]*floor(voi/constants[10])
    algebraic[2] = custom_piecewise([greater_equal(algebraic[1] , 0.00000) & less_equal(algebraic[1] , constants[11]*constants[10]), 1.00000-cos((3.14159*algebraic[1])/(constants[11]*constants[10])) , greater(algebraic[1] , constants[11]*constants[10]) & less_equal(algebraic[1] , constants[12]*constants[10]), 1.00000+cos((3.14159*(algebraic[1]-constants[11]*constants[10]))/((constants[12]-constants[11])*constants[10])) , greater(algebraic[1] , constants[12]*constants[10]) & less(algebraic[1] , constants[10]), 0.00000 , True, float('nan')])
    algebraic[3] = constants[9]+(algebraic[2]*(constants[8]-constants[9]))/2.00000
    algebraic[4] = constants[1]+algebraic[3]*(states[1]-constants[2])
    algebraic[5] = custom_piecewise([greater_equal(algebraic[4] , states[0]), 1.00000 , less(algebraic[4] , states[0]), 0.00000 , True, float('nan')])
    algebraic[6] = custom_piecewise([greater_equal(algebraic[4] , states[0]), constants[0]*algebraic[5]*(power(fabs(algebraic[4]-states[0]), 0.500000)) , less(algebraic[4] , states[0]), -1.00000*constants[0]*algebraic[5]*(power(fabs(states[0]-algebraic[4]), 0.500000)) , True, float('nan')])
    algebraic[15] = voi-constants[35]*floor(voi/constants[35])
    algebraic[18] = custom_piecewise([greater_equal(algebraic[15] , 0.00000) & less_equal(algebraic[15] , constants[34]*constants[35]), 0.500000*constants[32]*(1.00000-cos((2.00000*3.14159*algebraic[15])/(constants[34]*constants[35]))) , greater(algebraic[15] , constants[34]*constants[35]) & less(algebraic[15] , constants[35]), 0.00000 , True, float('nan')])
    algebraic[8] = voi-constants[10]*floor(voi/constants[10])
    algebraic[9] = custom_piecewise([greater_equal(algebraic[8] , 0.00000) & less_equal(algebraic[8] , ((constants[15]+constants[16])-1.00000)*constants[10]), 1.00000-cos((2.00000*3.14159*(algebraic[8]-(constants[15]-1.00000)*constants[10]))/(constants[16]*constants[10])) , greater(algebraic[8] , ((constants[15]+constants[16])-1.00000)*constants[10]) & less_equal(algebraic[8] , constants[15]*constants[10]), 0.00000 , greater(algebraic[8] , constants[15]*constants[10]) & less_equal(algebraic[8] , constants[10]), 1.00000-cos((2.00000*3.14159*(algebraic[8]-constants[15]*constants[10]))/(constants[16]*constants[10])) , True, float('nan')])
    algebraic[10] = constants[14]+(algebraic[9]*(constants[13]-constants[14]))/2.00000
    algebraic[11] = constants[5]+algebraic[10]*(states[2]-constants[6])
    algebraic[7] = algebraic[4]
    algebraic[12] = custom_piecewise([greater_equal(algebraic[11] , algebraic[7]), 1.00000 , less(algebraic[11] , algebraic[7]), 0.00000 , True, float('nan')])
    algebraic[13] = custom_piecewise([greater_equal(algebraic[11] , algebraic[7]), constants[4]*algebraic[12]*(power(fabs(algebraic[11]-algebraic[7]), 0.500000)) , less(algebraic[11] , algebraic[7]), -1.00000*constants[4]*algebraic[12]*(power(fabs(algebraic[7]-algebraic[11]), 0.500000)) , True, float('nan')])
    algebraic[21] = algebraic[13]+algebraic[18]
    algebraic[20] = (states[6]-algebraic[11])/constants[29]
    algebraic[17] = custom_piecewise([greater_equal(algebraic[15] , 0.00000) & less_equal(algebraic[15] , constants[34]*constants[35]), 0.00000 , greater(algebraic[15] , constants[34]*constants[35]) & less(algebraic[15] , constants[35]), 0.500000*constants[33]*(1.00000-cos((2.00000*3.14159*(algebraic[15]-constants[34]*constants[35]))/(constants[35]-constants[34]*constants[35]))) , True, float('nan')])
    algebraic[22] = states[3]+algebraic[17]
    algebraic[14] = states[5]
    algebraic[19] = states[6]+constants[28]*algebraic[14]
    algebraic[23] = algebraic[19]+constants[27]*states[5]
    algebraic[16] = algebraic[14]
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