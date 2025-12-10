# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 9
sizeConstants = 47
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "Cdc13T in component Cdc13T (dimensionless)"
    legend_constants[0] = "k1 in component Cdc13T (first_order_rate_constant)"
    legend_constants[1] = "k2_ in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[2] = "k2__ in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[3] = "k2___ in component kinetic_parameters (first_order_rate_constant)"
    legend_states[1] = "M in component M (dimensionless)"
    legend_states[2] = "Ste9 in component Ste9 (dimensionless)"
    legend_states[3] = "Slp1 in component Slp1 (dimensionless)"
    legend_states[4] = "preMPF in component preMPF (dimensionless)"
    legend_algebraic[3] = "kwee in component preMPF (first_order_rate_constant)"
    legend_constants[4] = "kwee_ in component preMPF (first_order_rate_constant)"
    legend_constants[5] = "kwee__ in component preMPF (first_order_rate_constant)"
    legend_constants[6] = "Vawee in component preMPF (first_order_rate_constant)"
    legend_constants[7] = "Viwee in component preMPF (first_order_rate_constant)"
    legend_constants[8] = "Jawee in component preMPF (dimensionless)"
    legend_constants[9] = "Jiwee in component preMPF (dimensionless)"
    legend_algebraic[5] = "k25 in component preMPF (first_order_rate_constant)"
    legend_constants[10] = "k25_ in component preMPF (first_order_rate_constant)"
    legend_constants[11] = "k25__ in component preMPF (first_order_rate_constant)"
    legend_constants[12] = "Va25 in component preMPF (first_order_rate_constant)"
    legend_constants[13] = "Vi25 in component preMPF (first_order_rate_constant)"
    legend_constants[14] = "Ja25 in component preMPF (dimensionless)"
    legend_constants[15] = "Ji25 in component preMPF (dimensionless)"
    legend_algebraic[2] = "MPF in component MPF (dimensionless)"
    legend_constants[16] = "k3_ in component Ste9 (first_order_rate_constant)"
    legend_constants[17] = "k3__ in component Ste9 (first_order_rate_constant)"
    legend_constants[18] = "k4 in component Ste9 (first_order_rate_constant)"
    legend_constants[19] = "k4_ in component Ste9 (first_order_rate_constant)"
    legend_constants[20] = "J3 in component Ste9 (dimensionless)"
    legend_constants[21] = "J4 in component Ste9 (dimensionless)"
    legend_states[5] = "SK in component SK (dimensionless)"
    legend_states[6] = "Slp1T in component Slp1T (dimensionless)"
    legend_constants[22] = "k5_ in component Slp1T (first_order_rate_constant)"
    legend_constants[23] = "k5__ in component Slp1T (first_order_rate_constant)"
    legend_constants[24] = "J5 in component Slp1T (dimensionless)"
    legend_constants[25] = "k6 in component kinetic_parameters (first_order_rate_constant)"
    legend_constants[26] = "k7 in component Slp1 (first_order_rate_constant)"
    legend_constants[27] = "k8 in component Slp1 (first_order_rate_constant)"
    legend_constants[28] = "J7 in component Slp1 (dimensionless)"
    legend_constants[29] = "J8 in component Slp1 (dimensionless)"
    legend_states[7] = "IEP in component IEP (dimensionless)"
    legend_constants[30] = "k9 in component IEP (first_order_rate_constant)"
    legend_constants[31] = "k10 in component IEP (first_order_rate_constant)"
    legend_constants[32] = "J9 in component IEP (dimensionless)"
    legend_constants[33] = "J10 in component IEP (dimensionless)"
    legend_states[8] = "Rum1T in component Rum1T (dimensionless)"
    legend_constants[34] = "k11 in component Rum1T (first_order_rate_constant)"
    legend_constants[35] = "k12 in component Rum1T (first_order_rate_constant)"
    legend_constants[36] = "k12_ in component Rum1T (first_order_rate_constant)"
    legend_constants[37] = "k12__ in component Rum1T (first_order_rate_constant)"
    legend_constants[38] = "k13 in component SK (first_order_rate_constant)"
    legend_constants[39] = "k14 in component SK (first_order_rate_constant)"
    legend_algebraic[4] = "TF in component TF (dimensionless)"
    legend_constants[40] = "mu in component M (first_order_rate_constant)"
    legend_algebraic[1] = "Trimer in component Trimer (dimensionless)"
    legend_algebraic[0] = "sum in component Trimer (dimensionless)"
    legend_constants[41] = "Kdiss in component Trimer (dimensionless)"
    legend_constants[42] = "k15 in component TF (first_order_rate_constant)"
    legend_constants[43] = "k16_ in component TF (first_order_rate_constant)"
    legend_constants[44] = "k16__ in component TF (first_order_rate_constant)"
    legend_constants[45] = "J15 in component TF (dimensionless)"
    legend_constants[46] = "J16 in component TF (dimensionless)"
    legend_rates[0] = "d/dt Cdc13T in component Cdc13T (dimensionless)"
    legend_rates[4] = "d/dt preMPF in component preMPF (dimensionless)"
    legend_rates[2] = "d/dt Ste9 in component Ste9 (dimensionless)"
    legend_rates[6] = "d/dt Slp1T in component Slp1T (dimensionless)"
    legend_rates[3] = "d/dt Slp1 in component Slp1 (dimensionless)"
    legend_rates[7] = "d/dt IEP in component IEP (dimensionless)"
    legend_rates[8] = "d/dt Rum1T in component Rum1T (dimensionless)"
    legend_rates[5] = "d/dt SK in component SK (dimensionless)"
    legend_rates[1] = "d/dt M in component M (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.2
    constants[0] = 0.03
    constants[1] = 0.03
    constants[2] = 1
    constants[3] = 0.1
    states[1] = 1
    states[2] = 1
    states[3] = 2.2
    states[4] = 0
    constants[4] = 0.15
    constants[5] = 1.3
    constants[6] = 0.25
    constants[7] = 1
    constants[8] = 0.01
    constants[9] = 0.01
    constants[10] = 0.05
    constants[11] = 5
    constants[12] = 1
    constants[13] = 0.25
    constants[14] = 0.01
    constants[15] = 0.01
    constants[16] = 1
    constants[17] = 10
    constants[18] = 35
    constants[19] = 2
    constants[20] = 0.01
    constants[21] = 0.01
    states[5] = 0
    states[6] = 0
    constants[22] = 0.005
    constants[23] = 0.3
    constants[24] = 0.3
    constants[25] = 0.1
    constants[26] = 1
    constants[27] = 0.25
    constants[28] = 0.001
    constants[29] = 0.001
    states[7] = 0
    constants[30] = 0.1
    constants[31] = 0.04
    constants[32] = 0.01
    constants[33] = 0.01
    states[8] = 0
    constants[34] = 0.1
    constants[35] = 0.01
    constants[36] = 1
    constants[37] = 3
    constants[38] = 0.1
    constants[39] = 0.1
    constants[40] = 0.005
    constants[41] = 0.001
    constants[42] = 1.5
    constants[43] = 1
    constants[44] = 2
    constants[45] = 0.01
    constants[46] = 0.01
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = constants[0]*states[1]-(constants[1]+constants[2]*states[2]+constants[3]*states[3])*states[0]
    rates[3] = (constants[26]*states[7]*(states[6]-states[3]))/((constants[28]+states[6])-states[3])-((constants[27]*states[3])/(constants[29]+states[3])+constants[25]*states[3])
    rates[1] = constants[40]*states[1]
    algebraic[0] = states[0]+states[8]+constants[41]
    algebraic[1] = (2.00000*states[0]*states[8])/(algebraic[0]+power(power(algebraic[0], 2.00000)-4.00000*states[0]*states[8], 1.0/2))
    algebraic[2] = ((states[0]-states[4])*(states[0]-algebraic[1]))/states[0]
    rates[2] = ((constants[16]+constants[17]*states[3])*(1.00000-states[2]))/((constants[20]+1.00000)-states[2])-((constants[19]*states[5]+constants[18]*algebraic[2])*states[2])/(constants[21]+states[2])
    rates[6] = (constants[22]+(constants[23]*(power(algebraic[2], 4.00000)))/(power(constants[24], 4.00000)+power(algebraic[2], 4.00000)))-constants[25]*states[6]
    rates[7] = (constants[30]*algebraic[2]*(1.00000-states[7]))/((constants[32]+1.00000)-states[7])-(constants[31]*states[7])/(constants[33]+states[7])
    rates[8] = constants[34]-(constants[35]+constants[36]*states[5]+constants[37]*algebraic[2])*states[8]
    algebraic[4] = (2.00000*constants[42]*states[1]*constants[46])/(((constants[43]+constants[44]*algebraic[2])-constants[42]*states[1])+(constants[43]+constants[44]*algebraic[2])*constants[45]+constants[42]*states[1]*constants[46]+power(power(((constants[43]+constants[44]*algebraic[2])-constants[42]*states[1])+(constants[43]+constants[44]*algebraic[2])*constants[45]+constants[42]*states[1]*constants[46], 2.00000)-4.00000*constants[42]*states[1]*constants[46]*((constants[43]+constants[44]*algebraic[2])-constants[42]*states[1]), 1.0/2))
    rates[5] = constants[38]*algebraic[4]-constants[39]*states[5]
    algebraic[3] = constants[4]+((constants[5]-constants[4])*2.00000*constants[6]*constants[9])/((constants[7]*algebraic[2]-constants[6])+constants[7]*algebraic[2]*constants[8]+constants[6]*constants[9]+power(power((constants[7]*algebraic[2]-constants[6])+constants[7]*algebraic[2]*constants[8]+constants[6]*constants[9], 2.00000)-4.00000*constants[6]*constants[9]*(constants[7]*algebraic[2]-constants[6]), 1.0/2))
    algebraic[5] = constants[10]+((constants[11]-constants[10])*2.00000*constants[12]*algebraic[2]*constants[15])/((constants[13]-constants[12]*algebraic[2])+constants[13]*constants[14]+constants[12]*algebraic[2]*constants[15]+power(power((constants[13]-constants[12]*algebraic[2])+constants[13]*constants[14]+constants[12]*algebraic[2]*constants[15], 2.00000)-4.00000*constants[12]*algebraic[2]*constants[15]*(constants[13]-constants[12]*algebraic[2]), 1.0/2))
    rates[4] = (algebraic[3]*(states[0]-states[4])-algebraic[5]*states[4])-(constants[1]+constants[2]*states[2]+constants[3]*states[3])*states[4]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = states[0]+states[8]+constants[41]
    algebraic[1] = (2.00000*states[0]*states[8])/(algebraic[0]+power(power(algebraic[0], 2.00000)-4.00000*states[0]*states[8], 1.0/2))
    algebraic[2] = ((states[0]-states[4])*(states[0]-algebraic[1]))/states[0]
    algebraic[4] = (2.00000*constants[42]*states[1]*constants[46])/(((constants[43]+constants[44]*algebraic[2])-constants[42]*states[1])+(constants[43]+constants[44]*algebraic[2])*constants[45]+constants[42]*states[1]*constants[46]+power(power(((constants[43]+constants[44]*algebraic[2])-constants[42]*states[1])+(constants[43]+constants[44]*algebraic[2])*constants[45]+constants[42]*states[1]*constants[46], 2.00000)-4.00000*constants[42]*states[1]*constants[46]*((constants[43]+constants[44]*algebraic[2])-constants[42]*states[1]), 1.0/2))
    algebraic[3] = constants[4]+((constants[5]-constants[4])*2.00000*constants[6]*constants[9])/((constants[7]*algebraic[2]-constants[6])+constants[7]*algebraic[2]*constants[8]+constants[6]*constants[9]+power(power((constants[7]*algebraic[2]-constants[6])+constants[7]*algebraic[2]*constants[8]+constants[6]*constants[9], 2.00000)-4.00000*constants[6]*constants[9]*(constants[7]*algebraic[2]-constants[6]), 1.0/2))
    algebraic[5] = constants[10]+((constants[11]-constants[10])*2.00000*constants[12]*algebraic[2]*constants[15])/((constants[13]-constants[12]*algebraic[2])+constants[13]*constants[14]+constants[12]*algebraic[2]*constants[15]+power(power((constants[13]-constants[12]*algebraic[2])+constants[13]*constants[14]+constants[12]*algebraic[2]*constants[15], 2.00000)-4.00000*constants[12]*algebraic[2]*constants[15]*(constants[13]-constants[12]*algebraic[2]), 1.0/2))
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