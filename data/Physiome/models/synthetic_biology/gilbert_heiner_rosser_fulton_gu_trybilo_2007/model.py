# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 4
sizeConstants = 17
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component temporalExistence (second)"
    legend_constants[0] = "alphaTF in component Device_TFgenerator (uM_per_second)"
    legend_constants[14] = "JTF in component Device_TFgenerator (uM_per_second)"
    legend_states[0] = "TFS in component BioEnv_TFSAssociation (uM)"
    legend_algebraic[5] = "Jgain_TF in component BioEnv_TFSAssociation_interface (uM_per_second)"
    legend_constants[15] = "Jgain_TFS in component BioEnv_TFSAssociation_interface (uM_per_second)"
    legend_algebraic[3] = "JTF in component Device_PhzMSgenerator (uM_per_second)"
    legend_states[1] = "TF in component BioEnv_TFSAssociation (uM)"
    legend_constants[1] = "s in component BioEnv_TFSAssociation (uM)"
    legend_constants[2] = "betaTFS in component BioEnv_TFSAssociation (per_uM_per_second)"
    legend_constants[3] = "kd in component BioEnv_TFSAssociation (per_second)"
    legend_constants[4] = "deltaTFS in component BioEnv_TFSAssociation (per_second)"
    legend_constants[5] = "deltaTF in component BioEnv_TFSAssociation (per_second)"
    legend_algebraic[0] = "Jassociation in component BioEnv_TFSAssociation (uM_per_second)"
    legend_constants[6] = "gammaPhzMS in component Device_PhzMSgenerator (uM)"
    legend_constants[7] = "betaPhzMS in component Device_PhzMSgenerator (uM_per_second)"
    legend_algebraic[1] = "JPhzMS in component Device_PhzMSgenerator (uM_per_second)"
    legend_constants[8] = "gammaTF in component Device_PhzMSgenerator (uM)"
    legend_constants[9] = "betaTF in component Device_PhzMSgenerator (uM_per_second)"
    legend_constants[10] = "feedbackOn in component Device_PhzMSgenerator (dimensionless)"
    legend_algebraic[4] = "Jgain_PhzMS in component BioEnv_PhzMStoPYO_interface (uM_per_second)"
    legend_constants[16] = "Jgain_PYO in component BioEnv_PhzMStoPYO_interface (uM_per_second)"
    legend_states[2] = "PhzMS in component BioEnv_PhzMStoPYO (uM)"
    legend_states[3] = "PYO in component BioEnv_PhzMStoPYO (uM)"
    legend_constants[11] = "deltaPhzMS in component BioEnv_PhzMStoPYO (per_second)"
    legend_constants[12] = "alphaPYO in component BioEnv_PhzMStoPYO (per_second)"
    legend_constants[13] = "deltaPYO in component BioEnv_PhzMStoPYO (per_second)"
    legend_algebraic[2] = "JPYOformation in component BioEnv_PhzMStoPYO (uM_per_second)"
    legend_rates[1] = "d/dt TF in component BioEnv_TFSAssociation (uM)"
    legend_rates[0] = "d/dt TFS in component BioEnv_TFSAssociation (uM)"
    legend_rates[2] = "d/dt PhzMS in component BioEnv_PhzMStoPYO (uM)"
    legend_rates[3] = "d/dt PYO in component BioEnv_PhzMStoPYO (uM)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 0.05
    states[0] = 0
    states[1] = 8.6207
    constants[1] = 5
    constants[2] = 1e6
    constants[3] = 4e6
    constants[4] = 3.851e-4
    constants[5] = 5.8e-3
    constants[6] = 5
    constants[7] = 0.1
    constants[8] = 4
    constants[9] = 0.07
    constants[10] = 1
    states[2] = 0
    states[3] = 0
    constants[11] = 8.0225e-6
    constants[12] = 1.3
    constants[13] = 5.8e-1
    constants[14] = constants[0]
    constants[15] = 0.00000
    constants[16] = 0.00000
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = constants[2]*constants[1]*states[1]-constants[3]*states[0]
    rates[0] = (algebraic[0]-constants[4]*states[0])+constants[15]
    algebraic[2] = constants[12]*states[2]
    rates[3] = (constants[16]-constants[13]*states[3])+algebraic[2]
    algebraic[1] = (constants[7]*states[0])/(constants[6]+states[0])
    algebraic[4] = algebraic[1]
    rates[2] = algebraic[4]-constants[11]*states[2]
    algebraic[3] = custom_piecewise([constants[10] != 0.00000, (constants[9]*states[0])/(constants[8]+states[0]) , True, 0.00000])
    algebraic[5] = constants[14]+algebraic[3]
    rates[1] = (-algebraic[0]+algebraic[5])-constants[5]*states[1]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = constants[2]*constants[1]*states[1]-constants[3]*states[0]
    algebraic[2] = constants[12]*states[2]
    algebraic[1] = (constants[7]*states[0])/(constants[6]+states[0])
    algebraic[4] = algebraic[1]
    algebraic[3] = custom_piecewise([constants[10] != 0.00000, (constants[9]*states[0])/(constants[8]+states[0]) , True, 0.00000])
    algebraic[5] = constants[14]+algebraic[3]
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