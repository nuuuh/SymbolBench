# Size of variable arrays:
sizeAlgebraic = 6
sizeStates = 2
sizeConstants = 15
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (minute)"
    legend_states[0] = "G in component G (mg_per_litre)"
    legend_algebraic[0] = "Gin in component G (mg_per_litre_minute)"
    legend_algebraic[2] = "f2_G in component f2_G (mg_per_minute)"
    legend_algebraic[3] = "f3_G in component f3_G (per_litre)"
    legend_algebraic[4] = "f4_I in component f4_I (mg_per_minute)"
    legend_algebraic[5] = "f5_I in component f5_I (mg_per_minute)"
    legend_states[1] = "I in component I (mU_per_litre)"
    legend_algebraic[1] = "Iin in component I (mU_per_litre_minute)"
    legend_constants[0] = "di in component model_parameters (first_order_rate_constant)"
    legend_constants[1] = "C2 in component model_parameters (mg_per_litre)"
    legend_constants[2] = "Ub in component model_parameters (mg_per_minute)"
    legend_constants[3] = "Vg in component model_parameters (litre)"
    legend_constants[4] = "C3 in component model_parameters (mg_per_litre)"
    legend_constants[5] = "C4 in component model_parameters (mU_per_litre)"
    legend_constants[6] = "Vi in component model_parameters (litre)"
    legend_constants[7] = "U0 in component model_parameters (mg_per_minute)"
    legend_constants[8] = "Um in component model_parameters (mg_per_minute)"
    legend_constants[9] = "beta in component model_parameters (dimensionless)"
    legend_constants[10] = "ti in component model_parameters (minute)"
    legend_constants[11] = "Vp in component model_parameters (litre)"
    legend_constants[12] = "C5 in component model_parameters (mU_per_litre)"
    legend_constants[13] = "Rg in component model_parameters (mg_per_minute)"
    legend_constants[14] = "alpha in component model_parameters (litre_per_mU)"
    legend_rates[0] = "d/dt G in component G (mg_per_litre)"
    legend_rates[1] = "d/dt I in component I (mU_per_litre)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 150.0
    states[1] = 5.0
    constants[0] = 0.0107
    constants[1] = 144
    constants[2] = 72.0
    constants[3] = 10.0
    constants[4] = 1000.0
    constants[5] = 80.0
    constants[6] = 11.0
    constants[7] = 40.0
    constants[8] = 940.0
    constants[9] = 1.77
    constants[10] = 100.0
    constants[11] = 3.0
    constants[12] = 26.0
    constants[13] = 180.0
    constants[14] = 0.29
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = custom_piecewise([less_equal(voi , 30.0000) & greater_equal(voi , 0.00000), 0.250000 , less(voi , 120.000) & greater_equal(voi , 30.0000), 0.250000+1.00000*(1.00000+(voi-120.000)/90.0000) , less(voi , 240.000) & greater_equal(voi , 120.000), 0.250000+1.00000*(1.00000-0.500000*((voi-120.000)/120.000)) , less_equal(voi , 480.000) & greater_equal(voi , 240.000), 0.250000+0.500000*(1.00000-(voi-240.000)/240.000) , less(voi , 530.000) & greater(voi , 480.000), 0.250000 , less(voi , 620.000) & greater_equal(voi , 530.000), 0.250000+1.00000*(1.00000+(voi-620.000)/90.0000) , less(voi , 740.000) & greater_equal(voi , 620.000), 0.250000+1.00000*(1.00000-0.500000*((voi-620.000)/120.000)) , less_equal(voi , 980.000) & greater_equal(voi , 740.000), 0.250000+0.500000*(1.00000-(voi-740.000)/240.000) , True, float('nan')])
    rates[1] = algebraic[1]-constants[0]*states[1]
    algebraic[0] = custom_piecewise([less(voi , 15.0000) & greater_equal(voi , 0.00000), 0.0500000+(5.00000/15.0000)*voi , less(voi , 45.0000) & greater_equal(voi , 15.0000), 0.0500000+5.00000*((45.0000-voi)/(45.0000-15.0000)) , less_equal(voi , 240.000) & greater_equal(voi , 45.0000), 0.0500000 , greater(voi , 240.000) & less(voi , 255.000), 0.0500000+(5.00000/255.000)*voi , less(voi , 285.000) & greater_equal(voi , 255.000), 0.0500000+5.00000*((285.000-voi)/(285.000-255.000)) , less_equal(voi , 480.000) & greater_equal(voi , 285.000), 0.0500000 , less(voi , 495.000) & greater(voi , 480.000), 0.0500000+(5.00000/495.000)*voi , less(voi , 525.000) & greater_equal(voi , 495.000), 0.0500000+5.00000*((525.000-voi)/(525.000-495.000)) , less_equal(voi , 720.000) & greater_equal(voi , 525.000), 0.0500000 , greater(voi , 720.000) & less(voi , 735.000), 0.0500000+(5.00000/735.000)*voi , less(voi , 765.000) & greater_equal(voi , 735.000), 0.0500000+5.00000*((765.000-voi)/(765.000-735.000)) , less_equal(voi , 960.000) & greater_equal(voi , 765.000), 0.0500000 , True, float('nan')])
    algebraic[2] = constants[2]*(1.00000-exp(-states[0]/(constants[1]*constants[3]*1.00000)))
    algebraic[3] = states[0]/(constants[4]*constants[3])
    algebraic[4] = constants[7]+(constants[8]-constants[7])/(1.00000+exp(-constants[9]*log(states[1]/(constants[5]*(1.00000/constants[6]+1.00000/(0.200000*constants[10]))))))
    algebraic[5] = constants[13]/(1.00000+exp(constants[14]*(states[1]/(1.00000*(1.00000*constants[11]-constants[12])))))
    rates[0] = (algebraic[0]+1.00000*algebraic[5])-(1.00000*algebraic[2]+algebraic[3]*algebraic[4])
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = custom_piecewise([less_equal(voi , 30.0000) & greater_equal(voi , 0.00000), 0.250000 , less(voi , 120.000) & greater_equal(voi , 30.0000), 0.250000+1.00000*(1.00000+(voi-120.000)/90.0000) , less(voi , 240.000) & greater_equal(voi , 120.000), 0.250000+1.00000*(1.00000-0.500000*((voi-120.000)/120.000)) , less_equal(voi , 480.000) & greater_equal(voi , 240.000), 0.250000+0.500000*(1.00000-(voi-240.000)/240.000) , less(voi , 530.000) & greater(voi , 480.000), 0.250000 , less(voi , 620.000) & greater_equal(voi , 530.000), 0.250000+1.00000*(1.00000+(voi-620.000)/90.0000) , less(voi , 740.000) & greater_equal(voi , 620.000), 0.250000+1.00000*(1.00000-0.500000*((voi-620.000)/120.000)) , less_equal(voi , 980.000) & greater_equal(voi , 740.000), 0.250000+0.500000*(1.00000-(voi-740.000)/240.000) , True, float('nan')])
    algebraic[0] = custom_piecewise([less(voi , 15.0000) & greater_equal(voi , 0.00000), 0.0500000+(5.00000/15.0000)*voi , less(voi , 45.0000) & greater_equal(voi , 15.0000), 0.0500000+5.00000*((45.0000-voi)/(45.0000-15.0000)) , less_equal(voi , 240.000) & greater_equal(voi , 45.0000), 0.0500000 , greater(voi , 240.000) & less(voi , 255.000), 0.0500000+(5.00000/255.000)*voi , less(voi , 285.000) & greater_equal(voi , 255.000), 0.0500000+5.00000*((285.000-voi)/(285.000-255.000)) , less_equal(voi , 480.000) & greater_equal(voi , 285.000), 0.0500000 , less(voi , 495.000) & greater(voi , 480.000), 0.0500000+(5.00000/495.000)*voi , less(voi , 525.000) & greater_equal(voi , 495.000), 0.0500000+5.00000*((525.000-voi)/(525.000-495.000)) , less_equal(voi , 720.000) & greater_equal(voi , 525.000), 0.0500000 , greater(voi , 720.000) & less(voi , 735.000), 0.0500000+(5.00000/735.000)*voi , less(voi , 765.000) & greater_equal(voi , 735.000), 0.0500000+5.00000*((765.000-voi)/(765.000-735.000)) , less_equal(voi , 960.000) & greater_equal(voi , 765.000), 0.0500000 , True, float('nan')])
    algebraic[2] = constants[2]*(1.00000-exp(-states[0]/(constants[1]*constants[3]*1.00000)))
    algebraic[3] = states[0]/(constants[4]*constants[3])
    algebraic[4] = constants[7]+(constants[8]-constants[7])/(1.00000+exp(-constants[9]*log(states[1]/(constants[5]*(1.00000/constants[6]+1.00000/(0.200000*constants[10]))))))
    algebraic[5] = constants[13]/(1.00000+exp(constants[14]*(states[1]/(1.00000*(1.00000*constants[11]-constants[12])))))
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