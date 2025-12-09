# Size of variable arrays:
sizeAlgebraic = 70
sizeStates = 7
sizeConstants = 43
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_states[0] = "V in component membrane (millivolt)"
    legend_algebraic[1] = "IK in component membrane (femtoA)"
    legend_algebraic[4] = "ICa in component membrane (femtoA)"
    legend_algebraic[2] = "IKCa in component membrane (femtoA)"
    legend_algebraic[29] = "IKATP in component KATP (femtoA)"
    legend_constants[0] = "Cm in component membrane (femtoF)"
    legend_constants[1] = "gK in component membrane (picoS)"
    legend_constants[2] = "VK in component KATP (millivolt)"
    legend_states[1] = "n in component membrane (dimensionless)"
    legend_constants[3] = "gKCa in component membrane (picoS)"
    legend_constants[4] = "kd in component membrane (micromolar)"
    legend_states[2] = "c in component calcium_handling (micromolar)"
    legend_constants[5] = "gCa in component membrane (picoS)"
    legend_algebraic[3] = "minf in component membrane (dimensionless)"
    legend_constants[6] = "VCa in component membrane (millivolt)"
    legend_constants[7] = "taun in component membrane (millisecond)"
    legend_algebraic[0] = "ninf in component membrane (dimensionless)"
    legend_states[3] = "cer in component calcium_handling (micromolar)"
    legend_constants[8] = "fcyt in component calcium_handling (dimensionless)"
    legend_algebraic[11] = "Jmem in component calcium_handling (flux)"
    legend_algebraic[9] = "Jer in component calcium_handling (flux)"
    legend_constants[9] = "fer in component calcium_handling (dimensionless)"
    legend_constants[10] = "sigmaV in component calcium_handling (dimensionless)"
    legend_constants[11] = "pleak in component calcium_handling (first_order_rate_constant)"
    legend_constants[12] = "Kserca in component calcium_handling (first_order_rate_constant)"
    legend_constants[13] = "lambdaer in component calcium_handling (dimensionless)"
    legend_constants[14] = "epser in component calcium_handling (dimensionless)"
    legend_constants[15] = "alpha in component calcium_handling (micromolar_per_femtoA_millisecond)"
    legend_constants[16] = "kpmca in component calcium_handling (first_order_rate_constant)"
    legend_algebraic[5] = "Jserca in component calcium_handling (flux)"
    legend_algebraic[7] = "Jleak in component calcium_handling (flux)"
    legend_algebraic[8] = "rgpdh in component glycolysis (flux)"
    legend_constants[17] = "Rgk in component glycolysis (per_second)"
    legend_constants[18] = "atot in component glycolysis (micromolar)"
    legend_constants[19] = "pfkbas in component glycolysis (dimensionless)"
    legend_algebraic[6] = "f6p in component glycolysis (micromolar)"
    legend_constants[20] = "lambda in component glycolysis (dimensionless)"
    legend_algebraic[69] = "pfk in component pfk (micromolar)"
    legend_states[4] = "g6p in component glycolysis (micromolar)"
    legend_states[5] = "fbp in component glycolysis (micromolar)"
    legend_constants[21] = "bottom1 in component pfk (dimensionless)"
    legend_constants[22] = "topa1 in component pfk (dimensionless)"
    legend_constants[23] = "k1 in component pfk (micromolar)"
    legend_constants[24] = "k2 in component pfk (micromolar)"
    legend_constants[25] = "k3 in component pfk (micromolar)"
    legend_constants[26] = "k4 in component pfk (micromolar)"
    legend_constants[27] = "cat in component pfk (dimensionless)"
    legend_algebraic[21] = "atp in component nucleotides (micromolar)"
    legend_algebraic[22] = "weight2 in component pfk (dimensionless)"
    legend_constants[42] = "topa2 in component pfk (dimensionless)"
    legend_algebraic[25] = "bottom2 in component pfk (dimensionless)"
    legend_algebraic[12] = "topa3 in component pfk (dimensionless)"
    legend_algebraic[10] = "weight3 in component pfk (dimensionless)"
    legend_algebraic[28] = "bottom3 in component pfk (dimensionless)"
    legend_constants[28] = "famp in component pfk (dimensionless)"
    legend_constants[29] = "fatp in component pfk (dimensionless)"
    legend_constants[30] = "ffbp in component pfk (dimensionless)"
    legend_constants[31] = "fbt in component pfk (dimensionless)"
    legend_constants[32] = "fmt in component pfk (dimensionless)"
    legend_algebraic[30] = "weight4 in component pfk (dimensionless)"
    legend_algebraic[31] = "topa4 in component pfk (dimensionless)"
    legend_algebraic[32] = "bottom4 in component pfk (dimensionless)"
    legend_algebraic[13] = "weight5 in component pfk (dimensionless)"
    legend_algebraic[33] = "topa5 in component pfk (dimensionless)"
    legend_algebraic[34] = "bottom5 in component pfk (dimensionless)"
    legend_algebraic[35] = "weight6 in component pfk (dimensionless)"
    legend_algebraic[36] = "topa6 in component pfk (dimensionless)"
    legend_algebraic[37] = "bottom6 in component pfk (dimensionless)"
    legend_algebraic[14] = "weight7 in component pfk (dimensionless)"
    legend_algebraic[38] = "topa7 in component pfk (dimensionless)"
    legend_algebraic[39] = "bottom7 in component pfk (dimensionless)"
    legend_algebraic[40] = "weight8 in component pfk (dimensionless)"
    legend_algebraic[41] = "topa8 in component pfk (dimensionless)"
    legend_algebraic[42] = "bottom8 in component pfk (dimensionless)"
    legend_algebraic[46] = "weight9 in component pfk (dimensionless)"
    legend_algebraic[43] = "topa9 in component pfk (dimensionless)"
    legend_algebraic[47] = "bottom9 in component pfk (dimensionless)"
    legend_algebraic[48] = "weight10 in component pfk (dimensionless)"
    legend_algebraic[44] = "topa10 in component pfk (dimensionless)"
    legend_algebraic[49] = "bottom10 in component pfk (dimensionless)"
    legend_algebraic[50] = "weight11 in component pfk (dimensionless)"
    legend_algebraic[51] = "topa11 in component pfk (dimensionless)"
    legend_algebraic[52] = "bottom11 in component pfk (dimensionless)"
    legend_algebraic[53] = "weight12 in component pfk (dimensionless)"
    legend_algebraic[54] = "topa12 in component pfk (dimensionless)"
    legend_algebraic[55] = "bottom12 in component pfk (dimensionless)"
    legend_algebraic[56] = "weight13 in component pfk (dimensionless)"
    legend_algebraic[57] = "topa13 in component pfk (dimensionless)"
    legend_algebraic[58] = "bottom13 in component pfk (dimensionless)"
    legend_algebraic[59] = "weight14 in component pfk (dimensionless)"
    legend_algebraic[60] = "topa14 in component pfk (dimensionless)"
    legend_algebraic[61] = "bottom14 in component pfk (dimensionless)"
    legend_algebraic[62] = "weight15 in component pfk (dimensionless)"
    legend_algebraic[63] = "topa15 in component pfk (dimensionless)"
    legend_algebraic[64] = "bottom15 in component pfk (dimensionless)"
    legend_algebraic[66] = "weight16 in component pfk (dimensionless)"
    legend_algebraic[67] = "topa16 in component pfk (dimensionless)"
    legend_algebraic[68] = "bottom16 in component pfk (dimensionless)"
    legend_algebraic[65] = "topb in component pfk (dimensionless)"
    legend_algebraic[45] = "amp in component nucleotides (micromolar)"
    legend_algebraic[15] = "mgadp in component KATP (micromolar)"
    legend_algebraic[16] = "adp3m in component KATP (micromolar)"
    legend_algebraic[23] = "atp4m in component KATP (micromolar)"
    legend_algebraic[17] = "topo in component KATP (dimensionless)"
    legend_algebraic[26] = "bottomo in component KATP (dimensionless)"
    legend_algebraic[27] = "katpo in component KATP (dimensionless)"
    legend_constants[33] = "gkatpbar in component KATP (picoS)"
    legend_states[6] = "adp in component nucleotides (micromolar)"
    legend_constants[34] = "kdd in component KATP (dimensionless)"
    legend_constants[35] = "ktd in component KATP (dimensionless)"
    legend_constants[36] = "ktt in component KATP (dimensionless)"
    legend_algebraic[19] = "fback in component nucleotides (dimensionless)"
    legend_constants[37] = "taua in component nucleotides (dimensionless)"
    legend_constants[38] = "r1 in component nucleotides (micromolar)"
    legend_constants[39] = "r in component nucleotides (dimensionless)"
    legend_algebraic[18] = "y in component nucleotides (dimensionless)"
    legend_constants[40] = "vg in component nucleotides (dimensionless)"
    legend_constants[41] = "kg in component nucleotides (flux)"
    legend_algebraic[20] = "rad in component nucleotides (dimensionless)"
    legend_algebraic[24] = "ratio in component nucleotides (dimensionless)"
    legend_rates[0] = "d/dt V in component membrane (millivolt)"
    legend_rates[1] = "d/dt n in component membrane (dimensionless)"
    legend_rates[2] = "d/dt c in component calcium_handling (micromolar)"
    legend_rates[3] = "d/dt cer in component calcium_handling (micromolar)"
    legend_rates[5] = "d/dt fbp in component glycolysis (micromolar)"
    legend_rates[4] = "d/dt g6p in component glycolysis (micromolar)"
    legend_rates[6] = "d/dt adp in component nucleotides (micromolar)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = -60
    constants[0] = 5300
    constants[1] = 2700
    constants[2] = -75
    states[1] = 0
    constants[3] = 600
    constants[4] = 0.5
    states[2] = 0.25
    constants[5] = 1000
    constants[6] = 25
    constants[7] = 20
    states[3] = 185
    constants[8] = 0.01
    constants[9] = 0.01
    constants[10] = 31
    constants[11] = 0.0002
    constants[12] = 0.4
    constants[13] = 1
    constants[14] = 1
    constants[15] = 0.00000450
    constants[16] = 0.2
    constants[17] = 0.2
    constants[18] = 3000
    constants[19] = 0.06
    constants[20] = 0.005
    states[4] = 200
    states[5] = 40
    constants[21] = 1
    constants[22] = 0
    constants[23] = 30
    constants[24] = 1
    constants[25] = 50000
    constants[26] = 1000
    constants[27] = 2
    constants[28] = 0.02
    constants[29] = 20
    constants[30] = 0.2
    constants[31] = 20
    constants[32] = 20
    constants[33] = 25000
    states[6] = 780
    constants[34] = 17
    constants[35] = 26
    constants[36] = 1
    constants[37] = 300000
    constants[38] = 0.35
    constants[39] = 1
    constants[40] = 2.2
    constants[41] = 10
    constants[42] = constants[22]
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = 1.00000/(1.00000+exp(-(16.0000+states[0]/1.00000)/5.00000))
    rates[1] = (algebraic[0]-states[1])/constants[7]
    algebraic[5] = constants[12]*states[2]
    algebraic[7] = constants[11]*(states[3]-states[2])
    algebraic[9] = (constants[14]*(algebraic[7]-algebraic[5]))/constants[13]
    rates[3] = -constants[9]*constants[10]*algebraic[9]
    algebraic[3] = 1.00000/(1.00000+exp(-(20.0000+states[0]/1.00000)/12.0000))
    algebraic[4] = constants[5]*algebraic[3]*(states[0]-constants[6])
    algebraic[11] = -(constants[15]*algebraic[4]+constants[16]*states[2])
    rates[2] = constants[8]*(algebraic[11]+algebraic[9])
    algebraic[20] = (power(fabs(power(states[6]-constants[18], 2.00000)-4.00000*(power(states[6], 2.00000))), 1.0/2))/1.00000
    algebraic[21] = 0.500000*((constants[18]-states[6])+algebraic[20]*1.00000)
    algebraic[8] = 0.200000*(power(fabs((states[5]*1.00000)/(power(1.00000, 2.00000))), 1.0/2))
    algebraic[18] = constants[40]*(algebraic[8]/(constants[41]+algebraic[8]))
    algebraic[19] = constants[39]+algebraic[18]
    rates[6] = (algebraic[21]-states[6]*exp(algebraic[19]*(1.00000-states[2]/constants[38])))/(constants[37]*1.00000)
    algebraic[1] = constants[1]*states[1]*(states[0]-constants[2])
    algebraic[2] = (constants[3]/(1.00000+power(constants[4]/states[2], 2.00000)))*(states[0]-constants[2])
    algebraic[15] = 0.165000*states[6]
    algebraic[17] = 0.0800000*(1.00000+(2.00000*algebraic[15])/(constants[34]*1.00000))+0.890000*(power(algebraic[15]/(constants[34]*1.00000), 2.00000))
    algebraic[16] = 0.135000*states[6]
    algebraic[23] = 0.0500000*algebraic[21]
    algebraic[26] = (power(1.00000+algebraic[15]/(constants[34]*1.00000), 2.00000))*(1.00000+algebraic[16]/(constants[35]*1.00000)+algebraic[23]/(constants[36]*1.00000))
    algebraic[27] = algebraic[17]/algebraic[26]
    algebraic[29] = constants[33]*algebraic[27]*(states[0]-constants[2])
    rates[0] = -(algebraic[1]+algebraic[4]+algebraic[2]+algebraic[29])/constants[0]
    algebraic[6] = 0.300000*states[4]
    algebraic[10] = (power(algebraic[6], 2.00000))/(constants[25]*1.00000)
    algebraic[12] = constants[42]+algebraic[10]
    algebraic[30] = (power(algebraic[6]*algebraic[21], 2.00000))/(constants[29]*constants[25]*constants[26]*(power(1.00000, 2.00000)))
    algebraic[31] = algebraic[12]+algebraic[30]
    algebraic[33] = algebraic[31]
    algebraic[36] = algebraic[33]
    algebraic[14] = (states[5]*(power(algebraic[6], 2.00000)))/(constants[24]*constants[25]*constants[30]*1.00000)
    algebraic[38] = algebraic[36]+algebraic[14]
    algebraic[40] = (states[5]*(power(algebraic[6], 2.00000))*(power(algebraic[21], 2.00000)))/(constants[24]*constants[25]*constants[26]*constants[30]*constants[31]*constants[29]*(power(1.00000, 2.00000)))
    algebraic[41] = algebraic[38]+algebraic[40]
    algebraic[43] = algebraic[41]
    algebraic[44] = algebraic[43]
    algebraic[45] = (states[6]*states[6])/algebraic[21]
    algebraic[50] = (algebraic[45]*(power(algebraic[6], 2.00000)))/(constants[23]*constants[25]*constants[28]*1.00000)
    algebraic[51] = algebraic[44]+algebraic[50]
    algebraic[53] = (algebraic[45]*(power(algebraic[6], 2.00000))*(power(algebraic[21], 2.00000)))/(constants[23]*constants[25]*constants[26]*constants[28]*constants[32]*constants[29]*(power(1.00000, 2.00000)))
    algebraic[54] = algebraic[51]+algebraic[53]
    algebraic[57] = algebraic[54]
    algebraic[60] = algebraic[57]
    algebraic[63] = algebraic[60]
    algebraic[66] = (algebraic[45]*states[5]*(power(algebraic[6], 2.00000))*(power(algebraic[21], 2.00000)))/(constants[23]*constants[24]*constants[25]*constants[26]*constants[30]*constants[28]*constants[31]*constants[32]*constants[29]*(power(1.00000, 2.00000)))
    algebraic[67] = algebraic[63]+algebraic[66]
    algebraic[22] = (power(algebraic[21], 2.00000))/(constants[26]*1.00000)
    algebraic[25] = constants[21]+algebraic[22]
    algebraic[28] = algebraic[25]+algebraic[10]
    algebraic[32] = algebraic[28]+algebraic[30]
    algebraic[13] = states[5]/constants[24]
    algebraic[34] = algebraic[32]+algebraic[13]
    algebraic[35] = (states[5]*(power(algebraic[21], 2.00000)))/(constants[24]*constants[26]*constants[31]*1.00000)
    algebraic[37] = algebraic[34]+algebraic[35]
    algebraic[39] = algebraic[37]+algebraic[14]
    algebraic[42] = algebraic[39]+algebraic[40]
    algebraic[46] = algebraic[45]/constants[23]
    algebraic[47] = algebraic[42]+algebraic[46]
    algebraic[48] = (algebraic[45]*(power(algebraic[21], 2.00000)))/(constants[23]*constants[26]*constants[32]*1.00000)
    algebraic[49] = algebraic[47]+algebraic[48]
    algebraic[52] = algebraic[49]+algebraic[50]
    algebraic[55] = algebraic[52]+algebraic[53]
    algebraic[56] = (algebraic[45]*states[5])/(constants[23]*constants[24])
    algebraic[58] = algebraic[55]+algebraic[56]
    algebraic[59] = (algebraic[45]*states[5]*(power(algebraic[21], 2.00000)))/(constants[23]*constants[24]*constants[26]*constants[31]*constants[32]*1.00000)
    algebraic[61] = algebraic[58]+algebraic[59]
    algebraic[62] = (algebraic[45]*states[5]*(power(algebraic[6], 2.00000)))/(constants[23]*constants[24]*constants[25]*constants[30]*constants[28]*1.00000)
    algebraic[64] = algebraic[61]+algebraic[62]
    algebraic[68] = algebraic[64]+algebraic[66]
    algebraic[65] = algebraic[62]
    algebraic[69] = 1.00000*((constants[19]*constants[27]*algebraic[67]+constants[27]*algebraic[65])/algebraic[68])
    rates[5] = constants[20]*(algebraic[69]/1.00000-0.500000*algebraic[8])
    rates[4] = constants[20]*(constants[17]*1.00000-algebraic[69]/1.00000)
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = 1.00000/(1.00000+exp(-(16.0000+states[0]/1.00000)/5.00000))
    algebraic[5] = constants[12]*states[2]
    algebraic[7] = constants[11]*(states[3]-states[2])
    algebraic[9] = (constants[14]*(algebraic[7]-algebraic[5]))/constants[13]
    algebraic[3] = 1.00000/(1.00000+exp(-(20.0000+states[0]/1.00000)/12.0000))
    algebraic[4] = constants[5]*algebraic[3]*(states[0]-constants[6])
    algebraic[11] = -(constants[15]*algebraic[4]+constants[16]*states[2])
    algebraic[20] = (power(fabs(power(states[6]-constants[18], 2.00000)-4.00000*(power(states[6], 2.00000))), 1.0/2))/1.00000
    algebraic[21] = 0.500000*((constants[18]-states[6])+algebraic[20]*1.00000)
    algebraic[8] = 0.200000*(power(fabs((states[5]*1.00000)/(power(1.00000, 2.00000))), 1.0/2))
    algebraic[18] = constants[40]*(algebraic[8]/(constants[41]+algebraic[8]))
    algebraic[19] = constants[39]+algebraic[18]
    algebraic[1] = constants[1]*states[1]*(states[0]-constants[2])
    algebraic[2] = (constants[3]/(1.00000+power(constants[4]/states[2], 2.00000)))*(states[0]-constants[2])
    algebraic[15] = 0.165000*states[6]
    algebraic[17] = 0.0800000*(1.00000+(2.00000*algebraic[15])/(constants[34]*1.00000))+0.890000*(power(algebraic[15]/(constants[34]*1.00000), 2.00000))
    algebraic[16] = 0.135000*states[6]
    algebraic[23] = 0.0500000*algebraic[21]
    algebraic[26] = (power(1.00000+algebraic[15]/(constants[34]*1.00000), 2.00000))*(1.00000+algebraic[16]/(constants[35]*1.00000)+algebraic[23]/(constants[36]*1.00000))
    algebraic[27] = algebraic[17]/algebraic[26]
    algebraic[29] = constants[33]*algebraic[27]*(states[0]-constants[2])
    algebraic[6] = 0.300000*states[4]
    algebraic[10] = (power(algebraic[6], 2.00000))/(constants[25]*1.00000)
    algebraic[12] = constants[42]+algebraic[10]
    algebraic[30] = (power(algebraic[6]*algebraic[21], 2.00000))/(constants[29]*constants[25]*constants[26]*(power(1.00000, 2.00000)))
    algebraic[31] = algebraic[12]+algebraic[30]
    algebraic[33] = algebraic[31]
    algebraic[36] = algebraic[33]
    algebraic[14] = (states[5]*(power(algebraic[6], 2.00000)))/(constants[24]*constants[25]*constants[30]*1.00000)
    algebraic[38] = algebraic[36]+algebraic[14]
    algebraic[40] = (states[5]*(power(algebraic[6], 2.00000))*(power(algebraic[21], 2.00000)))/(constants[24]*constants[25]*constants[26]*constants[30]*constants[31]*constants[29]*(power(1.00000, 2.00000)))
    algebraic[41] = algebraic[38]+algebraic[40]
    algebraic[43] = algebraic[41]
    algebraic[44] = algebraic[43]
    algebraic[45] = (states[6]*states[6])/algebraic[21]
    algebraic[50] = (algebraic[45]*(power(algebraic[6], 2.00000)))/(constants[23]*constants[25]*constants[28]*1.00000)
    algebraic[51] = algebraic[44]+algebraic[50]
    algebraic[53] = (algebraic[45]*(power(algebraic[6], 2.00000))*(power(algebraic[21], 2.00000)))/(constants[23]*constants[25]*constants[26]*constants[28]*constants[32]*constants[29]*(power(1.00000, 2.00000)))
    algebraic[54] = algebraic[51]+algebraic[53]
    algebraic[57] = algebraic[54]
    algebraic[60] = algebraic[57]
    algebraic[63] = algebraic[60]
    algebraic[66] = (algebraic[45]*states[5]*(power(algebraic[6], 2.00000))*(power(algebraic[21], 2.00000)))/(constants[23]*constants[24]*constants[25]*constants[26]*constants[30]*constants[28]*constants[31]*constants[32]*constants[29]*(power(1.00000, 2.00000)))
    algebraic[67] = algebraic[63]+algebraic[66]
    algebraic[22] = (power(algebraic[21], 2.00000))/(constants[26]*1.00000)
    algebraic[25] = constants[21]+algebraic[22]
    algebraic[28] = algebraic[25]+algebraic[10]
    algebraic[32] = algebraic[28]+algebraic[30]
    algebraic[13] = states[5]/constants[24]
    algebraic[34] = algebraic[32]+algebraic[13]
    algebraic[35] = (states[5]*(power(algebraic[21], 2.00000)))/(constants[24]*constants[26]*constants[31]*1.00000)
    algebraic[37] = algebraic[34]+algebraic[35]
    algebraic[39] = algebraic[37]+algebraic[14]
    algebraic[42] = algebraic[39]+algebraic[40]
    algebraic[46] = algebraic[45]/constants[23]
    algebraic[47] = algebraic[42]+algebraic[46]
    algebraic[48] = (algebraic[45]*(power(algebraic[21], 2.00000)))/(constants[23]*constants[26]*constants[32]*1.00000)
    algebraic[49] = algebraic[47]+algebraic[48]
    algebraic[52] = algebraic[49]+algebraic[50]
    algebraic[55] = algebraic[52]+algebraic[53]
    algebraic[56] = (algebraic[45]*states[5])/(constants[23]*constants[24])
    algebraic[58] = algebraic[55]+algebraic[56]
    algebraic[59] = (algebraic[45]*states[5]*(power(algebraic[21], 2.00000)))/(constants[23]*constants[24]*constants[26]*constants[31]*constants[32]*1.00000)
    algebraic[61] = algebraic[58]+algebraic[59]
    algebraic[62] = (algebraic[45]*states[5]*(power(algebraic[6], 2.00000)))/(constants[23]*constants[24]*constants[25]*constants[30]*constants[28]*1.00000)
    algebraic[64] = algebraic[61]+algebraic[62]
    algebraic[68] = algebraic[64]+algebraic[66]
    algebraic[65] = algebraic[62]
    algebraic[69] = 1.00000*((constants[19]*constants[27]*algebraic[67]+constants[27]*algebraic[65])/algebraic[68])
    algebraic[24] = algebraic[21]/states[6]
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