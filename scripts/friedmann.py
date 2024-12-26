#!/usr/bin/python2

#------------------------------
# Compute cosmic time in Gyrs
# for a LambdaCDM universe
# with given redshift
#------------------------------



#=================================================
def friedman(axp_min=1e-3, axp_max=1):
#=================================================
    """
    Integrate friedman equation to get table of
    expansion factors and times.
    Gives a in units of H0.
    See ramses/utils/f90/utils.f90/subroutine friedman
    """

    import numpy as np


    def dadtau(axp_tau):
        dadtau = axp_tau**3 * (omega_m + omega_l*axp_tau**3 + omega_k * axp_tau)
        return np.sqrt(dadtau)

    def dadt(axp_t):
        dadt = 1/axp_t * (omega_m + omega_l*axp_t**3 + omega_k*axp_t)
        return np.sqrt(dadt)


    omega_m     =  0.272
    omega_l     =  0.728
    omega_k     =  0.0

    epsilon = 1e-4

    axp_tau = 1   # expansion factor
    axp_t = 1     # expansion factor
    t = 0         # look-back time
    tau = 0       # conformal time



    #-----------------------------------
    # Integrate backwards in time
    #-----------------------------------
    a_past = np.zeros(1000000, dtype='float')
    t_past = np.zeros(1000000, dtype='float')
    #  tau_past = np.zeros(1000000, dtype='float')
    #  H_past = np.zeros(1000000, dtype='float')

    i = 0
    while axp_tau >= axp_min or axp_t >= axp_min:
        dtau = epsilon * axp_tau / dadtau(axp_tau)
        axp_tau_pre = axp_tau - dadtau(axp_tau)*dtau/2
        axp_tau = axp_tau - dadtau(axp_tau_pre)*dtau
        tau = tau - dtau

        dt = epsilon * axp_t / dadt(axp_t)
        axp_t_pre = axp_t - dadt(axp_t)*dt/2
        axp_t = axp_t - dadt(axp_t_pre)*dt
        t = t - dt

        if (abs((t - t_past[i])/t) >= 1e-5):
            a_past[i] = axp_t
            #  H_past[i] = dadtau(axp_tau)/axp_tau
            t_past[i] = t
            #  tau_past[i] = tau

            i += 1

    a_past[i] = axp_t
    t_past[i] = t
    #  tau_past[i] = tau
    #  H_past[i] = dadtau(axp_tau)/axp_tau
    i+=1
    pastind = i




    # Reset values

    axp_tau = 1   # expansion factor
    axp_t = 1     # expansion factor
    t = 0         # look-back time
    tau = 0       # conformal time



    #-----------------------------------
    # Integrate forwards in time
    #-----------------------------------
    a_fut = np.zeros(1000000, dtype='float')
    t_fut = np.zeros(1000000, dtype='float')
    #  tau_fut = np.zeros(1000000, dtype='float')
    #  H_fut = np.zeros(1000000, dtype='float')

    i = 0
    while axp_tau <= axp_max or axp_t <= axp_max:
        dtau = epsilon * axp_tau / dadtau(axp_tau)
        axp_tau_pre = axp_tau + dadtau(axp_tau)*dtau/2
        axp_tau = axp_tau + dadtau(axp_tau_pre)*dtau
        tau = tau + dtau

        dt = epsilon * axp_t / dadt(axp_t)
        axp_t_pre = axp_t + dadt(axp_t)*dt/2
        axp_t = axp_t + dadt(axp_t_pre)*dt
        t = t + dt

        if (abs((t - t_fut[i])/t) >= 1e-5):
            a_fut[i] = axp_t
            #  H_fut[i] = dadtau(axp_tau)/axp_tau
            t_fut[i] = t
            #  tau_fut[i] = tau

            i += 1

    a_fut[i] = axp_t
    t_fut[i] = t
    #  tau_fut[i] = tau
    #  H_fut[i] = dadtau(axp_tau)/axp_tau
    i+=1
    futind = i



    a_out = np.concatenate((a_fut[:futind][::-1], a_past[:pastind]))
    t_out = np.concatenate((t_fut[:futind][::-1], t_past[:pastind]))

    #  tau_out = np.concatenate((tau_fut[:futind], tau_past[:pastind]))
    #  H_out = np.concatenate((H_fut[:futind], H_past[:pastind]))


    return a_out, t_out



#========================================
def get_times(aexp):
#========================================
    """
    Compute times given the expansion factor.
    aexp: numpy array of a
    """
    import numpy as np


    H0          =  70.40
    times = np.zeros(aexp.shape[0], dtype='float')

    a_out, t_out = friedman(aexp.min(), aexp.max())

    i = 1
    for j,a in enumerate(aexp):
        while (a_out[i] > a) and (i < a_out.shape[0]-1):
            i += 1
        times[j] = t_out[i]*(a-a_out[i-1])/(a_out[i]-a_out[i-1]) + \
            t_out[i-1]*(a-a_out[i])/(a_out[i-1]-a_out[i])



    times *= 1.0/(H0*1e5/3.08e24)/(365*24*3600*1e9)


    return times




#=============================================
def get_aexp(times, axp_min, axp_max):
#=============================================
    """
    Compute aexp given times (in Gyrs)
    times: numpy array of times
    """
    import numpy as np

    H0          =  70.40
    a_exp = np.zeros(times.shape[0], dtype='float')

    # first convert times in units of H0^-1.
    # expected to be given in Gyrs
    times *= (H0*1e5/3.08e24)*(365*24*3600*1e9)

    a_out, t_out = friedman(axp_min, axp_max)

    i = 1
    for j,t in enumerate(times):
        while (t_out[i] > t) and (i < t_out.shape[0]-1):
            i += 1
        a_exp[j] = a_out[i]*(t-t_out[i-1])/(t_out[i]-t_out[i-1]) + \
            a_out[i-1]*(t-t_out[i])/(t_out[i-1]-t_out[i])

    return a_exp







#======================
def main():
#======================

    import numpy as np

    #  z = np.array([-4.07311850e-02,-2.23992776e-02,-4.09845370e-03, 2.01168586e-02,
    #    3.95292490e-02, 6.31204091e-02, 8.64874431e-02, 1.07127784e-01,
    #    1.34036305e-01, 1.58677825e-01, 1.87752345e-01, 2.17052984e-01,
    #    2.47434153e-01, 2.77904562e-01, 3.09750557e-01, 3.48462035e-01,
    #    3.82832009e-01, 4.25329504e-01, 4.68121478e-01, 5.11411521e-01,
    #    5.58585503e-01, 6.07055693e-01, 6.64086052e-01, 7.16900290e-01,
    #    7.81658371e-01, 8.42508308e-01, 9.14907978e-01, 9.90590865e-01,
    #    1.07629268e+00, 1.17064835e+00, 1.27116501e+00, 1.37560570e+00,
    #    1.49944556e+00, 1.61903904e+00, 1.77447897e+00, 1.92055716e+00,
    #    2.11234234e+00, 2.31882521e+00, 2.54526995e+00, 2.83223076e+00,
    #    3.12450873e+00, 3.51878224e+00, 3.94000857e+00, 4.24566403e+00,
    #    4.71069866e+00, 6.01239867e+00, 6.77034216e+00, 8.54055714e+00,
    #    1.07137145e+01, 1.49354445e+01, 2.30203750e+01, 4.82549864e+01,
    #    1.00000002e+02])

    #  a = 1/(1+z)

    #  times = get_times(a)


    times = np.array([0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7,  3., 3.3, 3.5,
    3.7,  3.9,  4.1,  4.3,  4.5,  4.7,  4.9,  5.1,  5.3,
    5.5,  5.7,  5.9,  6.1,  6.3,  6.5,  6.7,  6.9,  7.1,  7.3,  7.5,
    7.7,  7.9,  8.1,  8.3,  8.5,  8.7,  8.9,  9.1,  9.3,  9.5,  9.7,
    9.9, 10.1, 10.3, 10.5, 10.7, 10.9, 11.1, 11.3, 11.5, 11.7, 11.9,
    12.1, 12.3, 12.5, 12.7, 12.9, 13.1, 13.3, 13.5, 13.7, 13.9, 14.1,
    14.3])-13.8
    # WARNING: t=0 is expected to be today.
    # WARNING: array is expected to start with highest cosmological time, e.g. with a >= 1 and go towards a=0
    # => invert array
    times = times[::-1]


    aexp = get_aexp(times, axp_min=1e-4, axp_max=1.1)

    for a in aexp[::-1]:
        print str(round(a,3))+",",
    print

    print "nelements:",aexp.shape[0]

    return



if __name__ == '__main__':
    main()
