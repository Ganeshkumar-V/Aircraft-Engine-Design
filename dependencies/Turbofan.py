import math
import numpy as np
def TurboFan_without_AB(pi_c, pi_f, TL, M0, h, Fuel_name, alpha, d):
    ita_m = 0.99
    Ad = (math.pi / 4) * d ** 2

    if TL == 1:

        pi_dmax = 0.90
        e_c = 0.8
        e_f = 0.78
        pi_b = 0.9
        ita_b = 0.85
        e_t = 0.8
        pi_n = 0.93
        T_t4 = 1110



    elif TL == 2:

        pi_dmax = 0.95
        e_c = 0.84
        e_f = 0.82
        pi_b = 0.92
        ita_b = 0.91
        e_t = 0.85
        pi_n = 0.96
        T_t4 = 1390


    elif TL == 3:

        pi_dmax = 0.98
        e_c = 0.88
        e_f = 0.86
        pi_b = 0.94
        ita_b = 0.96
        e_t = 0.89
        pi_n = 0.97
        T_t4 = 1780


    else:

        pi_dmax = 0.995
        e_c = 0.90
        e_f = 0.89
        pi_b = 0.95
        ita_b = 0.99
        e_t = 0.9
        pi_n = 0.98
        T_t4 = 2000

    # input conditions
    pi_fn = pi_n
    if h == 0:
        T0 = 300
        P0 = 101325
    else:
        T0 = 216.66
        P0 = 22650 * math.exp(1.73 - 0.000157 * h * 1000)  # pa

    gamma_c = 1.4
    gamma_t = 1.3
    cp_c = 1.0048e3  # J/kg-K
    cp_t = 1.1556e3  # J/kg-K

    # Fuel
    if Fuel_name == "JET A":
        hR = 42.02e6  # MJ/kg
    elif Fuel_name == "JET A1":
        hR = 43.15e6  # MJ/kg
    elif Fuel_name == "AVGAS 100L":
        hR = 43.5e6  # MJ/kg
    else:
        hR = 43.5e6  # MJ/kg

    # engine
    R_c = ((gamma_c - 1) / (gamma_c)) * cp_c
    R_t = ((gamma_t - 1) / gamma_t) * cp_t
    a0 = math.sqrt(gamma_c * R_c * T0)
    V0 = a0 * M0
    print(" Velocity, V0 : ", V0)
    tau_r = 1 + ((gamma_c - 1) / 2) * M0 ** 2
    pi_r = tau_r ** (gamma_c / (gamma_c - 1))
    if M0 <= 1:
        ita_r = 1
    else:
        ita_r = 1 - 0.075 * (M0 - 1) ** 1.35
    pi_d = pi_dmax * ita_r
    tau_l = (cp_t * T_t4) / (cp_c * T0)
    tau_c = pi_c ** ((gamma_c - 1) / (gamma_c * e_c))
    ita_c = (pi_c ** ((gamma_c - 1) / gamma_c) - 1) / (tau_c - 1)
    tau_f = (pi_f) ** ((gamma_c - 1) / (gamma_c * e_f))
    ita_f = ((pi_f) ** ((gamma_c - 1) / gamma_c) - 1) / (tau_f - 1)
    f = (tau_l - tau_r * tau_c) / ((ita_b * hR) / (cp_c * T0) - tau_l)
    tau_t = 1 - (1 / (ita_m * (1 + f))) * (tau_r / tau_l) * (tau_c - 1 + alpha * (tau_f - 1))
    pi_t = tau_t ** (gamma_t / (gamma_t - 1) * e_c)
    ita_t = (1 - tau_t) / (1 - tau_t ** (1 / e_t))

    # Pr_t19=P_t19/P_19
    Pr_t19 = ((gamma_c + 1) / 2) ** (gamma_c / (gamma_c - 1))
    P0P19 = Pr_t19 / (pi_r * pi_d * pi_f * pi_fn)

    # Pr_t9=P_t9/P_9
    Pr_t9 = ((gamma_t + 1) / 2) ** (gamma_t / (gamma_t - 1))
    P0P9 = Pr_t9 / (pi_r * pi_d * pi_c * pi_b * pi_t * pi_n)

    M9 = math.sqrt((2 / (gamma_t - 1)) * (Pr_t9 ** ((gamma_t - 1) / gamma_t) - 1))

    T9T0 = (tau_l * tau_t * cp_c) / (Pr_t9 ** ((gamma_t - 1) / gamma_t) * cp_t)

    V9a0 = M9 * (((gamma_t * R_t * T9T0) / (gamma_c * R_c)) ** 0.5)

    M19 = math.sqrt((2 / (gamma_c - 1)) * (Pr_t19 ** ((gamma_c - 1) / gamma_c) - 1))

    T19T0 = (tau_r * tau_f) / (Pr_t19 ** ((gamma_c - 1) / gamma_c))

    V19a0 = M19 * math.sqrt(T19T0)
    rho_d = (P0 / (R_c * T0))
    print(" Density In : ", rho_d)
    m_o = rho_d * Ad * V0
    print(" mass flow rate: ", m_o)
    # SPECIFIC THRUST
    ST = (1 / (1 + alpha)) * a0 * (
                (1 + f) * V9a0 - M0 + (1 + f) * ((R_t * T9T0) / (R_c * V9a0)) * ((1 - P0P9) / gamma_c)) + (
                     alpha / (1 + alpha)) * a0 * (V19a0 - M0 + (T19T0 / V19a0) * ((1 - P0P19) / gamma_c))
    # TSFC
    S = f / ((1 + alpha) * ST)
    # Total thrust
    T = m_o * ST
    # efficiencies
    ita_th = (a0 ** 2 * ((1 + f) * (V9a0 ** 2) + alpha * (V19a0 ** 2) - (1 + alpha) * (M0 ** 2))) / (2 * f * hR)
    ita_p = (2 * M0 * ((1 + f) * V9a0 + alpha * V19a0 - (1 + alpha) * M0)) / (
    ((1 + f) * (V9a0 ** 2) + alpha * (V19a0 ** 2) - (1 + alpha) * (M0 ** 2)))

    # property values at every station
    Pt0 = P0 * pi_r
    Tt0 = T0 * tau_r
    Pt2 = pi_d * Pt0
    Tt2 = Tt0
    Pt3 = pi_c * Pt2
    Tt3 = tau_c * Tt2
    Pt4 = pi_b * Pt3  # doubt ask about tau_b
    Tt4 = T_t4
    Tt5 = tau_t * T_t4
    Pt5 = pi_t * Pt4
    Tt7 = Tt5
    Pt7 = Pt5
    P9 = P0 / (P0P9)
    Pt9 = Pr_t9 * P9
    Pt8 = Pt9 / pi_n
    T9 = T9T0 * T0
    Tt9 = (1 + ((gamma_t - 1) / 2) * M9 ** 2) * T9
    Tt8 = Tt9
    Pt13 = Pt2 * pi_f
    Tt13 = Tt2 * tau_f
    T19 = T19T0 * T0
    Tt19 = (1 + ((gamma_t - 1) / 2) * M19 ** 2) * T19
    P19 = P0 / (P0P19)
    Pt19 = Pr_t19 * P9
    Pt17 = Pt19 / pi_fn
    Tt17 = Tt13
    TSFC = S
    Ptl = [Pt0, Pt2, Pt3, Pt4, Pt5, Pt7, Pt8, Pt9, Pt13, Pt19]
    Ttl = [Tt0, Tt2, Tt3, Tt4, Tt5, Tt7, Tt8, Tt9, Tt13, Tt19]
    try:
        etal = [ita_th[0], ita_p[0], ita_th[0]*ita_p[0]]
    except:
        etal = [ita_th, ita_p, ita_th*ita_p]
    Thrustl = [ST, m_o, T]
    try:
        Fuell = [f[0], TSFC[0]]
    except:
        Fuell = [f, TSFC]
    Pt = np.array(Ptl, dtype=float)
    Tt = np.array(Ttl, dtype=float)
    eta = np.array(etal, dtype=float)
    Thrust = np.array(Thrustl, dtype=float)
    Fuel = np.array(Fuell, dtype=float)
    return Thrust, Fuel, eta, Pt, Tt

#
# pi_c = float(20)
# pi_f = float(2)
# TL = int(4)
# M0=float(0.89)
# h=float(12)
# Fuel_name="JET A"
# alpha = float(2)
# d = float(1.95)
# X=TurboFan_without_AB(pi_c,pi_f,TL,M0,h,Fuel_name,alpha,d)
# print("Specific Thrust", X[0])