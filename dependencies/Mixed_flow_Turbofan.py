# mixed flow  without afterburner
import math
import numpy as np

def mixed_flow_turbofan_without_afterburner(pi_c, pi_f, TL, M0, h, Fuel_name, d):
    ita_m = 0.99
    pi_M_max = 0.98
    M6 = 0.4
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
        T_t7 = 2500


    elif TL == 2:

        pi_dmax = 0.95
        e_c = 0.84
        e_f = 0.82
        pi_b = 0.92
        ita_b = 0.91
        e_t = 0.85
        pi_n = 0.96
        T_t4 = 1390
        T_t7 = 3000


    elif TL == 3:

        pi_dmax = 0.98
        e_c = 0.88
        e_f = 0.86
        pi_b = 0.94
        ita_b = 0.96
        e_t = 0.89
        pi_n = 0.97
        T_t4 = 1780
        T_t7 = 3600


    else:

        pi_dmax = 0.995
        e_c = 0.90
        e_f = 0.89
        pi_b = 0.95
        ita_b = 0.99
        e_t = 0.9
        pi_n = 0.98
        T_t4 = 2000
        T_t7 = 4000

    # input conditions
    pi_fn = pi_n
    if h == 0:
        T0 = 302
        P0 = 101325
    else:
        T0 = 216.66
        P0 = 22650 * math.exp(1.73 - 0.000157 * h * 1000)  # pa

    gamma_c = 1.4
    gamma_t = 1.3
    gamma_AB = gamma_t

    cp_c = 1.0048e3  # J/kg-K
    cp_t = 1.1556e3  # J/kg-K

    # Fuel
    if Fuel_name == "JET A":
        hR = 42.02e6  # MJ/kg
    elif Fuel_name == "JET A1":
        hR = 43.15e6  # MJ/kg
    else:
        hR = 43.5e6  # MJ/kg

    # engine
    R_c = ((gamma_c - 1) / gamma_c) * cp_c
    R_t = ((gamma_t - 1) / gamma_t) * cp_t
    cp_AB = cp_t
    R_AB = ((gamma_AB - 1) / gamma_AB) * cp_AB
    a0 = math.sqrt(gamma_c * R_c * T0)
    V0 = a0 * M0
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
    f = (tau_l - tau_r * tau_c) / ((ita_b * hR) / (cp_c * T0) - tau_l)
    tau_f = pi_f ** ((gamma_c - 1) / (gamma_c * e_f))
    ita_f = (pi_f ** ((gamma_c - 1) / gamma_c) - 1) / (tau_f - 1)
    alpha = ((ita_m) * (1 + f) * (tau_l / tau_r) * (1 - (pi_f / (pi_c * pi_b)) ** ((gamma_t - 1) * e_t / gamma_t)) - (
                tau_c - 1)) / (tau_f - 1)
    tau_t = 1 - (1 / (ita_m * (1 + f))) * (tau_r / tau_l) * (tau_c - 1 + alpha * (tau_f - 1))
    pi_t = tau_t ** (gamma_t / ((gamma_t - 1) * e_c))
    ita_t = (1 - tau_t) / (1 - tau_t ** (1 / e_t))
    tau_lAB = (cp_AB * T_t7) / (cp_c * T0)
    # Pr_t16=P_t16/P_t6
    Pr_t16 = pi_f / (pi_c * pi_b * pi_t)
    M16 = math.sqrt((2 / (gamma_c - 1)) * (
                (Pr_t16 * (1 + ((gamma_t - 1) / 2) * M6 ** 2) ** gamma_t / (gamma_t - 1)) ** (
                    (gamma_c - 1) / gamma_c) - 1))

    alpha_star = alpha / (1 + f)
    cp_6A = (cp_t + alpha_star * cp_c) / (1 + alpha_star)
    R_6A = (R_t + alpha_star * R_c) / (1 + alpha_star)
    gamma_6A = cp_6A / (cp_6A - R_6A)

    # T_t16 = T_t16/T_t6
    T_t16 = (T0 * tau_r * tau_f) / (T_t4 * tau_t)

    tau_M = (cp_t / cp_6A) * (1 + alpha_star * (cp_c / cp_t) * (T_t16)) / (1 + alpha_star)

    def PHI(M, gamma):
        phi_s = M ** 2 * (1 + ((gamma - 1) / 2) * M ** 2) / ((1 + gamma * M ** 2) ** 2)
        return phi_s

    phi_M6_gamma_6 = PHI(M6, gamma_t)
    phi_M16_gamma_16 = PHI(M16, gamma_c)
    phi = (((1 + alpha_star) / (1 / (math.sqrt(phi_M6_gamma_6)) + alpha_star * math.sqrt(
        ((R_c * gamma_t) / (R_t * gamma_t)) * (T_t16 / phi_M16_gamma_16)))) ** 2) * (R_6A / R_t) * (
                      gamma_t / gamma_6A) * tau_M
    M6A = ((2 * phi) / (1 - 2 * gamma_6A * phi + (1 - 2 * (gamma_6A + 1) * phi) ** 0.5)) ** 0.5
    # Area_ratio= A16/A6
    A16A6 = (alpha_star * math.sqrt(T_t16)) / ((M16 / M6) * math.sqrt(
        (gamma_c * R_t * (1 + ((gamma_c - 1) / 2) * M16 ** 2)) / (gamma_t * R_c * (1 + ((gamma_t - 1) / 2) * M6 ** 2))))

    def MFP(M, gamma, R):
        mfp = (M * math.sqrt(gamma / R)) / (1 + ((gamma - 1) / 2) * M ** 2) ** (gamma + 1) / (2 * (gamma - 1))
        return mfp

    MFP_6 = MFP(M6, gamma_t, R_t)
    MFP_6A = MFP(M6A, gamma_6A, R_6A)
    pi_M_ideal = (((1 + alpha_star) * math.sqrt(tau_M)) / (1 + A16A6)) * ((MFP_6) / (MFP_6A))
    pi_M = pi_M_max * pi_M_ideal
    # Pr_t9=P_t9/P_9
    Pr_t9 = ((gamma_t + 1) / 2) ** (gamma_t / (gamma_t - 1))
    P0P9 = (Pr_t9) / (pi_r * pi_d * pi_c * pi_b * pi_t * pi_M * pi_n)
    # Afterburner off
    cp_9 = cp_6A
    R_9 = R_6A
    gamma_9 = gamma_6A
    f_AB = 0
    T9T0 = ((T_t4 * tau_t * tau_M) / T0) / (Pr_t9 ** ((gamma_9 - 1) / gamma_9))
    # GENERAL
    M9 = math.sqrt((2 / (gamma_9 - 1)) * (Pr_t9 ** ((gamma_9 - 1) / gamma_9) - 1))
    V9a0 = M9 * (((gamma_9 * R_t * T9T0) / (gamma_c * R_c)) ** 0.5)
    f_O = (f / (1 + alpha)) + f_AB
    v9 = V9a0*a0
    # To find out overall mass flow rate
    rho_d = (P0 / (R_c * T0))
    m_o = rho_d * Ad * V0

    # Specific thrust
    ST = float(a0 * ((1 + f_O) * V9a0 - M0 + (1 + f_O) * ((R_9 * T9T0) / (R_c * V9a0)) * ((1 - P0P9) / gamma_c)))
    print(ST)
    # Total thrust
    T = m_o * ST
    # TSFC
    TSFC = f_O / ST
    # Thermal efficiency
    ita_th = (a0 ** 2 * ((1 + f_O) * (V9a0 ** 2) - M0 ** 2)) / (2 * f_O * hR)
    # Propulsive efficiency
    ita_p = 2*(V0/v9)/(1 + V0/v9)

    # Property values at different station
    Pt0 = P0 * pi_r
    Tt0 = T0 * tau_r
    Pt2 = pi_d * Pt0
    Tt2 = Tt0
    Pt3 = pi_c * Pt2
    Tt3 = tau_c * Tt2
    Pt4 = pi_b * Pt3
    Tt6 = tau_t * T_t4
    Pt6 = pi_t * Pt4
    # Mixer
    Tt6A = pi_M * Tt6
    Pt6A = tau_M * Pt6

    # As afterburner off
    Tt9 = T_t4 * tau_M * tau_t
    Pt9 = P0 * pi_r * pi_d * pi_c * pi_b * pi_t * pi_M * pi_n

    Pt7 = Pt9

    Tt8 = Tt9
    Pt8 = Pt9 / pi_n

    Tt16 = Tt6 * (T0 * tau_r * tau_f) / (T_t4 * tau_t)
    Pt16 = Pt6 * (pi_f / (pi_c * pi_b * pi_t))

    Pt13 = Pt2 * pi_f
    Tt13 = Tt2 * tau_f

    Thrustl = [ST, m_o, T, float(alpha)]
    Fuell = [float(f_O), float(TSFC)]
    Effl = [float(ita_th), float(ita_p), float(ita_th*ita_p)]
    Ptl = [Pt0, Pt2, Pt3, Pt4, Pt6, Pt6A, Pt7, Pt8, Pt9, Pt13, Pt16]
    Ttl = [Tt0, Tt2, Tt3, T_t4, Tt6, Tt6A, T_t7, Tt8, Tt9, Tt13, Tt16]
    Thrust = np.array(Thrustl)
    Fuel = np.array(Fuell)
    Eff = np.array(Effl)
    Pt = np.array(Ptl)
    Tt = np.array(Ttl)
    return Thrust, Fuel, Eff, Pt, Tt