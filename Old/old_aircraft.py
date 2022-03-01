import numpy as np

'''
Functions do not have to be edited, all other variables have to as specified.

Beware at the following variables:

- (1) aircraft.X_cg is measured from the aircraft nose (and thus negative);
- (2) all sweep input are measured at the half chord.

All parameters are abbreviated in accordance with the stability code guidelines.
'''


class aircraft():

# (1/6) general analysis definitions.

    # Set True for stability analysis, False for control analysis.
    stability = True

    # Set True for random generated gusts.
    gust = False

    # Set True in case of horizontal tail.
    horizontal = False

    # Set True for double (symmetric) vertical tail.
    double = True

    # Simulation time (in seconds).
    t_simulation = 30

# (2/6) aircraft properties.

    # (1/8) location dependant geometric properties.

        # (1/3) wing properties.

    def chord_wing(Y):

        '''
        Wing chord and aerodynamic center X-position as a function of wing span position.
        :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
        :return: wing chord and aerodynamic center X-position at the specified Y-position.
        '''

        c_w = aircraft.c_r_w - aircraft.c_r_w * (1 - aircraft.taper_w) * \
              ((abs(Y) * np.cos(aircraft.sweep_w)) / aircraft.b_w)

        X_ac_w = aircraft.X_LE_w - (aircraft.c_r_w / 2) - abs(Y) * np.sin(aircraft.sweep_w) \
                 + np.cos(aircraft.sweep_w) * c_w * (0.5 - aircraft.ac_c_w)

        return c_w, X_ac_w

    def moment_arms_wing(Y):

        '''
        X-, Y-, Z-moment arms from the cg as a function of wing span position.
        :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
        :return: X-, Y-, Z-moment arms at the specified Y-position for the wing.
        '''

        c_w, X_ac_w = aircraft.chord_wing(Y)

        l_X_w = aircraft.X_cg - X_ac_w
        l_Y_w = - (aircraft.Y_cg - np.cos(aircraft.sweep_w) * Y)
        l_Z_w = - (aircraft.Z_cg - np.sin(aircraft.dihedral_w) * Y * (Y / abs(Y)))

        return l_X_w, l_Y_w, l_Z_w

        # (2/3) vertical tail properties.

    def chord_vertical(Z):

        '''
        Vertical tail chord and aerodynamic center X-position as a function of vertical tail span position.
        :param Z: Z-position along the vertical tail (from 0 to - vertical tail span).
        :return: vertical tail chord and aerodynamic center X-position at the specified Z-position.
        '''

        c_v = aircraft.c_r_v - aircraft.c_r_v * (1 - aircraft.taper_v) * \
              ((abs(Z) * np.cos(aircraft.sweep_v)) / aircraft.b_v)

        X_ac_v = aircraft.X_LE_v - (aircraft.c_r_v / 2) + Z * np.sin(aircraft.sweep_v) + \
                 np.cos(aircraft.sweep_v) * c_v * (0.5 - aircraft.ac_c_v)

        return c_v, X_ac_v

    def moment_arms_vertical(Z):

        '''
        X-, Y-, Z-moment arms from the cg.
        :return: X-, Y-, Z-moment arms for the vertical tail.
        '''

        c_v, X_ac_v = aircraft.chord_vertical(Z)

        l_X_v = aircraft.X_cg - X_ac_v
        l_Y_v = 0
        l_Z_v = aircraft.Z_cg - Z * np.cos(aircraft.sweep_v)

        return l_X_v, l_Y_v, l_Z_v

        # (3/3) horizontal tail properties.

    def chord_horizontal(Y):

        '''
        Horizontal tail chord and aerodynamic center X-position as a function of horizontal tail span position.
        :param Y: Y-position along the horizontal tail (from - horizontal tail span / 2 to horizontal tail span / 2).
        :return: horizontal tail chord and aerodynamic center X-position at the specified Y-position.
        '''

        c_h = aircraft.c_r_h - aircraft.c_r_h * (1 - aircraft.taper_h) * \
              ((abs(Y) * np.cos(aircraft.sweep_h)) / aircraft.b_h)

        X_ac_h = aircraft.X_LE_h - (aircraft.c_r_h / 2) - abs(Y) * np.sin(aircraft.sweep_h) \
                 + np.cos(aircraft.sweep_h) * c_h * (0.5 - aircraft.ac_c_h)

        return c_h, X_ac_h

    def moment_arms_horizontal(Y):

        '''
        X-, Y-, Z-moment arms from the cg as a function of horizontal tail span position.
        :param Y: Y-position along the horizontal tail (from - horizontal tail span / 2 to horizontal tail span / 2).
        :return: X-, Y-, Z-moment arms at the specified Y-position for the horizontal tail.
        '''

        c_h, X_ac_h = aircraft.chord_horizontal(Y)

        l_X_h = aircraft.X_cg - X_ac_h
        l_Y_h = - (aircraft.Y_cg - np.cos(aircraft.sweep_h) * Y)
        l_Z_h = - (aircraft.Z_cg - np.sin(aircraft.dihedral_h) * Y * (Y / abs(Y)))

        return l_X_h, l_Y_h, l_Z_h

    # (2/8) location and angle dependant aerodynamic properties.

        # (1/3) wing properties.

    def coefficients_wing(alpha_w, Y):

        '''
        Lift and drag polars as a function of alpha and control input.
        :param alpha_w: wing alpha.
        :return: wing lift and drag coefficient.
        '''

        C_L_w = (alpha_w + aircraft.incidence_w - aircraft.alpha_0_w) * aircraft.C_L_w_alpha
        C_D_w = ((alpha_w + aircraft.incidence_w - aircraft.alpha_0_w) * aircraft.C_D_w_alpha) ** 2 \
                + aircraft.C_D_w_0
        C_M_ac_w = aircraft.C_M_ac_w_0

        if aircraft.Y_ailerons_start < abs(Y) < aircraft.Y_ailerons_end \
                and aircraft.t_ailerons_start < aircraft.t < aircraft.t_ailerons_end:

            C_L_w -= aircraft.C_L_w_ailerons * aircraft.ailerons * (Y / abs(Y))
            C_D_w -= aircraft.C_D_w_ailerons * aircraft.ailerons * (Y / abs(Y))

        if aircraft.Y_elevator_start < abs(Y) < aircraft.Y_elevator_end \
                and aircraft.t_elevator_start < aircraft.t < aircraft.t_elevator_end:

            C_L_w += aircraft.C_L_w_elevator * aircraft.elevator
            C_D_w += aircraft.C_D_w_elevator * aircraft.elevator

        if aircraft.t_elevator_start < aircraft.t < aircraft.t_elevator_end:

            C_M_ac_w += aircraft.C_M_ac_w_elevator * aircraft.elevator

        return C_L_w, C_D_w, C_M_ac_w

        # (2/3) vertical tail properties.

    def coefficients_vertical(beta_v, Z):

        '''
        Lift and drag polars as a function of beta and control input.
        :param beta_v: vertical tail beta.
        :return: vertical tail lift and drag coefficient.
        '''

        C_L_v = beta_v * aircraft.C_L_v_beta
        C_D_v = (beta_v * aircraft.C_D_v_beta) ** 2 + aircraft.C_D_v_0
        C_M_ac_v = aircraft.C_M_ac_v_0

        if aircraft.Z_rudder_start < abs(Z) < aircraft.Z_rudder_end \
                and aircraft.t_rudder_start < aircraft.t < aircraft.t_rudder_end:

            C_L_v += aircraft.C_L_v_rudder * aircraft.rudder
            C_D_v += aircraft.C_D_v_rudder * aircraft.rudder
            C_M_ac_v += aircraft.C_M_ac_v_rudder * aircraft.rudder

        return C_L_v, C_D_v, C_M_ac_v

        # (3/3) horizontal tail properties.

    def coefficients_horizontal(alpha_h, Y):

        '''
        Lift and drag polars as a function of alpha and control input.
        :param alpha_h: horizontal tail alpha.
        :return: horizontal tail lift and drag coefficient.
        '''

        C_L_h = (alpha_h * (1 - aircraft.downwash_h) + aircraft.incidence_h) * aircraft.C_L_h_alpha
        C_D_h = ((alpha_h * (1 - aircraft.downwash_h) + aircraft.incidence_h) * aircraft.C_D_h_alpha) ** 2 + \
                aircraft.C_D_h_0
        C_M_ac_h = aircraft.C_M_ac_h_0

        if aircraft.Y_elevator_start < abs(Y) < aircraft.Y_elevator_end \
                and aircraft.t_elevator_start < aircraft.t < aircraft.t_elevator_end:

            C_L_h += aircraft.C_L_h_elevator * aircraft.elevator
            C_D_h += aircraft.C_D_h_elevator * aircraft.elevator
            C_M_ac_h += aircraft.C_M_ac_h_elevator * aircraft.elevator

        return C_L_h, C_D_h, C_M_ac_h

    # (3/8) mass properties.

    m = 0.270
    g = 9.807
    W = m * g

    X_cg = - 0.04
    Y_cg = 0.0
    Z_cg = 0.0

    I_XX = 0.012
    I_YY = 0.002
    I_ZZ = 0.014

    # (4/8) thrust setting.

    T = 0

    # (5/8) aerodynamic center properties.

    ac_c_w = 0.25
    ac_c_v = 0.25
    ac_c_h = 0.25

    Y_ac_w = 0
    Y_ac_v = 0
    Y_ac_h = 0

    Z_ac_w = 0
    Z_ac_v = 0
    Z_ac_h = 0

    C_M_ac_w_0 = 0.055
    C_M_ac_v_0 = 0
    C_M_ac_h_0 = 0

    # (6/8) wing properties.

    incidence_w = 0 * (np.pi / 180)

    alpha_0_w = 0.7 * (np.pi / 180)

    C_L_w_alpha = 0.09 * (180 / np.pi)
    C_D_w_alpha = 0.02 * (180 / np.pi)
    C_D_w_0 = 0.0075

    X_LE_w = 0

    b_w = 1.2
    c_r_w = 0.2

    dihedral_w = 5 * (np.pi / 180)
    sweep_w = 5 * (np.pi / 180)
    taper_w = 0.5

    S_w = (b_w / 2) * c_r_w * (1 + taper_w)
    mac_w = (2 / 3) * c_r_w * ((1 + taper_w + taper_w ** 2) / (1 + taper_w))

    # (7/8) vertical tail properties.

    C_L_v_beta = 0.09 * (180 / np.pi)
    C_D_v_beta = 0.02 * (180 / np.pi)
    C_D_v_0 = 0.007

    X_LE_v = - 0.1

    b_v = 0.24
    c_r_v = 0.12

    sweep_v = 5 * (np.pi / 180)
    taper_v = 0.5

    S_v = (b_v / 2) * c_r_v * (1 + taper_v)
    mac_v = (2 / 3) * c_r_v * ((1 + taper_v + taper_v ** 2) / (1 + taper_v))

    # (8/8) horizontal tail properties.

    incidence_h = 0 * (np.pi / 180)
    downwash_h = 0 * (np.pi / 180)

    C_L_h_alpha = 0.0 * (180 / np.pi)
    C_D_h_alpha = 0.0 * (180 / np.pi)
    C_D_h_0 = 0.0

    X_LE_h = - 1.0

    b_h = 0.4
    c_r_h = 0.1

    dihedral_h = 0 * (np.pi / 180)
    sweep_h = 0 * (np.pi / 180)
    taper_h = 0.5

    S_h = (b_h / 2) * c_r_h * (1 + taper_h)
    mac_h = (2 / 3) * c_r_h * ((1 + taper_h + taper_h ** 2) / (1 + taper_h))

# (3/6) initial conditions.

    # (1/2) angle conditions.

    pitch = 0 * (np.pi / 180)
    pitch_dot = 0 * (np.pi / 180)

    roll = 10 * (np.pi / 180)
    roll_dot = 0 * (np.pi / 180)

    yaw = 0 * (np.pi / 180)
    yaw_dot = 0 * (np.pi / 180)

    # (2/2) dynamic conditions.

    rho = 1.225

    X_dot_e = 8
    Y_dot_e = 0
    Z_dot_e = 0

    velocity_vector_e = np.array([[X_dot_e],
                                  [Y_dot_e],
                                  [Z_dot_e]])

    velocity_vector_b = np.array([[0],
                                  [0],
                                  [0]])

    velocity_vector_a = np.array([[0],
                                  [0],
                                  [0]])

# (4/6) gust definitions.

    # Set gust maxima. This is a random generated number from the normal distribution (mean at zero) times this value.
    gust_intensity = 0.5

    # These are the gust values when random generated gusts are turned off.
    X_dot_gust = 0
    Y_dot_gust = 0
    Z_dot_gust = 0

# (5/6) control inputs.

    # (1/5) input input.

    ailerons = 0
    elevator = 0
    rudder = 0

    # (2/5) input time.

    t_ailerons_start = 0
    t_ailerons_end = 0

    t_elevator_start = 0
    t_elevator_end = 0

    t_rudder_start = 0
    t_rudder_end = 0

    # (3/5) ailerons properties.

    Y_ailerons_start = 0
    Y_ailerons_end = 0

    C_L_w_ailerons = 0
    C_D_w_ailerons = 0

    # (4/5) elevator properties.

    Y_elevator_start = 0
    Y_elevator_end = 0

    C_L_w_elevator = 0
    C_D_w_elevator = 0
    C_M_ac_w_elevator = 0

    C_L_h_elevator = 0
    C_D_h_elevator = 0
    C_M_ac_h_elevator = 0

    # (5/5) rudder properties.

    Z_rudder_start = 0
    Z_rudder_end = 0

    C_L_v_rudder = 0
    C_D_v_rudder = 0
    C_M_ac_v_rudder = 0

# (6/6) time step and Riemann sum step sizes.

    '''
    WARNING:
    
    Increasing time step and Riemann sum step sizes reduces runtime, but highly affects the solution of the program.
    It is advised to keep the steps as they are programmed.
    '''

    ds_w = b_w / 30
    ds_v = b_v / 10
    ds_h = b_h / 10

    t = 0
    dt = 0.01
