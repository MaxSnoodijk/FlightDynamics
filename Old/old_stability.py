import numpy as np
import matplotlib.pyplot as plt
import time

from old_aircraft import aircraft

'''
This code can be used for stability analysis of aircraft. 
All parameters in the Aircraft Python file have to be filled in for the to be tested aircraft.

The following limitations for the program apply:

- (1) lift is a linear function of alpha (or beta);
- (2) drag is quadratic function of alpha (or beta);
- (3) no twist for all lifting surfaces;
- (4) no dihedral for the vertical tail;
- (5) horizontal and vertical tail use a symmetric airfoil;
- (6) roll, pitch and yaw input do not exceed 90 degrees;
- (7) constant induced angle of attack along the spans.

Note: (1), (2), (3) and (7) could (relatively) easily be changed for the code by rewriting the coefficients
functions in the Aircraft Python file if needed, (4) is simply not necessary to implement,
(5) requires moment addition after Riemann sums in this Python file and 
(7) is due to Euler singularities in transformation matrices.

Assumptions:

- (1) no fuselage effects;
- (2) small angle approximation with rate derivatives;
- (3) pitch rate derivative for wing negligible.

Coordinate system:

Right hand Cartesian coordinates, X, Y, Z, M_X, M_Y, M_Z; 
X along longitudinal axis pointing away the aircraft nose (!);
Y along lateral axis pointing towards the right wing tip;
Z along normal axis pointing towards the Earth surface.

Variable naming:

Greek alphabet full name;
Specific variables full name;
General variables abbreviated.

e = inertial Earth magnitude_reference frame;
b = body magnitude_reference frame;
a = aerodynamic magnitude_reference frame.

dot = first derivative;
dot_dot = second derivative.
'''

begin = time.time()

# (1/7) rate derivatives per lifting surface.

    # (1/3) wing derivatives.

def roll_dot_derivative_wing(Y):

    '''
    Change in the velocity Z-component due to roll rate along the wingspan.
    :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
    :return: added Z-component to the wing aerodynamic velocity vector at the specified Y-position.
    '''

    if aircraft.stability:
        delta_Z_dot_a = 0
    else:
        delta_Z_dot_a = Y * aircraft.roll_dot

    return delta_Z_dot_a

def yaw_dot_derivative_wing(Y):

    '''
    Change in the velocity X-component due to yaw rate along the wingspan.
    :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
    :return: added X-component to the wing aerodynamic velocity vector at the specified Y-position.
    '''

    if aircraft.stability:
        delta_X_dot_a = 0
    else:
        delta_X_dot_a = - Y * aircraft.yaw_dot

    return delta_X_dot_a

    # (2/3) vertical tail derivatives.

def roll_dot_derivative_vertical(Z):

    '''
    Change in the velocity Y-component due to roll rate along the vertical tail span.
    :param Z: Z-position along the vertical tail (from 0 to - vertical tail span).
    :return: added Y-component to the vertical tail aerodynamic velocity vector.
    '''

    if aircraft.stability:
        delta_Y_dot_a = 0
    else:
        delta_Y_dot_a = - Z * aircraft.roll_dot

    return delta_Y_dot_a

def pitch_dot_derivative_vertical(Z):

    '''
    Change in the velocity Z-component due to pitch rate at the vertical tail.
    :return: added Z-component to the vertical tail aerodynamic velocity vector.
    '''

    if aircraft.stability:
        delta_Z_dot_a = 0

    else:
        c_v, X_ac_v = aircraft.chord_vertical(Z)
        delta_Z_dot_a = (aircraft.X_cg - X_ac_v) * aircraft.pitch_dot

    return delta_Z_dot_a

def yaw_dot_derivative_vertical(Z):

    '''
    Change in the velocity Y-component due to yaw rate at the vertical tail.
    :return: added Y-component to the vertical tail aerodynamic velocity vector.
    '''

    if aircraft.stability:
        delta_Y_dot_a = 0

    else:
        c_v, X_ac_v = aircraft.chord_vertical(Z)
        delta_Y_dot_a = - (aircraft.X_cg - X_ac_v) * aircraft.yaw_dot

    return delta_Y_dot_a

    # (3/3) horizontal tail derivatives.

def roll_dot_derivative_horizontal(Y):

    '''
    Change in the velocity Z-component due to roll rate along the wingspan.
    :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
    :return: added Z-component to the wing aerodynamic velocity vector at the specified Y-position.
    '''

    if aircraft.stability:
        delta_Z_dot_a = 0
    else:
        delta_Z_dot_a = Y * aircraft.roll_dot

    return delta_Z_dot_a

def pitch_dot_derivative_horizontal(Y):

    '''
    Change in the velocity Z-component due to pitch rate at the horizontal tail.
    :return: added Z-component to the horizontal tail aerodynamic velocity vector.
    '''

    if aircraft.stability:
        delta_Z_dot_a = 0

    else:
        c_h, X_ac_h = aircraft.chord_horizontal(Y)
        delta_Z_dot_a = (aircraft.X_cg - X_ac_h) * aircraft.pitch_dot

    return delta_Z_dot_a

def yaw_dot_derivative_horizontal(Y):

    '''
    Change in the velocity X-component due to yaw rate along the wingspan.
    :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
    :return: added X-component to the wing aerodynamic velocity vector at the specified Y-position.
    '''

    if aircraft.stability:
        delta_X_dot_a = 0
    else:
        delta_X_dot_a = - Y * aircraft.yaw_dot

    return delta_X_dot_a

# (2/7) output_body and moments per lifting surface.

def wing():

    '''
    1. Wing X-, Y-, Z-output_body and X-, Y-, Z-moments are calculated by Riemann sum along the wingspan.
    Reason for Riemann sum instead of integration is due to excessive runtime otherwise.
    2. Moment from the aerodynamic center is added to the Y-moment.
    :return: wing X-, Y-, Z-output_body and X-, Y-, Z-moments.
    '''

    def angles_velocities(Y):

        '''
        1. Wing aerodynamic velocity vector at the specified Y-position is calculated.
        2. Alpha, beta and the magnitude_velocity of the wing aerodynamic velocity vector are calculated.
        :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
        :return: alpha, beta and the magnitude_velocity of the wing aerodynamic velocity vector.
        '''

        # (1/2)

        X_dot_a_w = np.cos(aircraft.sweep_w) * aircraft.velocity_vector_a[0, 0] + yaw_dot_derivative_wing(Y) \
                    - (Y / abs(Y)) * np.sin(aircraft.sweep_w) * aircraft.velocity_vector_a[1, 0]
        Y_dot_a_w = np.cos(aircraft.sweep_w) * aircraft.velocity_vector_a[1, 0] \
                    + (Y / abs(Y)) * np.sin(aircraft.sweep_w) * aircraft.velocity_vector_a[0, 0]
        Z_dot_a_w = aircraft.velocity_vector_a[2, 0] + roll_dot_derivative_wing(np.cos(aircraft.sweep_w) * Y) \
                    + np.sin(aircraft.dihedral_w) * Y_dot_a_w * (Y / abs(Y))

        # (2/2)

        alpha_w = np.arctan(Z_dot_a_w / X_dot_a_w)
        beta_w = np.arctan(Y_dot_a_w / X_dot_a_w)
        V_w = np.sqrt(X_dot_a_w ** 2 + Y_dot_a_w ** 2 + Z_dot_a_w ** 2)

        return alpha_w, beta_w, V_w

    def forces_moments(Y):

        '''
        1. All required aerodynamic and geometric parameters are calculated using above formulas.
        2. Lift and drag components are calculated.
        3. Lift and drag contributions to X-, Y-, Z-output_body are calculated.
        4. Lift and drag contributions to X-, Y-, Z-moments are calculated.
        5. Local alpha and beta input are appended to local list to find (absolute) maximum per time step.
        :param Y: Y-position along the wing (from - wingspan / 2 to wingspan / 2).
        :return: wing X-, Y-, Z-output_body and X-, Y-, Z-moments at the specified Y-position.
        '''

        # (1/5)

        alpha_w, beta_w, V_w = angles_velocities(Y)
        c_w, X_ac_w = aircraft.chord_wing(Y)
        l_X_w, l_Y_w, l_Z_w = aircraft.moment_arms_wing(Y)
        C_L_w, C_D_w, C_M_ac_w = aircraft.coefficients_wing(alpha_w, Y)

        # (2/5)

        delta_L_w = 0.5 * aircraft.rho * (np.cos(beta_w) * V_w) ** 2 * C_L_w * c_w
        delta_D_w = 0.5 * aircraft.rho * (np.cos(beta_w) * V_w) ** 2 * C_D_w * c_w

        # (3/5)

        delta_F_X_w = np.cos(aircraft.sweep_w) * np.sin(alpha_w) * delta_L_w - np.cos(aircraft.sweep_w) * np.cos(alpha_w) * delta_D_w \
                      + np.sin(aircraft.dihedral_w) * np.sin(aircraft.sweep_w) * np.cos(alpha_w) * delta_L_w \
                      + np.sin(aircraft.dihedral_w) * np.sin(aircraft.sweep_w) * np.sin(alpha_w) * delta_D_w

        delta_F_Y_w = (Y / abs(Y)) * np.sin(aircraft.sweep_w) * np.sin(alpha_w) * delta_L_w \
                      - (Y / abs(Y)) * np.sin(aircraft.sweep_w) * np.cos(alpha_w) * delta_D_w \
                      - (Y / abs(Y)) * np.sin(aircraft.dihedral_w) * np.cos(aircraft.sweep_w) * np.cos(alpha_w) * delta_L_w \
                      - (Y / abs(Y)) * np.sin(aircraft.dihedral_w) * np.cos(aircraft.sweep_w) * np.sin(alpha_w) * delta_D_w

        delta_F_Z_w = - (np.cos(aircraft.dihedral_w) * np.cos(alpha_w) * delta_L_w + np.sin(alpha_w) * delta_D_w)

        # (4/5)

        delta_M_X_w = delta_F_Z_w * l_Y_w + delta_F_Y_w * l_Z_w
        delta_M_Y_w = delta_F_Z_w * l_X_w
        delta_M_Z_w = - (delta_F_X_w * l_Y_w)

        # (5/5)

        local_alpha_w_tab.append(alpha_w * (180 / np.pi))
        local_beta_w_tab.append(beta_w * (180 / np.pi))

        return delta_F_X_w, delta_F_Y_w, delta_F_Z_w, delta_M_X_w, delta_M_Y_w, delta_M_Z_w

    # These steps are required to ensure wing analysis remains symmetric after applying sweep.
    effective_span = aircraft.b_w / np.cos(aircraft.sweep_w)
    rounded_effective_span = round(effective_span / aircraft.ds_w) * aircraft.ds_w

    # (1/2) Riemann sum along the wingspan.

    steps = round(effective_span / aircraft.ds_w)
    F_X_w, F_Y_w, F_Z_w, M_X_w, M_Y_w, M_Z_w = 0, 0, 0, 0, 0, 0

    for i in range(steps):

        Y = (- rounded_effective_span / 2) + (aircraft.ds_w / 2) + i * aircraft.ds_w
        Y = round(Y, 3)

        # Wing cannot be evaluated at Y equal to zero after applying sweep due to asymmetry otherwise.
        if Y != 0:
            delta_F_X_w, delta_F_Y_w, delta_F_Z_w, delta_M_X_w, delta_M_Y_w, delta_M_Z_w = forces_moments(Y)
        else:
            delta_F_X_w, delta_F_Y_w, delta_F_Z_w, delta_M_X_w, delta_M_Y_w, delta_M_Z_w = 0, 0, 0, 0, 0, 0

        F_X_w += delta_F_X_w * aircraft.ds_w
        F_Y_w += delta_F_Y_w * aircraft.ds_w
        F_Z_w += delta_F_Z_w * aircraft.ds_w

        M_X_w += delta_M_X_w * aircraft.ds_w
        M_Y_w += delta_M_Y_w * aircraft.ds_w
        M_Z_w += delta_M_Z_w * aircraft.ds_w

    # (2/2) aerodynamic center moment addition.

    alpha_w, beta_w, V_w = angles_velocities(0.00001)
    C_L_w, C_D_w, C_M_ac_w = aircraft.coefficients_wing(alpha_w, 0)

    if aircraft.t > aircraft.dt:

        alpha_w_tab.append(alpha_w * 57.3)
        beta_w_tab.append(beta_w * 57.3)

    M_Y_ac_w = 0.5 * aircraft.rho * np.linalg.norm(aircraft.velocity_vector_b) ** 2 * \
               C_M_ac_w * aircraft.S_w * aircraft.mac_w
    M_Y_w += M_Y_ac_w

    return F_X_w, F_Y_w, F_Z_w, M_X_w, M_Y_w, M_Z_w

def vertical():

    '''
    Vertical tail X-, Y-, Z-output_body and X-, Y-, Z-moments are calculated by Riemann sum along the vertical tail span.
    Reason for Riemann sum instead of integration is due to excessive runtime otherwise.
    :return: vertical tail X-, Y-, Z-output_body and X-, Y-, Z-moments.
    '''

    def angles_velocities(Z):

        '''
        1. Vertical tail aerodynamic velocity vector at the specified Z-position is calculated.
        2. Alpha, beta and the magnitude_velocity of the vertical tail aerodynamic velocity vector is calculated.
        :param Z: Z-position along the vertical tail (from 0 to - vertical tail span).
        :return: alpha, beta and the magnitude_velocity of the vertical tail aerodynamic velocity vector.
        '''

        # (1/2)

        X_dot_a_v = np.cos(aircraft.sweep_v) * aircraft.velocity_vector_a[0, 0]
        Y_dot_a_v = aircraft.velocity_vector_a[1, 0] + roll_dot_derivative_vertical(Z) \
                    + yaw_dot_derivative_vertical(Z)
        Z_dot_a_v = aircraft.velocity_vector_a[2, 0] + pitch_dot_derivative_vertical(Z) \
                    + np.sin(aircraft.sweep_v) * aircraft.velocity_vector_a[0, 0]

        # (2/2)

        alpha_v = np.arctan(Z_dot_a_v / X_dot_a_v)
        beta_v = np.arctan(Y_dot_a_v / X_dot_a_v)
        V_v = np.sqrt(X_dot_a_v ** 2 + Y_dot_a_v ** 2 + Z_dot_a_v ** 2)

        return alpha_v, beta_v, V_v

    def forces_moments(Z):

        '''
        1. All required aerodynamic and geometric parameters are calculated using above formulas.
        2. Lift and drag components are calculated.
        3. Lift and drag contributions to X-, Y-, Z-output_body are calculated.
        4. Lift and drag contributions to X-, Y-, Z-moments are calculated.
        5. Local alpha and beta input are appended to local list to find (absolute) maximum per time step.
        :param Z: Z-position along the vertical tail (from 0 to - vertical tail span).
        :return: vertical tail X-, Y-, Z-output_body and X-, Y-, Z-moments at the specified Z-position.
        '''

        # (1/5)

        alpha_v, beta_v, V_v = angles_velocities(Z)
        c_v, X_ac_v = aircraft.chord_vertical(Z)
        l_X_v, l_Y_v, l_Z_v = aircraft.moment_arms_vertical(Z)
        C_L_v, C_D_v, C_M_ac_v = aircraft.coefficients_vertical(beta_v, Z)

        # (2/5)

        delta_L_v = 0.5 * aircraft.rho * np.cos(alpha_v) * V_v ** 2 * C_L_v * c_v
        delta_D_v = 0.5 * aircraft.rho * np.cos(alpha_v) * V_v ** 2 * C_D_v * c_v

        # (3/5)

        delta_F_X_v = np.sin(beta_v) * delta_L_v - np.cos(beta_v) * np.cos(aircraft.sweep_v) * delta_D_v
        delta_F_Y_v = - (np.cos(beta_v) * delta_L_v + np.sin(beta_v) * np.cos(aircraft.sweep_v) * delta_D_v)
        delta_F_Z_v = np.sin(aircraft.sweep_v) * delta_D_v

        # (4/5)

        delta_M_X_v = delta_F_Y_v * l_Z_v
        delta_M_Y_v = - (delta_F_X_v * l_Z_v) - Z * delta_F_Z_v
        delta_M_Z_v = - (delta_F_Y_v * l_X_v)

        # (5/5)

        local_alpha_v_tab.append(alpha_v * (180 / np.pi))
        local_beta_v_tab.append(beta_v * (180 / np.pi))

        return delta_F_X_v, delta_F_Y_v, delta_F_Z_v, delta_M_X_v, delta_M_Y_v, delta_M_Z_v

    # Riemann sum along the vertical tail span.

    steps = round(aircraft.b_v / aircraft.ds_v)
    F_X_v, F_Y_v, F_Z_v, M_X_v, M_Y_v, M_Z_v = 0, 0, 0, 0, 0, 0

    for i in range(steps):

        Z = (- aircraft.ds_v / 2) - i * aircraft.ds_v
        delta_F_X_v, delta_F_Y_v, delta_F_Z_v, delta_M_X_v, delta_M_Y_v, delta_M_Z_v = forces_moments(Z)

        F_X_v += delta_F_X_v * aircraft.ds_v
        F_Y_v += delta_F_Y_v * aircraft.ds_v
        F_Z_v += delta_F_Z_v * aircraft.ds_v

        M_X_v += delta_M_X_v * aircraft.ds_v
        M_Y_v += delta_M_Y_v * aircraft.ds_v
        M_Z_v += delta_M_Z_v * aircraft.ds_v

    if aircraft.double:
        return 2 * F_X_v, 2 * F_Y_v, 2 * F_Z_v, 2 * M_X_v, 2 * M_Y_v, 2 * M_Z_v
    else:
        return F_X_v, F_Y_v, F_Z_v, M_X_v, M_Y_v, M_Z_v

def horizontal():

    '''
    Horizontal tail X-, Y-, Z-output_body and X-, Y-, Z-moments are calculated by Riemann sum along the horizontal tail span.
    Reason for Riemann sum instead of integration is due to excessive runtime otherwise.
    :return: horizontal tail X-, Y-, Z-output_body and X-, Y-, Z-moments.
    '''

    def angles_velocities(Y):

        '''
        1. Horizontal tail aerodynamic velocity vector at the specified Y-position is calculated.
        2. Alpha, beta and the magnitude_velocity of the horizontal tail aerodynamic velocity vector are calculated.
        :param Y: Y-position along the horizontal tail (from - horizontal tail span / 2 to horizontal tail span / 2).
        :return: alpha, beta and the magnitude_velocity of the horizontal tail aerodynamic velocity vector.
        '''

        # (1/2)

        X_dot_a_h = np.cos(aircraft.sweep_h) * aircraft.velocity_vector_a[0, 0] + yaw_dot_derivative_horizontal(Y) \
                    + (Y / abs(Y)) * np.sin(aircraft.sweep_h) * aircraft.velocity_vector_a[1, 0]
        Y_dot_a_h = np.cos(aircraft.sweep_h) * aircraft.velocity_vector_a[1, 0] \
                    - (Y / abs(Y)) * np.sin(aircraft.sweep_h) * aircraft.velocity_vector_a[0, 0]
        Z_dot_a_h = aircraft.velocity_vector_a[2, 0] + roll_dot_derivative_horizontal(np.cos(aircraft.sweep_h) * Y) \
                    + np.sin(aircraft.dihedral_h) * Y_dot_a_h * (Y / abs(Y)) + pitch_dot_derivative_horizontal(Y)

        # (2/2)

        alpha_h = np.arctan(Z_dot_a_h / X_dot_a_h)
        beta_h = np.arctan(Y_dot_a_h / X_dot_a_h)
        V_h = np.sqrt(X_dot_a_h ** 2 + Y_dot_a_h ** 2 + Z_dot_a_h ** 2)

        return alpha_h, beta_h, V_h

    def forces_moments(Y):

        '''
        1. All required aerodynamic and geometric parameters are calculated using above formulas.
        2. Lift and drag components are calculated.
        3. Lift and drag contributions to X-, Y-, Z-output_body are calculated.
        4. Lift and drag contributions to X-, Y-, Z-moments are calculated.
        5. Local alpha and beta input are appended to local list to find (absolute) maximum per time step.
        :param Y: Y-position along the horizontal tail (from - horizontal tail span / 2 to horizontal tail span / 2).
        :return: horizontal tail X-, Y-, Z-output_body and X-, Y-, Z-moments at the specified Y-position.
        '''

        # (1/5)

        alpha_h, beta_h, V_h = angles_velocities(Y)
        c_h, X_ac_h = aircraft.chord_horizontal(Y)
        l_X_h, l_Y_h, l_Z_h = aircraft.moment_arms_horizontal(Y)
        C_L_h, C_D_h, C_M_ac_h = aircraft.coefficients_horizontal(alpha_h, Y)

        # (2/5)

        delta_L_h = 0.5 * aircraft.rho * np.cos(beta_h) * V_h ** 2 * C_L_h * c_h
        delta_D_h = 0.5 * aircraft.rho * np.cos(beta_h) * V_h ** 2 * C_D_h * c_h

        # (3/5)

        delta_F_X_h = np.cos(aircraft.sweep_h) * np.sin(alpha_h) * delta_L_h - np.cos(aircraft.sweep_h) * np.cos(alpha_h) * delta_D_h \
                      + np.sin(aircraft.dihedral_h) * np.sin(aircraft.sweep_h) * np.cos(alpha_h) * delta_L_h \
                      + np.sin(aircraft.dihedral_h) * np.sin(aircraft.sweep_h) * np.sin(alpha_h) * delta_D_h

        delta_F_Y_h = (Y / abs(Y)) * np.sin(aircraft.sweep_h) * np.sin(alpha_h) * delta_L_h \
                      - (Y / abs(Y)) * np.sin(aircraft.sweep_h) * np.cos(alpha_h) * delta_D_h \
                      - (Y / abs(Y)) * np.sin(aircraft.dihedral_h) * np.cos(aircraft.sweep_h) * np.cos(alpha_h) * delta_L_h \
                      - (Y / abs(Y)) * np.sin(aircraft.dihedral_h) * np.cos(aircraft.sweep_h) * np.sin(alpha_h) * delta_D_h

        delta_F_Z_h = - (np.cos(aircraft.dihedral_h) * np.cos(alpha_h) * delta_L_h + np.sin(alpha_h) * delta_D_h)

        # (4/5)

        delta_M_X_h = delta_F_Z_h * l_Y_h + delta_F_Y_h * l_Z_h
        delta_M_Y_h = delta_F_Z_h * l_X_h
        delta_M_Z_h = - (delta_F_X_h * l_Y_h)

        # (5/5)

        local_alpha_h_tab.append(alpha_h * (180 / np.pi))
        local_beta_h_tab.append(beta_h * (180 / np.pi))

        return delta_F_X_h, delta_F_Y_h, delta_F_Z_h, delta_M_X_h, delta_M_Y_h, delta_M_Z_h

    # These steps are required to ensure horizontal tail analysis remains symmetric after applying sweep.
    effective_span = aircraft.b_h / np.cos(aircraft.sweep_h)
    rounded_effective_span = round(effective_span / aircraft.ds_h) * aircraft.ds_h

    # Riemann sum along the wingspan.

    steps = round(effective_span / aircraft.ds_h)
    F_X_h, F_Y_h, F_Z_h, M_X_h, M_Y_h, M_Z_h = 0, 0, 0, 0, 0, 0

    for i in range(steps):

        Y = (- rounded_effective_span / 2) + (aircraft.ds_h / 2) + i * aircraft.ds_h
        Y = round(Y, 3)

        # Horizontal tail cannot be evaluated at Y equal to zero after applying sweep due to asymmetry otherwise.
        if Y != 0:
            delta_F_X_h, delta_F_Y_h, delta_F_Z_h, delta_M_X_h, delta_M_Y_h, delta_M_Z_h = forces_moments(Y)
        else:
            delta_F_X_h, delta_F_Y_h, delta_F_Z_h, delta_M_X_h, delta_M_Y_h, delta_M_Z_h = 0, 0, 0, 0, 0, 0

        F_X_h += delta_F_X_h * aircraft.ds_w
        F_Y_h += delta_F_Y_h * aircraft.ds_w
        F_Z_h += delta_F_Z_h * aircraft.ds_w

        M_X_h += delta_M_X_h * aircraft.ds_w
        M_Y_h += delta_M_Y_h * aircraft.ds_w
        M_Z_h += delta_M_Z_h * aircraft.ds_w

    return F_X_h, F_Y_h, F_Z_h, M_X_h, M_Y_h, M_Z_h

# (3/7) equations of motion (6 DOF) for the body magnitude_reference frame.

def equations_of_motion():

    '''
    1. Forces and moments from wing and vertical tail are calculated.
    2. Total output_body and moments are calculated.
    3. Translational accelerations are calculated.
    4. Angular accelerations are calculated.
    :return: translational and angular accelerations.
    '''

    # (1/4)

    F_X_w, F_Y_w, F_Z_w, M_X_w, M_Y_w, M_Z_w = wing()
    F_X_v, F_Y_v, F_Z_v, M_X_v, M_Y_v, M_Z_v = vertical()

    if aircraft.horizontal:
        F_X_h, F_Y_h, F_Z_h, M_X_h, M_Y_h, M_Z_h = horizontal()
    else:
        F_X_h, F_Y_h, F_Z_h, M_X_h, M_Y_h, M_Z_h = 0, 0, 0, 0, 0, 0

    # (2/4)

    F_X = F_X_w + F_X_v + F_X_h - np.sin(aircraft.pitch) * aircraft.W
    F_Y = F_Y_w + F_Y_v + F_Y_h + np.sin(aircraft.roll) * aircraft.W
    F_Z = F_Z_w + F_Z_v + F_Z_h + np.cos(aircraft.pitch) * aircraft.W

    M_X = M_X_w + M_X_v + M_X_h
    M_Y = M_Y_w + M_Y_v + M_Y_h
    M_Z = M_Z_w + M_Z_v + M_Z_h

    # (3/4)

    X_dot_dot_b = F_X / aircraft.m
    Y_dot_dot_b = F_Y / aircraft.m
    Z_dot_dot_b = F_Z / aircraft.m

    # (4/4)

    roll_dot_dot = M_X / aircraft.I_XX
    pitch_dot_dot = M_Y / aircraft.I_YY
    yaw_dot_dot = M_Z / aircraft.I_ZZ

    return X_dot_dot_b, Y_dot_dot_b, Z_dot_dot_b, roll_dot_dot, pitch_dot_dot, yaw_dot_dot

# (4/7) transformation matrices for the Earth and body magnitude_reference frames.

def T_X(roll):

    '''
    :param roll: roll angle of the aircraft.
    :return: transformation matrix around the X-axis.
    '''

    T_x = np.array([[1, 0, 0],
                    [0, np.cos(roll), np.sin(roll)],
                    [0, - np.sin(roll), np.cos(roll)]])
    return T_x

def T_Y(pitch):

    '''
    :param pitch: pitch angle of the aircraft.
    :return: transformation matrix around the Y-axis.
    '''

    T_y = np.array([[np.cos(pitch), 0, - np.sin(pitch)],
                    [0, 1, 0],
                    [np.sin(pitch), 0, np.cos(pitch)]])
    return T_y

def T_Z(yaw):

    '''
    :param yaw: yaw angle of the aircraft.
    :return: transformation matrix around the Z-axis.
    '''

    T_z = np.array([[np.cos(yaw), np.sin(yaw), 0],
                    [- np.sin(yaw), np.cos(yaw), 0],
                    [0, 0, 1]])
    return T_z

def T_be():

    '''
    Transformation matrix from the Earth magnitude_reference frame to the body magnitude_reference frame.
    :return: None (aircraft class modified inside function).
    '''

    velocity_e_1 = np.dot(T_X(aircraft.roll), aircraft.velocity_vector_e)
    velocity_e_2 = np.dot(T_Y(aircraft.pitch), velocity_e_1)
    aircraft.velocity_vector_b = np.dot(T_Z(aircraft.yaw), velocity_e_2)

    return

def T_eb():

    '''
    Transformation matrix from the body magnitude_reference frame to the Earth magnitude_reference frame.
    :return: None (aircraft class modified inside function).
    '''

    velocity_b_1 = np.dot(T_Z(- aircraft.yaw), aircraft.velocity_vector_b)
    velocity_b_2 = np.dot(T_Y(- aircraft.pitch), velocity_b_1)
    aircraft.velocity_vector_e = np.dot(T_X(- aircraft.roll), velocity_b_2)

    return

# (5/7) aerodynamic velocity vector function.

def gusts():

    '''
    Aerodynamic velocity vector is calculated using the body velocity vector and (randomly generated) gust velocity vector.
    :return: None (aircraft class modified inside function).
    '''

    aircraft.t = round(aircraft.t, 3)

    if aircraft.t % 1 == 0 and aircraft.gust == True:

        aircraft.X_dot_gust = np.random.normal() * aircraft.gust_intensity
        aircraft.Y_dot_gust = np.random.normal() * aircraft.gust_intensity
        aircraft.Z_dot_gust = np.random.normal() * aircraft.gust_intensity

    aircraft.gust_velocity_vector = np.array([[aircraft.X_dot_gust],
                                              [aircraft.Y_dot_gust],
                                              [aircraft.Z_dot_gust]])

    aircraft.velocity_vector_a = aircraft.velocity_vector_b + aircraft.gust_velocity_vector

    return

# (6/7) simulation using data lists.

t_tab = []

alpha_w_tab = []
alpha_v_tab = []
alpha_h_tab = []

beta_w_tab = []
beta_v_tab = []
beta_h_tab = []

roll_tab = []
pitch_tab = []
yaw_tab = []

X_dot_e_tab = []
Y_dot_e_tab = []
Z_dot_e_tab = []

while aircraft.t < aircraft.t_simulation:

    '''
    1. Start with Earth magnitude_reference frame vector and transform to body magnitude_reference frame vector.
    2. Check if gusts are applicable and reset local lists.
    3. Apply translational accelerations.
    4. Transform back to Earth magnitude_reference frame.
    5. Apply angular accelerations.
    This order ensures that translational and angular accelerations are independent.
    6. Append data to lists used for plotting (first step skipped due to irregularities at start)
    Repeat for the simulation duration.
    '''

    # (1/6)

    T_be()

    # (2/6)

    gusts()

    local_alpha_w_tab = []
    local_alpha_v_tab = []
    local_alpha_h_tab = []

    local_beta_w_tab = []
    local_beta_v_tab = []
    local_beta_h_tab = []

    # (3/6)

    aircraft.X_dot_dot_b, aircraft.Y_dot_dot_b, aircraft.Z_dot_dot_b, aircraft.roll_dot_dot, \
    aircraft.pitch_dot_dot, aircraft.yaw_dot_dot = equations_of_motion()

    aircraft.velocity_vector_b[0] += aircraft.X_dot_dot_b * aircraft.dt
    aircraft.velocity_vector_b[1] += aircraft.Y_dot_dot_b * aircraft.dt
    aircraft.velocity_vector_b[2] += aircraft.Z_dot_dot_b * aircraft.dt

    aircraft.roll_dot += aircraft.roll_dot_dot * aircraft.dt
    aircraft.pitch_dot += aircraft.pitch_dot_dot * aircraft.dt
    aircraft.yaw_dot += aircraft.yaw_dot_dot * aircraft.dt

    # (4/6)

    T_eb()

    # (5/6)

    aircraft.roll += aircraft.roll_dot * aircraft.dt
    aircraft.pitch += aircraft.pitch_dot * aircraft.dt
    aircraft.yaw += aircraft.yaw_dot * aircraft.dt

    # (6/6)

    if aircraft.t > aircraft.dt:

        t_tab.append(aircraft.t)

        #alpha_w_tab.append(np.max(np.abs(local_alpha_w_tab)))
        alpha_v_tab.append(np.max(np.abs(local_alpha_v_tab)))

        #beta_w_tab.append(np.max(np.abs(local_beta_w_tab)))
        beta_v_tab.append(np.max(np.abs(local_beta_v_tab)))

        if aircraft.horizontal:

            alpha_h_tab.append(np.max(np.abs(local_alpha_h_tab)))
            beta_h_tab.append(np.max(np.abs(local_beta_h_tab)))

        roll_tab.append(aircraft.roll * (180 / np.pi))
        pitch_tab.append(aircraft.pitch * (180 / np.pi))
        yaw_tab.append(aircraft.yaw * (180 / np.pi))

        X_dot_e_tab.append(aircraft.velocity_vector_e[0])
        Y_dot_e_tab.append(aircraft.velocity_vector_e[1])
        Z_dot_e_tab.append(aircraft.velocity_vector_e[2])

    aircraft.t += aircraft.dt

end = time.time()

# (7/7) print (absolute) maxima and plot results.

print("\n", "Max. alpha wing: ", np.round(np.max(alpha_w_tab), 1))
print("Max. beta wing: ", np.round(np.max(beta_w_tab), 1))

print("\n", "Max. alpha vertical tail: ", np.round(np.max(alpha_v_tab), 1))
print("Max. beta vertical tail: ", np.round(np.max(beta_v_tab), 1))

if aircraft.horizontal:

    print("\n", "Max. alpha horizontal tail: ", np.round(np.max(alpha_h_tab), 1))
    print("Max. beta horizontal tail: ", np.round(np.max(beta_h_tab), 1))

print("\n", "Max. pitch: ", np.round(np.max(pitch_tab), 1))
print("Max. roll: ", np.round(np.max(roll_tab), 1))
print("Max. yaw: ", np.round(np.max(yaw_tab), 1))

print("\n", "Calculation time: ", round(end - begin, 2), "seconds.")

plt.subplot(1, 3, 1)
plt.xlabel("Time [s]")
plt.ylabel("Angle [deg]")
line1, = plt.plot(t_tab, alpha_w_tab, label = 'Alpha wing')
line2, = plt.plot(t_tab, beta_w_tab, label = 'Beta wing')

if aircraft.horizontal:

    line3, = plt.plot(t_tab, alpha_h_tab, label = 'Alpha horizontal')
    plt.legend(handles = [line1, line2, line3], loc = 'lower right')

else:
    plt.legend(handles=[line1, line2], loc='lower right')


plt.subplot(1, 3, 2)
plt.xlabel("Time [s]")
plt.ylabel("Angle [deg]")
line5, = plt.plot(t_tab, roll_tab, label = 'Roll')
line6, = plt.plot(t_tab, pitch_tab, label = 'Pitch')
line7, = plt.plot(t_tab, yaw_tab, label = 'Yaw')
plt.legend(handles = [line5, line6, line7], loc = 'lower right')

plt.subplot(1, 3, 3)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
line8, = plt.plot(t_tab, X_dot_e_tab, label = 'X dot')
line9, = plt.plot(t_tab, Y_dot_e_tab, label = 'Y dot')
line10, = plt.plot(t_tab, Z_dot_e_tab, label = 'Z dot')
plt.legend(handles = [line8, line9, line10], loc = 'lower right')

plt.suptitle("Stability analysis")
plt.show()
