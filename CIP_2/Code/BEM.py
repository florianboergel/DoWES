# Only BEM

# Defining initial conditions

# Rotational speed of the rotor is constant and given by design
RotationalSpeed = Ud *lambdad / R

# Initializing vector with induction factors
a   = np.zeros(Ne)
ap  = np.zeros(Ne)
print a

# Defining tolerance for convergence
a_tolerance = 0.001

# Iterating over each ring
#for i,BladeElementRadius in enumerate(BladeElementRadii):

# Let us do it for ring number 3
RingNumber = 3
BladeElementRadius = BladeElementRadii[RingNumber-1]


# 1- Initialize induction factors
a_tmp1  = 0
ap_tmp1 = 0
iter = 1
MaxIter = 10
HasConverged = False

alpha_pitch = 0

while ( (HasConverged == False) & (iter <= MaxIter) ):

    print 'Iteration: %d' % iter
    print 'Induction value (guess): %0.3f' % a_tmp1

    # 2- Evaluate inflow angle

    # Local tip-speed ratio
    lambdal = (RotationalSpeed * BladeElementRadius)/Ud
    print 'Local tip-speed ratio: %0.3f' % lambdal

    # Slig%rhtly reformulated equation as in the lecture 04
    alpha_tmp = arctan( 1.0/lambdal * (1 - a_tmp1)/(1 - ap_tmp1) )
    print 'Inflow angle: %0.3f deg' % (alpha_tmp * 180/pi)

    # 3- Evaluating angle of attack
    alpha_A_tmp = alpha_tmp - (alpha_twist[i] + alpha_pitch)
    print 'Angle of attack: %0.3f deg' % (alpha_A_tmp * 180/pi)

    # 4- Evaluate the lift coefficient
    cl_tmp = GetLiftCoefficient(alpha_A_tmp, profiledata)
    print 'Lift coefficient: %0.3f' % (cl_tmp)

    # 5- Compute Cthrust
    ct_tmp = cl_tmp * cos(alpha_tmp)
    print 'Thrust coefficient: %0.3f' % (ct_tmp)

    # 6- Evaluate induction value
    sigma = Nblades*chords[RingNumber-1]/(2.0*pi*BladeElementRadius)
    K = 4.0 * square(sin(alpha_tmp)) / (sigma * ct_tmp)
    print 'Spera factor K: %0.3f' % K


    a_tmp2 = 1 / (K +1)
    print 'Induction value: %0.3f' % a_tmp2

    # Applying Spera correction if needed
    ac = 0.2
    if a_tmp2 >= ac :
        a_tmp2_prev = a_tmp2   # Storing value for printing log
        a_tmp2 = 0.5 * (2 + K * (1 - 2*ac) - sqrt(square((K*(1-2*ac))+2)+4*(K*square(ac)-1)))
        print 'Correcting induction %0.3f -> %0.3f' % (a_tmp2_prev, a_tmp2)

    # 7- Checking convergence
    delta_a = a_tmp2 - a_tmp1
    if abs(delta_a) < a_tolerance:
        print 'Induction has converged with delta : %0.3f' % delta_a
        HasConverged = True
        a[RingNumber-1] = a_tmp2
        break
    else:
        # We follow the iteration with the new induction factor
        a_tmp1 = a_tmp2


    print '---------------'
