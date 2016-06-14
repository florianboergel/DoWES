# Loading needed modules
import numpy as np
from pylab import *


# This line configures matplotlib to show figures embedded in the notebook,
# instead of opening a new window for each figure.
# %matplotlib inline

# Enabling zooming in notebook
#import mpld3
#mpld3.enable_notebook()

rc('text', usetex = True)
rc('font', size = 16)


# Setting precision for printing arrays
np.set_printoptions(precision=3)

D = 62.18*2           # RotorDiameter
N = 3                # Number of blades
RotorDiameter = 0.95 # Electrical conversion efficiency
Ud = 10.0            # Design wind speed in [m/s]
lambdad = 8.20        # Design tip-speed ratio
cld = 1.345          # Design lift coefficient selected for max. gliding factor
rho = 1.225          # Air density [kg/m^3] @ standard conditions

Ne = 10              # Number of elements
R = D / 2.0          # Calculate the radius of the turbine

# Create an array of radii for each of the elements
# 0 ---be1---|---be2---|---be3---| ... |---be10--->

# Calculating the size of each element assuming that the
# first blade element begins at the hub (r = 0) and the last
# ends at the tip of the blade (r = R)

BladeRootOffset = 1.25  #  Defining radial position at which the first element begins

BladeElementLength = (R - BladeRootOffset) / 8

# Calculating the radial position of the center of each blade element
BladeElementRadii = np.linspace(BladeRootOffset + BladeElementLength / 2,
                                R - BladeElementLength / 2,
                                Ne)

# We calculate first the angles
alpha_Ad_deg = 10.0
alpha_Ad    = alpha_Ad_deg * np.pi/180

# Iterating over all blade elements
# and storing results in arrays:
# - alpha_twist
# - alpha_twist_deg
alpha_twist = np.array([])     # Initializing empty array to store angular data
alpha_twist_deg = np.array([]) # Initializing empty array to store angular data
alpha_1 = np.array([])     # Initializing empty array to store angular data

for BladeElementRadius in BladeElementRadii:
    alpha_1i = np.arctan(R/(lambdad * BladeElementRadius))
    alpha_i = 2.0/3.0 * alpha_1i
    alpha_twisti = alpha_i - alpha_Ad
    # Storing results in arrays
    alpha_1     = np.append(alpha_1, alpha_1i)
    #alpha       = np.append(alpha, alpha_i)
    alpha_twist = np.append(alpha_twist, alpha_twisti)
    alpha_twist_deg = np.append(alpha_twist_deg, alpha_twisti*180 / np.pi)

print 'Twist angles : '
print alpha_twist_deg

# Plotting alpha twist
plot(BladeElementRadii, alpha_twist_deg, '-o')
xlabel('Blade span ([m])')
ylabel(r'Twist angle ([deg]) $\alpha_{twist}$')

# Constructing array of chords

Nblades = 3  # Number of blades


# Iterating over all blade elements
# and storing results in array
# - chords
chords = np.array([])     # Initializing empty array to store chord data

for i,BladeElementRadius in enumerate(BladeElementRadii):

    chord_i = 16.0 * np.pi * BladeElementRadius / (Nblades * cld) * np.square(np.sin(1.0/3*alpha_1[i]))

    # Storing results in array
    chords = np.append(chords, chord_i)
    # Printing
print 'Blade element radius [m]\tChords [m]\t\t\tTwist[]'
for i,BladeElementRadius in enumerate(BladeElementRadii):
    print ' %0.2f' % BladeElementRadius + '\t\t\t\t %0.2f' % chords[i] + '\t\t\t\t %0.2f' % alpha_twist_deg[i]

# Plotting chord
plot(BladeElementRadii, chords, '-o')
xlabel('Blade span ([m])')
ylabel(r'Twist angle ([deg]) $\alpha_{twist}$')

# Overwriting the chord with Davide's results
#chords = np.array([13.59, 8.22, 5.46, 4.03, 3.18, 2.62, 2.22, 1.93, 1.71, 1.53])
#alpha_twist = np.array([52.1, 22.2, 10.7, 5.0, 1.8, -0.4, -1.8, -2.9, -3.8, -4.5]) * pi /180

# Plotting chord
plot(BladeElementRadii, chords, '-o')
xlabel('Blade span ([m])')
ylabel(r'Twist angle ([deg]) $\alpha_{twist}$')

# Loading data
profiledata = np.loadtxt(open("NACA-64-415.dat","rb"),delimiter="   ")

# Plotting profile
#plot(profiledata[:,0], profiledata[:,1])
#xlabel('Angle of attack [deg]')
#ylabel(r'Lift coefficient $c_l$')

from scipy.interpolate import interp1d

def GetLiftCoefficient(alpha, profiledata) :
    '''

        Summary:

        This function interpolates the lift coefficient
        for a given angle of attack. The profile data
        are provided as a 2D numpy array. The interpolation
        is performed with a cubic spline.


        Input:


        alpha

        Angle at which the interpolation should be done

        profiledata

        This function takes a table of blade profile
        data stored in a numpy 2D array with columns
        defined as follows:

        |Angle of attack [deg]| Lift coefficient | ...


    '''

    # Converting to degrees since our table is given in degrees
    alpha_deg = alpha * 180/pi

    x = profiledata[:,0]
    y = profiledata[:,1]

    f = interp1d(x, y, kind='cubic')

    return f(alpha_deg)


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
