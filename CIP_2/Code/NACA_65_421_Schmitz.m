%% Blade Design according to Schmitz NACA 65-421


D = 62.18*2          % RotorDiameter
N = 3                % Number of blades
RotorDiameter = 0.95 % Electrical conversion efficiency
Ud = 10.0            % Design wind speed in [m/s]
lambdad = 8.20       % Design tip-speed ratio
cld = 1.255         % Design lift coefficient selected for max. gliding factor
rho = 1.225        % Air density [kg/m^3] @ standard conditions
alpha_Ad_deg = 11.0  % Angle of Attack Design
alpha_Ad    = alpha_Ad_deg * pi/180

Ne = 10              % Number of elements
R = D / 2.0          % Calculate the radius of the turbine

RotationalSpeed = Ud * lambdad / R
BladeRootOffset = 1.25  %  Defining radial position at which the first element begins
BladeElementLength = (R - BladeRootOffset) / 8

% Calculate BladeElementRadii
BladeElementRadii(1) = BladeRootOffset + BladeElementLength / 2
for i = 2:Ne-1
   BladeElementRadii(i) = BladeElementRadii(i-1) + BladeElementLength
end
BladeElementRadii(Ne-1) = BladeElementRadii(Ne-2) + BladeElementLength/2

% Calcualte Twist-Angle and Chordlength
count = 1
for BladeElementRadius = BladeElementRadii
    alpha_li    =   atan(R/(lambdad * BladeElementRadius))
    alpha_i     =   2.0/3.0 * alpha_li
    alpha_twisti=   alpha_i - alpha_Ad
    % Storing results in arrays
    alpha_1(count) = alpha_li
    %alpha = np.append(alpha, alpha_i)
    alpha_twist(count) = alpha_twisti
    alpha_twist_deg(count) =  alpha_twisti*180 / pi
    chord_i = 16.0 * pi * BladeElementRadius / (3 * cld) * (sin(1.0/3*alpha_li))^2
    chord(count) = chord_i
    count = count +1;
end

% BEM Algorithm

alpha_pitch = 0
count_indice = 1

for r = BladeElementRadii
    %step 1
    a = 0;
    a_old = 1000;
    a_Dash = 0;
    while (count==0 || abs(a-a_old) > 0.001)
        a_old = a;
        %step 2, step 3
        angleOfAttack = 10;
        %alpha = atan((v1*(1-a))/(omega*r*(1+a_Dash)))*180/pi
        alpha = atan(Ud*(1-a)/(RotationalSpeed*r));
        % evaluate Angle of Attack
        alpha_A_tmp = alpha - (alpha_twist(count_indice) + alpha_pitch)

       % disp(alpha*180/pi);

        %step 4
        interpolationTable = xlsread('Interpol.xlsx','NACA 65-421','B6:D81');
        C_L = interp1(interpolationTable(:,1),interpolationTable(:,2),alpha_A_tmp*180/pi,'spline');
        C_D = interp1(interpolationTable(:,1),interpolationTable(:,3),alpha_A_tmp*180/pi,'spline');

        %step 5
        C_T = C_L * cos(alpha) % + C_D *sin(alpha);
        C_Q = C_L * cos(alpha) - C_D *sin(alpha);

        %step 6
        %c = 16*pi*r/(N*C_L)*(sin(1/3*alpha))^2;
        sigma = N*chord(count_indice)/(2*pi*r);
        a_Dash_new = 1/((4*sin(alpha) * cos(alpha)/(sigma*C_Q)-1)) ;
        a = 1/(4*sin(alpha) * sin(alpha)/(sigma*C_T)+1);

        %step 7
        a_c = 0.2;
        K = 4*sin(alpha)^2/(sigma*C_T);
        a=0.5*(2+K*(1-2*a_c)-sqrt((K*(1-2*a_c)+2)^2+4*(K*a_c^2-1)));
        count = count+1;
    end;
    count_indice = count_indice + 1;
end;






