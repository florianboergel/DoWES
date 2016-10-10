%% Task 1 Wind Speed vs hours per year
cut_in_ws = 3.5
cut_out_ws = 25
rated_ws = 11
k = 2
vAve=10
A = vAve/sqrt(pi/2)
data = wblpdf(0:0.01:30,9,k)
data = data*8760/100 %must be fitted to bin size 1

figure;
plot(0:0.01:30,data*100)
y1=get(gca,'ylim')
hold on
plot([cut_in_ws cut_in_ws],y1)
text(cut_in_ws,450,'Cut in')
plot([rated_ws rated_ws],y1)
text(rated_ws,450,'Rated wind speed')
plot([cut_out_ws cut_out_ws],y1)
text(cut_out_ws,450,'Cut Out')
title('wind speed vs. hours in a year')
xlabel('Wind speed in [m/s]')
ylabel('Time in [h/a]')
saveas(gcf,'windspeed_vs_hours_year.png')
hold off
%% Power Curve
rho = 1.225
total_eff = 0.4705
R = 54.0

power = @(V) 0.5*rho*total_eff*pi*R^2*V^3;
j=1;
for i= 0:0.01:30
    data_power(j) = power(i);
    if data_power(j) > 3500000
        data_power(j) = 3500000;
    end
    if i<cut_in_ws || i>= cut_out_ws
        data_power(j) = 0;
    end

    j = j +1;
end

figure;
plot(0:0.01:30, data_power)
y1 = [0 3800000];
hold on
plot([cut_in_ws cut_in_ws],y1)
text(cut_in_ws,2000000,'Cut in')
plot([rated_ws rated_ws],y1)
text(rated_ws,2000000,'Rated wind speed')
plot([cut_out_ws cut_out_ws],y1)
text(cut_out_ws,2000000,'Cut Out')
axis([0,30,0,3800000])
title('Power Curve')
xlabel('Wind speed in [m/s]')
ylabel('Power in [W]')
saveas(gcf,'power_curve.png')
hold off

%% Energy-Yield
for i= 1:length(data_power)
    energy_yield(i) = data_power(i)/1000*data(i);
end
figure;
plot(0:0.01:30, energy_yield/10)
title('Energy yield')
xlabel('Wind speed in [m/s]')
ylabel('Energy in [MWh/a]')
saveas(gcf, 'energy_yield.png')

total_energy_yield = sum(energy_yield);
disp('Sum Energy yield in [GWh/a]')
disp(total_energy_yield)
disp('Money 100 % available 0.08 € / kWh')
disp(total_energy_yield*0.08)
disp('Money 95 % available 0.08 € / kWh')
disp(total_energy_yield*0.08*0.95)
disp('Difference in Money')
disp(total_energy_yield*0.08-total_energy_yield*0.08*0.95)

%% Calculate Turbulence Intensity for different wind speeds
I_ref_b = 0.14;

% According to NTM, Slide 6
count = 1;
for i = 1:25
    sigma(count) = I_ref_b*(0.75*i+5.6);
    turbulence_intensity(count) = sigma(count)/i;
    count = count + 1 ;
end

disp('Turbulence Intensity')
disp(turbulence_intensity)
%% Wake Conditions
count = 1;
p_w = 0.06;
N = 1;
di_4 = 4;
di_8 = 8;

%According to Frandsen, Slide 7
for i = 1:25
    sigma_t(count,1) = sqrt(0.9*i^2/(1.5+0.3*di_4*sqrt(i))^2+sigma(count));
    sigma_t(count,2) = sqrt(0.9*i^2/(1.5+0.3*di_8*sqrt(i))^2+sigma(count));
    I_eff_4(count) = 1/i*((1-N*p_w)*sigma(count)^4+p_w*sigma_t(count,1)^4)^(1/4);
    I_eff_8(count) = 1/i*((1-N*p_w)*sigma(count)^4+p_w*sigma_t(count,2)^4)^(1/4);

    count = count + 1 ;
end

figure;
plot(I_eff_4,'-or');
hold on;
plot(I_eff_8,'-ob');
plot(turbulence_intensity,'-og');
xlabel('Windspeed in [m/s]');
ylabel('Turbulence Intensity in [%]');
title('Wake Effects with neighbouring wind turbines');
legend('Wake with D=4', 'Wake with D=8', 'Free Stream');
saveas(gcf, 'wake_effects_5.png')

disp('I_eff_4, for ws 5,10,15,25');
disp(I_eff_4(5));
disp(I_eff_4(10));
disp(I_eff_4(15));
disp(I_eff_4(25));

disp('I_eff_8, for ws 5,10,15,25');
disp(I_eff_8(5));
disp(I_eff_8(10));
disp(I_eff_8(15));
disp(I_eff_8(25));
