%% Task 1 Wind Speed vs hours per year
cut_in_ws = 3.5
cut_out_ws = 25
rated_ws = 11
data = wblpdf(1:0.5:30,9,2)
data = data*8760

figure;
plot(1:0.5:30,data)
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
for i= 1:0.5:30
    data_power(j) = power(i);
    if data_power(j) > 3500000
        data_power(j) = 3500000;
    end
    j = j +1;
end

figure;
plot(1:0.5:30, data_power)
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
ylabel('Power in [kW]')
saveas(gcf,'power_curve.png')
hold off

%% Energy-Yield
for i= 1:length(data_power)
    energy_yield(i) = data_power(i)/1000*data(i);
end
figure;
plot(1:0.5:30, energy_yield)
title('Energy yield')
xlabel('Wind speed in [m/s]')
ylabel('Energy in [kWh/a]')
saveas(gcf, 'energy_yield.png')

total_energy_yield = sum(energy_yield);
disp('Sum Energy yield in [kWh/a]')
disp(total_energy_yield)
disp('Money 100 % available 0.08 € / kWh')
disp(total_energy_yield*0.08)
disp('Money 95 % available 0.08 € / kWh')
disp(total_energy_yield*0.08*0.95)
disp('Difference in Money')
disp(total_energy_yield*0.08-total_energy_yield*0.08*0.95)

%% Calculate Turbulence Intensity for different wind speeds
I_ref_b = 0.14;

count = 1;
for i = [5 10 15 25]
    sigma(count) = I_ref_b*(0.75*i+5.6);
    turbulence_intensity(count) = sigma(count)/i;
    count = count + 1 ;
end

disp('Turbulence Intensity')
disp(turbulence_intensity)
