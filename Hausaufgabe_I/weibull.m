clear all;
clear all;
%% Weibull Distribution 7.5 Mean
clear all;
disp(' Weibull Distribution 7.5 Mean')
subplot(221)
sample = wblrnd(8.51, 2, 8760, 1);
mean = mean(sample);
disp(mean)
pd = histfit(sample,25, 'wbl');
x_values = (0:1:25)';
hist = histogram(sample,26);
y_values = hist.Values';
subplot(222)
plot(x_values,y_values,'LineWidth',2)
% Data vesgas
x_vestas = (0:25)';
y_vestas = [0, 0,0,0,91,200,362,588,889,1255,1604,1769,1798,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800]';
subplot(223)
plot(x_vestas, y_vestas);
x_vestas;
for i = 1:length(y_values)
    y_vestas_yield(i,1) = y_values(i)*y_vestas(i);
end
subplot(224)
plot(x_vestas, y_vestas_yield);
netyield_vestas = sum(y_vestas_yield)*0.97*0.98*0.9;
clear sample;
clear mean;
%% Weibull Distribution 8.5
disp(' Weibull Distribution 8.5 Mean')
subplot(221)
sample = wblrnd(9.55, 2, 8760, 1);
mean = mean(sample);
disp(mean);
pd = histfit(sample,24, 'wbl');
x_values = (0:1:24)';
hist = histogram(sample,25);
y_values = hist.Values;
subplot(222)
plot(x_values,y_values,'LineWidth',2)

x_nordex = (0:24)';
y_nordex = [0,0,0,0,80,206,378,596,894,1250,1674,2064,2374,2454,2500,2511,2523,2534,2534,2534,2534,2523,2523,2500,2500]';
subplot(223)
plot(x_nordex, y_nordex);
for i = 1:length(y_values)
    y_nordex_yield(i,1) = y_values(i)*y_nordex(i);
end
subplot(224)
plot(x_nordex, y_nordex_yield);
netyield_nordex =sum(y_nordex_yield)*0.97*0.98*0.9;
clear sample;
clear mean;
%% Enercon; Weibull 9.5
disp(' Weibull Distribution 9.5 Mean')
subplot(221)
sample = wblrnd(11.3, 2, 8760, 1);
mean = mean(sample);
disp(mean)
pd = histfit(sample,25, 'wbl');
x_values = (0:1:25)';
hist = histogram(sample,26);
y_values = hist.Values;
subplot(222)
plot(x_values,y_values,'LineWidth',2)

x_enercon = (0:25)';
y_enercon = [0,0,3,25,82,174,321,532,815,1180,1580,1900,2200,2400,2480,2700,2850,2950,3020,3020,3020,3020,3020,3020,3020,3020]';
subplot(223)
plot(x_enercon, y_enercon);
for i = 1:length(y_values)
    y_enercon_yield(i,1) = y_values(i)*y_enercon(i);
end
subplot(224)
plot(x_enercon, y_enercon_yield);
netyield_enercon =sum(y_enercon_yield)*0.97*0.98*0.9;
%% Calculation
full_load_hours_vestas = netyield_vestas/1800;
full_load_hours_nordex = netyield_nordex/2500;
full_load_hours_enercon = netyield_enercon/3000;

capacity_vestas = netyield_vestas/(1800*8760);
capacity_nordex = netyield_nordex/(2500*8760);
capacity_enercon = netyield_enercon/(3000*8760);