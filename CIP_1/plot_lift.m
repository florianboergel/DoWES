% Load Excelsheet
NACA_65_415 = xlsread('CIP_01_DoWES.Rotor_design_profile_data.SS2016.01.xlsx','NACA 65-415','B6:F81');
NACA_65_421 = xlsread('CIP_01_DoWES.Rotor_design_profile_data.SS2016.01.xlsx','NACA 65-421','B6:F81');
% Plot
figure;
plot(NACA_65_415(19:62,1), NACA_65_415(19:62,5),'-or')
hold on
plot([NACA_65_415(25,1) NACA_65_415(25,1)],[0 floor(max(NACA_65_415(:,5)))], '-r')
text(NACA_65_415(26,1),50,'Maximum')
xlabel('Angle of Attack in [°]')
ylabel('Lift-to-drag Ratio')
legend('NACA\_65\_415')
hold off
saveas(gcf,'lift_to_drag_ratio_415.png')

figure;
plot(NACA_65_421(19:62,1), NACA_65_421(19:62,5),'-ob')
hold on
plot([NACA_65_421(29,1) NACA_65_421(29,1)],[0 floor(max(NACA_65_421(:,5)))])
text(NACA_65_421(30,1),50,'Maximum')
xlabel('Angle of Attack in [°]')
ylabel('Lift-to-drag Ratio')
legend('NACA\_65\_421')
saveas(gcf,'lift_to_drag_ratio_421.png')