% Script: Plot Fast Timeseries
% Created to visualize aeroelastic simulations created with FAST and used  
% in DOWES class (SS2014)
% by Luis, rev 0.0

% Matlab functions required:
% 'select_file_from_folder', 'fast2mat', 'fast2plot' & 'dequiv'

%% Main code

% add path with matlab functions
addpath('C:/Fastextreme/MatlabFunctions');

% select FAST output file
FAST_output_file_name  = select_file_from_folder('C:\Fastextreme\*.out');

% create structure in workspace with results
FAST_results = fast2mat(FAST_output_file_name);

% create folder with timeseries per each variable
create_figures = fast2plot(FAST_results);

% % calculate damage equivalent loads
% Damage_Eq_Loads = struct;
% % blade root flapwise bending moment (Whöler coefficient = 10)
% [Damage_Eq_Loads.RootMFlp_Req R]= dequiv(FAST_results.Time, FAST_results.RootMFlp2,10);
% Damage_Eq_Loads.RootMFlp_m = 10;
% % blade root edgewise bending moment (Whöler coefficient = 10)
% [Damage_Eq_Loads.RootMedg_Req R]= dequiv(FAST_results.Time, FAST_results.RootMedg2,10);
% Damage_Eq_Loads.RootMedg_m = 10;
% % tower base longitudinal (fore-aft) bending moment (Whöler coefficient = 4)
% [Damage_Eq_Loads.TwrBsMyt_Req R]= dequiv(FAST_results.Time, FAST_results.TwrBsMyt,4);
% Damage_Eq_Loads.TwrBsMyt_m = 4;
% % tower base transversal (side-side) bending moment (Whöler coefficient = 4)
% [Damage_Eq_Loads.TwrBsMxt_Req R]= dequiv(FAST_results.Time, FAST_results.TwrBsMxt,4);
% Damage_Eq_Loads.TwrBsMxt_m = 4;
% clear R

% %% power spectrum density
% % uncomment this section to analyse PSD
% 
% % sensor (load value) to analyse
% data = FAST_results.RootMFlp2; % blade root flapwise bending moment
% %data = FAST_results.TwrBsMyt; % tower bottom fore-aft bending moment 
% 
% % plot power spectrum density
% figure1 = figure;
% axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
%     'YMinorGrid','on','YGrid','on','XScale','log','XMinorTick','on',...
%     'XMinorGrid','on','XGrid','on');
% 
% % select the power spectrum
% Hs = spectrum.periodogram;
% % sampling frequency
% sampling = 100000; 
% psd(Hs,data,'Fs',sampling) 
% % update graph for better display
% hline = findobj(gcf, 'type', 'line');
% set(hline,'Marker','.','LineStyle','none')
% 
% % axis in log format
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% 
% % scale axis
% xlim(axes1,[0 20]);
% ylim(axes1,[0.1 100]);
% 
% % labels and title
% xlabel('Frequency (Hz)');
% ylabel('Power/frequency (dB/Hz)');
% title('Power Spectral Density');