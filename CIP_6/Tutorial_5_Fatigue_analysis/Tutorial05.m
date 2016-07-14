%% Convert signal to peaks

% Loading data 
data = load('LoadsTimeSeries.csv');
FlapwiseBendingMoment = data(:,3);



% Taking original time series and simplifying it
% to a time series of direction changes
[sigpeaks_FlapwiseBendingMoment, exttime] = sig2ext(FlapwiseBendingMoment);

figure;
plot(FlapwiseBendingMoment, '-b');
hold on
plot(exttime, sigpeaks_FlapwiseBendingMoment, '.r');



%% Let us make our rain flow counting

rfc_FlapwiseBendingMoment = rainflow(sigpeaks_FlapwiseBendingMoment);


%% Let us show off with the Markov matrix

figure;
rfmatrix(rfc_FlapwiseBendingMoment);
[markovmat, mx, my] = rfmatrix(rfc_FlapwiseBendingMoment);