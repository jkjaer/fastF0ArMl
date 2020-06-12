clear
clc
close all

load('results/timingVsPitchOrderIterative.mat');

timing2plot = squeeze(min(timingMtx,[],2));

semilogy(maxPitchOrderList', timing2plot, 'lineWidth',2)
xlabel('Max. pitch order');
ylabel('Computation time [s]');
grid on;
legend(methods)

% store data in text file so that it is ready for plotting using tikz
comment = '';
legend_str = ['order ', strjoin(methods)];
data = [maxPitchOrderList', timing2plot];
mtx_to_tbl_pgfplots('results/timingVsPitchOrderIterative',comment,legend_str,data);
