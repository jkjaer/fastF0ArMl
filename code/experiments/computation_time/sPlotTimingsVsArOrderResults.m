clear
clc
close all

load('results/timingVsArOrder.mat');

timing2plot = squeeze(min(timingMtx,[],2));

semilogy(maxArOrderList', timing2plot, 'lineWidth',2)
xlabel('Max. AR order');
ylabel('Computation time [s]');
grid on;
legend(methods)

% store data in text file so that it is ready for plotting using tikz
comment = '';
legend_str = ['order ', strjoin(methods)];
data = [maxArOrderList', timing2plot];
mtx_to_tbl_pgfplots('results/timingVsArOrder',comment,legend_str,data);
