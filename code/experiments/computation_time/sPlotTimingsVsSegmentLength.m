clear
clc
close all

load('results/timingVsSegmentLength.mat');

timing2plot = squeeze(min(timingMtx,[],2));

semilogy(nDataList', timing2plot, 'lineWidth',2)
xlabel('Segment length [samples]');
ylabel('Computation time [s]');
grid on;
legend(methods)

% store data in text file so that it is ready for plotting using tikz
comment = '';
legend_str = ['order ', strjoin(methods)];
data = [nDataList', timing2plot];
mtx_to_tbl_pgfplots('results/timingVsSegmentLength',comment,legend_str,data);
