clear
clc
close all

addpath('../../lib/');

load('results/f0VsSnrIterative.mat');

[fpeMtx, gpeMtx] = ...
    computeF0ErrorMetrics(f0Mtx(:,:,1), f0Mtx(:,:,2:end), 0.2);
fpeMtx = [sqrt(mean(crlbMtx,2)), fpeMtx];

figure(1)
semilogy(snrDbList', fpeMtx, '-o', 'lineWidth',2)
xlabel('SNR [dB]');
ylabel('RMSE [cycles/sample]');
grid on;
legend('CRLB', methods{:})

figure(2)
plot(snrDbList', gpeMtx, '-o', 'lineWidth',2)
xlabel('SNR [dB]');
ylabel('GPE [%]');
grid on;
legend(methods{:})

% store data in text file so that it is ready for plotting using tikz
comment = '';
legend_str = ['snr crlb ', strjoin(methods)];
data = [snrDbList', fpeMtx];
mtx_to_tbl_pgfplots('results/fpeVsSnrIterative',comment,legend_str,data,8);

comment = '';
legend_str = ['snr ', strjoin(methods)];
data = [snrDbList', gpeMtx];
mtx_to_tbl_pgfplots('results/gpeVsSnrIterative',comment,legend_str,data);
