clear
clc
close all

addpath('../../lib/');

load('results/f0VsSnr.mat');

[fpeMtx, gpeMtx] = ...
    computeF0ErrorMetrics(f0Mtx(:,:,1), f0Mtx(:,:,2:end), 0.2);
fpeMtx = [fpeMtx, sqrt(mean(crlbMtx,2))];

figure(1)
semilogy(snrDbList', fpeMtx, '-o', 'lineWidth',2)
xlabel('SNR [dB]');
ylabel('FPE [cycles/sample]');
grid on;
legend(methods{:}, 'CRLB')

figure(2)
plot(snrDbList', gpeMtx, '-o', 'lineWidth',2)
xlabel('SNR [dB]');
ylabel('GPE [%]');
grid on;
legend(methods{:})

% store data in text files so that it is ready for plotting using tikz
comment = '';
legend_str = ['snr ', strjoin(methods), ' crlb'];
data = [snrDbList', fpeMtx];
mtx_to_tbl_pgfplots('results/fpeVsSnr',comment,legend_str,data, 8);

comment = '';
legend_str = ['snr ', strjoin(methods)];
data = [snrDbList', gpeMtx];
mtx_to_tbl_pgfplots('results/gpeVsSnr',comment,legend_str,data);
