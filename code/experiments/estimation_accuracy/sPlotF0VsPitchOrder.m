clear
clc
close all

addpath('../../lib/');

load('results/f0VsPitchOrder.mat');

[fpeMtx, gpeMtx] = ...
    computeF0ErrorMetrics(f0Mtx(:,:,1), f0Mtx(:,:,2:end), 0.2);
fpeMtx = [sqrt(mean(crlbMtx,2)), fpeMtx];

figure(1)
semilogy(maxPitchOrderList', fpeMtx, '-o', 'lineWidth',2)
xlabel('Pitch order');
ylabel('FPE [cycles/sample]');
grid on;
legend('CRLB', methods{:})

figure(2)
plot(nDataList', gpeMtx, '-o', 'lineWidth',2)
xlabel('Segment length [samples]');
ylabel('GPE [%]');
grid on;
legend(methods{:})

% store data in text file so that it is ready for plotting using tikz
comment = '';
legend_str = ['pitchOrder crlb ', strjoin(methods)];
data = [maxPitchOrderList', rmseMtx];
mtx_to_tbl_pgfplots('results/fpeVsPitchOrder',comment,legend_str,data, 8);

comment = '';
legend_str = ['pitchOrder ', strjoin(methods)];
data = [maxPitchOrderList', gpeMtx];
mtx_to_tbl_pgfplots('results/gpeVsPitchOrder',comment,legend_str,data);
