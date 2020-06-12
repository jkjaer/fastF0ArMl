clear
clc
close all

addpath('../../lib/');

load('results/f0VsSegmentLength.mat');

[fpeMtx, gpeMtx] = ...
    computeF0ErrorMetrics(f0Mtx(:,:,1), f0Mtx(:,:,2:end), 0.2);
fpeMtx = [sqrt(mean(crlbMtx,2)), fpeMtx];

figure(1)
semilogy(nDataList', fpeMtx, '-o', 'lineWidth',2)
xlabel('Segment length [samples]');
ylabel('FPE [cycles/sample]');
grid on;
legend('CRLB', methods{:})

figure(2)
plot(nDataList', gpeMtx, '-o', 'lineWidth',2)
xlabel('Segment length [samples]');
ylabel('GPE [%]');
grid on;
legend(methods{:})

% store data in text files so that it is ready for plotting using tikz
comment = '';
legend_str = ['nData crlb ', strjoin(methods)];
data = [nDataList', fpeMtx];
mtx_to_tbl_pgfplots('results/fpeVsSegmentLength',comment,legend_str,data, 8);

comment = '';
legend_str = ['nData ', strjoin(methods)];
data = [nDataList', gpeMtx];
mtx_to_tbl_pgfplots('results/gpeVsSegmentLength',comment,legend_str,data);