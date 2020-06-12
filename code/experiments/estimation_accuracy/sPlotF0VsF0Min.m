clear
clc
close all

addpath('../../lib/');

load('results/f0VsF0Min.mat');

nSnrs = length(snrDbList);
for rr = 1:nSnrs
    [fpeMtx, gpeMtx] = ...
        computeF0ErrorMetrics(f0Mtx(:,:,1,rr), f0Mtx(:,:,2:end,rr), 0.2);
    fpeMtx = [fpeMtx, sqrt(mean(crlbMtx(:,:,rr),2))];
    figure(rr)
    subplot(2,1,1)
    semilogy(minPitchList'*nData, fpeMtx, '-o', 'lineWidth',2)
    title(['SNR = ', num2str(snrDbList(rr)), ' dB']);
    xlabel('F0 min. [cycles/segment]');
    ylabel('FPE [cycles/sample]');
    grid on;
    legend(methods{:}, 'CRLB')

    subplot(2,1,2)
    plot(minPitchList'*nData, gpeMtx, '-o', 'lineWidth',2)
    xlabel('F0 min. [cycles/segment]');
    ylabel('GPE [%]');
    grid on;
    legend(methods{:})
    
    % store data in text files so that it is ready for plotting using tikz
    comment = '';
    legend_str = ['f0Min ', strjoin(methods), ' crlb'];
    data = [minPitchList'*nData, fpeMtx];
    mtx_to_tbl_pgfplots(...
        ['results/fpeVsF0Min_',num2str(snrDbList(rr)),'dB'],...
        comment,legend_str,data, 8);

    comment = '';
    legend_str = ['f0Min ', strjoin(methods)];
    data = [minPitchList'*nData, gpeMtx];
    mtx_to_tbl_pgfplots(...
        ['results/gpeVsF0Min_',num2str(snrDbList(rr)),'dB'],...
        comment,legend_str,data);
end
