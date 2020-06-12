clear
close all
clc

addpath('../../estimator/');
addpath('../../lib/');
addpath('../../reference_methods/iterative/');

% setup
rng(10)
maxArOrder = 3;
maxPitchOrder = 6;
nData = 512;
nMc = 10000;
nIter = 10;
refinementTol = 1e-5;
pitchBounds = [2/nData, 0.4];
samplingFreq = 1;
snrDbList = -10:10; % dB
nSnrs = length(snrDbList);
methods = {'F0-AR-ML-E', 'Iterative'};
nMethods = length(methods);
f0Mtx = nan(nSnrs, nMc, nMethods+1);
crlbMtx = nan(nSnrs, nMc);

for ii = 1:nSnrs
    ii
    snrDb = snrDbList(ii);
    for jj = 1:nMc
        [dataVector, truePitch, sinusAmps, arParams, exVar] = ...
            generateData(nData, snrDb, maxPitchOrder, maxArOrder, ...
            pitchBounds);
        freqs = truePitch*(1:length(sinusAmps))';
        crlbMtx(ii,jj) = computeAsympCrlb(nData, freqs, sinusAmps, ...
            arParams, exVar, samplingFreq);
        f0Mtx(ii,jj,1) = truePitch;
        for kk = 1:nMethods
            switch methods{kk}
                case 'F0-AR-ML-E'
                    EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                        pitchBounds, 1, 'E');
                case 'Iterative'
                    EstObj = IterArF0(nData, maxPitchOrder, maxArOrder, ...
                        pitchBounds, 1, 'E');
            end
            % make sure that a pitch is estimated
            EstObj.setValidPitchOrders(1:maxPitchOrder);
            % estimate the pitch
            switch methods{kk}
                case 'F0-AR-ML-E'
                    estPitch = EstObj.estimate(dataVector, refinementTol);
                case 'Iterative'
                    estPitch = EstObj.estimate(dataVector, nIter, ...
                        refinementTol);
            end
            f0Mtx(ii,jj,kk+1) = estPitch;
        end
    end
end

save('results/f0VsSnrIterative.mat', 'f0Mtx', 'crlbMtx', 'snrDbList', ...
    'nData', 'maxArOrder', 'maxPitchOrder', 'methods', 'nMc');
