clear
close all
clc

addpath('../../estimator/');
addpath('../../lib/');

% setup
rng(10)
maxArOrder = 10;
nData = 512;
pitchBounds = [2/nData, 0.4];
nMc = 500;
samplingFreq = 1;
snrDb = 5; % dB
maxPitchOrderList = 1:10;
nPitchOrders = length(maxPitchOrderList);
methods = {'F0-AR-ML-E', 'F0-AR-ML-A', 'F0-AR-ML-A1'};
nMethods = length(methods);
f0Mtx = nan(nPitchOrders, nMc, nMethods+1);
crlbMtx = nan(nPitchOrders, nMc);

for ii = 1:nPitchOrders
    ii
    maxPitchOrder = maxPitchOrderList(ii);
    for jj = 1:nMc
        [dataVector, truePitch, sinusAmps, arParams, exVar] = ...
            generateData(nData, snrDb, maxPitchOrder);
        freqs = truePitch*(1:length(sinusAmps))';
        crlbMtx(ii,jj) = computeAsympCrlb(nData, freqs, sinusAmps, ...
            arParams, exVar, samplingFreq);
        f0Mtx(ii,jj,1) = truePitch;
        for kk = 1:nMethods
            switch methods{kk}
                case 'F0-AR-ML-E'
                    EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                        pitchBounds, 1, 'E');
                case 'F0-AR-ML-A'
                    EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                        pitchBounds, 1, 'A');
                case 'F0-AR-ML-A1'
                    EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                        pitchBounds, 1, 'A1');
            end
            % make sure that a pitch is estimated
            EstObj.setValidPitchOrders(1:maxPitchOrder);
            % estimate the pitch
            estPitch = EstObj.estimate(dataVector, 1e-6);
            f0Mtx(ii,jj,kk+1) = estPitch;
        end
    end
end

save('results/f0VsPitchOrder.mat', 'f0Mtx', 'crlbMtx', 'snrDb', ...
    'nData', 'maxArOrder', 'maxPitchOrderList', 'methods', 'nMc');
