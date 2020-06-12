clear
close all
clc

addpath('../../estimator/');
addpath('../../lib/');

% setup
rng(10)
maxArOrder = 3;
maxPitchOrder = 6;
nData = 512;
nMc = 10000;
samplingFreq = 1;
snrDbList = [0,5,10]; % dB
nSnrs = length(snrDbList);
maxPitch = 0.4;
minPitchList = (1:0.25:8)/nData;
nPitchBounds = length(minPitchList);
methods = {'F0-AR-ML-E', 'F0-AR-ML-A', 'F0-AR-ML-A1'};
nMethods = length(methods);
f0Mtx = nan(nPitchBounds, nMc, nMethods+1, nSnrs);
crlbMtx = nan(nPitchBounds, nMc, nSnrs);

for rr = 1:nSnrs
    snrDb = snrDbList(rr);
    for ii = 1:nPitchBounds
        [rr,ii]
        pitchBounds = [minPitchList(ii), maxPitch];
        for jj = 1:nMc
            [dataVector, truePitch, sinusAmps, arParams, exVar] = ...
                generateData(nData, snrDb, maxPitchOrder, maxArOrder, ...
                pitchBounds);
            freqs = truePitch*(1:length(sinusAmps))';
            crlbMtx(ii,jj,rr) = computeAsympCrlb(nData, freqs, sinusAmps, ...
                arParams, exVar, samplingFreq);
            f0Mtx(ii,jj,1,rr) = truePitch;
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
                f0Mtx(ii,jj,kk+1,rr) = estPitch;
            end
        end
    end
end

save('results/f0VsF0Min.mat', 'f0Mtx', 'crlbMtx', 'minPitchList', ...
    'snrDbList', 'nData', 'maxArOrder', 'maxPitchOrder', 'methods', 'nMc');
