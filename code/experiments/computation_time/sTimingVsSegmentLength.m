clear
close all
clc

addpath('../../estimator/');

% setup
maxPitchOrder = 3;
maxArOrder = 3;
nMc = 50;
samplingFreq = 1;
nDataList = 2.^(6:12);
nDataLengths = length(nDataList);
methods = {'F0-AR-ML-E', 'F0-AR-ML-A', 'F0-AR-ML-A1', 'F0-AR-ML-A2', ...
    'Naive'};
nMethods = length(methods);
timingMtx = nan(nDataLengths, nMc, nMethods);

for ii = 1:nDataLengths
    nData = nDataList(ii);
    pitchBounds = [1.5/nData, 0.4];
	for jj = 1:nMethods
        switch methods{jj}
            case 'F0-AR-ML-E'
                EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'E', true);
            case 'F0-AR-ML-A'
                EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'A', true);
            case 'F0-AR-ML-A1'
                EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'A1', true);
            case 'F0-AR-ML-A2'
                EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'A2', true);
            case 'Naive'
                EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'N', true);
        end
        for kk = 1:nMc
            [ii,jj,kk]
            dataVector = randn(nData,1);
            EstObj.estimate(dataVector);
            timingMtx(ii,kk,jj) = EstObj.compTime;
        end
    end
end

save('results/timingVsSegmentLength.mat', 'timingMtx', 'maxArOrder', ...
    'nDataList', 'maxPitchOrder', 'methods', 'nMc');
