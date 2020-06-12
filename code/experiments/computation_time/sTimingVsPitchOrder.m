clear
close all
clc

addpath('../../estimator/');

% setup
maxArOrder = 3;
nData = 512;
nMc = 50;
pitchBounds = [1.5/nData, 0.4];
samplingFreq = 1;
maxPitchOrderList = 1:15;
nPitchOrders = length(maxPitchOrderList);
methods = {'F0-AR-ML-E', 'F0-AR-ML-A', 'F0-AR-ML-A1', 'F0-AR-ML-A2', ...
    'Naive'};
nMethods = length(methods);
timingMtx = nan(nPitchOrders, nMc, nMethods);

for ii = 1:nPitchOrders
    maxPitchOrder = maxPitchOrderList(ii);
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

save('results/timingVsPitchOrder.mat', 'timingMtx', 'maxPitchOrderList', ...
    'nData', 'maxArOrder', 'methods', 'nMc');
