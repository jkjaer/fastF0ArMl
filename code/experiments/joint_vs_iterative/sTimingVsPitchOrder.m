clear
close all
clc

addpath('../../estimator/');
addpath('../../reference_methods/iterative/');

% setup
maxArOrder = 3;
nData = 512;
nMc = 50;
nIter = 10;
pitchBounds = [1.5/nData, 0.4];
samplingFreq = 1;
maxPitchOrderList = 1:15;
nPitchOrders = length(maxPitchOrderList);
methods = {'F0-AR-ML-E', 'Iterative'};
nMethods = length(methods);
timingMtx = nan(nPitchOrders, nMc, nMethods);

for ii = 1:nPitchOrders
    maxPitchOrder = maxPitchOrderList(ii);
	for jj = 1:nMethods
        switch methods{jj}
            case 'F0-AR-ML-E'
                EstObj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'E', true);
            case 'Iterative'
                EstObj = IterArF0(nData, maxPitchOrder, maxArOrder, ...
                    pitchBounds, 1, 'E', true);
        end
        for kk = 1:nMc
            [ii,jj,kk]
            dataVector = randn(nData,1);
            switch methods{jj}
                case 'F0-AR-ML-E'
                    EstObj.estimate(dataVector);
                case 'Iterative'
                    EstObj.estimate(dataVector, nIter);
            end
            timingMtx(ii,kk,jj) = EstObj.compTime;
        end
    end
end

save('results/timingVsPitchOrderIterative.mat', 'timingMtx', ...
    'maxPitchOrderList', 'nData', 'maxArOrder', 'methods', 'nMc');
