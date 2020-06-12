function arDataMatrix = computeArDataMatrix(dataVector, arOrder, ...
        method)
    if arOrder > 0 && (nargin < 3 || strcmp(method, 'autocorrelation'))
        dataVector = [dataVector; zeros(arOrder,1)];
        method = 'prewindow';
    end
    nData = length(dataVector);
    if arOrder > 0
        arDataMatrix = nan(nData, arOrder);
        if strcmp(method, 'circulant')
            for ii = 1:arOrder
                shiftIdx = mod((0:nData-1)'-ii,nData)+1;
                arDataMatrix(:, ii) = dataVector(shiftIdx);
            end
        elseif strcmp(method, 'prewindow')
            for ii = 1:arOrder
                arDataMatrix(:, ii) = [...
                    zeros(ii,1); ...
                    dataVector(1:end-ii)];
            end
        else
            error('Method is not implemented yet!');
        end
    else
        arDataMatrix = zeros(nData, 0);
    end
end
