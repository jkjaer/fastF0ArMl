function estPitch = estimate(Obj, dataVector, refinementTol)
    % store the data
    Obj.dataVector = dataVector;
    % zero-pad the data vector (we use the autocorrelation method)
    dataVector = [dataVector; zeros(Obj.maxArOrder, 1)];
%     if length(dataVector) ~= Obj.nData
%         error('Data vector has the wrong length.');
%     end
    if nargin < 3
        refinementTol = Obj.samplingFreq/Obj.nDft;
    end
    % Compute the objective on a grid
    [Obj.objValNoPitch, Obj.objValPitch] = ...
        Obj.computeObjOnGrid(dataVector);

    % Estimate the model order
    [Obj.estArOrder, Obj.estPitchOrder] = Obj.estimateModelOrders;

    % estimate the pitch, harmonic amplitudes and phases, 
    % the AR parameters, and the excitation variance
    [Obj.estPitch, Obj.estSinusAmps, Obj.estSinusPhases, ...
        Obj.estArParams, Obj.estExVar] = ...
        Obj.estimatePitchArParameters(dataVector, ...
        Obj.estPitchOrder, Obj.estArOrder, refinementTol);

    % return the pitch
    estPitch = Obj.estPitch;
end
