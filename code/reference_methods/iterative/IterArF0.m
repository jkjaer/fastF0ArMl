classdef IterArF0 < handle
    properties (SetAccess=immutable)
        nData;
        maxPitchOrder;
        maxArOrder;
        samplingFreq = 1; % Hz
    end
    
    properties (SetAccess=private)
        FastF0Ml;
        estArOrder;
        estPitchOrder;
        estPitch; % Hz
        estExVar; 
        estSinusAmps;
        estSinusPhases;
        estArParams;
        benchmarkAlg = false;
        compTime;
    end
    
    properties (Hidden = true)
        timeIndices;
    end
    
    methods (Access = public)
        % constructor
        function Obj = IterArF0(nData, maxPitchOrder, maxArOrder, ...
                pitchBounds, samplingFreq, f0Alg, benchmarkAlg)
            Obj.FastF0Ml = FastF0ArMl(nData, maxPitchOrder, 0, ...
                pitchBounds, samplingFreq, f0Alg);
            % we are using the autocorrelation method so we pad the data
            % with maxArOrder zeros
            Obj.nData = nData;
            Obj.maxPitchOrder = maxPitchOrder;
            Obj.maxArOrder = maxArOrder;
            Obj.samplingFreq = samplingFreq;
            Obj.timeIndices = (0:nData-1)'-(nData-1)/2;
            if nargin > 6
                Obj.benchmarkAlg = benchmarkAlg;
            end
        end
        
        function setValidPitchOrders(Obj, validPitchOrders)
            Obj.FastF0Ml.setValidPitchOrders(validPitchOrders);
        end

        function estPitch = estimate(Obj, dataVector, nIter, refinementTol)
            if length(dataVector) ~= Obj.nData
                error('Data vector has the wrong length.');
            end
            if nargin < 4
                refinementTol = 1/Obj.FastF0Ml.nDft;
            end
            estPitch = Obj.estArFirst(dataVector, nIter, refinementTol);
            if Obj.benchmarkAlg
                f = @() Obj.estArFirst(dataVector, nIter, refinementTol);
                Obj.compTime = timeit(f);
            end
        end
            
    end
    methods (Access = private)
        function estPitch = estArFirst(Obj, dataVector, nIter, refinementTol)
            arSignal = dataVector;
            for ii = 1:nIter
                % compute AR parameters
                [predErrVar, arParamMtx] = myLpc(arSignal,Obj.maxArOrder);
                % estimate model order using BIC
                [~, estIdx] = max(-Obj.nData*log(predErrVar)/2-...
                    ((0:Obj.maxArOrder)'/2)*log(Obj.nData));
                Obj.estArOrder = estIdx-1;
                Obj.estExVar = predErrVar(Obj.estArOrder+1);
                % remove predicted AR signal from the observation
                if Obj.estArOrder > 0
                    Obj.estArParams = ...
                        arParamMtx(1:Obj.estArOrder,Obj.estArOrder);
                    modelledArSignal = ...
                        filter([0; Obj.estArParams],1,arSignal);
                    perSignal = dataVector-modelledArSignal;
                else
                    Obj.estArParams = [];
                    perSignal = dataVector;
                end
                % compute pitch parameters
                Obj.estPitch = ...
                    Obj.FastF0Ml.estimate(perSignal, refinementTol);
                Obj.estPitchOrder = Obj.FastF0Ml.estPitchOrder;
                if Obj.estPitchOrder > 0
                    Obj.estSinusAmps = Obj.FastF0Ml.estSinusAmps;
                    Obj.estSinusPhases = Obj.FastF0Ml.estSinusPhases;
                    cisMtx = exp(1i*2*pi*Obj.estPitch*Obj.timeIndices*...
                        (1:Obj.estPitchOrder)/Obj.samplingFreq);
                    modelledPerSignal = real(cisMtx*...
                        (Obj.estSinusAmps.*exp(1i*Obj.estSinusPhases)));
                    arSignal = dataVector-modelledPerSignal;
                else
                    Obj.estSinusAmps = [];
                    Obj.estSinusPhases = [];
                    arSignal = dataVector;
                end
            end
            estPitch = Obj.estPitch;
        end
    end
end
