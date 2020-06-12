classdef FastF0ArMl < handle
%
% Estimates the pitch and AR parameters and the model orders of a 
% periodic signal in AR noise.
% The signal model is
%     y(n) = s(n) + e(n)
% where y(n) is the observed data, s(n) is the periodic signal, and e(n) is
% the autoregressive noise.
%
% The estimation is performed using one of the following methods
% E) F0-AR-ML-E (fast, exact implementation)
% A) F0-AR-ML-A (fast, approximate implementation)
% A1) F0-AR-ML-A1 (fast, approximate implementation)
% A2) F0-AR-ML-A2 (fast, approximate implementation)
% N) F0-AR-ML-N (naïve, exact implementation)
% Method E) and N) give the exact same results, but E) is much faster than
% N). Method A1) and A2) also give the exact same results, but are only 
% approximate implementations of method E) and N). Finally, method A) is
% also an approximate implementation, but is based on a different
% approximation than method A1) and A2).
% Note that exact here is when the autocorrelation window is used for a 
% data segment.
%
% To use this object, you first have to contruct it. This is done by
% running, e.g., 
% <name of object> = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
%     pitchBounds, samplingFreq, alg, benchmarkAlg)
% where:
%
%    nData: Signal length (integer, scalar)    
%    
%    maxPitchOrder: maximum number of harmonic components (integer, scalar).
%
%    maxArOrder: maximum number of AR parameters (integer, scalar).
%
%    pitchBounds: Lower and upper bound on the pitch/fundamental frequency
%        in Hz. The lower bound should not be lower than the frequency
%        correponding to 1 period/segment and cannot be higher than half
%        the sampling frequency. (2D vector)
%
%    samplingFreq: Sampling frequency in Hz. Set to 1 if you want the unit
%        to be cycles/sample. (scaler)
%
%    alg: Set to either 'E', 'A', 'A1', 'A2', and 'N'. 'E' is default.
%
%    benchmarkAlg: Boolean variable indicating if the algorithm should
%    be timed. Default is false.
%    
% The contructor returns a FastF0ArMl handle.
%
% Once the object has been contructed, you can invoke a number of methods.
% These are invoked by running
% [output variables] = <name of object>.<name of method>(input variables)
% listed below.
%
% Methods:
%
%   setValidPitchOrders(validPitchOrders): set the valid sinusoidal orders
%
%     Input: 
%    
%         validPitchOrders: a vector containing the valid sinusoidal
%         orders. The list can only include integers from the range 
%         [0,maxPitchOrder]. The default value is all orders, i.e., 
%         validPitchOrders = 0:maxPitchOrder.
% 
%   estimate(dataVector): Estimates all model parameters and orders
%
%     Input: 
%    
%         dataVector: the data vector (vector, length nData)    
%
%     Output:
%
%         estPitch: the estimated pitch in Hz. (scalar)
%
%   If values of, e.g., the estimated harmonic amplitudes, the
%   AR-parameteres, or the orders can be retrieved by looking at the
%   properties of the object. For example
%   estPitchOrder = <name of object>.estPitchOrder;
%   gives the estimated number of harmonic components. To see a list of all
%   properties, write
%   <name of object>
%   in the command prompt.
%
%   plotObjective(pitchOrder, arOrder): Plot the objective function for a
%   specific pitch and AR order.
%
%   Input:
%
%       pitchOrder: (optional value is the estimated order) The number of 
%           harmonics.
%
%       arOrder: (optional value is the estimated order) The number of 
%           AR parameters.
%
%   computeModelledSpectrum: Compute the modelled power spectral density.
%
%   Input: 
%
%       nDft: Number of frequency points. The default is simply the class
%       property of the same name.
%
%   Output:
%
%       modelledSpectrum: The modelled PSD in dB/Hz.
%
%       arPsd: The modelled AR PSD in dB/Hz.
%
%       pitchPsd: The modelled harmonic PSD in dB/Hz.
%
%       freqVector: A vector of frequencies corresponding to the point in
%           the computed PSDs
%       
%   plotModelledSpectrum: Plots the modelled PSDs on top of the periodogram
%       estimate.
%
%   Input:
%   
%       plotRangeDb: (optional with 100 dB as default) limit the plotting 
%           range.
%
%   computePrewhitenedSpectrum: Compute the pre-whitened spectral data
%       on top of the periodogram of the original data
%
%   Input: 
%
%       nDft: Number of frequency points. The default is simply the class
%       property of the same name.
%
%   Output:
%
%       preWhitenedSpectrum: The pre-whitened spectrum in dB/Hz.
%
%
%       freqVector: A vector of frequencies corresponding to the point in
%           the computed PSDs
%
%   plotPrewhitenedSpectrum: Plots the periodogram of the pre-whitened data
%       on top of the periodogram of the original data
%
%   Input:
%   
%       plotRangeDb: (optional with 100 dB as default) limit the plotting 
%           range.
%
%   reportProblem: A report a discovered problem.
%   
%
% Current version: 2.0.0 (2020-06-04)
%
%
% Authors: J. K. Nielsen, jkn@create.aau.dk
% 

    % Can only be set in the constructor
    properties (SetAccess=immutable)
        nData;
        maxPitchOrder;
        maxArOrder;
        pitchBounds; % Hz
        samplingFreq = 1; % Hz
        nDft;
        pitchResolution; % Hz
        alg = 'E';
        startIndex;
    end
    
    properties (SetAccess=private)
        dataVector;
        objValNoPitch;
        objValPitch;
        fullPitchGrid; % Hz
        dftRange;
        bayesFactor; % probabilities
        estArOrder;
        estPitchOrder;
        validPitchOrders;
        estPitch; % Hz
        estExVar; 
        estSinusAmps;
        estSinusPhases;
        estArParams;
        benchmarkAlg = false;
        compTime;
    end
    
    properties (Hidden=true)
        % only if alg. A is used
        gammaTpH=nan;
        gammaTmH=nan;
    end
    
    methods (Access = public)
        % constructor
        function Obj = FastF0ArMl(nData, maxPitchOrder, maxArOrder, ...
                pitchBounds, samplingFreq, alg, benchmarkAlg)
            if nargin > 4
                Obj.samplingFreq = samplingFreq;
            end
            if nargin > 5
                if strcmp(alg, 'E')
                    Obj.alg = 'E';
                elseif strcmp(alg, 'A')
                    Obj.alg = 'A';
                elseif strcmp(alg, 'A1')
                    Obj.alg = 'A1';
                elseif strcmp(alg, 'A2')
                    Obj.alg = 'A2';
                elseif strcmp(alg, 'N')
                    Obj.alg = 'N';
                else
                    error(['You have not selected a valid '...
                        'algorithm (either E, A, A1, A2, or N).']);
                end
            end
            % we are using the autocorrelation method so we pad the data
            % with maxArOrder zeros
            minNDft = 5*nData*maxPitchOrder;
            Obj.nDft = 2^ceil(log2(minNDft));
            Obj.pitchResolution = samplingFreq/Obj.nDft;
            Obj.pitchBounds = pitchBounds;
            Obj.maxPitchOrder = maxPitchOrder;
            Obj.validPitchOrders = 0:maxPitchOrder;
            Obj.maxArOrder = maxArOrder;
            Obj.nData = nData;
            Obj.startIndex = -(Obj.nData+maxArOrder-1)/2;
            if nargin > 6
                Obj.benchmarkAlg = benchmarkAlg;
            end
            if pitchBounds(1) < samplingFreq/Obj.nData
                warning(['The lower pitch bound is set lower than', ...
                    ' 1 periods/segment. Inaccurate results will ', ...
                    'most likely be produced!']);
            end
            if maxPitchOrder == 0
                error(['The maximum pitch order is set to 0,', ...
                    ' and this is not supported.', ...
                    ' Please use the lpc function instead.']);
            end
            % pre-compute some quantaties
            [Obj.fullPitchGrid, Obj.dftRange] = ...
                computePitchGrid(Obj.nDft, pitchBounds, maxPitchOrder, ...
                samplingFreq);
            % for F0-AR-ML-E (alg. A) only
            if strcmp(Obj.alg, 'E') || strcmp(Obj.alg, 'A')
                [Obj.gammaTpH, Obj.gammaTmH] = ...
                    Obj.dataIndependentStep(Obj.nData+maxArOrder);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % other public methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the valid pitch orders
        setValidPitchOrders(Obj, validPitchOrders);
        % estimate the parameters and orders
        estPitch = estimate(Obj, dataVector, refinementTol);
        % plot the objective function for specified orders
        plotObjective(Obj, pitchOrder, arOrder);
        % compute the modelled spectrum
        [modelledSpectrum, arPsd, pitchPsd, freqVector] = ...
            computeModelledSpectrum(Obj, nDft);
        % plot the computed spectrum
        plotModelledSpectrum(Obj, plotRangeDb);
        % compute the pre-whitened spectrum
        [preWhitenedSpectrum, freqVector] = ...
            computePrewhitenedSpectrum(Obj, nDft);
        % plot the pre-whitened spectrum
        plotPrewhitenedSpectrum(Obj, plotRangeDb);
        % report a problem
        reportProblem(Obj);
    end
    methods (Access = private)
        % compute the objective on a grid of frequencies
        [objValNoPitch, objValPitch, pitchGrids] = ...
        	computeObjOnGrid(Obj, dataVector, samplingFreq, alg);
        % compute the objective function using the F0-AR-ML-E method
        [objValNoPitch, objValPitch] = ...
            computeArPitchFastExactObj(Obj, dataVector);
        % compute the objective function using the F0-AR-ML-A method
        [objValNoPitch, objValPitch] = ...
            computeArPitchFastApproxObj(Obj, dataVector);
        % compute the objective function using the F0-AR-ML-A1 method
        [objValNoPitch, objValPitch] = ...
            computeArPitchFastApprox1Obj(Obj, dataVector);
        % compute the objective function using the F0-AR-ML-A2 method
        [objValNoPitch, objValPitch] = ...
            computeArPitchFastApprox2Obj(Obj, dataVector);
        % compute the objective function using the naive method
        [objValNoPitch, objValPitch] = ...
            computeArPitchExactObj(Obj, dataVector);
        % estimate the model order
        [estArOrder, estPitchOrder] = estimateModelOrders(Obj, ...
                dataVector, samplingFreq);
        % naïve computation of all parameters
        [estPitch, estSinusAmps, estSinusPhases, estArParams, estExVar] ...
            = estimatePitchArParameters(Obj, dataVector, ...
            pitchOrder, arOrder, refinementTol);
        % compute the estimated linear parameters for the sinuoids
        [estSinusAmps, estSinusPhases] = ...
            estimatePitchAmpsAndPhases(Obj, estArParams, ...
            shapedPitchLinearParameters, estPitch);
        % compute the data independent step
        [gammaTpH, gammaTmH] = dataIndependentStep(Obj, nData);
    end
end
