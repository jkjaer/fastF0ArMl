function bayesFactor = pitchArBicModelComparison(nData, objValNoPitch, ...
        objValPitch, validPitchOrders)
    maxArOrder = length(objValNoPitch)-1;
    maxPitchOrder = size(objValPitch,3);
    logBayesFactor = -inf(maxPitchOrder+1, maxArOrder+1);
    for ii = validPitchOrders
        for jj = 0:maxArOrder
            % estimate the excitation variance
            if ii == 0
                estExVar = min(objValNoPitch(jj+1));
            else
                estExVar = min(objValPitch(jj+1,:,ii));
            end
            if estExVar < 0
                warning(['The estimated excitation variance is ' ...
                    'negative. Expect unreliable results.']);
                estExVar = inf;
            end
            logBayesFactor(ii+1,jj+1) = ...
                -nData*log(estExVar)/2-(ii+jj/2)*log(nData);
        end
    end
    % normalise the Bayes' factors
    bayesFactor = exp(logBayesFactor-max(max(logBayesFactor)));
    bayesFactor = bayesFactor/sum(sum(bayesFactor));
end