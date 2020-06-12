function bayesFactor = pitchArModelComparison(nData, objValNoPitch, ...
        objValPitch, dftRange, nDft, delta, validPitchOrders)
    [arOrderPlus1, ~, pitchOrder] = size(objValPitch);
    arOrder = arOrderPlus1-1;
    logBayesFactor = -inf(pitchOrder+1, arOrder+1);
    deltaFreq = 1/nDft;
    for ii = validPitchOrders
        if ii > 0
            nPitches = dftRange(ii,2)-dftRange(ii,1)+1;
        end
        for jj = 0:arOrder
            if ii == 0
                % compute the coefficient of determination for no pitch
                cod = 1-objValNoPitch(jj+1)/objValNoPitch(1);
            else
                % compute the coefficient of determination for a pitch
                cod = 1-objValPitch(jj+1,1:nPitches,ii)/objValNoPitch(1);
            end
            if any(cod > 1)
                warning(['The estimated excitation variance is ' ...
                    'negative. Expect unreliable results.']);
                cod(cod > 1) = 0;
            end
            [gHat, tauVar] = computeLaplaceParameters(cod, 1, ...
                (nData-(2*ii+jj)-delta)/2, nData/2);
            logMarginalLikelihood = log(gHat*(delta-2)/2)+...
                (nData-(2*ii+jj)-delta)*log(1+gHat)/2-...
                nData*log(1+gHat.*(1-cod))/2+log(2*pi*tauVar)/2;
            if ii == 0
                logBayesFactor(ii+1,jj+1) = logMarginalLikelihood;
            else
                % Avoid numerical problem when doing the integration
                maxLogMarginalLikelihood = max(logMarginalLikelihood);
                logBayesFactor(ii+1,jj+1) = maxLogMarginalLikelihood+...
                    log(deltaFreq*trapz(exp(...
                    logMarginalLikelihood-maxLogMarginalLikelihood)));
            end
        end
    end
    % normalise the Bayes' factors
    bayesFactor = exp(logBayesFactor-max(max(logBayesFactor)));
    bayesFactor = bayesFactor/sum(sum(bayesFactor));
end

function [gHat, tauVar] = computeLaplaceParameters(cod, v, w, u)
    a = (v+w-u)*(1-cod);
    b = (u-v)*cod+2*v+w-u;
    gHat = (b+sqrt(b.^2-4.*a.*v))./(-2*a);
    tauVar = 1./(gHat.*(1-cod)*u./(1+gHat.*(1-cod)).^2-...
        gHat.*w./(1+gHat).^2);
end
