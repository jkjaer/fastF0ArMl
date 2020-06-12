function [objValNoPitch, objValPitch] = ...
        computeArPitchFastExactObj(Obj, dataVector)
    % no periodic signal is present
    covVct = computeCovVct(dataVector, Obj.maxArOrder, Obj.nData);
    [objValNoPitch, cholVct] = genSchurToeplitzAlg(covVct);
    % periodic signal is present
    dftData = fft(dataVector, Obj.nDft);
    % lower triangular matrix of cholesky indices
    cholIdx = tril((-(Obj.maxArOrder-1):0)'*ones(1,Obj.maxArOrder)+...
        ones(Obj.maxArOrder,1)*(cumsum(Obj.maxArOrder:-1:1)));
    for q = 1:Obj.maxPitchOrder
        qShiftedDftIdx = q*(Obj.dftRange(q,1):Obj.dftRange(q,2))-...
            Obj.dftRange(q,1)+1;
        nPitches = length(qShiftedDftIdx);
        f0Idx = 1:nPitches;
        if q == 1
            % initialise array
            objValPitch = ...
                nan(Obj.maxArOrder+1, nPitches, Obj.maxPitchOrder);
            % the +1 is to compensate for MATLAB's indexing
            neededDftBins = Obj.dftRange(q,1):max(...
                (1:Obj.maxPitchOrder)'.*Obj.dftRange(:,2));
            Yhat = exp(-1i*((0:Obj.maxArOrder)'+Obj.startIndex)*...
                2*pi*neededDftBins/Obj.nDft).*...
                ((ones(Obj.maxArOrder+1,1)*dftData(neededDftBins+1).'));
            dftDataTpH = nan(Obj.maxPitchOrder, nPitches, Obj.maxArOrder+1);
            dftDataTmH = nan(Obj.maxPitchOrder, nPitches, Obj.maxArOrder+1);
            if Obj.maxArOrder > 0
                qRho = covVct(2:end)*ones(1,nPitches);
                % Cholesky factor with stacked columns
                cholVct = cholVct*ones(1,nPitches);
            end
        else
            % remove unneeded entries
            if Obj.maxArOrder > 0
                qRho = qRho(:,f0Idx);
                cholVct = cholVct(:,f0Idx);
            end
        end
        lambdaTpH = nan(Obj.maxArOrder+1, nPitches);
        lambdaTmH = nan(Obj.maxArOrder+1, nPitches);
        for p = 0:Obj.maxArOrder
            dftDataTpH(q,f0Idx,p+1) = real(Yhat(p+1,qShiftedDftIdx));
            dftDataTmH(q,f0Idx,p+1) = -imag(Yhat(p+1,qShiftedDftIdx));
            lambdaTpH(p+1,:) = sum(Obj.gammaTpH{q}(1:q, f0Idx).*...
                    dftDataTpH(1:q,f0Idx,p+1),1)/sqrt(Obj.nData);
            lambdaTmH(p+1,:) = sum(Obj.gammaTmH{q}(1:q, f0Idx).*...
                    dftDataTmH(1:q,f0Idx,p+1),1)/sqrt(Obj.nData);
            if p == 0
                if q == 1
                    objValPitch(1,:,1) = objValNoPitch(1)-...
                        lambdaTpH(1,:).^2-lambdaTmH(1,:).^2;
                else
                    objValPitch(1,f0Idx,q) = objValPitch(1,f0Idx,q-1)-...
                        lambdaTpH(1,:).^2-lambdaTmH(1,:).^2;
                end
            else % p > 0
                qRho(p,:) = qRho(p,:)-...
                    lambdaTpH(1,:).*lambdaTpH(p+1,:)-...
                    lambdaTmH(1,:).*lambdaTmH(p+1,:);
            end
        end
        if Obj.maxArOrder > 0
            % update Cholesky decomposition using two rank one downdating
            % steps
            cholVct = cholDowndating(cholVct, lambdaTpH(2:end,:));
            cholVct = cholDowndating(cholVct, lambdaTmH(2:end,:));
            % recursively compute the objective using forward substitution
            objValPitch(2:p+1,f0Idx,q) = ...
                ones(Obj.maxArOrder,1)*objValPitch(1,f0Idx,q)-...
                orderRecursiveForwardSubstitution(cholVct, qRho);
        end
    end
end