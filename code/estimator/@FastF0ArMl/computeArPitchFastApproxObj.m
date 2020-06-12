function [objValNoPitch, objValPitch] = ...
        computeArPitchFastApproxObj(Obj, dataVector)
    % no periodic signal is present
    [objValNoPitch, rho] = ldr(dataVector, Obj.maxArOrder, Obj.nData);  
    % periodic signal is present
    dftData = fft(dataVector, Obj.nDft);
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
            qRho = rho*ones(1,nPitches);
        else
            % remove unneeded entries
            qRho = qRho(:,f0Idx);
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
            else % p >= 1
                % we run the recursion in the reversed version of beta to 
                % avoid one flipud
                qRho(p,:) = qRho(p,:)-...
                    lambdaTpH(1,:).*lambdaTpH(p+1,:)-...
                    lambdaTmH(1,:).*lambdaTmH(p+1,:);
                if p == 1
                    betaRev = qRho(1,:)./objValPitch(1,f0Idx,q);
                    objValPitch(2,f0Idx,q) = objValPitch(1,f0Idx,q).*...
                        (1-betaRev.^2);
                else % p >= 2
                    nu = (qRho(p,:)-sum(betaRev.*qRho(1:p-1,:),1))./...
                    	objValPitch(p,f0Idx,q);
                    objValPitch(p+1,f0Idx,q) = ...
                        objValPitch(p,f0Idx,q).*(1-nu.^2);
                    if p < Obj.maxArOrder
                        betaRev = ...
                            [nu; betaRev-flipud(betaRev).*(ones(p-1,1)*nu)];
                    end
                end
            end
        end
    end
end