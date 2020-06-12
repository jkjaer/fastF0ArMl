function [objValNoPitch, objValPitch] = ...
        computeArPitchFastExactObj(Obj, dataVector)
    % no periodic signal is present
    [objValNoPitch, rho] = ldr(dataVector, Obj.maxArOrder);  
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
            y = objValNoPitch(1)*ones(Obj.maxArOrder,nPitches);
            if Obj.maxArOrder > 1
                eta12 = rho(1)*ones(1,nPitches);
                if Obj.maxArOrder > 2
                    eta13 = rho(2)*ones(1,nPitches);
                    eta23 = rho(1)*ones(1,nPitches);
                end
            end
        else
            % remove unneeded entries
            qRho = qRho(:,f0Idx);
            y = y(:,f0Idx);
            if Obj.maxArOrder > 1
                eta12 = eta12(f0Idx);
                if Obj.maxArOrder > 2
                    eta13 = eta13(f0Idx);
                    eta23 = eta23(f0Idx);
                end
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
            else % p >= 1
                if p == 2
                    eta12 = eta12-lambdaTpH(p,:).*lambdaTpH(p+1,:)-...
                        lambdaTmH(p,:).*lambdaTmH(p+1,:);
                    chi = eta12./y(1,:);
                    eta = eta12;
                elseif p == 3
                    eta13 = eta13-lambdaTpH(p-1,:).*lambdaTpH(p+1,:)-...
                        lambdaTmH(p-1,:).*lambdaTmH(p+1,:);
                    eta23 = eta23-lambdaTpH(p,:).*lambdaTpH(p+1,:)-...
                        lambdaTmH(p,:).*lambdaTmH(p+1,:);
                    chi2 = (eta23-chi.*eta13)./c;
                    chi = [eta13./y(1,:)-chi2.*chi;chi2];
                    eta = [eta13;eta23];
                end
                qRho(p,:) = qRho(p,:)-...
                    lambdaTpH(1,:).*lambdaTpH(p+1,:)-...
                    lambdaTmH(1,:).*lambdaTmH(p+1,:);
                y(p,:) = y(p,:)-lambdaTpH(p+1,:).^2-...
                    lambdaTmH(p+1,:).^2;
                if p == 1
                    objValPitch(p+1,f0Idx,q) = objValPitch(p,f0Idx,q)-...
                        qRho(p,:).^2./y(p,:);
                else % p >= 2
                    Lambda = qRho(p,:)-sum(chi.*qRho(1:p-1,:),1);
                    c = y(p,:)-sum(eta.*chi,1);
                    objValPitch(p+1,f0Idx,q) = objValPitch(p,f0Idx,q)-...
                        Lambda.^2./c;
                end
            end
        end
    end
end