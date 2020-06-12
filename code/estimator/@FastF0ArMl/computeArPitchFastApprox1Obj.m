function [objValNoPitch, objValPitch] = ...
        computeArPitchFastApprox1Obj(Obj, dataVector)
    % no periodic signal is present
    [objValNoPitch, rho] = ldr(dataVector, Obj.maxArOrder, Obj.nData);  
    % periodic signal is present
    phi = 2*abs((fft(dataVector,Obj.nDft))).^2/...
        (Obj.nData*(Obj.nData+Obj.maxArOrder));
    for q = 1:Obj.maxPitchOrder
        qDftIdx = q*(Obj.dftRange(q,1):Obj.dftRange(q,2))+1;
        qShiftedDftIdx = qDftIdx-Obj.dftRange(q,1);
        nPitches = length(qDftIdx);
        f0Idx = 1:nPitches;
        if q == 1
            % initialise array
            objValPitch = ...
                nan(Obj.maxArOrder+1, nPitches, Obj.maxPitchOrder);
            % the +1 is to compensate for MATLAB's indexing
            objValPitch(1,:,1) = objValNoPitch(1)-phi(qDftIdx)';
            neededDftBins = Obj.dftRange(q,1):max(...
                (1:Obj.maxPitchOrder)'.*Obj.dftRange(:,2));
            Phi = cos((1:Obj.maxArOrder)'*2*pi*neededDftBins/Obj.nDft).*...
                (ones(Obj.maxArOrder,1)*phi(neededDftBins+1)');
            rho = rho*ones(1,nPitches)-Phi(:,qShiftedDftIdx);
        else
            % the +1 is to compensate for MATLAB's indexing
            objValPitch(1,f0Idx,q) = ...
                objValPitch(1,f0Idx,q-1)-phi(qDftIdx)';
            rho = rho(:,f0Idx)-Phi(:,qShiftedDftIdx);
        end
        % we run the recursion in the reversed version beta to avoid one flipud
        if Obj.maxArOrder > 0
            betaRev = rho(1,:)./objValPitch(1,f0Idx,q);
            objValPitch(2,f0Idx,q) = ...
                objValPitch(1,f0Idx,q).*(1-betaRev.^2);
            for p = 2:Obj.maxArOrder
                nu = (rho(p,:)-sum(betaRev.*rho(1:p-1,:),1))./...
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