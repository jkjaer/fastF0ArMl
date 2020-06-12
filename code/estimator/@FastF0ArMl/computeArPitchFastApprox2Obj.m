function [objValNoPitch, objValPitch] = ...
        computeArPitchFastApprox2Obj(Obj, dataVector)
    % no periodic signal is present
    [objValNoPitch, ~, nu] = ldr(dataVector, Obj.maxArOrder, Obj.nData);    
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
            neededDftBins = Obj.dftRange(q,1):max(...
                (1:Obj.maxPitchOrder)'.*Obj.dftRange(:,2));
            neededDftData = ...
                exp(-1i*Obj.startIndex*2*pi*neededDftBins/Obj.nDft).*...
                dftData(neededDftBins+1).'/Obj.nData;
            shiftVct = exp(-1i*2*pi*neededDftBins/Obj.nDft);
            neededShiftedDftData = shiftVct.*neededDftData;
            % the +1 is to compensate for MATLAB's indexing
            g = nan(Obj.maxPitchOrder,nPitches,Obj.maxArOrder);
            h = nan(Obj.maxPitchOrder,nPitches,Obj.maxArOrder);
            g(1,:,1) = neededDftData(1,qShiftedDftIdx);
            h(1,:,1) = neededShiftedDftData(1,qShiftedDftIdx);
            W = shiftVct(qShiftedDftIdx);
            objValPitch(1,:,1) = objValNoPitch(1)-...
                2*abs(g(1,:,1)).^2*Obj.nData/(Obj.nData+Obj.maxArOrder);
        else
            g(q,f0Idx,1) = neededDftData(1,qShiftedDftIdx);
            h(q,f0Idx,1) = neededShiftedDftData(1,qShiftedDftIdx);
            W = [W(:,f0Idx);shiftVct(qShiftedDftIdx);];
            % the +1 is to compensate for MATLAB's indexing
            objValPitch(1,f0Idx,q) = objValPitch(1,f0Idx,q-1)-...
                2*abs(g(q,f0Idx,1)).^2*Obj.nData/(Obj.nData+Obj.maxArOrder);
        end
        d = 2*g(1:q,f0Idx,1)*Obj.nData/(Obj.nData+Obj.maxArOrder);
        mu = 2*h(1:q,f0Idx,1)*Obj.nData/(Obj.nData+Obj.maxArOrder);
        for p = 1:Obj.maxArOrder
            h_iOm_g = sum(real(h(1:q,f0Idx,p)).*real(d)+...
                imag(h(1:q,f0Idx,p)).*imag(d),1);
            kappa = ...
                (nu(p)*objValNoPitch(p)-h_iOm_g)./objValPitch(p,f0Idx,q);
            objValPitch(p+1,f0Idx,q) = objValPitch(p,f0Idx,q).*...
                (1-kappa.^2);
            if p < Obj.maxArOrder
                g(q,f0Idx,p+1) = g(q,f0Idx,p) - nu(p)*h(q,f0Idx,p);
                h(q,f0Idx,p+1) = ...
                    W(q,:).*(h(q,f0Idx,p)-nu(p)*g(q,f0Idx,p));
                dOld = d;
                d = dOld-(ones(q,1)*kappa).*mu;
                if p < Obj.maxArOrder-1
                    mu = W.*(mu-(ones(q,1)*kappa).*dOld);
                end
            end
        end
    end
end