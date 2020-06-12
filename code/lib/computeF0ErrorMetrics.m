function [fpeMtx, gpeMtx] = computeF0ErrorMetrics(f0True, f0Est, threshold)
    [nCases, nMc, nMethods] = size(f0Est);
    fpeMtx = nan(nCases, nMc, nMethods);
    gpeMtx = zeros(nCases, nMc, nMethods);
    for ii = 1:nCases
        for jj = 1:nMc
            jjF0True = f0True(ii,jj);
            for kk = 1:nMethods
                f0Error = jjF0True-f0Est(ii,jj,kk);
                if abs(f0Error) > threshold*jjF0True
                    % GPE
                    gpeMtx(ii, jj, kk) = 1;
                else
                    % FPE
                    fpeMtx(ii, jj, kk) = f0Error;
                end
            end
        end
    end
    gpeMtx = 100*squeeze(mean(gpeMtx,2));
    fpeMtx = sqrt(squeeze(nanmean(fpeMtx.^2,2)));
end
