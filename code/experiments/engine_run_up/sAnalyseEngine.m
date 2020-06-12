clear
clc
close all

% add paths to the data and the estimator
addpath('../../estimator/');
addpath('../../lib/');
addpath data

% load the data
filename = 'carEngineRunUp';
[rawData, samplingFreq] = audioread([filename, '.wav']);
tacho_raw = extractTachoSignalFromWav(filename, 'wav');

% resample the data and extract the left channel
resamplingFreq = samplingFreq/20; % Hz
resampledData = resample(rawData(:,1), resamplingFreq, samplingFreq);
resampledTacho = resample(tacho_raw(:,2), resamplingFreq, samplingFreq);


% set up the segment-by-segment processing (no overlap)
nData = length(resampledData);
segmentTime = 0.25; % seconds
segmentLength = round(segmentTime*resamplingFreq);
nSegments = floor(nData/segmentLength);

% set up the analysis
maxArOrder = 5;
maxPitchOrder = 25;
pitchBounds = [10, 100]; % Hz

% do the analysis
methods = {'Clean', 'F0-AR-ML-E', 'F0-ML-E'};
nMethods = length(methods);
estArOrder = nan(nMethods,nSegments);
estPitchOrder = nan(nMethods,nSegments);
estPitchHz = nan(nMethods,nSegments);
exampleSegmentNo = 40;
exampleNDft = 2048;
for ii = 1:nMethods
    disp(['Processing method ', num2str(ii), ' of ', ...
            num2str(nMethods)]);
    switch methods{ii}
        case 'Clean'
            Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
                maxArOrder, pitchBounds, resamplingFreq, 'E');
            obsData = resampledTacho;
        case 'F0-AR-ML-E'
            Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
                maxArOrder, pitchBounds, resamplingFreq, 'E');
            obsData = resampledData;
        case 'F0-AR-ML-A'
            Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
                maxArOrder, pitchBounds, resamplingFreq, 'A');
            obsData = resampledData;
        case 'F0-AR-ML-A1'
            Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
                maxArOrder, pitchBounds, resamplingFreq, 'A1');
            obsData = resampledData;
        case 'F0-ML-E'
            Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
                0, pitchBounds, resamplingFreq, 'E');
            obsData = resampledData;
    end
    Est.setValidPitchOrders(20:maxPitchOrder);
    idx = 1:segmentLength;
    for jj = 1:nSegments
        disp(['Processing segment ', num2str(jj), ' of ', ...
            num2str(nSegments)]);
        dataSegment = obsData(idx);
        % estimate parameters and orders
        estPitchHz(ii,jj) = Est.estimate(dataSegment);
        estPitchOrder(ii,jj) = Est.estPitchOrder;
        estArOrder(ii,jj) = Est.estArOrder;
        if strcmp(methods{ii},'F0-AR-ML-E')  && jj == exampleSegmentNo
            [modelledSpectrum, ~, ~, ~] = ...
                Est.computeModelledSpectrum(exampleNDft);
            [preWhitenedSpectrum, freqVector] = ...
                Est.computePrewhitenedSpectrum(exampleNDft);
            perEst = abs(fft(dataSegment, exampleNDft)).^2/...
                (Est.samplingFreq*segmentLength);
        end
        % update the data indices
        idx = idx + segmentLength;
    end
end
%% plot the results
% spectrogram settings
segmentLengthRawData = round(segmentTime*samplingFreq);
win = hann(segmentLengthRawData);
nOverlap = round(segmentLengthRawData*3/4);
nFft = 16*segmentLengthRawData;

% do the plotting
timeMidSegments = (1:nSegments)*segmentTime-segmentTime/2;
[S,F,T] = spectrogram(rawData(:,1), win, nOverlap, nFft, samplingFreq);
figure(1)
spec = 10*log10(dynamicRangeLimiting(...
    abs(S).^2/(samplingFreq*segmentLength),60));
imagesc(T,F,spec)
set(gca,'YDir','normal')
hold on
% plot(timeMidSegments, estPitchHz(2,:),'*r','LineWidth',2)
hold off
% legend('F0-AR Est.')
xlabel('Time [s]')
ylabel('Freq. [Hz]')
ylim([0,resamplingFreq/2])
colormap(gray)
axis off
print -dpng results/spectrogram
axis on
colorbar;

figure(2)
plot(timeMidSegments, estPitchOrder, 'o')
legend(methods)
xlabel('Time [s]')

figure(3)
plot(timeMidSegments, estArOrder, 'o')
legend(methods)
xlabel('Time [s]')

figure(4)
plot(timeMidSegments, estPitchHz(1,:), '-', 'LineWidth',2)
hold on
plot(timeMidSegments, estPitchHz(2:end,:), 'o', 'LineWidth',2)
hold off
legend(methods)
xlabel('Time [s]')
ylabel('Freq. [Hz]')

figure(5)
plot(freqVector, 10*log10([perEst, preWhitenedSpectrum, ...
    modelledSpectrum]), 'LineWidth',2)
legend('Periodogram', 'Prewhitened periodogram', 'Modelled spectrum')
xlabel('Freq. [Hz]')
ylabel('PSD [dB/Hz]')
xlim([0, resamplingFreq/2])
ylim([-80,-20])

%% store data in files
comment = '';
legend_str = 'time tacho F0-AR-ML-E F0-ML-E';
data = [timeMidSegments', estPitchHz'];
mtx_to_tbl_pgfplots('results/engineRunUpF0Est',comment,legend_str,data);

idx = (1:exampleNDft/2+1);
comment = '';
legend_str = 'freq per preWhitenedPer modelledPsd';
data = [freqVector(idx), 10*log10([perEst(idx), preWhitenedSpectrum(idx), ...
    modelledSpectrum(idx)])];
mtx_to_tbl_pgfplots('results/engineRunUpExample',comment,legend_str,data);
