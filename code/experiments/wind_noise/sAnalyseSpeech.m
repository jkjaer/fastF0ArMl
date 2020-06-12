clear
clc
close all

% add paths to the data and the estimator
addpath('../../estimator/');
addpath('../../lib/');
addpath data

% load the data
[rawData, samplingFreq] = audioread('speechInWind_5dB.wav');
[rawCleanData, cleanSamplingFreq] = audioread('speechInWind_100dB.wav');

% resample the data and extract the left channel
resamplingFreq = 16000; % Hz
resampledData = resample(rawData(:,1), resamplingFreq, samplingFreq);
resampledCleanData = resample(rawCleanData(:,1), resamplingFreq, ...
    cleanSamplingFreq);

% set up the segment-by-segment processing (no overlap)
nData = length(resampledData);
segmentTime = 0.025; % seconds
segmentLength = round(segmentTime*resamplingFreq);
nSegments = floor(nData/segmentLength);

% set up the analysis
maxArOrder = 10;
maxPitchOrder = 15;
pitchBounds = [60, 400]; % Hz

% do the analysis
methods = {'Clean', 'F0-AR-ML-E', 'F0-ML-E'};
nMethods = length(methods);
estArOrder = nan(nMethods,nSegments);
estPitchOrder = nan(nMethods,nSegments);
estPitchHz = nan(nMethods,nSegments);
exampleSegmentNo = 55;
exampleNDft = 2048;
for ii = 1:nMethods
    disp(['Processing method ', num2str(ii), ' of ', ...
            num2str(nMethods)]);
    switch methods{ii}
        case 'Clean'
            Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
                maxArOrder, pitchBounds, resamplingFreq, 'E');
            obsData = resampledCleanData;
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
    Est.setValidPitchOrders([0,3:maxPitchOrder]);
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
        idx = idx + segmentLength;
    end
end

%% plot the results
% spectrogram settings
segmentLengthRawData = round(segmentTime*samplingFreq);
win = hann(segmentLengthRawData);
nOverlap = round(segmentLengthRawData*3/4);
nFft = 4*segmentLengthRawData;

% do the plotting
timeMidSegments = (1:nSegments)*segmentTime-segmentTime/2;
[S,F,T] = spectrogram(rawData(:,1), win, nOverlap, nFft, samplingFreq);
figure(1)
spec = 10*log10(dynamicRangeLimiting(...
    abs(S).^2/(samplingFreq*segmentLength),80));
imagesc(T,F,spec)
set(gca,'YDir','normal')
xlabel('Time [s]')
ylabel('Freq. [Hz]')
ylim([0,resamplingFreq/2])
colormap(gray)
axis off
print -dpng results/spectrogram
axis on
colorbar;

figure(2)
plot(timeMidSegments, estPitchOrder, 'o', 'LineWidth',2)
legend(methods)
xlabel('Time [s]')

figure(3)
plot(timeMidSegments, estArOrder, 'o', 'LineWidth',2)
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
ylim([-100,-40])

%% store data in files
comment = '';
legend_str = 'time clean F0-AR-ML-E F0-ML-E';
data = [timeMidSegments', estPitchHz'];
mtx_to_tbl_pgfplots('results/speechInWind5DbF0Est',comment,legend_str,data);

idx = (1:exampleNDft/2+1);
comment = '';
legend_str = 'freq per preWhitenedPer modelledPsd';
data = [freqVector(idx), 10*log10([perEst(idx), preWhitenedSpectrum(idx), ...
    modelledSpectrum(idx)])];
mtx_to_tbl_pgfplots('results/speechInWind5DbExample',comment,legend_str,data);
