clear
clc
close all

% add paths to the data and the estimator
addpath('../../estimator/');
addpath('../../lib/');

% load the data
filename = 'carEngineRunUp';
[rawData, samplingFreq] = audioread([filename, '.wav']);

% resample the data and extract the left channel
resamplingFreq = samplingFreq/20; % Hz
resampledData = resample(rawData(:,1), resamplingFreq, samplingFreq);


% set up the segment-by-segment processing (no overlap)
nData = length(resampledData);
segmentTime = 0.25; % seconds
segmentLength = round(segmentTime*resamplingFreq);
nSegments = floor(nData/segmentLength);

% spectrogram settings
segmentLengthRawData = round(segmentTime*samplingFreq);
win = hann(segmentLengthRawData);
nOverlap = round(segmentLengthRawData*3/4);
nFft = 16*segmentLengthRawData;
[S,F,T] = spectrogram(rawData(:,1), win, nOverlap, nFft, samplingFreq);
spec = 10*log10(dynamicRangeLimiting(...
    abs(S).^2/(samplingFreq*segmentLength),60));

% set up the analysis
maxArOrder = 5;
maxPitchOrder = 25;
pitchBounds = [10, 100]; % Hz

% do the analysis
estArOrder = nan(1,nSegments);
estPitchOrder = nan(1,nSegments);
estPitchHz = nan(1,nSegments);
estRpm = nan(1,nSegments);

Est = FastF0ArMl(segmentLength, maxPitchOrder, ...
    maxArOrder, pitchBounds, resamplingFreq, 'E');
Est.setValidPitchOrders(20:maxPitchOrder);

    
% prepare for plotting
timeMidSegments = (1:nSegments)*segmentTime-segmentTime/2;
freqVector = resamplingFreq*(0:Est.nDft-1)'/Est.nDft;
figure(1);
subplot(2,2,1)
imagesc(T,F,spec)
set(gca,'YDir','normal')
hold on
h11 = plot(timeMidSegments, estPitchHz,'*r','LineWidth',2);
h11.XDataSource = 'timeMidSegments';
h11.YDataSource = 'estPitchHz';
hold off
xlabel('Time [s]')
ylabel('Freq. [Hz]')
ylim([0,resamplingFreq/2])
colormap(gray)
colorbar;

subplot(2,2,3)
h21 = plot(timeMidSegments, estRpm ,'*');
h21.XDataSource = 'timeMidSegments';
h21.YDataSource = 'estRpm';
xlim([0,nData/resamplingFreq])
ylim(60*pitchBounds)
xlabel('Time [s]')
ylabel('RPM [r/min]')

subplot(2,2,2)
h12 = plot(freqVector, nan(length(freqVector),2), 'LineWidth',2);
h12(1).XDataSource = 'freqVector';
h12(2).XDataSource = 'freqVector';
per = nan(Est.nDft,1);
h12(1).YDataSource = 'per';
estPsd = nan(Est.nDft,1);
h12(2).YDataSource = 'estPsd';
legend('Periodogram', 'Modelled spectrum')
xlabel('Freq. [Hz]')
ylabel('PSD [dB/Hz]')
xlim([0, resamplingFreq/2])
ylim([-80,-20])

subplot(2,2,4)
h22 = plot(freqVector, nan(length(freqVector),2), 'LineWidth',2);
h22(1).XDataSource = 'freqVector';
h22(2).XDataSource = 'freqVector';
h22(1).YDataSource = 'per';
prewhitenedPer = nan(Est.nDft,1);
h22(2).YDataSource = 'prewhitenedPer';
legend('Periodogram', 'Pre-whitened spectrum')
xlabel('Freq. [Hz]')
ylabel('PSD [dB/Hz]')
xlim([0, resamplingFreq/2])
ylim([-80,-20])

idx = 1:segmentLength;
for ii = 1:nSegments
    disp(['Processing segment ', num2str(ii), ' of ', ...
        num2str(nSegments)]);
    dataSegment = resampledData(idx);
    % estimate parameters and orders
    estPitchHz(ii) = Est.estimate(dataSegment);
    estPitchOrder(ii) = Est.estPitchOrder;
    estArOrder(ii) = Est.estArOrder;
    estRpm(ii) = estPitchHz(ii)*60;
    % compute spectra
    per = 10*log10(abs(fft(dataSegment, Est.nDft)).^2/...
        (resamplingFreq*segmentLength));
    estPsd = 10*log10(Est.computeModelledSpectrum);
    prewhitenedPer = 10*log10(Est.computePrewhitenedSpectrum);
    % update the plots
    refreshdata
    drawnow
    % update the data indices
    idx = idx + segmentLength;
end
