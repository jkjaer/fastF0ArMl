clear
clc
close all

resamplingFreq = 16000;
% load speech and resample left channel
[rawSpeech, speechSamplingFreq] = audioread('FB07_01.wav');
speech = resample(rawSpeech(:,1), resamplingFreq, speechSamplingFreq);
nData = length(speech);
speechPower = speech'*speech;
% load noise and resample left channel
[rawNoise, noiseSamplingFreq] = audioread('mic2_10cm.wav');
noiseLong = resample(rawNoise(:,1), ...
    resamplingFreq, noiseSamplingFreq);
% I asssume that the length of the noise signal is longer than
% that of the speech signal
noise = noiseLong(1:nData);
noisePower = noise'*noise;
% desired SNRs
snrListDb = [-10:5:25, 100];
nSnrs = length(snrListDb);
% generate the audio files
for ii = 1:nSnrs
    % generate noisy speech with the desired SNR
    noiseGain = sqrt(10^(-snrListDb(ii)/10)*...
        (speechPower/noisePower));
    noisySpeech = speech+noiseGain*noise;
    % ensure that the biggest value is smaller than 1 (otherwise clipping
    % is happening in the wav-file)
    if max(abs(noisySpeech)) >= 1
        noisySpeech = noisySpeech/(max(abs(noisySpeech))+eps);
    end
    % store the file
    filename = ['speechInWind_',num2str(snrListDb(ii)),'dB.wav'];
    audiowrite(filename, noisySpeech, resamplingFreq);
end
