[Tokyo,fs] = audioread('Tokyo_vocal.wav');
y = Tokyo(28.5*fs:34*fs);
n = length(y);

figure
window = 1024;
[S,fc,t1] = melSpectrogram(y',fs, ...
                   "Window",hann(window,'periodic'),...
                   "OverlapLength",window/2,...
                   'NumBands',128, ...
                   'FrequencyRange',[0,8000]);
S = 20*log10(S+0.0001);
S = S/max(max(abs(S)))+1;
imagesc(t1, fc, S);
set(gca, 'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
title('Spectrogram')
colorbar

[row,col] = size(S);
SCORE = zeros(row,col-1);
score = zeros(1,col-1);
for i = 1:col-1
    SCORE(:,i) = S(:,i).*S(:,i+1)/(norm(S(:,i)) * norm(S(:,i+1)));
    score(i) = length(find(SCORE(:,i) > 0.001));
end
% score = abs(diff(score));
score = score/max(score);

% S = diff(S,1,2);
% score = zeros(1,col-1);
% for i = 1:col-1
%     score(i) = norm(S(:,i));
% end

figure
tn = [0:n-1]/fs;
plot(tn,y)
hold on
t = t1(1:end-1)+diff(t1)/2;
% t = t1(2:end-1);
plot(t,score,'Marker','o')

% samp = y(0.93*fs:1.09*fs);
% n = length(samp);
% fn = (0:1/n:1-1/n)*fs;
% figure
% S = 20*log10(abs(fft(samp)));
% S = S - max(S);
% plot(fn,S)
% xlim([0 5000])
% 
% samp_w = y(fs:fs+1024);
% n = length(samp_w);
% fn = (0:1/n:1-1/n)*fs;
% figure
% S = 20*log10(abs(fft(samp_w)));
% S = S - max(S);
% plot(fn,S)
% xlim([0 5000])