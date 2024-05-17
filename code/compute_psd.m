% method: welch, fft
% data 
% fs
% window_sec
% overlap_percent
% freqBandHz: [low high]
% normalize_flg: create pdf
function [P, f] = compute_psd(method, data, fs, window_sec, overlap_percent, freqBandHz, normalize_flg)

N = window_sec*fs;
specWindow = ones(N,1);
% specWindow = hamming(N,'periodic'); %hann hamming blackman
specOverlap = round(N * overlap_percent);

[nElec,L] = size(data);
if L<N
    error('data too short!');
end

f = (0:floor(N/2))/window_sec;
switch method
    case 'welch'
        P = pwelch(data',specWindow,specOverlap,f,fs);
        if nElec > 1
            P = P';
        end
    case 'fft'
        for iElec = 1:nElec
            [buffData,~] = buffer(data(iElec,:),N,specOverlap,"nodelay" );
            buffDataPower = (abs(fft(buffData.*specWindow))/N).^2;
            P(iElec,:) = mean(buffDataPower(1:floor(N/2)+1,:),2)';
        end
end

P = double(P(:,f>=freqBandHz(1) & f<=freqBandHz(2)));
f = f(f>=freqBandHz(1) & f<=freqBandHz(2));

if normalize_flg
    P = P./sum(P,2);
end
