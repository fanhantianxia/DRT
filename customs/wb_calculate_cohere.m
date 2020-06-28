function Cohere = wb_calculate_cohere(data,srate,FreqBand,NFFT)
% Calculate Magnitude squared coherence using Welch's averaged modified 
% periodogram method.
% Input:
%      data: The EEG potentials,channels X time points. e.g.
%            62 channels X 5000 time points.
%      srate: sampling rate of data.
%      FreqBand: frequency bands. e.g. [1,4] or [1,4;4,8]; If FreqBand is
%            empty, coherence of fullband will be calculated.
%      NFFT: FFT length which determines the frequencies at which the
%            coherence is estimated.Maximum of 256 or the next power of 2 
%            greater than the length of data.
%            Default is NFFT = 2^nextpow2(length of data);
% Output:
%      Cohere: connection matrix of coherence
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Li_dong729@163.com)
% $ 2018.1.22
% -------------------------------------------------------------------------

if nargin < 2
    error ('Two inputs are reqiured at least!!!!!');
elseif nargin == 2
    FreqBand = [];
    NFFT = 2^nextpow2(size(data,2));
elseif nargin == 3
    NFFT = 2^nextpow2(size(data,2));
end
% ----------------------
[N_channs,~] = size(data);
if ~isempty(FreqBand)
    N_Freq = size(FreqBand,1);
    Cohere = zeros(N_channs,N_channs,N_Freq);% connection matrix
else
    Cohere = zeros(N_channs,N_channs);% connection matrix
end

if isempty(NFFT)
    NFFT = 2^nextpow2(size(data,2));
end

for k1 = 1: N_channs-1
    for k2 = k1+1:N_channs
        [temp_coher,temp_Freq] = mscohere(data(k1,:),data(k2,:),[],[],NFFT,srate); % default Hamming window, default noverlap to obtain 50% overlap.
        if ~isempty(FreqBand)
            for k3 = 1:N_Freq
                temp1 = temp_coher(temp_Freq >= FreqBand(k3,1) & temp_Freq < FreqBand(k3,2));
                Cohere(k1,k2,k3) = mean(temp1(~isnan(temp1)));
            end
        else
            Cohere(k1,k2) = mean(temp_coher(~isnan(temp_coher)));
        end
    end
end