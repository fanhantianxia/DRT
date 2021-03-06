function [results] = wb_EEG_QA(EEG,WindowSeconds,HighPassband,selechanns,badWindowThreshold,...
          robustDeviationThreshold,PowerFrequency,FrequencyNoiseThreshold,flagNotchFilter,correlationThreshold,...
          ransacCorrelationThreshold,ransacChannelFraction,ransacSampleSize)
% Quality Assessmenting of a continuous EEG data, automatically. The bad data in small windows 
% of each channel could be detected by kinds of 4 methods, and a number of 
% indices related to the data quality will be calculated. Meanwhile, 
% the overall data quality rating will be also provided, including levels 
% of A, B, C, D (corresponding to perfect, good, poor, bad).

% Method 1:Detect constant or NaN/Inf singals in each window.
%
% Method 2: Detecting unusually high or low amplitude using robust standard 
%           deviation across time points in each window (Method 2). If the z 
%           score of robust time deviation falls below ��robustDeviationThreshold�� 
%           or the absolute amplitude exceeds 200 microvolts (��V), the small 
%           window is considered to be bad.
% Method 3: Detecting high or power frequency noises in each window by 
%           calculating the noise-to-signal ratio based on Christian Kothe's 
%           method (Method 3). If the z score of estimate of signal above 40 Hz 
%           (power frequency - 10Hz) to that below 40 Hz above ��highFrequencyNoiseThreshold��
%           or absolute NSR exceeds 0.5, the small window is considered to be bad.
%           Noting that if the sampling rate is below 2*power frequency, 
%           this step will be skipped.
% Method 4a: Detecting low correlations with other channels in each window using 
%           Pearson correlation (Method 4a). For Pearson correlation, If the
%           maximum correlation of the window of a channel to the other channels 
%           falls below ��correlationThreshold��, the window is considered bad.
% Method 4b: Detecting low correlations with other channels in each window using 
%            RANSAC correlation (Method 4b).For RANSAC correlation, each 
%            window of a channel is predicted using RANSAC interpolation based 
%            on a RANSAC fraction of the channels. If the correlation of the 
%            prediction to the actual behavior falls below��ransacCorrelationThreshold��
%            or calculation is too long, the window is marked as bad. The 
%            time cost of this method is high, and the channel locations are 
%            required. The RANSAC correlation is optional and default is not performed.
%
% Assumptions:
%  - The signal is a structure of continuous data with data, srate at least
%  - No segments of the EEG data have been removed.

% Methods 2 and 4b are adapted from code by Christian Kothe and Methods 3
% and 4a are adapted from code by Nima Bigdely-Shamlo
% -------------------------------------------------------------------------
% Input:
%      EEG:  EEG structure loaded by EEGLAB (EEG data is raw data).
%            EEG.data and EEG.srate are required at least.
%      WindowSeconds : window size in seconds (default = 1 sec)
%      HighPassband :  lower edge of the frequency for high pass filtering,
%                      Hz.Default is 1.
%      selechanns : number with indices of the selected channels 
%                      (e.g. [1:4,7:30] or 'all').Default is 'all';
%      badWindowThreshold : cutoff fraction of bad windows (default = 0.4)
%                      for detecting bad channels.
%
%      robustDeviationThreshold : Z-score cutoff for robust channel deviation (default = 5)
%
%      PowerFrequency : power frequency. default is 50 Hz (in Chinese). 
%                      Noting that in USA, power frequency is 60Hz
%      flagNotchFilter : flagNotchFilter = 1: remove 0.5* power frequency 
%                      noise using notch filtering. Default is off (flagNotchFilter = 0).
%      FrequencyNoiseThreshold = : Z-score cutoff for NSR (signal above 40 Hz).
%                                Default is 3;
%
%      correlationThreshold : maximal correlation below which window is bad (range is (0,1), default = 0.6)
% 
%      ransacSampleSize : samples for computing ransac (default = 50)
%      ransacChannelFraction : fraction of channels for robust reconstruction (default = 0.3)
%      ransacCorrelationThreshold : cutoff correlation for abnormal wrt 
%                           neighbors(default = [] | --> not
%                           performed).Default is 0.6.

% Output: a structure array of QA results.

% Qulaity Measures: 

% results.ONS   % Overall ratio of No Signal windows
% results.OHA   % Overall ratio of windows of High Amplitudes
% results.OFN   % Overall ratio of windows of high Frequency and power frequency Noise
% results.OLC   % Overall ratio of windows of Low Correlation
% results.OLRC  % Overall ratio of windows of Low RANSAC Correlation (optional)
% results.badChannels   % Bad channels based on overall bad windows
% results.NBC           % No. of bad channels
% results.OBC           % Overall ratio of Bad Channels
% results.OBClus        % Overall ratio of Bad Clusters 
% results.ODQ           % Overall Data Quality: overall ratio of windows of good data
% results.DataQualityRating  % Overall Data Quality Rating.
%                            % A: ODQ >= 90
%                            % B: ODQ >= 80 && ODQ < 90
%                            % C: ODQ >= 60 && ODQ < 80
%                            % D: ODQ < 60
% 
% results.allMAV      % mean absolute value of all windows
% results.badMAV      % mean absolute value of bad windows
% results.goodMAV     % mean absolute value of good windows
%  
% results.NoSignalMask                 % mask of windows with no signals
% results.AmpliChannelMask             % mask of windows with high amplitudes
% results.FrequencyNoiseMask           % mask of windows with high frequency (and power frequency, if applicable) noise
% results.LowCorrelationMask           % mask of windows with low correlations
% results.RansacBadWindowMask          % mask of windows with RANSAC low correlations
% results.OverallBadMask               % mask of windows with overall bad signals
% results.fractionBadWindows           % fractions of bad windows for each channel
% results.badChannelsFromAll           % bad channels from all methods

% Paramerters:
% results.parameters.srate                    % sampling rate
% results.parameters.WindowSeconds            % window size in seconds (default = 1 sec)
% results.parameters.HighPassband             % lower edge of the frequency for high pass filtering, Hz
% results.parameters.selechanns               % number with indices of the selected channels (e.g. [1:4,7:30] or 'all').Default is 'all';
% results.parameters.badWindowThreshold       % cutoff fraction of bad windows (default = 0.4)
% results.parameters.PowerFrequency           % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
% 
% results.parameters.robustDeviationThreshold     % Z-score cutoff for robust channel deviation (default = 5)
% results.parameters.FrequencyNoiseThreshold      % Z-score cutoff for NSR (signal above 40 Hz).
% results.parameters.correlationThreshold         % maximal correlation below which window is bad (range is (0,1), default = 0.4)
% 
% results.parameters.chanlocsflag                % flag of channel locations.flag == 1: have channel locations.
% results.parameters.chanlocsXYZ                 % xyz coordinates of selected channels
% results.parameters.chanlocs                    % channel locations of selected channels
% results.parameters.ransacSampleSize            % samples for computing ransac (default = 50)
% results.parameters.ransacChannelFraction       % fraction of channels for robust reconstruction (default = 0.3)
% results.parameters.ransacCorrelationThreshold  % cutoff correlation for abnormal wrt neighbors(default = [] | --> not performed)
% -------------------------------------------------------------------------
% usage example: 
% 
% WindowSeconds = 1;       % window size in seconds (default = 1 sec)
% HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
% selechanns = 'all'; %[1:31,33:62];%'all';      % number with indices of the selected channels (e.g. [1:4,7:30] or 'all').Default is 'all';
% badWindowThreshold = 0.4;         % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
% PowerFrequency = 50;              % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
% flagNotchFilter = 0;              % flagNotchFilter = 1: remove 0.5* power frequency noise using notch filtering.
% 
% robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
% FrequencyNoiseThreshold = 3;      % Z-score cutoff for NSR (signal above 40 Hz).Default is 3.
% correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
% 
% ransacSampleSize = 50;            % samples for computing ransac (default = 50)
% ransacChannelFraction = 0.3;      % fraction of channels for robust reconstruction (default = 0.3)
% ransacCorrelationThreshold = [];  % cutoff correlation for abnormal wrt neighbors(default = [] | --> not performed)
% 
%  [results] = wb_EEG_QA(EEG,WindowSeconds,HighPassband,selechanns,badWindowThreshold,...
%          robustDeviationThreshold,PowerFrequency,FrequencyNoiseThreshold,flagNotchFilter,correlationThreshold,...
%          ransacCorrelationThreshold,ransacChannelFraction,ransacSampleSize)
% -------------------------------------------------------------------------
% Written by Li Dong,UESTC (Lidong@uestc.edu.cn)
% $ 2019.10.2
% -------------------------------------------------------------------------
if nargin < 1
    error ('One input is reqiured at least!!!!!');
elseif nargin == 1
    WindowSeconds = 1;       
    HighPassband = 1;       
    selechanns = 'all';      
    badWindowThreshold = 0.4;         
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
    ransacCorrelationThreshold = [];         
    ransacChannelFraction = 0.3; 
    ransacSampleSize = 50;
elseif nargin == 2
    HighPassband = 1;       
    selechanns = 'all';      
    badWindowThreshold = 0.4;         
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
    ransacCorrelationThreshold = [];         
    ransacChannelFraction = 0.3; 
    ransacSampleSize = 50;
elseif nargin == 3
    selechanns = 'all';      
    badWindowThreshold = 0.4;         
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
    ransacCorrelationThreshold = [];         
    ransacChannelFraction = 0.3; 
    ransacSampleSize = 50;
elseif nargin == 4
    badWindowThreshold = 0.4;         
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
    ransacCorrelationThreshold = [];         
    ransacChannelFraction = 0.3; 
    ransacSampleSize = 50;
elseif nargin == 5
    robustDeviationThreshold = 5;
    PowerFrequency = 50;
    FrequencyNoiseThreshold = 3;
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
    ransacCorrelationThreshold = [];
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 6
    PowerFrequency = 50;
    FrequencyNoiseThreshold = 3;
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
    ransacCorrelationThreshold = [];
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 7
    FrequencyNoiseThreshold = 3;
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
    ransacCorrelationThreshold = [];
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 8
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
    ransacCorrelationThreshold = [];
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 9
    correlationThreshold = 0.6;
    ransacCorrelationThreshold = [];
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 10
    ransacCorrelationThreshold = [];
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 11
    ransacChannelFraction = 0.3;
    ransacSampleSize = 50;
elseif nargin == 12
    ransacSampleSize = 50;
end
% -------------------------------------------------------------------------
% checking all inputs
if isempty(WindowSeconds)
    WindowSeconds = 1;
end
if isempty(HighPassband)
    HighPassband = 1;
end
if isempty(selechanns)
    selechanns= 'all';
end
if isempty(badWindowThreshold)
    badWindowThreshold = 0.4;
end
if isempty(robustDeviationThreshold)
    robustDeviationThreshold = 5;
end
if isempty(PowerFrequency)
    PowerFrequency = 50;
end
if isempty(FrequencyNoiseThreshold)
    FrequencyNoiseThreshold = 3;
end
if isempty(flagNotchFilter)
    flagNotchFilter = 0;
end
if isempty(correlationThreshold)
    correlationThreshold = 0.6;
end
if isempty(ransacChannelFraction)
    ransacChannelFraction = 0.3;
end
if isempty(ransacSampleSize)
    ransacSampleSize = 50;
end
% check sampling rate
try 
    srate = EEG.srate; % sampling rate
    if isfinite(srate)
        disp(['sampling rate = ',num2str(srate)]);
    else
        disp('sampling rate is invalide');
        error('sampling rate is invalide');
    end
catch
    disp('sampling rate is not found');
    error('sampling rate is not found');
end

% check EEG data
try
    if isempty(EEG.data)
        error('EEG.data is empty!!!!!');
    end
catch
    error('EEG.data is not exist!!!!');
end

% check channs
if isequal(selechanns,'all')
    selechanns = 1:size(EEG.data,1);
end
Nchan = length(selechanns); % No. of selected channs
disp(['No. of selected channels: ',num2str(Nchan)]);

% check channel locations
chanlocsflag = 0;
try
    if isfield(EEG,'chanlocs')
        if ~(isfield(EEG.chanlocs,'X') && isfield(EEG.chanlocs,'Y') && isfield(EEG.chanlocs,'Z') && all([length([EEG.chanlocs.X]),length([EEG.chanlocs.Y]),length([EEG.chanlocs.Z])] > length(EEG.chanlocs)*0.5))
            chanlocsflag = 0;
            warning('because most of your channels do not have X,Y,Z location measurements, the RANSAC method is ignored');
        else
            % get the matrix of selected channel locations [3xN]
            [x,y,z] = deal({EEG.chanlocs.X},{EEG.chanlocs.Y},{EEG.chanlocs.Z});
            x = x(selechanns);
            y = y(selechanns);
            z = z(selechanns);
            usable_channels = find(~cellfun('isempty',x) & ~cellfun('isempty',y) & ~cellfun('isempty',z));
            chanlocsXYZ = [cell2mat(x(usable_channels));cell2mat(y(usable_channels));cell2mat(z(usable_channels))];
            chanlocsflag = 1;
            if size(chanlocsXYZ,2) ~= Nchan
                chanlocsflag = 0;
                warning('some selected channels do not have xyz coordinates, the RANSAC method is ignored');
            end
            
        end
    else
        chanlocsflag = 0;
        warning('channel location is not found, the RANSAC method is ignored');
    end
catch
     warning('there are some unknow errors of channel locations, the RANSAC method is ignored');
end
% check other parameters

% =========================================================================
% Detecting bad segements using several methods
disp('------------');
disp('Detecting bad windows using several methods...');
disp('------------');
% --------------------------------------
% high pass filtering
disp('High Pass filtering...');
EEG = pop_eegfiltnew(EEG,HighPassband,[],[],0);
% --------------------------------------
% get the data required
Nt = size(EEG.data,2);
winlenth = srate * WindowSeconds; % window length (time points)
Nwin = floor(Nt/winlenth); % No. of windows 
selectdata = EEG.data(selechanns,1:winlenth*Nwin);
reshapedata = reshape(selectdata,Nchan,winlenth,Nwin);
disp(['No. of time points: ',num2str(Nt)]);
disp(['No. of windows: ',num2str(Nwin)]);
% --------------------------------------
% Initialization
% NoSignalMask = zeros(Nchan,Nwin);
% HighAmpliMask = zeros(Nchan,Nwin);
% 
% LowCorrelationMask = zeros(Nchan,Nwin);
% --------------------------------------
% Method 1:
% Detect constant or NaN/Inf signals in each window
disp('------------');
disp('Detecting constant, Inf, or NaN channels......');
median1 = reshape(mad(reshapedata, 1, 2),Nchan,Nwin);
std1 = reshape(std(reshapedata, 1, 2),Nchan,Nwin);
NanSignalMask = reshape(sum(~isfinite(reshapedata),2),Nchan,Nwin);
NoSignalMask = double( median1 < 10e-10 | std1 < 10e-10) + NanSignalMask;
% ---------------------------------------
% Method 2:
% Detect unusually high or low amplitude using robust STD
disp('------------');
disp('Detecting unusually high or low amplitude using robust STD......');
disp(['Robust deviation threshold: ',num2str(robustDeviationThreshold)]);

index1 = abs(reshapedata) > 200;  % absolute amplitude > 200 ��V
high1 = reshape(sum(index1,2),Nchan,Nwin) > 0;

channelDeviation = reshape(0.7413 *iqr(reshapedata,2),Nchan,Nwin); % Robust estimate of SD
channelDeviationSD = 0.7413*iqr(channelDeviation(:));
channelDeviationMedian = nanmedian(channelDeviation(:),1);
robustChannelDeviation = (channelDeviation - channelDeviationMedian) / channelDeviationSD;
HighAmpliMask = abs(robustChannelDeviation) > robustDeviationThreshold | isnan(robustChannelDeviation) | high1;
% ----------------------------------------
% Method 3: Compute the noise-to-signal ratio (based on Christian Kothe's clean_channels)
% Note: RANSAC and global correaltion uses the filtered values X of the data
disp('------------');
if flagNotchFilter == 1
    disp('Detecting high frequency noise and power frequency noise using noise-to-signal ratio......');
else
    disp('Detecting high frequency noise using noise-to-signal ratio......');
end

disp(['Frequency noise threshold: ',num2str(FrequencyNoiseThreshold)]);
% FrequencyNoiseMask = zeros(Nchan,Nwin);
if srate > 2*PowerFrequency
    % Remove signal content above 40Hz/50Hz and below 1 Hz
    disp('low Pass filtering...');
    EEG1 = pop_eegfiltnew(EEG,[],PowerFrequency - 10,[],0); % In Chinese, the power frequency is 50Hz, so set the high pass frequency as 50-10=40.
    if flagNotchFilter == 1
        disp('Notch filtering for 0.5*power frequency...');
        EEG1 = pop_eegfiltnew(EEG1,0.5*PowerFrequency - 5, 0.5*PowerFrequency + 5,[],1); % notch filtering for 0.5* power frequency
    end
    % checking high frequency noise 
    % -------------- 
    X = EEG1.data(selechanns,1:winlenth*Nwin);
    X = reshape(X,Nchan,winlenth,Nwin);
    % Determine z-scored level of EM noise-to-signal ratio for each window
    noisiness = mad(reshapedata - X, 1, 2)./mad(X, 1, 2);
    noisiness = reshape(noisiness,Nchan,Nwin);
        
    noisinessMedian = nanmedian(noisiness(:));
    noisinessSD = mad(noisiness(:), 1)*1.4826; % median absolute deviation
    zscoreFreNoiseTemp = (noisiness - noisinessMedian) ./ noisinessSD;
        
%     noisinessMedian = nanmedian(noisiness);
%     noisinessSD = mad(noisiness, 1)*1.4826;
%     zscoreFreNoiseTemp = bsxfun ( @minus, noisiness, noisinessMedian);
%     zscoreFreNoiseTemp = bsxfun ( @rdivide, zscoreFreNoiseTemp,noisinessSD);
    
    FrequencyNoiseMask = (abs(zscoreFreNoiseTemp) > FrequencyNoiseThreshold) | isnan(zscoreFreNoiseTemp) | (abs(noisiness) > 0.5); % or the absolute noise-to-signal ratio > 0.5
    FrequencyNoiseMask = FrequencyNoiseMask .* (abs(noisiness) > 0.0075) == 1;   % the error between signals of twice low pass fitering is about 0.0075, so the noisiness < 0.0075 may be indistinguishable.
                                                                                          
else
    warning('The sampling rate is below 2*PowerFrequency (too low), detecting high frequency noise is skipped');
    X = EEG.data(selechanns,1:winlenth*Nwin);
    X = reshape(X,Nchan,winlenth,Nwin);
    FrequencyNoiseMask = zeros(Nchan,Nwin);
end
    
% -----------------------------------------
% Method 4: Global correlation criteria in time domain (from Nima Bigdely-Shamlo)
disp('------------');
disp('Detecting low correlation with other channels......');
disp(['correlation threshold: ',num2str(correlationThreshold)]);

channelCorrelations = zeros(Nchan,Nwin);
for k1 = 1:Nwin 
    eegPortion = squeeze(X(:, :, k1))'; % using filtered data X
    windowCorrelation = corrcoef(eegPortion);
    abs_corr = abs(windowCorrelation - eye(Nchan,Nchan));
    channelCorrelations(:,k1)  = quantile(abs_corr, 0.98); % approximate maximal correlation: quantile of 98%
end;

dropOuts = isnan(channelCorrelations) | isnan(noisiness);
channelCorrelations(dropOuts) = 0;
% noisiness(dropOuts) = 0;
LowCorrelationMask = channelCorrelations < correlationThreshold;

% -------------------------------------------------
% Method 5: Detecting low correlation using RANSAC correlation (may not be performed if channel location is empty)
if isempty(ransacCorrelationThreshold) || ~isfinite(ransacCorrelationThreshold)
  RansacBadWindowMask = zeros(Nchan,Nwin);
else
  if chanlocsflag == 1 % if have channel location
    disp('------------');
    disp('Detecting low correlation using RANSAC method......');
    disp(['RANSAC Correlation Threshold:',num2str(ransacCorrelationThreshold)]);
    
    subset_size = round(ransacChannelFraction * Nchan);
    
    % caculate all-channel reconstruction matrices from random channel subsets
    P = hlp_microcache('cleanchans',@calc_projector,chanlocsXYZ,ransacSampleSize,subset_size);
    RansacCorrelation = zeros(Nchan,Nwin);
    
    % calculate each channel's correlation to its RANSAC reconstruction for each window
    % timePassedList = zeros(Nwin,1);
    for iw = 1:Nwin
      % tic; % makoto
      XX = X(:,:,iw)';
      YY = sort(reshape(XX*P,winlenth,Nchan,ransacSampleSize),3);
      YY = YY(:,:,round(size(YY,3)/2));
      RansacCorrelation(:,iw) = sum(XX.*YY)./(sqrt(sum(XX.^2)).*sqrt(sum(YY.^2)));
      %         timePassedList(iw) = toc; % makoto
      %         medianTimePassed = median(timePassedList(1:iw));
      %         disp(sprintf('clean_channel: %3.0d/%d, %.1f minutes remaining.', iw,Nwin, medianTimePassed*(Nwin-iw)/60)); % makoto
    end
    RansacBadWindowMask = RansacCorrelation < ransacCorrelationThreshold | isnan(RansacCorrelation);
  else
    RansacBadWindowMask = zeros(Nchan,Nwin);
  end
end

% =========================================================================
% quality assessment
disp('---------------');
disp('Quality Assessment...');
disp('Detecting bad channels...');
disp(['Bad Window Threshold:',num2str(badWindowThreshold)]);
OverallBadMask = NoSignalMask + HighAmpliMask + FrequencyNoiseMask + LowCorrelationMask + RansacBadWindowMask; % considering all methods
OverallBadMask = OverallBadMask > 0;
fractionBadWindows = sum(OverallBadMask,2)/Nwin;
badChannelsFromAll = fractionBadWindows > badWindowThreshold; % bad channels

% Calculate the quality measures
N_all = Nchan*Nwin;
% overall ratio of no signal windows
ONS = nansum(NoSignalMask(:))/N_all; 

% overall ratio of windows of high amplitude
OHA = nansum(HighAmpliMask(:))/N_all;

% overall ratio of windows of high frequency and power frequency noise
OFN = nansum(FrequencyNoiseMask(:))/N_all;

% overall ratio of windows of low correlation
OLC = nansum(LowCorrelationMask(:))/N_all;

% overall ratio of windows of low RANSAC correlation (optional)
if chanlocsflag == 1 && ~isempty(ransacCorrelationThreshold) && isfinite(ransacCorrelationThreshold)
    OLRC = nansum(RansacBadWindowMask(:))/N_all;
else
    OLRC = [];
end

% bad channels based on overall bad windows
badChannels = selechanns(badChannelsFromAll);

% No. of bad channels
NBC = nansum(badChannelsFromAll);

% overall ratio of bad channels 
OBC = NBC/Nchan;

% Overall data quality: overall ratio of windows of good data  
ODQ = 100*(1 - (nansum(OverallBadMask(:))/N_all));

% Overall ratio of bad clusters 
L_temp1 = bwlabel(OverallBadMask);
OBClus = max(L_temp1(:))/(round(size(L_temp1,1)/2)*(round(size(L_temp1,2)/2)-1));

% unthresholded mean absolute voltage
M1 = reshape(nanmean(abs(reshapedata),2),Nchan,Nwin);
allMAV = nanmean(abs(reshapedata(:)));     % mean absolute value of all windows
badMAV = nanmean(M1(OverallBadMask==1));   % mean absolute value of bad windows
goodMAV = nanmean(M1(OverallBadMask==0));  % mean absolute value of good windows
% --------------------------------
% data quality rating
DataQualityRating = [];
if ODQ < 60
    DataQualityRating = 'D';
elseif ODQ >= 60 && ODQ < 80
    DataQualityRating = 'C';
elseif ODQ >= 80 && ODQ < 90
    DataQualityRating = 'B';
elseif ODQ >= 90
    DataQualityRating = 'A';
end

% % plot check
% figure;
% subplot(3,3,1);imagesc(NoSignalMask);title('NoSignalMask');xlabel('windows');ylabel('channels');
% subplot(3,3,2);imagesc(HighAmpliMask);title('HighAmplitudeMask');xlabel('windows');ylabel('channels');
% subplot(3,3,3);imagesc(FrequencyNoiseMask);title('FrequencyNoiseMask');xlabel('windows');ylabel('channels');
% subplot(3,3,4);imagesc(LowCorrelationMask);title('LowCorrelationMask');xlabel('windows');ylabel('channels');
% subplot(3,3,5);imagesc(RansacBadWindowMask);title('RansacBadWindowMask');xlabel('windows');ylabel('channels');
% subplot(3,3,6);imagesc(OverallBadMask);title('OverallBadMask');xlabel('windows');ylabel('channels');
% 
% subplot(3,3,7);bar(fractionBadWindows);grid on;title('FractionOfBadWindows');xlabel('channels');ylabel('raito');
% hold on; plot(0:length(fractionBadWindows)+2,ones(length(fractionBadWindows)+3)*badWindowThreshold,'--r');hold off;
% subplot(3,3,8);bar(double(badChannelsFromAll));grid on;title('OverallBadChannels');xlabel('channels');ylabel('values');
% subplot(3,3,9);bar([allMAV,badMAV,goodMAV]);grid on;title('Mean Absolute Values');ylabel('values');set(gca,'XTickLabel',{'allMAV','badMAV','goodMAV'});
% -------------------------------------------------------------------------
% save quality measures and parameters


results.ONS = ONS;  % Overall ratio of No Signal windows
results.OHA = OHA;  % Overall ratio of windows of High Amplitudes
results.OFN = OFN;  % Overall ratio of windows of high Frequency and power frequency Noise
results.OLC = OLC;  % Overall ratio of windows of Low Correlation
results.OLRC = OLRC; % Overall ratio of windows of Low RANSAC Correlation (optional)
results.badChannels = badChannels;  % Bad channels based on overall bad windows
results.NBC = NBC;                  % No. of bad channels
results.OBC = OBC;                  % Overall ratio of Bad Channels 
results.ODQ = ODQ;                  % Overall Data Quality: overall ratio of windows of good data
results.OBClus = OBClus;            % Overall ratio of bad clusters 
results.DataQualityRating = DataQualityRating; % Overall Data Quality Rating.
                                               % A: ODQ >= 90
                                               % B: ODQ >= 80 && ODQ < 90
                                               % C: ODQ >= 60 && ODQ < 80
                                               % D: ODQ < 60

results.allMAV = allMAV;     % mean value of all windows
results.badMAV = badMAV;     % mean value of bad windows
results.goodMAV = goodMAV;   % mean value of good windows

results.NoSignalMask = NoSignalMask;                % mask of windows with no signals
results.HighAmpliMask = HighAmpliMask;              % mask of windows with high amplitudes
results.FrequencyNoiseMask = FrequencyNoiseMask;    % mask of windows with high frequency noise
results.LowCorrelationMask = LowCorrelationMask;    % mask of windows with low correlations
results.RansacBadWindowMask = RansacBadWindowMask;  % mask of windows with RANSAC low correlations
results.OverallBadMask = OverallBadMask;            % mask of windows with overall bad signals
results.fractionBadWindows = fractionBadWindows;    % fractions of bad windows for each channel
results.badChannelsFromAll = badChannelsFromAll;    % bad channels from all methods.

% paramerters
results.parameters.srate = srate;                       % sampling rate
results.parameters.WindowSeconds = WindowSeconds;       % window size in seconds (default = 1 sec)
results.parameters.HighPassband = HighPassband;         % lower edge of the frequency for high pass filtering, Hz
results.parameters.selechanns = selechanns;             % number with indices of the selected channels (e.g. [1:4,7:30] or 'all').Default is 'all';
results.parameters.badWindowThreshold = badWindowThreshold;  % cutoff fraction of bad windows (default = 0.4)
results.parameters.PowerFrequency = PowerFrequency;          % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz

results.parameters.robustDeviationThreshold = robustDeviationThreshold;    % Z-score cutoff for robust channel deviation (default = 4)
results.parameters.FrequencyNoiseThreshold = FrequencyNoiseThreshold;      % Z-score cutoff for NSR (signal above 40 Hz).
results.parameters.correlationThreshold = correlationThreshold;            % correlation below which window is bad (range is (0,1), default = 0.4)

results.parameters.chanlocsflag = chanlocsflag;
try results.parameters.chanlocsXYZ = chanlocsXYZ; catch;end;
results.parameters.chanlocs = EEG.chanlocs(1,selechanns);          % channel locations of selected channels
results.parameters.ransacSampleSize = ransacSampleSize;            % samples for computing ransac (default = 50)
results.parameters.ransacChannelFraction = ransacChannelFraction;      % fraction of channels for robust reconstruction (default = 0.25)
results.parameters.ransacCorrelationThreshold = ransacCorrelationThreshold; % cutoff correlation for abnormal wrt neighbors(default = [] | --> not performed)
