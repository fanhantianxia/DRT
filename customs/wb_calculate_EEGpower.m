function EEG_results= wb_calculate_EEGpower(EEG,epochLenth,eventlabel,bandLimit,bandName,selechanns,proportion)
% Calculate power indices using time-frequency analysis of EEGALB.
% Calculating power indices consists of
% [1] Specific event data are first extracted according to the input
%     'eventlabel'. If the input 'eventlabel' is empty, all data will be used.
%     If applicable, EEG segments in bad block (label 9999, marked by wb_pipeline_EEG_Mark)
%     will also be rejected automatically, and NOT used to calculate power indices.
% [2] Specific event EEG signals will be divided into small epochs.
% [3] EEG data of each epoch (default is 5s epoch) was subjected to
%     time-frequency analysis with Fast-Fourier Transform (FFT) to obtain
%     the absolute EEG band power at each electrode in the specific bands.

% Noting that data in bad blocks (labeled as 9999) will not be used to
% calculate power indices !!!
% Power in specific band = 2*abs(complex number calculated by FFT).^2/length(epoch);
% Default list of 25 indices:
%    delta: mean power acorss 1 - 4 Hz
%    theta: mean power acorss 4 - 8 Hz
%    alpha1: mean power acorss 8 - 10.5 Hz
%    alpha2: mean power acorss 10.5 - 12.5 Hz
%    beta1: mean power acorss 12.5 - 18.5 Hz
%    beta2: mean power acorss 18.5 - 21 Hz
%    beta3: mean power acorss 21 - 30 Hz
%    gamma1: mean power acorss 30 - 40 Hz
%    gamma2: mean power acorss 40 - 60 Hz
%    fullband: mean power acorss 1-60 Hz
%    relative power in delta/theta/alpha1/alpha2/beta1/beta2/beta3/gamma1/gamma2 band
%       = power of specific band/total power across fullband.
%    r1 = theta/(alpha1 + alpha2 + beta1);
%    r2 = (delta + theta)/(alpha1 + alpha2 + beta1 + beta2);
%    r3 = theta/alpha =  theta/(alpha1+alpha2);
%    r4 = theta/beta = theta/(beta1 + beta2 + beta3);
%    r5 = delta/theta;
%    r6 = alpha/beta = (alpha1 + alpha2)/(beta1 + beta2 + beta3);
%    PAF (peak of alpha frequency) = max power in alpha(alpha1+alpha2) band.
% References:
%   Nuwer, M. R., et al. (1994). "IFCN Guidelines for Topographic and Frequency-Analysis of EEGs and EPs - Report of an IFCN Committee." Electroencephalography and Clinical Neurophysiology 91(1): 1-5.
%   Thatcher, R. W., et al. (2005). "EEG and intelligence: Relations between EEG coherence, EEG phase delay and power." Clinical Neurophysiology 116(9): 2129-2141.
%   Chen, A. C., et al. (2008). "EEG default mode network in the human brain: spectral regional field powers." Neuroimage 41(2): 561-574.
%   Snaedal, J., et al. (2010). "The use of EEG in Alzheimer's disease, with and without scopolamine - A pilot study." Clinical Neurophysiology 121(6): 836-841.
%   Malver, L. P., et al. (2014). "Electroencephalography and analgesics." Br J Clin Pharmacol 77(1): 72-95.
%   Kashefpoor, M., et al. (2016). "Automatic Diagnosis of Mild Cognitive Impairment Using Electroencephalogram Spectral Features." J Med Signals Sens 6(1): 25-32.
% -------------------------------------------------------------------------
% Input:
%    EEG: EEG structure imported using EEGLAB. EEG.data should be channels
%         X time points OR channels X time points X epoches.
%    epochLenth: length of small epochs to calculate power (no overlapped).
%         unit is second.Default is 5s. If epochLenth is negative, it means
%         that if possible, data before event labels (eventlabel) will be used
%         (no overlapped).
%    eventlabel: Event label which means specific event data. Default is empty.
%        If it is empty, all data will be used to calculate indices. If
%        eventlabel is not found in events, NO data will be epoched and
%        calculated. If structure event (eventlabel) doesn¡¯t include duration,
%        the duration will be equal to epochLenth.
%    bandLimit: a cell array with specific frequency bands. Default is
%        [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21.0;21,30;30,40;40,60;1,60;];
%    bandName: a cell array with names of frequency bands. bandName
%        must be corresponding to bandLimit. Names of frequency
%        bands are required and included in {'delta';'theta';'alpha1';
%        'alpha2';'beta1';'beta2';'beta3'} for r1-r6 and PAF,
%        or included in {'delta';'theta';'alpha','beta'} for r3-r6 and PAF.
%        For relative power indices, 'fullband' must be included in
%        bandName({'fullband'}).
%    selechanns: number with indices of the selected channels
%                   (e.g. [1:4,7:30] or 'all').Default is 'all';
%    proportion: overlapped percentage for each segments/sliding windows. It
%           should be [0,1). Default is 0 (no overlapped).
%  Output:
%    EEG_results: results including power indices, mean indices across epoches
%        and parameters.
%         EEG_results.type: type of results;
%         EEG_results.spectrum: frequency spectrum for each epochs
%         EEG_results.Power: power acorss frequency bands;
%         EEG_results.Power_relative : relative power acorss frequency bands;
%         EEG_results.PAF : max power in alpha(alpha1+alpha2) band;
%         EEG_results.R1 : r1;
%         EEG_results.R2 : r2;
%         EEG_results.R3 : r3;
%         EEG_results.R4 : r4;
%         EEG_results.R5 : r5;
%         EEG_results.R6 : r6;
%         EEG_results.Power_mean : mean power acorss epoches;
%         EEG_results.Power_relative_mean : mean relative power acorss epoches;
%         EEG_results.R1_mean : mean r1 across epochs;
%         EEG_results.R2_mean : mean r2 across epochs;
%         EEG_results.R3_mean : mean r3 across epochs;
%         EEG_results.R4_mean : mean r4 across epochs;
%         EEG_results.R5_mean : mean r5 across epochs;
%         EEG_results.R6_mean : mean r6 across epochs;
%         EEG_results.PAF_mean : mean PAF across epochs;
%         EEG_results.spectrum: frequency spectrum across each epochs
%         EEG_results.spectrum_mean: mean frequency spectrum across each epochs

%         EEG_results.parameter.bandLimit: An array with specific frequency bands;
%         EEG_results.parameter.bandName: A cell array with band names;
%         EEG_results.parameter.eventlabel: An eventlabel which means good quality data;
%         EEG_results.parameter.selechanns: An array with selected channels;
%         EEG_results.parameter.epochLenth: Length of small epochs. Unit is time point.
%         EEG_results.parameter.srate: Sampling rate of EEG data.
%         EEG_results.parameter.proportion: overlapped percentage for each segments/sliding windows.;
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Li_dong729@163.com)
% $ 2018.5.25
% -------------------------------------------------------------------------
% exmaple:
% epochLenth = 5;
% eventlabel = '2001';
% bandLimit = [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21;21,30;30,40;40,60;1,60];
% bandName = {'delta';'theta';'alpha1';'alpha2';'beta1';'beta2';...
%     'beta3';'gamma1';'gamma2';'fullband'};% different frequency band's names
% selechanns = [1:10,12:15];
% EEG_results= wb_caculate_EEGpower(EEG);
% EEG_results= wb_caculate_EEGpower(EEG,[],eventlabel,[],[],selechanns);
% EEG_results= wb_caculate_EEGpower(EEG,epochLenth,eventlabel,bandLimit,bandName,selechanns)
% -------------------------------------------------------------------------
if nargin < 1
    error ('One input is reqiured at least!!!!!');
elseif nargin == 1
    epochLenth = 5;
    eventlabel = [];
    bandLimit = [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21;21,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha1';'alpha2';'beta1';'beta2';...
        'beta3';'gamma1';'gamma2';'fullband'};% different frequency band's names
    selechanns = 'all';
    proportion = 0;
elseif nargin == 2
    if isempty(epochLenth)
        epochLenth = 5;
    end
    eventlabel = [];
    bandLimit = [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21;21,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha1';'alpha2';'beta1';'beta2';...
        'beta3';'gamma1';'gamma2';'fullband'};% different frequency band's names
    selechanns = 'all';
    proportion = 0;
elseif nargin == 3
    
    if isempty(epochLenth)
        epochLenth = 5;
    end
    
    bandLimit = [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21;21,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha1';'alpha2';'beta1';'beta2';...
        'beta3';'gamma1';'gamma2';'fullband'};%different freqband's names
    selechanns = 'all';
    proportion = 0;
elseif nargin == 4
    if isempty(epochLenth)
        epochLenth = 5;
    end
    
    if isempty(bandLimit)
        bandLimit = [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21;21,30;30,40;40,60;1,60];
    end
    bandName = {'delta';'theta';'alpha1';'alpha2';'beta1';'beta2';...
        'beta3';'gamma1';'gamma2';'fullband'};% different frequency band's names
    selechanns = 'all';
    proportion = 0;
elseif nargin == 5
    if isempty(epochLenth)
        epochLenth = 5;
    end
    
    if isempty(bandLimit)
        bandLimit = [1,4;4,8;8,10.5;10.5,12.5;12.5,18.5;18.5,21;21,30;30,40;40,60;1,60];
    end
    
    if isempty(bandName)
        bandName = {'delta';'theta';'alpha1';'alpha2';'beta1';'beta2';...
            'beta3';'gamma1';'gamma2';'fullband'};% different frequency band's names
    end
    selechanns = 'all';
    proportion = 0;
elseif nargin == 6
    proportion = 0;
end
% ----------------
% check inputs
try
    srate = EEG.srate;
    if isfinite(srate) && length(size(EEG.data)) == 2
        disp(['sampling rate = ',num2str(srate)])
        if epochLenth >= 0
            disp(['epoch length = ',num2str(epochLenth),' s']);
            epochLenth1 = epochLenth * srate;
            sign_epochLenth = 1;
        else
            disp('epochLenth is negative, set to positive.');
            disp('If possible, data before event labels (eventlabel) will be used (no overlapped)');
            disp(['epoch length = ',num2str(abs(epochLenth)),' s']);
            epochLenth1 = abs(epochLenth * srate); 
            sign_epochLenth = 0;
        end
    end
catch
    disp('sampling rate is not found in EEG.');
    error('sampling rate is not found in EEG.');
end

% No. of timepoints
Nt = size(EEG.data,2);
disp(['EEG data/epoch legnth:',num2str(Nt)]);

% channs
if isequal(selechanns,'all')
    selechanns = 1:size(EEG.data,1);
end
N_channs = length(selechanns); % No. of selected channs
disp(['no. of selected channels: ',num2str(N_channs)]);

% frequency bands
N_freq = size(bandLimit,1); % number of frequency bands
N_freq_name = length(bandName); % number of frequency band names
if N_freq == 0
    disp('number of frequency bands is 0');
    error('number of frequency bands is 0');
end

if  N_freq ~= N_freq_name
    warning('N_freq is not euqal to N_freq_name, frequency-band are renamed as band#. Only calculate power in these bands.');
    bandName = [];
    for i = 1:N_freq
        bandName{i} = ['band-',num2str(i)];
    end
end
disp('frequency bands: ');
for i = 1:N_freq
    disp([bandName{i},': ',num2str(bandLimit(i,1)),'-',num2str(bandLimit(i,2)),' Hz']);
end

% proportion
try
    if proportion < 0 || proportion >= 1 || ~isfinite(proportion) || (sign_epochLenth == 0 && ~isempty(eventlabel))
        proportion = 0;
        disp('Input proportion is invalid, set proportion = 0');
    end
catch
    proportion = 0;
    disp('Input proportion may be invalid, set proportion = 0');
end


disp(['overlapped proportion for each segments/sliding windows: ',num2str(proportion)]);

% ---------------
val_1 = cellfun(@(x) isequal(x,'delta'), bandName);
temp_band_1 = bandLimit(val_1,:);
val_2 = cellfun(@(x) isequal(x,'theta'), bandName);
temp_band_2 = bandLimit(val_2,:);
val_3 = cellfun(@(x) isequal(x,'alpha1'), bandName);
temp_band_3 = bandLimit(val_3,:);
val_4 = cellfun(@(x) isequal(x,'alpha2'), bandName);
temp_band_4 = bandLimit(val_4,:);
val_5 = cellfun(@(x) isequal(x,'beta1'), bandName);
temp_band_5 = bandLimit(val_5,:);
val_6 = cellfun(@(x) isequal(x,'beta2'), bandName);
temp_band_6 = bandLimit(val_6,:);
val_7 = cellfun(@(x) isequal(x,'beta3'), bandName);
temp_band_7 = bandLimit(val_7,:);
val_8 = cellfun(@(x) isequal(x,'alpha'), bandName);
temp_band_8 = bandLimit(val_8,:);
val_9 = cellfun(@(x) isequal(x,'beta'), bandName);
temp_band_9 = bandLimit(val_9,:);
% -------------------------------------------------------------------------
% events
IndBadBlock = []; % Index of bad blocks
if ~isempty(eventlabel)
    % find eventlabel and bad block
    if isfield(EEG,'event')
        allevents = EEG.event;
        if ~isempty(allevents)            
            
            IndexInd = wb_findevent(eventlabel,allevents);% find eventlabel
            IndexInd_badblock = wb_findevent('9999',allevents); % find bad block
            if isempty(IndexInd_badblock)
                IndexInd_badblock = wb_findevent(9999,allevents);
            end
            
            if ~isempty(IndexInd)
                EEG_Blocks = allevents(IndexInd.index);
                flag1 = 1; % flag1 = 1: event label is found in EEG events;
                disp(['use the data labeled as: ',num2str(eventlabel)]);
            else
                flag1 = 2; % flag1 = 2: eventlabel is not found in events.
                warning(['eventlabel(', num2str(eventlabel),')is not found in events,No data is used to calculate power indices']);
            end
            
            if ~isempty(IndexInd_badblock)
                BadBlocks = allevents(IndexInd_badblock.index);
                IndBadBlock = zeros(1,Nt);
                for k1 = 1:length(BadBlocks)
                    IndBadBlock(BadBlocks(1,k1).latency:BadBlocks(1,k1).latency+BadBlocks(1,k1).duration)=1;
                end
                disp('Bad blocks (label 9999) were found in EEG events, data in bad block will not be used');
            else
                disp('No bad blocks (label 9999) were found');
            end
        else
            warning('events are empty in EEG, all data will be used');
            flag1 = 0;
        end
    else
        warning('events are not found in EEG, all data will be used');
        flag1 = 0;
    end
else
    warning('eventlabel is empty, all data will be used');
    flag1 = 0;
    % find bad block only
    if isfield(EEG,'event')
        allevents = EEG.event;
        if ~isempty(allevents)
            % ----------
            IndexInd_badblock = wb_findevent('9999',allevents); % find bad block
            if isempty(IndexInd_badblock)
                IndexInd_badblock = wb_findevent(9999,allevents);
            end

            if ~isempty(IndexInd_badblock)
                BadBlocks = allevents(IndexInd_badblock.index);
                IndBadBlock = zeros(1,Nt);
                for k1 = 1:length(BadBlocks)
                    IndBadBlock(BadBlocks(1,k1).latency:BadBlocks(1,k1).latency+BadBlocks(1,k1).duration)=1;
                end
                disp('Bad blocks (label 9999) were found in EEG events, data in bad block will not be used');
            else
                disp('No bad blocks (label 9999) were found');
            end
        end
    end
end

if length(size(EEG.data)) == 3
    warning('EEG.data has been epoched, all data will be used and bad blocks (label 9999, if exist) are not considered');
    flag1 = 0;
    IndBadBlock = [];
end
% -------------------------------------------------------------------------
% Default parameters of time-frequency analysis
WaveletCycles = 0;
WaveletMethod = 'dftfilt3';
TaperingFunction = 'hanning';
DetrendStr = 'on';
disp('---------')
disp('Settings of time-frequency analysis: ');
disp(['WaveletCycles = ',num2str(WaveletCycles)]);
disp(['WaveletMethod = ',WaveletMethod]);
disp(['TaperingFunction = ',TaperingFunction]);
disp(['DetrendStr = ',DetrendStr]);
% -------------------------------------------------------------------------
% caculate power indices
% indices
Power = [];
Power_relative=[];
spectrum = [];
R1 = [];
R2 = [];
R3 = [];
R4 = [];
R5 = [];
R6 = [];
PAF = [];
freqs1 = [];
if flag1 == 1 % flag1 = 1: eventlabel is found in EEG events
    N_Blocks = length(EEG_Blocks);
    for i1 = 1:N_Blocks
        if ~isfield(EEG_Blocks(i1),'duration')
            EEG_Blocks(i1).duration = epochLenth1; % if event doesn't include duration, set duration as epochLenth.
        else
            if isempty(EEG_Blocks(i1).duration)
                EEG_Blocks(i1).duration = epochLenth1;
            end
        end
        block_length = EEG_Blocks(i1).duration;
        
        % make sure the EEG blocks does not contain bad blocks.
        if sign_epochLenth == 1 % sign of epochLenth is positive
            if EEG_Blocks(i1).latency + block_length <= Nt && ~isempty(IndBadBlock)
                Ind_badblock = IndBadBlock(1,EEG_Blocks(i1).latency:EEG_Blocks(i1).latency + block_length);
                temp_badblock = sum(Ind_badblock(:));
            else
                temp_badblock = 0;
            end
        elseif sign_epochLenth == 0 % sign of epochLenth is negative
            if EEG_Blocks(i1).latency - epochLenth1 > 0 && ~isempty(IndBadBlock)
                Ind_badblock = IndBadBlock(1,EEG_Blocks(i1).latency - epochLenth1 : EEG_Blocks(i1).latency);
                temp_badblock = sum(Ind_badblock(:));
            else
                temp_badblock = 0;
            end
        end
        
        if sign_epochLenth == 1 % sign of epochLenth is positive
            
            if EEG_Blocks(i1).latency + block_length <= Nt && temp_badblock == 0
                % N_epoches = fix(block_length/epochLenth);
                temp_data = EEG.data(selechanns,EEG_Blocks(i1).latency:EEG_Blocks(i1).latency + block_length);
                temp_windata = wb_slidewin(temp_data,epochLenth1,proportion);
                N_epoches = length(temp_windata);
                
                if N_epoches >= 1
                    temp_power_epoch = zeros(N_channs,N_freq,N_epoches);
                    temp_power_epoch_relative=zeros(N_channs,N_freq,N_epoches);
                    temp_alpha_max = zeros(N_channs,N_epoches);
                    temp_r1 = zeros(N_channs,N_epoches);
                    temp_r2 = zeros(N_channs,N_epoches);
                    temp_r3 = zeros(N_channs,N_epoches);
                    temp_r4 = zeros(N_channs,N_epoches);
                    temp_r5 = zeros(N_channs,N_epoches);
                    temp_r6 = zeros(N_channs,N_epoches);
                    for i2 = 1:N_epoches
                        % extract epoch data
                        % t1 = EEG_Blocks(i1).latency + epochLenth * (i2-1) + 1 ;
                        % t2 = EEG_Blocks(i1).latency + epochLenth * i2;
                        % temp_data1 = EEG.data(selechanns,t1:t2);
                        
                        temp_data1 = temp_windata{1,i2};
                        
                        for i3 = 1:N_channs
                            temp_data2 = temp_data1(i3,:);
                            % time-frequency analysis
                            [tf, freqs1, ~] = timefreq(temp_data2,...
                                srate,'cycles',WaveletCycles,'wletmethod',WaveletMethod,...
                                'ffttaper',TaperingFunction, 'detrend',DetrendStr);
                            % calculate power
                            temp_power1 = 2*abs(tf).^2/length(temp_data1);
                            temp_MeanPower = mean(temp_power1,2);
                            temp_MeanPower2(i3,:,i2) = temp_MeanPower;
                            % power and relative power
                            for j = 1:N_freq
                                Hz1 = bandLimit(j,1);
                                Hz2 = bandLimit(j,2);
                                temp_power_epoch(i3,j,i2) = mean(temp_MeanPower(freqs1>Hz1 & freqs1<=Hz2,:));
                                % relative power
                                val = cellfun(@(x) isequal(x,'fullband'), bandName);
                                temp_band = bandLimit(val,:);
                                if ~isempty(temp_band)
                                    temp_power_epoch_relative(i3,j,i2) = sum(temp_MeanPower(freqs1>Hz1 & freqs1<=Hz2,:))...
                                        /sum(temp_MeanPower(freqs1>temp_band(1)&freqs1<=temp_band(2),:));
                                end
                            end
                            %    r1 = theta/(alpha1 + alpha2 + beta1);
                            %    r2 = (delta + theta)/(alpha1 + alpha2 + beta1 + beta2);
                            %    r3 = theta/alpha =  theta/(alpha1+alpha2);
                            %    r4 = theta/beta = theta/(beta1 + beta2 + beta3);
                            %    r5 = delta/theta;
                            %    r6 = alpha/beta = (alpha1 + alpha2)/(beta1 + beta2 + beta3);
                            %    PAF (peak of alpha frequency) = max power in alpha(alpha1+alpha2) band.
                            if ~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5)
                                temp_r1(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:)));
                                %                     else
                                %                         disp('r1 bands(theta,alpha1,alpha2,beta1) not found');
                            end
                            
                            if ~isempty(temp_band_1)&&~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5) && ~isempty(temp_band_6)
                                temp_r2(i3,i2) = (sum(temp_MeanPower(freqs1>temp_band_1(1)&freqs1<=temp_band_1(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:)))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:)));
                                %                     else
                                %                         disp('r2 bands(delta,theta,alpha1,alpha2,beta1,beta2) not found');
                            end
                            
                            if  ~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4)
                                temp_r3(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)));
                            elseif ~isempty(temp_band_2) && ~isempty(temp_band_8)
                                temp_r3(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    sum(temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:));
                                %                     else
                                %                         disp('r3 bands([theta,alpha1,alpha2] or [theta,alpha]) not found');
                            end
                            
                            if ~isempty(temp_band_2) && ~isempty(temp_band_5) && ~isempty(temp_band_6) && ~isempty(temp_band_7)
                                temp_r4(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_7(1)&freqs1<=temp_band_7(2),:)));
                            elseif ~isempty(temp_band_2) && ~isempty(temp_band_9)
                                temp_r4(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_9(1)&freqs1<=temp_band_9(2),:)));
                                %                     else
                                %                         disp('r4 bands([theta,beta1,beta2,beta3] or [theta,beta]) not found');
                            end
                            
                            if ~isempty(temp_band_1) && ~isempty(temp_band_2)
                                temp_r5(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_1(1)&freqs1<=temp_band_1(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:)));
                                %                     else
                                %                         disp('r5 bands(delta,theta) not found');
                            end
                            
                            if ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5) && ~isempty(temp_band_6) && ~isempty(temp_band_7)
                                temp_r6(i3,i2) = (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_7(1)&freqs1<=temp_band_7(2),:)));
                            elseif ~isempty(temp_band_8) && ~isempty(temp_band_9)
                                temp_r6(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_9(1)&freqs1<=temp_band_9(2),:)));
                                %                     else
                                %                         disp('r6 bands([alpha1,alpha2,beta1,beta2,beta3] or [alpha,beta]) not found');
                            end
                            
                            if ~isempty(temp_band_3) && ~isempty(temp_band_4)
                                temp_alpha_max(i3,i2) = max([0;temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:);...
                                    temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)]);
                            elseif ~isempty(temp_band_8)
                                temp_alpha_max(i3,i2) = max([0;temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:)]);
                                %                     else
                                %                         disp('PAF bands([alpha1,alpha2] or [alpha]) not found');
                            end
                        end
                    end
                    
                    Power = cat(3,Power,temp_power_epoch); % concatenate arrays
                    Power_relative = cat(3,Power_relative,temp_power_epoch_relative);
                    spectrum = cat(3,spectrum,temp_MeanPower2);
                    
                    
                    if any(temp_alpha_max)
                        PAF = [PAF,temp_alpha_max];
                    end;
                    if any(temp_r1)
                        R1 = [R1,temp_r1];
                    end
                    if any(temp_r2)
                        R2 = [R2,temp_r2];
                    end
                    if any(temp_r3)
                        R3 = [R3,temp_r3];
                    end
                    if any(temp_r4)
                        R4 = [R4,temp_r4];
                    end;
                    if any(temp_r5)
                        R5 = [R5,temp_r5];
                    end
                    if any(temp_r6)
                        R6 = [R6,temp_r6];
                    end
                end
            end
        elseif sign_epochLenth == 0 % sign of epochLenth is negative
            if EEG_Blocks(i1).latency - epochLenth1 > 0 && temp_badblock == 0
                % N_epoches = fix(block_length/epochLenth);
                temp_windata{1,1} = EEG.data(selechanns,EEG_Blocks(i1).latency - epochLenth1 : EEG_Blocks(i1).latency);
                % temp_windata = wb_slidewin(temp_data,epochLenth1,proportion);
                N_epoches = 1;
                
                if N_epoches >= 1
                    temp_power_epoch = zeros(N_channs,N_freq,N_epoches);
                    temp_power_epoch_relative=zeros(N_channs,N_freq,N_epoches);
                    temp_alpha_max = zeros(N_channs,N_epoches);
                    temp_r1 = zeros(N_channs,N_epoches);
                    temp_r2 = zeros(N_channs,N_epoches);
                    temp_r3 = zeros(N_channs,N_epoches);
                    temp_r4 = zeros(N_channs,N_epoches);
                    temp_r5 = zeros(N_channs,N_epoches);
                    temp_r6 = zeros(N_channs,N_epoches);
                    for i2 = 1:N_epoches
                        % extract epoch data
                        % t1 = EEG_Blocks(i1).latency + epochLenth * (i2-1) + 1 ;
                        % t2 = EEG_Blocks(i1).latency + epochLenth * i2;
                        % temp_data1 = EEG.data(selechanns,t1:t2);
                        
                        temp_data1 = temp_windata{1,i2};
                        
                        for i3 = 1:N_channs
                            temp_data2 = temp_data1(i3,:);
                            % time-frequency analysis
                            [tf, freqs1, ~] = timefreq(temp_data2,...
                                srate,'cycles',WaveletCycles,'wletmethod',WaveletMethod,...
                                'ffttaper',TaperingFunction, 'detrend',DetrendStr);
                            % calculate power
                            temp_power1 = 2*abs(tf).^2/length(temp_data1);
                            temp_MeanPower = mean(temp_power1,2);
                            temp_MeanPower2(i3,:,i2) = temp_MeanPower;
                            % power and relative power
                            for j = 1:N_freq
                                Hz1 = bandLimit(j,1);
                                Hz2 = bandLimit(j,2);
                                temp_power_epoch(i3,j,i2) = mean(temp_MeanPower(freqs1>Hz1 & freqs1<=Hz2,:));
                                % relative power
                                val = cellfun(@(x) isequal(x,'fullband'), bandName);
                                temp_band = bandLimit(val,:);
                                if ~isempty(temp_band)
                                    temp_power_epoch_relative(i3,j,i2) = sum(temp_MeanPower(freqs1>Hz1 & freqs1<=Hz2,:))...
                                        /sum(temp_MeanPower(freqs1>temp_band(1)&freqs1<=temp_band(2),:));
                                end
                            end
                            %    r1 = theta/(alpha1 + alpha2 + beta1);
                            %    r2 = (delta + theta)/(alpha1 + alpha2 + beta1 + beta2);
                            %    r3 = theta/alpha =  theta/(alpha1+alpha2);
                            %    r4 = theta/beta = theta/(beta1 + beta2 + beta3);
                            %    r5 = delta/theta;
                            %    r6 = alpha/beta = (alpha1 + alpha2)/(beta1 + beta2 + beta3);
                            %    PAF (peak of alpha frequency) = max power in alpha(alpha1+alpha2) band.
                            if ~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5)
                                temp_r1(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:)));
                                %                     else
                                %                         disp('r1 bands(theta,alpha1,alpha2,beta1) not found');
                            end
                            
                            if ~isempty(temp_band_1)&&~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5) && ~isempty(temp_band_6)
                                temp_r2(i3,i2) = (sum(temp_MeanPower(freqs1>temp_band_1(1)&freqs1<=temp_band_1(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:)))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:)));
                                %                     else
                                %                         disp('r2 bands(delta,theta,alpha1,alpha2,beta1,beta2) not found');
                            end
                            
                            if  ~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4)
                                temp_r3(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)));
                            elseif ~isempty(temp_band_2) && ~isempty(temp_band_8)
                                temp_r3(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    sum(temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:));
                                %                     else
                                %                         disp('r3 bands([theta,alpha1,alpha2] or [theta,alpha]) not found');
                            end
                            
                            if ~isempty(temp_band_2) && ~isempty(temp_band_5) && ~isempty(temp_band_6) && ~isempty(temp_band_7)
                                temp_r4(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_7(1)&freqs1<=temp_band_7(2),:)));
                            elseif ~isempty(temp_band_2) && ~isempty(temp_band_9)
                                temp_r4(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_9(1)&freqs1<=temp_band_9(2),:)));
                                %                     else
                                %                         disp('r4 bands([theta,beta1,beta2,beta3] or [theta,beta]) not found');
                            end
                            
                            if ~isempty(temp_band_1) && ~isempty(temp_band_2)
                                temp_r5(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_1(1)&freqs1<=temp_band_1(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:)));
                                %                     else
                                %                         disp('r5 bands(delta,theta) not found');
                            end
                            
                            if ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5) && ~isempty(temp_band_6) && ~isempty(temp_band_7)
                                temp_r6(i3,i2) = (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:))+...
                                    sum(temp_MeanPower(freqs1>temp_band_7(1)&freqs1<=temp_band_7(2),:)));
                            elseif ~isempty(temp_band_8) && ~isempty(temp_band_9)
                                temp_r6(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:))/...
                                    (sum(temp_MeanPower(freqs1>temp_band_9(1)&freqs1<=temp_band_9(2),:)));
                                %                     else
                                %                         disp('r6 bands([alpha1,alpha2,beta1,beta2,beta3] or [alpha,beta]) not found');
                            end
                            
                            if ~isempty(temp_band_3) && ~isempty(temp_band_4)
                                temp_alpha_max(i3,i2) = max([0;temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:);...
                                    temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)]);
                            elseif ~isempty(temp_band_8)
                                temp_alpha_max(i3,i2) = max([0;temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:)]);
                                %                     else
                                %                         disp('PAF bands([alpha1,alpha2] or [alpha]) not found');
                            end
                        end
                    end
                    
                    Power = cat(3,Power,temp_power_epoch); % concatenate arrays
                    Power_relative = cat(3,Power_relative,temp_power_epoch_relative);
                    spectrum = cat(3,spectrum,temp_MeanPower2);
                    
                    
                    if any(temp_alpha_max)
                        PAF = [PAF,temp_alpha_max];
                    end;
                    if any(temp_r1)
                        R1 = [R1,temp_r1];
                    end
                    if any(temp_r2)
                        R2 = [R2,temp_r2];
                    end
                    if any(temp_r3)
                        R3 = [R3,temp_r3];
                    end
                    if any(temp_r4)
                        R4 = [R4,temp_r4];
                    end;
                    if any(temp_r5)
                        R5 = [R5,temp_r5];
                    end
                    if any(temp_r6)
                        R6 = [R6,temp_r6];
                    end
                end
            end
        end
    end
elseif flag1 == 0  % flag1 = 0: events are not contained in EEG or eventlabel is empty,all data will be used.
    
    if length(size(EEG.data)) == 2
        % block_length = size(EEG.data,2);
        % N_epoches = fix(block_length/epochLenth);
        temp_data = EEG.data(selechanns,:);
        temp_windata = wb_slidewin(temp_data,epochLenth1,proportion);
        
        % make sure the EEG epoches does not contain bad blocks.
        if ~isempty(IndBadBlock)
            temp_winbadblock = wb_slidewin(IndBadBlock,epochLenth1,proportion);
            temp_badblock3 = zeros(length(temp_winbadblock),1);
            for k2 = 1:length(temp_winbadblock)
                temp_badblock2 = temp_winbadblock{1,k2};
                temp_badblock3(k2) = sum(temp_badblock2(:));
            end
            index1 = find(temp_badblock3==0);
            N_epoches = length(index1);
        else
            N_epoches = length(temp_windata);
            index1 = 1:N_epoches;
        end
    elseif length(size(EEG.data)) == 3
        N_epoches = size(EEG.data,3);
    end
    
    if N_epoches >= 1
        temp_power_epoch = zeros(N_channs,N_freq,N_epoches);
        temp_power_epoch_relative=zeros(N_channs,N_freq,N_epoches);
        temp_alpha_max = zeros(N_channs,N_epoches);
        temp_r1 = zeros(N_channs,N_epoches);
        temp_r2 = zeros(N_channs,N_epoches);
        temp_r3 = zeros(N_channs,N_epoches);
        temp_r4 = zeros(N_channs,N_epoches);
        temp_r5 = zeros(N_channs,N_epoches);
        temp_r6 = zeros(N_channs,N_epoches);
        for i2 = 1:N_epoches
            % extract epoch data
            if length(size(EEG.data)) == 2
                % t1 = epochLenth * (i2-1)+1;
                % t2 = epochLenth * i2;
                % temp_data1 = EEG.data(selechanns,t1:min(t2,Nt));
                
                temp_data1 = temp_windata{1,index1(i2)};
                
            elseif length(size(EEG.data)) == 3
                temp_data1 = EEG.data(selechanns,:,i2);
            end
            if ~isempty(temp_data1)
                % -----
                for i3 = 1:N_channs
                    temp_data2 = temp_data1(i3,:);
                    % time-frequency analysis
                    [tf, freqs1, ~] = timefreq(temp_data2,...
                        srate,'cycles',WaveletCycles,'wletmethod',WaveletMethod,...
                        'ffttaper',TaperingFunction, 'detrend',DetrendStr);
                    % calculate power
                    temp_power1 = 2*abs(tf).^2/length(temp_data1);
                    temp_MeanPower = mean(temp_power1,2);
                    temp_MeanPower2(i3,:,i2) = temp_MeanPower;
                    % power and relative power
                    for j = 1:N_freq
                        Hz1 = bandLimit(j,1);
                        Hz2 = bandLimit(j,2);
                        temp_power_epoch(i3,j,i2) = mean(temp_MeanPower(freqs1>Hz1 & freqs1<=Hz2,:));
                        % relative power
                        val = cellfun(@(x) isequal(x,'fullband'), bandName);
                        temp_band = bandLimit(val,:);
                        if ~isempty(temp_band)
                            temp_power_epoch_relative(i3,j,i2) = sum(temp_MeanPower(freqs1>Hz1 & freqs1<=Hz2,:))...
                                /sum(temp_MeanPower(freqs1>temp_band(1)&freqs1<=temp_band(2),:));
                        end
                    end
                    %    r1 = theta/(alpha1 + alpha2 + beta1);
                    %    r2 = (delta + theta)/(alpha1 + alpha2 + beta1 + beta2);
                    %    r3 = theta/alpha =  theta/(alpha1+alpha2);
                    %    r4 = theta/beta = theta/(beta1 + beta2 + beta3);
                    %    r5 = delta/theta;
                    %    r6 = alpha/beta = (alpha1 + alpha2)/(beta1 + beta2 + beta3);
                    %    PAF (peak of alpha frequency) = max power in alpha(alpha1+alpha2) band.
                    if ~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5)
                        temp_r1(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                            (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:)));
                        %                 else
                        %                     disp('r1 bands(theta,alpha1,alpha2,beta1) not found');
                    end
                    
                    if ~isempty(temp_band_1)&&~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5) && ~isempty(temp_band_6)
                        temp_r2(i3,i2) = (sum(temp_MeanPower(freqs1>temp_band_1(1)&freqs1<=temp_band_1(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:)))/...
                            (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:)));
                        %                 else
                        %                     disp('r2 bands(delta,theta,alpha1,alpha2,beta1,beta2) not found');
                    end
                    
                    if  ~isempty(temp_band_2) && ~isempty(temp_band_3) && ~isempty(temp_band_4)
                        temp_r3(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                            (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)));
                    elseif ~isempty(temp_band_2) && ~isempty(temp_band_8)
                        temp_r3(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                            sum(temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:));
                        %                 else
                        %                     disp('r3 bands([theta,alpha1,alpha2] or [theta,alpha]) not found');
                    end
                    
                    if ~isempty(temp_band_2) && ~isempty(temp_band_5) && ~isempty(temp_band_6) && ~isempty(temp_band_7)
                        temp_r4(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                            (sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_7(1)&freqs1<=temp_band_7(2),:)));
                    elseif ~isempty(temp_band_2) && ~isempty(temp_band_9)
                        temp_r4(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:))/...
                            (sum(temp_MeanPower(freqs1>temp_band_9(1)&freqs1<=temp_band_9(2),:)));
                        %                 else
                        %                     disp('r4 bands([theta,beta1,beta2,beta3] or [theta,beta]) not found');
                    end
                    
                    if ~isempty(temp_band_1) && ~isempty(temp_band_2)
                        temp_r5(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_1(1)&freqs1<=temp_band_1(2),:))/...
                            (sum(temp_MeanPower(freqs1>temp_band_2(1)&freqs1<=temp_band_2(2),:)));
                        %                 else
                        %                     disp('r5 bands(delta,theta) not found');
                    end
                    
                    if ~isempty(temp_band_3) && ~isempty(temp_band_4) && ~isempty(temp_band_5) && ~isempty(temp_band_6) && ~isempty(temp_band_7)
                        temp_r6(i3,i2) = (sum(temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)))/...
                            (sum(temp_MeanPower(freqs1>temp_band_5(1)&freqs1<=temp_band_5(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_6(1)&freqs1<=temp_band_6(2),:))+...
                            sum(temp_MeanPower(freqs1>temp_band_7(1)&freqs1<=temp_band_7(2),:)));
                    elseif ~isempty(temp_band_8) && ~isempty(temp_band_9)
                        temp_r6(i3,i2) = sum(temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:))/...
                            (sum(temp_MeanPower(freqs1>temp_band_9(1)&freqs1<=temp_band_9(2),:)));
                        %                 else
                        %                     disp('r6 bands([alpha1,alpha2,beta1,beta2,beta3] or [alpha,beta]) not found');
                    end
                    
                    if ~isempty(temp_band_3) && ~isempty(temp_band_4)
                        temp_alpha_max(i3,i2) = max([0;temp_MeanPower(freqs1>temp_band_3(1)&freqs1<=temp_band_3(2),:);...
                            temp_MeanPower(freqs1>temp_band_4(1)&freqs1<=temp_band_4(2),:)]);
                    elseif ~isempty(temp_band_8)
                        temp_alpha_max(i3,i2) = max([0;temp_MeanPower(freqs1>temp_band_8(1)&freqs1<=temp_band_8(2),:)]);
                        %                 else
                        %                     disp('PAF bands([alpha1,alpha2] or [alpha]) not found');
                    end
                end
            end
        end
        
        Power = cat(3,Power,temp_power_epoch); % concatenate arrays
        Power_relative = cat(3,Power_relative,temp_power_epoch_relative);
        spectrum = cat(3,spectrum,temp_MeanPower2);
        
        if any(temp_alpha_max)
            PAF = [PAF,temp_alpha_max];
        end;
        if any(temp_r1)
            R1 = [R1,temp_r1];
        end
        if any(temp_r2)
            R2 = [R2,temp_r2];
        end
        if any(temp_r3)
            R3 = [R3,temp_r3];
        end
        if any(temp_r4)
            R4 = [R4,temp_r4];
        end;
        if any(temp_r5)
            R5 = [R5,temp_r5];
        end
        if any(temp_r6)
            R6 = [R6,temp_r6];
        end
    end
end

% save results
EEG_results = [];
EEG_results.type = 'power';
try EEG_results.filename = EEG.filename;catch;end;
EEG_results.Power = Power;
EEG_results.spectrum = spectrum;
EEG_results.freqs = freqs1;
EEG_results.Power_relative = Power_relative;
EEG_results.PAF = PAF;
EEG_results.R1 = R1;
EEG_results.R2 = R2;
EEG_results.R3 = R3;
EEG_results.R4 = R4;
EEG_results.R5 = R5;
EEG_results.R6 = R6;
if flag1 == 2
    EEG_results.Block_percentage = 0;
elseif length(size(EEG.data)) == 3
    EEG_results.Block_percentage = 1;
else
    Ngd = size(Power,3);
    EEG_results.Block_percentage = (Ngd * epochLenth1)/size(EEG.data,2);
end
% save parameters
EEG_results.parameter.WaveletCycles=0;
EEG_results.parameter.WaveletMethod='dftfilt3';
EEG_results.parameter.TaperingFunction='hanning';
EEG_results.parameter.DetrendStr = 'on';
EEG_results.parameter.bandLimit = bandLimit;
EEG_results.parameter.bandName = bandName;
EEG_results.parameter.eventlabel = eventlabel;
EEG_results.parameter.selechanns = selechanns;
EEG_results.parameter.epochLenth = epochLenth;
EEG_results.parameter.srate = EEG.srate;
EEG_results.parameter.proportion = proportion;
EEG_results.parameter.chanlocs = EEG.chanlocs;
EEG_results.parameter.ref = EEG.ref;

% calulate mean value acorss epoches
EEG_results.Power_mean = nanmean(Power,3);
EEG_results.spectrum_mean = nanmean(spectrum,3);
EEG_results.Power_relative_mean = nanmean(Power_relative,3);
EEG_results.R1_mean = mean(R1,2);
EEG_results.R2_mean = mean(R2,2);
EEG_results.R3_mean = mean(R3,2);
EEG_results.R4_mean = mean(R4,2);
EEG_results.R5_mean = mean(R5,2);
EEG_results.R6_mean = mean(R6,2);
EEG_results.PAF_mean = mean(PAF,2);