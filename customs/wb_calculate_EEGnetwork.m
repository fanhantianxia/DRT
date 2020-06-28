function EEG_results = wb_calculate_EEGnetwork(EEG,epochLenth,eventlabel,bandLimit,selechanns,method,proportion)
% Calculate EEG network between EEG channels.
% Calculating EEG netowork consists of
% [1] Specific event data can be extracted according to the input 
%    ¡®eventlabel¡¯. If the input ¡®eventlabel¡¯ is empty, all data will be 
%     used. If applicable, EEG segments in bad block (label 9999, marked by 
%     wb_pipeline_EEG_Mark) will be rejected automatically, and NOT used to 
%     calculate network.
% [2] Specific EEG signals will be divided into small epochs.
% [3] EEG data of each epoch (default is 5-s epoch) was subjected to 
%     calculate correlation/coherence/PSI/PLI to obtain the EEG network 
%     across electrodes in the specific bands.


% Noting that data in bad blocks (labeled as 9999) will not be used to
% calculate power indices !!!

% EEG networks including:
% Corr = C(i,j)/sqrt(C(i,i)C(j,j));C is the covariance matrix.
% Coherence = |Pxy(f)|.^2/(Pxx(f)Pyy(f));  Magnitude squared coherence.
% Phase Synchronization Index (PSI) = sqrt(<cos(xy_theta)>.^2 + <sin(xy_theta)>.^2);
% Phase Locking Index (PLI) = abs(sum(exp(i*xy_theta)))/length of data;
%
% Default frequency bands:
%    delta: 1 - 4 Hz
%    theta: 4 - 8 Hz
%    alpha: 8 - 12.5 Hz
%    beta: 12.5 - 25 Hz
%    high beta: 25 - 30 Hz
%    gamma1: mean power acorss 30 - 40 Hz
%    gamma2: mean power acorss 40 - 60 Hz
%    fullband: 1 - 60 Hz

% References:
% Frequency bands:
%   Nuwer, M. R., et al. (1994). "IFCN Guidelines for Topographic and Frequency-Analysis of EEGs and EPs - Report of an IFCN Committee." Electroencephalography and Clinical Neurophysiology 91(1): 1-5.
%   Thatcher, R. W., et al. (2005). "EEG and intelligence: Relations between EEG coherence, EEG phase delay and power." Clinical Neurophysiology 116(9): 2129-2141.
%   Chen, A. C., et al. (2008). "EEG default mode network in the human brain: spectral regional field powers." Neuroimage 41(2): 561-574.
%   Malver, L. P., et al. (2014). "Electroencephalography and analgesics." Br J Clin Pharmacol 77(1): 72-95.
% Networks:
%   Bob, P., et al. (2008). "EEG phase synchronization in patients with paranoid schizophrenia." Neurosci Lett 447(1): 73-77.
%   Lee, Y. Y. and S. Hsieh (2014). "Classifying different emotional states by means of EEG-based functional connectivity patterns." PLoS One 9(4): e95415.
%   Xu, P., et al. (2014). "Differentiating between psychogenic nonepileptic seizures and epilepsy based on common spatial pattern of weighted EEG resting networks." IEEE Trans Biomed Eng 61(6): 1747-1755.
%   Edagawa, K. and M. Kawasaki (2017). "Beta phase synchronization in the frontal-temporal-cerebellar network during auditory-to-motor rhythm learning." Scientific Reports 7.
% -------------------------------------------------------------------------
% Input:
%    EEG: EEG structure imported using EEGLAB. EEG.data should be channels
%         X time points OR channels X time points X epoches.
%    epochLenth: length of small epochs to calculate EEG network.
%         unit is second.Default is 5s.If epochLenth is negative, it means
%         that if possible, data before event labels (eventlabel) will be used
%         (no overlapped).
%    eventlabel: Event label which means specific event data. Default is empty. 
%        If it is empty, all data will be used. If eventlabel is not found, 
%        NO data will be epoched and calculated. If structure event (eventlabel)
%        doesn't include duration, the duration will be equal to epochLenth.
%    bandLimit: A cell array with specific frequency bands. Default is
%        [1,4;4,8;8,12.5;12.5,25;25,30;30,40;40,60;1,60];
%    selechanns: number with indices of the selected channels
%                   (e.g. [1:4,7:30] or 'all').Default is 'all';
%    method: Method used to calcualte EEG network. Default is 'psi'.
%            'corr': Pearson's correlation
%            'cohere': Magnitude squared coherence
%            'psi':   Phase Synchronization Index
%            'plv':   Phase Locking Value
%   proportion: overlapped percentage for each segments/sliding windows. It
%           should be [0,1). Default is 0 (no overlapped).
%  Output:
%    EEG_results: results including connection matrices, mean connection 
%        matrices across epoches, z-score connection matrices and parameters.

%        EEG_results.type: type of results, i.e. ¡®network¡¯.
%        EEG_results.nettype: type of network.
%             ¡®BU¡¯: binary undirected network;
%             ¡®BD¡¯: binary directed network;
%             ¡®WU¡¯: weighted undirected network;
%             ¡®WD¡¯: weighted directed network.
%         EEG_results.M: connection matrix (symmetric) across frequency bands;
%               Size of M is channels X channels X epoches X frequency bands.
%         EEG_results.M_zscore : Fisher's z-score connection matrix across 
%               frequency bands;Size of M_zscore is channels X channels X epoches X frequency bands.
%         EEG_results.M_mean : mean connection matrix across epoches;
%               Size of M_mean is channels X channels X frequency bands
%         EEG_results.M_zscore_mean : mean Fisher's z-score connection matrix across epoches;
%               Size of M_zscore_mean is channels X channels X frequency bands
%         EEG_results.Block_percentage: percentage of EEG data used to calculate network.
%         EEG_results.filename: filename of EEG data.

%         EEG_results.parameter.bandLimit: An array with specific frequency bands;
%         EEG_results.parameter.bandName: A cell array with band names;
%         EEG_results.parameter.eventlabel: An eventlabel which means good quality data;
%         EEG_results.parameter.selechanns: An array with selected channels;
%         EEG_results.parameter.epochLenth: Length of small epochs. Unit is time point.
%         EEG_results.parameter.srate: Sampling rate of EEG data.
%         EEG_results.parameter.method: Method used to calcualte EEG network.
%         EEG_results.parameter.chanlocs = EEG.chanlocs;
%         EEG_results.parameter.ref = EEG.ref;
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Li_dong729@163.com)
% $ 2018.5.29
% -------------------------------------------------------------------------
% exmaple:
% epochLenth = 5;
% eventlabel = [];
% bandLimit = [1,4;4,8;8,12.5;12.5,25;25,30;30,40;40,60;1,60];
% bandLimit = [1,4;4,8];
% bandName = {'delta';'theta';'alpha';'beta';'high beta';'gamma1';'gamma2';'fullband'};
% selechanns = [1:10];
% method = 'psi';
% EEG_network = wb_calculate_EEGnetwork(EEG);
% EEG_network = wb_calculate_EEGnetwork(EEG,[],eventlabel,[],[],selechanns);
% EEG_network = wb_calculate_EEGnetwork(EEG,epochLenth,eventlabel,bandLimit,selechanns,method);
% -------------------------------------------------------------------------
if nargin < 1
    error ('One input is reqiured at least!!!!!');
elseif nargin == 1
    epochLenth = 5;
    eventlabel = [];
    bandLimit = [1,4;4,8;8,12.5;12.5,25;25,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha';'beta';'high beta';'gamma1';'gamma2';'fullband'};% different frequency band's names
    selechanns = 'all';
    method = 'psi';
    proportion = 0;
elseif nargin == 2
    eventlabel = [];
    bandLimit = [1,4;4,8;8,12.5;12.5,25;25,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha';'beta';'high beta';'gamma1';'gamma2';'fullband'};% different frequency band's names
    selechanns = 'all';
    method = 'psi';
    proportion = 0;
elseif nargin == 3
    bandLimit = [1,4;4,8;8,12.5;12.5,25;25,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha';'beta';'high beta';'gamma1';'gamma2';'fullband'};% different frequency band's names
    selechanns = 'all';
    method = 'psi';
    proportion = 0;
elseif nargin == 4
    if ~isempty(bandLimit)
        N_freq = size(bandLimit,1); % number of frequency bands
        for j = 1:N_freq
            bandName{j,1} = ['band-',num2str(j)];
        end
    end
    selechanns = 'all';
    method = 'psi';
    proportion = 0;
elseif nargin == 5
    if ~isempty(bandLimit)
        N_freq = size(bandLimit,1); % number of frequency bands
        for j = 1:N_freq
            bandName{j,1} = ['band-',num2str(j)];
        end
    end
    method = 'psi';
    proportion = 0;
elseif nargin == 6
    if ~isempty(bandLimit)
        N_freq = size(bandLimit,1); % number of frequency bands
        for j = 1:N_freq
            bandName{j,1} = ['band-',num2str(j)];
        end
    end
    proportion = 0;
end
% ----------------
% check inputs: 

% epochLenth
if isempty(epochLenth)
    epochLenth = 5;
end
% sampling rate
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
disp(['No. of selected channels: ',num2str(N_channs)]);

disp(['connection method: ',method]);
% frequency bands
if isempty(bandLimit) || size(bandLimit,2) == 1
    bandLimit = [1,4;4,8;8,12.5;12.5,25;25,30;30,40;40,60;1,60];
    bandName = {'delta';'theta';'alpha';'beta';'high beta';'gamma1';'gamma2';'fullband'};% different frequency band's names
end
N_freq = size(bandLimit,1); % number of frequency bands

% proportion
if proportion < 0 || proportion >= 1 || ~isfinite(proportion) || (sign_epochLenth == 0 && ~isempty(eventlabel))
    proportion = 0;
    disp('Input proportion is invalid, set proportion = 0');
end
disp(['overlapped proportion for each segments/sliding windows: ',num2str(proportion)]);
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
disp('----------');
if flag1 == 2
    disp('calcualtion skipped ...');
else
    disp('calculating connection matries ......');
end
% connection matrix
M = []; 
if flag1 == 1 % flag1 = 1: eventlabel is found in EEG events
    for j = 1:N_freq
        disp(['calculating frequency band: ',num2str(bandLimit(j,1)),'-',num2str(bandLimit(j,2)),' Hz']);
        EEG_temp = pop_eegfiltnew(EEG,bandLimit(j,1),bandLimit(j,2),[],0); % filtering
        N_Blocks = length(EEG_Blocks);
        M_temp1 = [];
        for i1 = 1:N_Blocks
            if ~isfield(EEG_Blocks(i1),'duration') % if event doesn't include duration, set duration as epochLenth.
                EEG_Blocks(i1).duration = epochLenth1; 
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
                    
                    temp_data = EEG_temp.data(selechanns,EEG_Blocks(i1).latency:EEG_Blocks(i1).latency + block_length);
                    temp_windata = wb_slidewin(temp_data,epochLenth1,proportion);
                    N_epoches = length(temp_windata);
                    if N_epoches >= 1
                        M_temp = zeros(N_channs,N_channs,N_epoches);
                        for i2 = 1:N_epoches
                            % extract epoch data
                            % t1 = EEG_Blocks(i1).latency + epochLenth * (i2-1) + 1 ;
                            % t2 = EEG_Blocks(i1).latency + epochLenth * i2;
                            % temp_data1 = EEG_temp.data(selechanns,t1:t2);
                            temp_data1 = temp_windata{1,i2};
                            % -----
                            switch method
                                case 'corr'
                                    % correlation
                                    temp1 = corr(temp_data1');
                                    M_temp(:,:,i2) = triu(temp1,1) + triu(temp1,1)';
                                case 'cohere'
                                    % coherence
                                    temp2 = wb_calculate_cohere(temp_data1,srate,[bandLimit(j,1),bandLimit(j,2)]);
                                    M_temp(:,:,i2) = temp2 + temp2';
                                case 'psi'
                                    % PSI
                                    [temp3,~] = wb_calculate_PLV(temp_data1);
                                    M_temp(:,:,i2) = temp3 + temp3';
                                case 'plv'
                                    % PLV
                                    [~,temp3] = wb_calculate_PLV(temp_data1);
                                    M_temp(:,:,i2) = temp3 + temp3';
                            end
                        end
                        M_temp1 = cat(3,M_temp1,M_temp); % concatenate arrays
                    end
                end
            elseif sign_epochLenth == 0 % sign of epochLenth is negative
                if EEG_Blocks(i1).latency - epochLenth1 > 0 && temp_badblock == 0
                    temp_windata{1,1} = EEG_temp.data(selechanns,EEG_Blocks(i1).latency - epochLenth1 : EEG_Blocks(i1).latency);
                    % temp_windata = wb_slidewin(temp_data,epochLenth1,proportion);
                    N_epoches = 1; 
                    if N_epoches >= 1
                        M_temp = zeros(N_channs,N_channs,N_epoches);
                        for i2 = 1:N_epoches
                            % extract epoch data
                            % t1 = EEG_Blocks(i1).latency + epochLenth * (i2-1) + 1 ;
                            % t2 = EEG_Blocks(i1).latency + epochLenth * i2;
                            % temp_data1 = EEG_temp.data(selechanns,t1:t2);
                            temp_data1 = temp_windata{1,i2};
                            % -----
                            switch method
                                case 'corr'
                                    % correlation
                                    temp1 = corr(temp_data1');
                                    M_temp(:,:,i2) = triu(temp1,1) + triu(temp1,1)';
                                case 'cohere'
                                    % coherence
                                    temp2 = wb_calculate_cohere(temp_data1,srate,[bandLimit(j,1),bandLimit(j,2)]);
                                    M_temp(:,:,i2) = temp2 + temp2';
                                case 'psi'
                                    % PSI
                                    [temp3,~] = wb_calculate_PLV(temp_data1);
                                    M_temp(:,:,i2) = temp3 + temp3';
                                case 'plv'
                                    % PLV
                                    [~,temp3] = wb_calculate_PLV(temp_data1);
                                    M_temp(:,:,i2) = temp3 + temp3';
                            end
                        end
                        M_temp1 = cat(3,M_temp1,M_temp); % concatenate arrays
                    end
                end
            end
        end
        M(:,:,:,j) = M_temp1;
    end
elseif flag1 == 0 % flag1 = 0: events are not contained in EEG or eventlabel is empty,all data will be used.
    for j = 1:N_freq
        disp(['calculating frequency band: ',num2str(bandLimit(j,1)),'-',num2str(bandLimit(j,2)),' Hz']);
        EEG_temp = pop_eegfiltnew(EEG,bandLimit(j,1),bandLimit(j,2),[],0); % filtering
        
        if length(size(EEG.data)) == 2
            % block_length = size(EEG.data,2);
            % N_epoches = fix(block_length/epochLenth);
            
            temp_data = EEG_temp.data(selechanns,:);
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
            M_temp = zeros(N_channs,N_channs,N_epoches);
            for i2 = 1:N_epoches
                % extract epoch data
                if length(size(EEG.data)) == 2
                    % t1 = epochLenth * (i2-1)+1;
                    % t2 = epochLenth * i2;
                    % temp_data1 = EEG_temp.data(selechanns,t1:min(t2,Nt));
                    temp_data1 = temp_windata{1,index1(i2)};
                elseif length(size(EEG.data)) == 3
                    temp_data1 = EEG_temp.data(selechanns,:,i2);
                end
                % -----
                switch method
                    case 'corr'
                        % correlation
                        temp1 = corr(temp_data1');
                        M_temp(:,:,i2) = triu(temp1,1) + triu(temp1,1)';
                    case 'cohere'
                        % coherence
                        temp2 = wb_calculate_cohere(temp_data1,srate,[bandLimit(j,1),bandLimit(j,2)]);
                        M_temp(:,:,i2) = temp2 + temp2';
                    case 'psi'
                        % PSI
                        [temp3,~] = wb_calculate_PLV(temp_data1);
                        M_temp(:,:,i2) = temp3 + temp3';
                    case 'plv'
                        % PLV
                        [~,temp3] = wb_calculate_PLV(temp_data1);
                        M_temp(:,:,i2) = temp3 + temp3';
                end
            end
            M(:,:,:,j) = M_temp;
        end
    end
end
% fisher z-tranforming
M_zscore = wb_fisherZ(M);
% -------------------------------------------------------------------------
% save results
EEG_results = [];
EEG_results.type = 'network';
EEG_results.nettype = 'WU';
try EEG_results.filename = EEG.filename;catch;end;
EEG_results.M = M;
EEG_results.M_zscore = M_zscore;

if flag1 == 2
    EEG_results.Block_percentage = 0;
elseif length(size(EEG.data)) == 3
    EEG_results.Block_percentage = 1;
else
    Ngd = size(M,3);
    EEG_results.Block_percentage = (Ngd * epochLenth1)/size(EEG.data,2);
end
% save parameters

try EEG_results.parameter.bandName = bandName;catch;end;
EEG_results.parameter.bandLimit = bandLimit;
EEG_results.parameter.eventlabel = eventlabel;
EEG_results.parameter.selechanns = selechanns;
EEG_results.parameter.epochLenth = epochLenth;
EEG_results.parameter.srate = EEG.srate;
EEG_results.parameter.method = method;
EEG_results.parameter.proportion = proportion;
EEG_results.parameter.chanlocs = EEG.chanlocs;
EEG_results.parameter.ref = EEG.ref;

% calulate mean value acorss epoches
if ~isempty(M)
    EEG_results.M_mean = reshape(nanmean(M,3),N_channs,N_channs,N_freq);
    EEG_results.M_zscore_mean = reshape(nanmean(M_zscore,3),N_channs,N_channs,N_freq);
end

