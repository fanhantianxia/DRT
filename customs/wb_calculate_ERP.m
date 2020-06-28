function EEG = wb_calculate_ERP(EEG,eventl,epochlimits,valuelim_1,valuelim_2,valuelim_3,selechanns,marker1,t1)
% Creating averaged event related potential (ERP)
% Consisting of
%  [1] Extract epochs and Baseline correstion: Convert a continuous EEG 
%      dataset to epoched data by extracting data epochs time locked to 
%      specified event types or event indices. If applicable, time locked 
%      events corresponding to correct-response marker will be extracted (i.e. event2).
%  [2] Artifact rejection in epoched data using simple voltage threshold.
%      Three criterions including amplitude, gradient and max-min criterions
%      were used.If applicable, events in bad block (label 9999, marked by
%      wb_pipeline_EEG_Mark) will also be rejected automatically.
%  [3] Epochs (default is [-0.2 0.8] s) will be generated and averaged.

% References:

% -------------------------------------------------------------------------
% Input:
%    EEG: EEG structure imported using EEGLAB. EEG.data should be channels
%         X time points OR channels X time points X epoches. sampling rate
%         must be contained in EEG as EEG.srate.
%    epochlimits: epoch latency range [start, end] in seconds relative
%         to the time-locking events. Default is [-0.2, 0.8].
%    event1: specified event types or event indices (e.g. event label)
%        If event label is not found, NO data will be epoched and calculated.
%    valuelim_1: Amplitude criterion: Lower and upper bound latencies for
%        trial data.If one positive value is given,the opposite value is
%        used for lower bound.For example, use [-100,100] to remove artifactual
%        epoch.Default is [-100,100].
%    valuelim_2: Gradient criterion: maximum allowed voltage step/sampling
%        point.Default is 50.
%    valuelim_3: Max-min criterion: maximum allowed absolute difference in
%       the segment/epoch.Default is 150.
%    selechanns: number with indices of the selected channels
%        (e.g. [1:4,7:30] or 'all').Default is 'all';
%    marker1: correct-response marker corresponding to the specified event
%        (e.g. event1). Default is [].
%    t1: duration (in seconds) before correct-response marker. Default is 2s.
%  Output:
%    EEG: Input dataset. Data may already be epoched; in this case,
%         extract (shorter) subepochs time locked to epoch events.

%         EEG.eventlist: list of accepted events;
%         EEG.erp: averaged event-related potentials for each channel.
%         EEG.trials: No. of trials;
%         EEG.xmin: Epoch latency limits [start] in seconds;
%         EEG.xmax: Epoch latency limits [end] in seconds;
%         EEG.epoch: filling with values of other events in the same epochs;
%         
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Li_dong729@163.com)
% $ 2018.4.20
% Fix bug: try EEG.chanlocs = EEG.chanlocs(1,selechanns);catch;end;  Revised by Li Dong $ 2018.8.8
% -------------------------------------------------------------------------
% exmaple:
% eventl = 'S 22';
% epochlimits = [-0.2,0.8];
% valuelim_1 = 100;
% selechanns = 'all';
% marker1 = 'S222';
% t1 = 2;
% valuelim_2 = 50;
% valuelim_3 = 150;
% EEG = wb_calculate_ERP(EEG,eventl,epochlimits,valuelim_1,valuelim_2,valuelim_3,selechanns,marker1,t1);
% EEG = wb_calculate_ERP(EEG,eventl);
% EEG = wb_calculate_ERP(EEG,eventl,[],valuelim_1,valuelim_2,valuelim_3,[],marker1,t1)
% % -------------------------------------------------------------------------
if nargin < 2
    error ('2 inputs are reqiured at least!!!!!');
elseif nargin == 2
    epochlimits = [-0.2,0.8];
    valuelim_1 = 100;
    valuelim_2 = 50;
    valuelim_3 = 150;
    selechanns = 'all';
    marker1 = [];
    t1 = 2;
elseif nargin == 3
    valuelim_1 = 100;
    valuelim_2 = 50;
    valuelim_3 = 150;
    selechanns = 'all';
    marker1 = [];
    t1 = 2;
elseif nargin == 4
    valuelim_2 = 50;
    valuelim_3 = 150;
    selechanns = 'all';
    marker1 = [];
    t1 = 2;
elseif nargin == 5
    valuelim_3 = 150;
    selechanns = 'all';
    marker1 = [];
    t1 = 2;
elseif nargin == 6
    selechanns = 'all';
    marker1 = [];
    t1 = 2;
elseif nargin == 7
    marker1 = [];
    t1 = 2;
elseif nargin == 8
    t1 = 2;
end
% % ----------------

if length(size(EEG.data)) == 3
    disp('EEG.data has been epoched,averaging ERP ONLY');
    DIM = size(EEG.data);
    disp(['No. of channels:',num2str(DIM(1))]);
    disp(['No. of trials:',num2str(DIM(3))]);
    disp(['EEG epoch legnth:',num2str(DIM(2))]);
    EEG.erp = nanmean(EEG.data,3);
    EEG.eventlist = event;
else
    % ---------------------
    % check inputs
    try
        srate = EEG.srate;
        if isfinite(srate) && length(size(EEG.data)) == 2
            disp(['sampling rate = ',num2str(srate)])
            if isempty(epochlimits)
                disp('use default epoch range')
                epochlimits = [-0.2,0.8];
            end
            disp(['epoch range = [',num2str(epochlimits(1)),',',num2str(epochlimits(2)),'] s']);
            epochLenth1 = round(epochlimits * srate);
        end
    catch
        disp('sampling rate is not found in EEG.');
        error('sampling rate is not found in EEG.');
    end
    
    % No. of timepoints
    Nt = size(EEG.data,2);
    disp(['EEG data length:',num2str(Nt)]);
    
    % valuelim_1
    if isempty(valuelim_1)
        valuelim_1 = [-100,100];
    else
        if length(valuelim_1) == 1
            valuelim_1 = [-abs(valuelim_1),abs(valuelim_1)];
        end
    end
    disp(['Amplitude criterion [min,max]:  [',num2str(valuelim_1(1)),',',num2str(valuelim_1(2)),']']);
    
    % valuelim_2
    if isempty(valuelim_2)
        valuelim_2 = 50;
    end
    disp(['Gradient criterion [maximum allowed voltage step/sampling point]: ',num2str(valuelim_2)]);
    
    % valuelim_3
    if isempty(valuelim_3)
        valuelim_3 = 150;
    end
    
    disp(['Max-Min criterion [maximum allowed absolute difference]: ',num2str(valuelim_3)]);
    
    
    % channs
    if isequal(selechanns,'all') || isempty(selechanns)
        selechanns = 1:size(EEG.data,1);
    end
    N_channs = length(selechanns); % No. of selected channs
    disp(['No. of selected channels: ',num2str(N_channs)]);
    
    % marker1 and t1
    if isempty(marker1)
        flag1 = 0;
    else
        flag1 = 1;
        if isempty(t1) || t1 < 0
            disp('t1 is empty or < 0, use default value 2s');
            t1 = 2;
        end
    end
    % ---------------------------------------------------------------------
    % generate event list
    disp('Generate event list');
    % events
    IndBadBlock = []; % Index of bad blocks
    Indmarker1 = [];  % Index of marker1
    if flag1 == 0 % marker1 is empty
        if ~isempty(eventl)
            % find eventlabel and bad block
            if isfield(EEG,'event')
                allevents = EEG.event;
                if ~isempty(allevents)
                    IndexInd = wb_findevent(eventl,allevents);
                    IndexInd_badblock = wb_findevent('9999',allevents); % find bad block
                    if isempty(IndexInd_badblock)
                        IndexInd_badblock = wb_findevent(9999,allevents);
                    end
                    % ----------                    
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
                    
                    if ~isempty(IndexInd)
                        EEG_epochs = allevents(IndexInd.index);
                        flag2 = 1; % flag1 = 1: event type is found in EEG events;
                        disp(['event type: ',num2str(eventl),' is found in events']);
                    else
                        flag2 = 0; % flag1 = 0: event type is not found in events.
                        warning(['event type (', num2str(eventl),')is not found in events,No data is used to calculate ERP']);
                        EEG = [];
                        % error(['event type (', num2str(eventl),')is not found in events,No data is used to calculate ERP']);
                    end
                else
                    disp('events are empty in EEG');
                    error('events are empty in EEG');
                end
            else
                disp('events are not found in EEG');
                error('events are not found in EEG');
            end
        else
            disp('event1 is empty');
            error('event1 is empty');
        end
        %-----------------------------------------------
        trials = [];
        k2 = 1;
        C1 = [];
        C2 = [];
        C3 = [];
        if flag2 == 1 % flag1 = 1: event1 is found in EEG events
            % extract epochs
            disp('Extract epochs and baseline correction...');
            N_epochs = length(EEG_epochs);
            for i1 = 1:N_epochs
                % make sure the EEG blocks does not contain bad blocks.
                if  ~isempty(IndBadBlock)
                    Ind_badblock = IndBadBlock(1,max(1,EEG_epochs(i1).latency + epochLenth1(1)):min(Nt,EEG_epochs(i1).latency + epochLenth1(2)));
                    temp_badblock = sum(Ind_badblock(:));
                else
                    temp_badblock = 0;
                end
                
                if temp_badblock == 0
                    temp_data = EEG.data(selechanns,max(1,EEG_epochs(i1).latency + epochLenth1(1)):min(Nt,EEG_epochs(i1).latency + epochLenth1(2)));
                    temp_data = temp_data - repmat(mean(temp_data(:,1:abs(epochLenth1(1))),2),1,size(temp_data,2)); % baseline correction
                    
                    trials(:,:,k2) = temp_data;
                    
                    % artifact detection
                    for j1 = 1:size(temp_data,1) % channels
                        temp1 = temp_data(j1,:);
                        
                        % Amplitude Criterion: maximum and minimum
                        % amplitude
                        value1 = [min(temp1),max(temp1)];
                        C1(j1,k2) = value1(1) >= valuelim_1(1) && value1(2) <= valuelim_1(2);
                        
                        % Gradient Criterion: maximum allowed voltage
                        % step/sampling point
                        value2 = max(abs(diff(temp1)));
                        C2(j1,k2) = value2 <= valuelim_2;
                        
                        % Max-Min Criterion: maximum allowed absolute
                        % difference in the segment
                        value3 = max(temp1) - min(temp1);
                        C3(j1,k2) = value3 <= valuelim_3;
                        
                    end
                    k2 = k2 + 1;
                end
            end
            % artifact rejection and averaging
            disp('Artifact rejection...');
            all_C = C1 + C2 + C3;
            Logic_1 = all_C == 3; % satisfying 3 Criterion at same time
            Logic_2 = sum(Logic_1) == size(trials,1); % clean trials for each channel
            trials = trials(:,:,Logic_2);
            
            disp('Averaging...')
            % temp_ERP.ntrial = size(trials,3);
            temp_ERP = mean(trials,3);
            
            %             for j1 = 1:size(temp_data,1) % channels
            %                 temp_clean_trial = trials(j1,:,Logic_1(j1,:));
            %
            %                 temp_ERP.ntrial(j1,1) = size(temp_clean_trial,3); % retained trials for each channel
            %                 temp_ERP.erp(j1,:) = mean(temp_clean_trial,3);  % averaging;
            %             end
            % -----------------------
            % save in EEG
            EEG.erp = temp_ERP;
            EEG.xmin = epochlimits(1);
            EEG.xmax = epochlimits(2);
            EEG.data = trials;
            EEG.trials = size(trials,3);
            EEG.pnts = size(trials,2);
            EEG.icaact = [];
            EEG.eventlist = EEG_epochs(1,Logic_2);
            try EEG.chanlocs = EEG.chanlocs(1,selechanns);catch;end;
            EEG.nbchan = length(selechanns);
            
            EEG.event = EEG_epochs(1,Logic_2);
            for i1 = 1:size(trials,3)
                EEG.event(i1).latency = 1 - epochLenth1(1) + EEG.pnts*(i1-1);
                EEG.event(i1).epoch = i1;
            end
            EEG.epoch = [];
            EEG.times = 1:size(trials,2);
            EEG = eeg_checkset(EEG, 'eventconsistency');
            
            disp(['No. of trials:',num2str(EEG.trials)]);
            disp(['EEG epoch legnth:',num2str(EEG.pnts)]);
        end
        
        
    elseif flag1 == 1 % marker1 is used
        
        if ~isempty(eventl)
            % find eventlabel and bad block
            if isfield(EEG,'event')
                allevents = EEG.event;
                if ~isempty(allevents)
                    
                    IndexInd = wb_findevent(eventl,allevents); % find eventlabel
                    IndexInd_marker1 = wb_findevent(marker1,allevents); % find correct-response marker
                    IndexInd_badblock = wb_findevent('9999',allevents); % find bad block
                    if isempty(IndexInd_badblock)
                        IndexInd_badblock = wb_findevent(9999,allevents);
                    end
                    % --------------
                    
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
                    
                    
                    if ~isempty(IndexInd_marker1)
                        Marker1Blocks = allevents(IndexInd_marker1.index);
                        Indmarker1 = zeros(1,Nt);
                        for k1 = 1:length(Marker1Blocks)
                            Indmarker1(max(Marker1Blocks(1,k1).latency - round(t1*srate),1):Marker1Blocks(1,k1).latency) = 1;
                        end
                        disp(['correct-response marker (', num2str(marker1),'), were found in EEG events']);
                    else
                        disp(['correct-response marker (', num2str(marker1),') is not found in events']);
                        disp(['If possible, all events of ', num2str(eventl), ' will be used to calculate ERP']);
                    end
                    
                    if ~isempty(IndexInd)
                        EEG_epochs = allevents(IndexInd.index);
                        flag2 = 1; % flag2 = 1: event1 is found in EEG events;
                        disp(['event type: ',num2str(eventl),' is found in events']);
                    else
                        flag2 = 0; % flag2 = 0: event1 is not found in events.
                        warning(['event1(', num2str(eventl),')is not found in events, No data is used to calculate ERP']);
                        EEG = [];
                        % error(['eventlabel(', num2str(eventl),')is not found in events, No data is used to calculate ERP']);
                    end
                    
                    
                else
                    disp('events are empty in EEG');
                    error('events are empty in EEG');
                end
            else
                disp('events are not found in EEG');
                error('events are not found in EEG');
            end
        else
            disp('event1 is empty');
            error('event1 is empty');
        end
        % -----------------------------------------------------------------
        trials = [];
        k2 = 1;
        C1 = [];
        C2 = [];
        C3 = [];
        if flag2 == 1 % flag1 = 1: event1 is found in EEG events
            % extract epochs
            disp('Extract epochs and baseline correction');
            N_epochs = length(EEG_epochs);
            for i1 = 1:N_epochs
                % make sure the EEG epochs does not contain bad blocks.
                if  ~isempty(IndBadBlock)
                    Ind_badblock = IndBadBlock(1,max(1,EEG_epochs(i1).latency + epochLenth1(1)):min(Nt,EEG_epochs(i1).latency + epochLenth1(2)));
                    temp_badblock = sum(Ind_badblock(:));
                else
                    temp_badblock = 0;
                end
                
                % make sure the EEG epochs contain marker1  
                if  ~isempty(Indmarker1)
                    Ind_marker1 = Indmarker1(1,max(1,EEG_epochs(i1).latency + epochLenth1(1)):min(Nt,EEG_epochs(i1).latency + epochLenth1(2)));
                    temp_marker1 = any(Ind_marker1);
                else
                    temp_marker1 = 1;
                end
                
                if temp_badblock == 0 && temp_marker1 == 1
                    temp_data = EEG.data(selechanns,max(1,EEG_epochs(i1).latency + epochLenth1(1)):min(Nt,EEG_epochs(i1).latency + epochLenth1(2)));
                    temp_data = temp_data - repmat(mean(temp_data(:,1:abs(epochLenth1(1))),2),1,size(temp_data,2)); % baseline correction
                    
                    trials(:,:,k2) = temp_data;
                    
                    % artifact detection
                    for j1 = 1:size(temp_data,1) % channels
                        temp1 = temp_data(j1,:);
                        
                        % Amplitude Criterion: maximum and minimum
                        % amplitude
                        value1 = [min(temp1),max(temp1)];
                        C1(j1,k2) = value1(1) >= valuelim_1(1) && value1(2) <= valuelim_1(2);
                        
                        % Gradient Criterion: maximum allowed voltage
                        % step/sampling point
                        value2 = max(abs(diff(temp1)));
                        C2(j1,k2) = value2 <= valuelim_2;
                        
                        % Max-Min Criterion: maximum allowed absolute
                        % difference in the segment
                        value3 = max(temp1) - min(temp1);
                        C3(j1,k2) = value3 <= valuelim_3;
                        
                    end
                    k2 = k2 + 1;
                end
            end
            % artifact rejection and averaging
            disp('Artifact rejection...');
            all_C = C1 + C2 + C3;
            Logic_1 = all_C == 3; % satisfying 3 Criterion at same time
            Logic_2 = sum(Logic_1) == size(trials,1); % clean trials for each channel
            trials = trials(:,:,Logic_2);
            
            disp('Averaging...')
            % temp_ERP.ntrial = size(trials,3);
            temp_ERP = mean(trials,3);
            
            %             for j1 = 1:size(temp_data,1) % channels
            %                 temp_clean_trial = trials(j1,:,Logic_1(j1,:));
            %
            %                 temp_ERP.ntrial(j1,1) = size(temp_clean_trial,3); % retained trials for each channel
            %                 temp_ERP.erp(j1,:) = mean(temp_clean_trial,3);  % averaging;
            %             end
            % -----------------------
            % save in EEG
            EEG.erp = temp_ERP;
            EEG.xmin = epochlimits(1);
            EEG.xmax = epochlimits(2);
            EEG.data = trials;
            EEG.trials = size(trials,3);
            EEG.pnts = size(trials,2);
            EEG.icaact = [];
            EEG.eventlist = EEG_epochs(1,Logic_2);
            try EEG.chanlocs = EEG.chanlocs(1,selechanns);catch;end;
            EEG.nbchan = length(selechanns);
            
            EEG.event = EEG_epochs(1,Logic_2);
            for i1 = 1:size(trials,3)
                EEG.event(i1).latency = 1 - epochLenth1(1) + EEG.pnts*(i1-1);
                EEG.event(i1).epoch = i1;
            end
            EEG.epoch = [];
            EEG.times = 1:size(trials,2);
            EEG = eeg_checkset(EEG, 'eventconsistency');
                   
            disp(['No. of trials:',num2str(EEG.trials)]);
            disp(['EEG epoch legnth:',num2str(EEG.pnts)]);
        end
    end
    
end