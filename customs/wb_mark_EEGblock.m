function [EEG_marked,percent,data_Z] = wb_mark_EEGblock(EEG,flag1,flag2,Thre,WinLenth,k1)
% Automatically mark the bad block/good quality EEG data based on
% thresholding z-scores.
% Detecting arifacts consists of
% (1) Z-transforming the EEG data/calculating global field power. 
%     Per channel/electrode every timepoint is z-normalized 
%     (mean subtracted and divided by standard deviation). Or the standard 
%     deviation (global field power,GFP) of the signal at all electrodes are calculated.
% (2) Averaging z-values/using global field power over channels/electrodes 
%     allows evidence for an artifact to accumulateand averaging it over channels.
% (3) Threshold the accumulated z-score/global field power for each
%     epcoh/window.The bad block data will be NOT marked as good qulity
%     data.
% -------------------------------------------------------------------------
% Input:
%      EEG: EEG structure loaded by EEGLAB (EEG data is fileterd).
%           EEG.data and EEG.srate are required at least.
%      flag1: flag1 = 0: mark bad blocks (Default).
%             flag1 = 1: mark good quality data.
%      flag2: flag2 = 0: global field power 
%             flag2 = 1: z-tranforming (Default)
%      Thre: Threshold of z-score/global field power. Default is 3.
%      WinLenth: Length of the window. unit is second.Default is 1 sec.
%      k1: Index (positive integer) of EEG channels (row number X 1, 
%            e.g. 1:32 or 32). Default is all channels (1:end).
% Output:
%      EEG_marked: EEG structure with marker. Bad blocks are labeled by
%             '9999'(EEG.event.type is '9999'). Good quality data were
%             labeled by '2001'(EEG.event.type is '2001').
%      percent: Percentage of marked event (duration).
%      data_Z: z-score or GFP of data.
%
% -------------------------------------------------------------------------
% Written by Li Dong (Li_dong729@163.com)
% $ 2017.9.29
% revised by Li Dong $ 2018.8.1
% update by Li Dong % 2019.11.11: 1)use robust estimate of SD (global filed
%                                   power) line 167
%                                 2) add isnan (lines 246 and 264) 
%                                 3) add max(1,5% time points) (line 247 and 265)
% -------------------------------------------------------------------------
if nargin < 1
    error ('One input is reqiured at least!!!!!');
elseif nargin == 1
    flag1 = 0;
    flag2 = 1;
    Thre = 3;
    WinLenth = 1;
    k1 = 1:size(EEG.data,1);
elseif nargin == 2
    if isempty(flag1)
       flag1 = 0; 
    end
    flag2 = 1;
    Thre = 3;
    WinLenth = 1;
    k1 = 1:size(EEG.data,1);
elseif nargin == 3
    if isempty(flag1)
        flag1 = 0;
    end
    if isempty(flag2)
        flag2 = 1;
    end
    Thre = 3;
    WinLenth = 1;
    k1 = 1:size(EEG.data,1);
elseif nargin == 4
    if isempty(flag1)
        flag1 = 0;
    end
    if isempty(flag2)
        flag2 = 1;
    end
    if isempty(Thre)
        Thre = 3;
    end
    WinLenth = 1;
    k1 = 1:size(EEG.data,1);
elseif nargin == 5
    if isempty(flag1)
        flag1 = 0;
    end
    if isempty(flag2)
        flag2 = 1;
    end
    if isempty(Thre)
        Thre = 3;
    end
    if isempty(WinLenth)
        WinLenth = 1;
    end
    k1 = 1:size(EEG.data,1);
elseif nargin == 6
    if isempty(flag1)
        flag1 = 0;
    end
    if isempty(flag2)
        flag2 = 1;
    end
    if isempty(Thre)
        Thre = 3;
    end
    if isempty(WinLenth)
        WinLenth = 1;
    end
end
% -------------------------------------------------------------------------
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

try
    if isempty(EEG.data)
        error('EEG.data is empty!!!!!');
    end
catch
    error('EEG.data is not exist!!!!');
end

% find label 9999 (bad block) in events
bad_index = zeros(size(EEG.data,2),1);
eventtypeflag = 0; % default event type is number

if isfield(EEG,'event')
    allevents = EEG.event;
    if ~isempty(allevents)
        % ----------
        % find bad block with label '9999'
        IndexInd = wb_findevent('9999',allevents); % find bad block
        if isempty(IndexInd)
            IndexInd = wb_findevent(9999,allevents);
        end
        if ~isempty(IndexInd)
            EEG_BadBlocks = allevents(IndexInd.index);
            disp('Find bad blocks (label 9999) in events');
            disp('Data in bad blocks will be not marked as good qulity data');
            for j2 = 1:length(EEG_BadBlocks)
                bad_index(EEG_BadBlocks(j2).latency : (EEG_BadBlocks(j2).latency + EEG_BadBlocks(j2).duration - 1)) = 1;
            end
        end
        % -------------
        % event type is character or number?
        if isfield(allevents(1),'type')
            if ischar(allevents(1).type)
                eventtypeflag = 1; % is character
            end
        end
        % -----------
    end
end

% --------------------------------------
if flag2 == 1
    data_Z = zscore(EEG.data(k1,:),0,2); % z-score
    data_Z = mean(data_Z,1);
else
    if size(EEG.data(k1,:),1)>=2
        data_Z = zscore(std(EEG.data(k1,:),0,1));    % global field power (zscored)
        
%         channelDeviation1 = 0.7413 *iqr(EEG.data(k1,:),1); % Robust estimate of SD
%         channelDeviationSD = 0.7413 * iqr(channelDeviation1(:));
%         channelDeviationMedian = nanmedian(channelDeviation1(:),1);
%         data_Z = (channelDeviation1 - channelDeviationMedian) / channelDeviationSD;
    else
        data_Z = EEG.data(k1,:); 
    end
end

if ~isempty(k1)
    [Latency] = MarkBadBlock(data_Z,Thre,srate,flag1,WinLenth,bad_index);
    Nepoch = size(Latency,1);
    percent = (sum(Latency(:,2)-Latency(:,1)))./size(EEG.data,2);
    % --------------
    % Add event in EEG structure
    if ~isfield(EEG,'event')
        k2 = 0;
    else
        k2 = length(EEG.event);
    end
    if flag1 == 0 % mark bad blocks
        if eventtypeflag == 0 % event type is number
            for i = 1:Nepoch
                EEG.event(1,k2+1).type = 9999;
                EEG.event(1,k2+1).latency = Latency(i,1);
                EEG.event(1,k2+1).duration = Latency(i,2) - Latency(i,1)+1;
                EEG.event(1,k2+1).urevent = [];
                k2 = k2 + 1;
            end
        else % event type is character
            for i = 1:Nepoch
                EEG.event(1,k2+1).type = '9999';
                EEG.event(1,k2+1).latency = Latency(i,1);
                EEG.event(1,k2+1).duration = Latency(i,2) - Latency(i,1)+1;
                EEG.event(1,k2+1).urevent = [];
                k2 = k2 + 1;
            end
        end
    else  
        % mark good quality data
        if eventtypeflag == 0 % event type is number
            for i = 1:Nepoch
                EEG.event(1,k2+1).type = 2001;
                EEG.event(1,k2+1).latency = Latency(i,1);
                EEG.event(1,k2+1).duration = Latency(i,2) - Latency(i,1)+1;
                EEG.event(1,k2+1).urevent = [];
                k2 = k2 + 1;
            end
        else % event type is character
            for i = 1:Nepoch
                EEG.event(1,k2+1).type = '2001';
                EEG.event(1,k2+1).latency = Latency(i,1);
                EEG.event(1,k2+1).duration = Latency(i,2) - Latency(i,1)+1;
                EEG.event(1,k2+1).urevent = [];
                k2 = k2 + 1;
            end
        end
    end
end
EEG_marked = EEG;
% =========================================================================
% subfunction
function [latency] = MarkBadBlock(data,Thre1,srate,flag,epoLenth,bad_index1)
    % Input:
    %      data: EEG data (1 X time points).
    %      Thre1: Threhold.
    %      flag: flag = 0: mark bad blocks.
    %            flag = 1: mark good quality data.
    %      epoLenth: length of the window. unit is second.
    %      bad_index1£ºN X 1 vector in which 1 means bad data points.
    % Output:
    %      latency: latencies with onset and end (N X 2).
    % ------------------------------------------
    Nt1 = length(data);
    L1 = srate * epoLenth;
    Nw = floor(Nt1/L1); % No. of windows
    
    if flag == 0
        index1 = reshape(abs(data(1:L1*Nw)) >= Thre1 | isnan(data(1:L1*Nw)),L1,Nw); % |z-score| >= threhold
        temp1 = repmat(double(sum(index1,1)>=fix(max(1,0.01*size(index1,1)))),L1,1);  % mark as a bad epoch, if bad time points > 1% in a epcoh.
%         if any(bad_index1)
%             temp2 = (temp1(:) + bad_index1(1:L1*Nw))>0;
%         else
%             temp2 = temp1(:);
%         end
        temp2 = temp1(:);
        L_temp1 = bwlabel(temp2); % find the connected components (i.e. label of continuous data).
        Labels_temp = unique(L_temp1);
        Labels_temp(Labels_temp == 0) = []; % take out 0.
        latency = zeros(length(Labels_temp),2);
        for j = 1:length(Labels_temp) % No. of all epochs
            Index_temp2 = find(L_temp1 == Labels_temp(j)); % index of epoch data
            latency(j,1) = Index_temp2(1);
            latency(j,2) = Index_temp2(end);
        end
    else
        index1 = reshape(abs(data(1:L1*Nw)) >= Thre1 | isnan(data(1:L1*Nw)),L1,Nw); % |z-score| >= threhold 
        temp1 = repmat(double(sum(index1,1)<=fix(max(1,0.05*size(index1,1)))),L1,1); % mark as a good epoch, if bad time points < 5% in a epoch.
        if any(bad_index1)
            temp2 = (temp1(:) - bad_index1(1:L1*Nw))>0;
        else
            temp2 = temp1(:);
        end
       
        L_temp1 = bwlabel(temp2); % find the connected components (i.e. label of continuous data).
        Labels_temp = unique(L_temp1);
        Labels_temp(Labels_temp == 0) = []; % take out 0.
        latency = zeros(length(Labels_temp),2);
        for j = 1:length(Labels_temp) % No. of all epochs
            Index_temp2 = find(L_temp1 == Labels_temp(j)); % index of epoch data
            latency(j,1) = Index_temp2(1);
            latency(j,2) = Index_temp2(end);
        end
    end
end
end