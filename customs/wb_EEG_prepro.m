function EEGprepro = wb_EEG_prepro(EEG,selechanns,EOGchanns,thre_ODQ,passband,PowerFrequency,...
                keepUnselectChannsFlag,badChannelInterploateFlag,residualArtifactRemovalFlag,...
                MARA_thre,WindowSeconds,HighPassband,badWindowThreshold,robustDeviationThreshold,...
                FrequencyNoiseThreshold,correlationThreshold)
% -------------------------------------------------------------------------
% Instruction of EEG preprocessing pipeline
% step 1£©: Quality Assessment of EEG data. Noting that QA will not change the EEG data
% step 2£©: Passband and notch filtering
% step 3£©: Artifact removal: EOG regression
% step 4) : Artifact removal: residual artifact removal.
% step 5£©: Bad channel interpolation and re-referencing
% step 6) : Quality Assessment of EEG data after artifact removing
% step 7) : Marking residual bad block with unusually high or low amplitude using zscored STD.

% References:
% Pedroni, A., et al. (2019). "Automagic: Standardized preprocessing of big EEG data." Neuroimage 200: 460-473.
% Bigdely-Shamlo, N., et al. (2015). "The PREP pipeline: standardized preprocessing for large-scale EEG analysis." Front Neuroinform 9.
% Cowley, B. U., et al. (2017). "Computational testing for automated preprocessing: a Matlab toolbox to enable large scale electroencephalography data processing." PeerJ Computer Science 3: e108.
% Mognon, A., Jovicich, J., Bruzzone, L., & Buiatti, M. (2011). ADJUST: An automatic EEG artifact detector based on the joint use of spatial and temporal features. Psychophysiology, 48(2), 229-240.
% Winkler, I., et al. (2011). "Automatic classification of artifactual ICA-components for artifact removal in EEG signals." Behav Brain Funct 7: 30.
% -------------------------------------------------------------------------
% Input:
%      EEG:  EEG structure loaded by EEGLAB (EEG data is raw data).
%            EEG.data and EEG.srate are required at least. EEG.channlocs is
%            stronly recommended to be contained in EEG.
%      selechanns : number with indices of the selected channels 
%                      (e.g. [1:4,7:30] or 'all').Default is 'all';
%      EOGchanns:  number with indices of the EOG channels.Default is [].
%      thre_ODQ: threshold of Overall Data Quality. Default is if ODQ > 80, 
%                the preprocessing could be continue.
%      passband: pass bands of filtering. default is [1,40].
%      PowerFrequency:         power frequency. Default is 50 Hz (in
%                              Chinese). Noting that in USA, power frequency is 60Hz.
%      keepUnselectChannsFlag: keepUnselectChannsFlag = 0: do not keep unselected channels (default); 
%                              keepUnselectChannsFlag = 1: keep all channels;
%      badChannelInterploateFlag: badChannelInterploateFlag = 0: do NOT interpolate, and if have 
%                                      channel locations in EEG.chanlocs, then re-referencing to REST; 
%                                 badChannelInterploateFlag = 1: interpolate the bad channels rows of EEG.data 
%                                      using reference electrode standardization interpolation technique (RESIT); 
%                                      default is using RESIT (The bad channels will be interpolated with REST reference); 
%                                 badChannelInterploateFlag = 2: interpolate the bad channels rows of EEG.data using 
%                                      spherical spline interpolation (SSI), and then re-referencing to REST.
%      residualArtifactRemovalFlag: residualArtifactRemovalFlag = 0: not removal;
%                                   residualArtifactRemovalFlag = 1: ICA based MARA;
%                                   residualArtifactRemovalFlag = 2: ICA based ADJUST;
%                                   residualArtifactRemovalFlag = 3: rPCA;
%                                   residualArtifactRemovalFlag = 4: ASR method;
%      MARA_thre:      cuttoff posterior probability for each IC of being
%                      an artefact while using MARA.
%      WindowSeconds : window size in seconds (default = 1 sec)
%      HighPassband :  lower edge of the frequency for high pass filtering,
%                      Hz.Default is 1 Hz.
%      badWindowThreshold : cutoff fraction of bad windows (default = 0.4)
%                      for detecting bad channels.

%      robustDeviationThreshold : Z-score cutoff for robust channel deviation (default = 5)
%      FrequencyNoiseThreshold: Z-score cutoff for SNR (signal above 40 Hz).Default is 3;
%      correlationThreshold : maximal correlation below which window is bad (range is (0,1), default = 0.6)

% Output: a structure of EEG.
%  Output:
%    EEG.preprocessed: details of each steps including parameters and results.
% -------------------------------------------------------------------------
% Usage of function
% clear all;clc;
% 
% EEG = pop_loadset('F:\Works\work1\WeBrain_Platform\WeBrain_pipelines\EEG\TEST_data\ERP-data\s1_avg_filtered1-30Hz.set');
% % EEG = pop_importdata('data','F:\Works\work1\WeBrain_Platform\WeBrain_pipelines\EEG\TEST_data\Matlab-DAT\','dataformat','matlab');
% 
% selechanns = [1:31,33:62];    % number with indices of the selected channels (e.g. [1:4,7:30] or 'all').Default is 'all';
% EOGchanns = [32,63];       % number with indices of the EOG channels.
% thre_ODQ = 80;          % threshold of Overall Data Quality. Default is if ODQ > 80, the preprocessing could be continue.
% 
% passband = [1,40];      % pass band 
% PowerFrequency = 50;              % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
% 
% keepUnselectChannsFlag = 1;   % keepUnselectChannsFlag = 0: do not keep unselected channels (default); 
%                               % keepUnselectChannsFlag = 1: keep all channels;
% badChannelInterploateFlag = 1; % badChannelInterploateFlag = 0: do NOT interpolate, and if have 
                                 % channel locations in EEG.chanlocs, then re-referencing to REST; 
                                 % badChannelInterploateFlag = 1: interpolate the bad channels rows of EEG.data 
                                 % using reference electrode standardization interpolation technique (RESIT); 
                                 % default is using RESIT (The bad channels will be interpolated with REST reference); 
                                 % badChannelInterploateFlag = 2: interpolate the bad channels rows of EEG.data using 
                                 % spherical spline interpolation (SSI), and then re-referencing to REST.
% residualArtifactRemovalFlag = 1;  % residualArtifactRemovalFlag = 0: not removal (default);
%                                   % residualArtifactRemovalFlag = 1: MARA;
%                                   % residualArtifactRemovalFlag = 2: ADJUST;
%                                   % residualArtifactRemovalFlag = 3: rPCA;
%                                   % residualArtifactRemovalFlag = 4: cleandata;
% MARA_thre = 0.8;         % cuttoff posterior probability for each IC of being an artefact using MARA
% 
% WindowSeconds = 1;       % window size in seconds (default = 1 sec)
% HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
% badWindowThreshold = 0.35;         % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
% robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
% FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
% correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
% 
% EEG = wb_EEG_prepro(EEG,selechanns,EOGchanns,thre_ODQ,passband,PowerFrequency,...
%                 keepUnselectChannsFlag,badChannelInterploateFlag,residualArtifactRemovalFlag,...
%                 MARA_thre,WindowSeconds,HighPassband,badWindowThreshold,robustDeviationThreshold,...
%                 FrequencyNoiseThreshold,correlationThreshold);
% -------------------------------------------------------------------------
% Code Summary for working in School of Life Science and Technology,UESTC.
% Author: Li Dong, e-mail: Lidong@uestc.edu.cn
% This template is for non commercial use only.
% It is freeware but not in the public domain.

% Written by Li Dong (Lidong@uestc.edu.cn), UESTC
% $ 2019.11.5
% -------------------------------------------------------------------------
% check all inputs
if nargin < 1
    error ('One input is reqiured at least!!!!!');
elseif nargin == 1
    selechanns = 'all';     % number with indices of the selected channels (e.g. [1:4,7:30] or 'all').Default is 'all';
    EOGchanns = [];         % number with indices of the EOG channels.
    thre_ODQ = 80;          % threshold of Overall Data Quality. Default is if ODQ > 80, the preprocessing could be continue.
    passband = [1,40];      % pass band
    PowerFrequency = 50;    % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default); 
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 2
    EOGchanns = [];         % number with indices of the EOG channels.
    thre_ODQ = 80;          % threshold of Overall Data Quality. Default is if ODQ > 80, the preprocessing could be continue.
    passband = [1,40];      % pass band
    PowerFrequency = 50;    % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default);
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 3
    thre_ODQ = 80;          % threshold of Overall Data Quality. Default is if ODQ > 80, the preprocessing could be continue.
    passband = [1,40];      % pass band
    PowerFrequency = 50;    % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default);
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 4
    passband = [1,40];      % pass band
    PowerFrequency = 50;    % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default);
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 5
    PowerFrequency = 50;    % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default);
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 6
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default);
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 7
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 8
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 9
    MARA_thre = 0.70;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.7.
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 10
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 11
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 12
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 13
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 14
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
elseif nargin == 15
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
end
% -------------------------------------------------------------------------
% checking all inputs
if isempty(selechanns)
    selechanns = 'all';     % number with indices of the selected channels (e.g. [1:4,7:30] or 'all').Default is 'all';
end
if any(~isfinite(EOGchanns))
    EOGchanns = [];         % number with indices of the EOG channels.
end
if isempty(thre_ODQ)
    thre_ODQ = 80;          % threshold of Overall Data Quality. Default is if ODQ > 80, the preprocessing could be continue.
end
if isempty(passband)
    passband = [1,40];      % pass band
end
if isempty(PowerFrequency)
    PowerFrequency = 50;    % power frequency. default is 50 Hz (in Chinese). Noting that in USA, power frequency is 60Hz
end

if isempty(keepUnselectChannsFlag) || keepUnselectChannsFlag < 0
    keepUnselectChannsFlag = 0;      % keepUnselectChannsFlag = 0: do not keep unselected channels (default);
end
if isempty(badChannelInterploateFlag)
    badChannelInterploateFlag = 1;   % badChannelInterploateFlag = 1: REST interpolate and referencing to REST;
end
if isempty(residualArtifactRemovalFlag)
    residualArtifactRemovalFlag = 0; % residualArtifactRemovalFlag = 0: not removal (default);
end
if isempty(MARA_thre)
    MARA_thre = 0.7;         % cuttoff posterior probability for each IC of being an artefact using MARA.default is 0.8.
end
if isempty(WindowSeconds)
    WindowSeconds = 1;       % window size in seconds (default = 1 sec)
end
if isempty(HighPassband)
    HighPassband = 1;        % lower edge of the frequency for high pass filtering, Hz. default is 1 Hz.
end
if isempty(badWindowThreshold)
    badWindowThreshold = 0.4;        % cutoff fraction of bad windows (default = 0.4) for detecting bad channels.
end
if isempty(robustDeviationThreshold)
    robustDeviationThreshold = 5;     % Z-score cutoff for robust channel deviation (default = 5)
end
if isempty(FrequencyNoiseThreshold)
    FrequencyNoiseThreshold = 3;      % Z-score cutoff for SNR (signal above 40 Hz).Default is 3.
end
if isempty(correlationThreshold)
    correlationThreshold = 0.6;       % correlation below which window is bad (range is (0,1), default = 0.6)
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
            warning('because most of your channels do not have X,Y,Z location measurements, some methods (inlcuding the bad channel interpolation, RANSAC method, ADJUST and MARA) are ignored');
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
                warning('some selected channels do not have xyz coordinates,some methods (inlcuding the bad channel interpolation, RANSAC method, ADJUST and MARA) are ignored');
            end
            
        end
    else
        chanlocsflag = 0;
        warning('channel location is not found, some methods (inlcuding the bad channel interpolation, RANSAC method, ADJUST and MARA) are ignored');
    end
catch
     warning('there are some unknow errors of channel locations, some methods (inlcuding the bad channel interpolation, RANSAC method, ADJUST and MARA) are ignored');
end
% check other parameters

if ~isempty(PowerFrequency) && (PowerFrequency > 5)
    NotchBand = [PowerFrequency - 5, PowerFrequency + 5];
else
    NotchBand = [];
end

% =========================================================================
disp('-------------------------------------------------------------');
disp('Quality Assessment of EEG data before artifact removing...');
% -------------------------------------------------------------------------
% [1] Quality Assessment of EEG data

ransacSampleSize = 50;            % samples for computing ransac (default = 50)
ransacChannelFraction = 0.3;      % fraction of channels for robust reconstruction (default = 0.3)
ransacCorrelationThreshold = [];  % cutoff correlation for abnormal wrt neighbors(default = [] | --> not performed)
flagNotchFilter = 0;              % flagNotchFilter = 1: remove 0.5* power frequency noise using notch filtering.Default is off (flagNotchFilter = 0).

[results_QA] = wb_EEG_QA(EEG,WindowSeconds,HighPassband,selechanns,badWindowThreshold,...
    robustDeviationThreshold,PowerFrequency,FrequencyNoiseThreshold,flagNotchFilter,correlationThreshold,...
    ransacCorrelationThreshold,ransacChannelFraction,ransacSampleSize);
% EEG.preprocessed.QA.check = 'yes';
disp(['Overall Data Quality (raw data): ',num2str(results_QA.ODQ)]);
disp(['Data Quality Rating (raw data): ',results_QA.DataQualityRating]);
disp('-------------------------------------------------------------');
if results_QA.ODQ >= thre_ODQ
    % ---------------------------------------------------------------------
    % [2] filtering
    % passband filtering
    if ~isempty(passband) && isnumeric(passband) && all(isfinite(selechanns))
        disp('Passband filter data using Hamming windowed sinc FIR filter');
        disp(['Passband:',num2str(passband(1)),'-',num2str(passband(2)),' Hz']);
        EEG = pop_eegfiltnew(EEG,passband(1),passband(2),[],0);
        EEG.preprocessed.PassbandFilter.check = 'yes';
        EEG.preprocessed.PassbandFilter.passband = passband;
        EEG.preprocessed.PassbandFilter.comments = 'Hamming windowed sinc FIR filter';
    else
        disp('No passband filtering')
        EEG.preprocessed.PassbandFilter.check = 'no';
    end
    %
    % notch filtering
    if (~isempty(passband) && max(passband)>=min(NotchBand)) || (~isempty(NotchBand)&&isempty(passband))
        disp('Notch filtering...');
        disp(['Notch band: [',num2str(NotchBand),']']);
        EEG = pop_eegfiltnew(EEG,NotchBand(1),NotchBand(2),[],1);
        EEG.preprocessed.NotchFilter.check = 'yes';
        EEG.preprocessed.NotchFilter.notchband = NotchBand;
    else
        disp('Skip Notch filtering: max(passband)<min(NotchBand) or NotchBand is uncorrect or empty');
        EEG.preprocessed.NotchFilter.check = 'no';
    end
    % ---------------------------------------------------------------------
    % [3] artifact removal: EOG regression
    if ~isempty(EOGchanns) && all(isfinite(EOGchanns))
        disp('Artifact Removal: EOG regression...');
        try badchanns = results_QA.badChannels;catch; badchanns = [];end;
        goodchanns = setdiff(selechanns,badchanns);
        eegY = EEG.data(goodchanns,:)';
        eogX = [zscore(EEG.data(EOGchanns,:)'),ones(size(eegY,1),1)];
        beta1 = eogX \ eegY;
        eegclean =  eegY - eogX(:,1:size(beta1,1)-1) * beta1(1:size(beta1,1)-1,:);
        EEG.data(goodchanns,:) = eegclean';
        % Write back what has happened
        EEG.preprocessed.EOGregression.check = 'yes';
        EEG.preprocessed.EOGregression.EOGchanns = EOGchanns;
    else
        EEG.preprocessed.EOGregression.check = 'no';
        EEG.preprocessed.EOGregression.EOGchanns = [];
        disp('Skip EOG regression: EOG channel is empty or unfinite...');
    end
    % ---------------------------------------------------------------------
    % [4] residual artifact removal
    if residualArtifactRemovalFlag >= 0 && residualArtifactRemovalFlag <= 4
        switch residualArtifactRemovalFlag
            case 0
                EEG.preprocessed.residualArtifactRemoval.check = 'no';
                disp('Skip residual artifact removal: unused...');
            case 1  % ICA based MARA             
                if chanlocsflag == 1
                    disp('Artifact Removal: removing residual artifact using ICA based MARA...');
                    try badchanns = results_QA.badChannels;catch; badchanns = [];end;
                    goodchanns = setdiff(selechanns,badchanns);
                    % goodchanns = selechanns; % ICA on all selected channels
                    temp_EEG = EEG;
                    % unionchanns = union(selechanns,EOGchanns);
                    temp_EEG.data = EEG.data(goodchanns,:);
                    temp_EEG.chanlocs = EEG.chanlocs(1,goodchanns);
                    temp_EEG.nbchan = length(goodchanns);
                    if  isfield(temp_EEG,'icaweights') && isempty(temp_EEG.icaweights)
                        ICs = 'default';      % number of ICA components to compute (default -> chans or 'pca' arg)
                        Ntrain = 0;           % perform tanh() "extended-ICA" with sign estimation N training blocks.
                        PCs = size(temp_EEG.data,1)-1;      % decompose a principal component (default = 0 -> off).
                        % PCs = 'default';
                        stop = 1e-6;      % stop training when weight-change < this
                        MaxSteps = 512;   % max number of ICA training steps
                        sphering = 'on';  % ['on'/'off'] flag sphering of data
                        temp_EEG = wb_runICA(temp_EEG,'all',ICs,Ntrain,PCs,stop,MaxSteps,sphering);
                        
                        EEG.preprocessed.residualArtifactRemoval.MARA.ICs = ICs;
                        EEG.preprocessed.residualArtifactRemoval.MARA.ICANtrain = Ntrain;
                        EEG.preprocessed.residualArtifactRemoval.MARA.ICAstop = stop;
                        EEG.preprocessed.residualArtifactRemoval.MARA.ICAMaxSteps = MaxSteps;
                        EEG.preprocessed.residualArtifactRemoval.MARA.ICAsphering = sphering;
                        
                    else
                        disp('Skip ICA decomposition: ICA has been conducted on EEG data...');
                    end
                    [~, MARAinfo] = MARA(temp_EEG); % MARA
                    try artcomps = find(MARAinfo.posterior_artefactprob > MARA_thre);catch; artcomps = [];end; % posterior probability for each IC of being an artefact > Thre
                    
                    if isempty(artcomps)
                        disp('''artcomps'' is empty, there may be NO artifact components detected by MARA....');
                    end
                    temp_EEG = pop_subcomp(temp_EEG,artcomps,0,0); % remove ICs from EEG data
                    disp(['removed ICs: ',num2str(artcomps)]);
                    
                    EEG.data(goodchanns,:) = temp_EEG.data;
                    
                    EEG.preprocessed.residualArtifactRemoval.MARA.check = 'yes';
                    EEG.preprocessed.residualArtifactRemoval.MARA.artcomps = artcomps;
                    EEG.preprocessed.residualArtifactRemoval.MARA.MARAinfo = MARAinfo;
                    EEG.preprocessed.residualArtifactRemoval.MARA.MARA_thre = MARA_thre;
                else
                    warning('Skip residual artifact removal: channel location is uncorrect or empty');
                    disp('Skip residual artifact removal: unused...');
                    EEG.preprocessed.residualArtifactRemoval.check = 'no';
                end
                
            case 2  % ICA based ADJUST
                if chanlocsflag == 1
                    disp('Artifact Removal: removing residual artifact using ICA based ADJUST...');
                    try badchanns = results_QA.badChannels;catch; badchanns = [];end;
                    goodchanns = setdiff(selechanns,badchanns);
                    % goodchanns = selechanns; % ICA on all selected channels
                    temp_EEG = EEG;
                    % unionchanns = union(selechanns,EOGchanns);
                    temp_EEG.data = EEG.data(goodchanns,:);
                    temp_EEG.chanlocs = EEG.chanlocs(1,goodchanns);
                    temp_EEG.nbchan = length(goodchanns);
                    if  isfield(temp_EEG,'icaweights') && isempty(temp_EEG.icaweights)
                        ICs = 'default';      % number of ICA components to compute (default -> chans or 'pca' arg)
                        Ntrain = 0;           % perform tanh() "extended-ICA" with sign estimation N training blocks.
                        PCs = size(temp_EEG.data,1)-1;      % decompose a principal component (default = 0 -> off).
                        stop = 1e-6;      % stop training when weight-change < this
                        MaxSteps = 512;   % max number of ICA training steps
                        sphering = 'on';  % ['on'/'off'] flag sphering of data
                        temp_EEG = wb_runICA(temp_EEG,'all',ICs,Ntrain,PCs,stop,MaxSteps,sphering);
                        
                        EEG.preprocessed.residualArtifactRemoval.ADJUST.ICs = ICs;
                        EEG.preprocessed.residualArtifactRemoval.ADJUST.ICANtrain = Ntrain;
                        EEG.preprocessed.residualArtifactRemoval.ADJUST.ICAstop = stop;
                        EEG.preprocessed.residualArtifactRemoval.ADJUST.ICAMaxSteps = MaxSteps;
                        EEG.preprocessed.residualArtifactRemoval.ADJUST.ICAsphering = sphering;
                    else
                        disp('Skip ICA decomposition: ICA has been conducted on EEG data...');
                    end
                    
                    try [artcomps, horizcomps, vertcomps, blinkcomps, disccomps,soglia_DV, diff_var, soglia_K, ...
                            ~, meanK, soglia_SED, ~, SED, soglia_SAD, ~, SAD, ...
                            soglia_GDSF, ~, GDSF, soglia_V, ~, nuovaV, ~, ~] = ADJUST (temp_EEG,'adjust.txt');  % ADJUST
                    catch
                        artcomps = [];
                        horizcomps = [];
                        vertcomps = [];
                        blinkcomps = [];
                        disccomps = [];
                        soglia_DV = [];
                        diff_var = [];
                        soglia_K = [];
                        % med2_K = [];
                        meanK = [];
                        soglia_SED = [];
                        % med2_SED = [];
                        SED = [];
                        soglia_SAD = [];
                        % med2_SAD = [];
                        SAD = [];
                        soglia_GDSF = [];
                        % med2_GDSF = [];
                        GDSF = [];
                        soglia_V = [];
                        % med2_V = [];
                        nuovaV = [];
                    end
                    temp_EEG = pop_subcomp(temp_EEG,artcomps,0, 0); % remove ICs from EEG data
                    disp(['removed ICs: ',num2str(artcomps)]);
                    
                    EEG.data(goodchanns,:) = temp_EEG.data;
                    %   artcomps        - List of artifacted ICs
                    %   horizcomps      - List of HEM ICs
                    %   vertcomps       - List of VEM ICs
                    %   blinkcomps      - List of EB ICs
                    %   disccomps       - List of GD ICs
                    %   soglia_DV  - SVD threshold
                    %   diff_var   - SVD feature values
                    %   soglia_K   - TK threshold
                    %   meanK      - TK feature values
                    %   soglia_SED - SED threshold
                    %   SED        - SED feature values
                    %   soglia_SAD - SAD threshold
                    %   SAD        - SAD feature values
                    %   soglia_GDSF- GDSF threshold
                    %   GDSF       - GDSF feature values
                    %   soglia_V   - MEV threshold
                    %   nuovaV     - MEV feature values
                    
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.check = 'yes';
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.artcomps = artcomps;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.horizcomps = horizcomps;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.vertcomps = vertcomps;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.blinkcomps = blinkcomps;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.disccomps = disccomps;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.soglia_DV = soglia_DV;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.diff_var = diff_var;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.soglia_K = soglia_K;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.meanK = meanK;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.soglia_SED = soglia_SED;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.SED = SED;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.soglia_SAD = soglia_SAD;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.SAD = SAD;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.soglia_GDSF = soglia_GDSF;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.GDSF = GDSF;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.soglia_V = soglia_V;
                    EEG.preprocessed.residualArtifactRemoval.ADJUST.nuovaV = nuovaV;
                else
                    warning('Skip residual artifact removal: channel location is uncorrect or empty');
                    disp('Skip residual artifact removal: unused...');
                    EEG.preprocessed.residualArtifactRemoval.check = 'no';
                end
            case 3  % Robust PCA
                disp('Artifact Removal: removing residual artifact using rPCA...');
                try badchanns = results_QA.badChannels;catch; badchanns = [];end;
                goodchanns = setdiff(selechanns,badchanns);
                % goodchanns = selechanns; % ICA on all selected channels
                Nt = size(EEG.data,2);
                lambda = 1 / sqrt(Nt); % weight on sparse error term in the cost function
                tol = 1e-7;     % tolerance for stopping criterion.
                maxIter = 1000; % maximum number of iterations
                [A_hat, ~,~] = inexact_alm_rpca(EEG.data(goodchanns,:)', lambda, tol, maxIter);
                EEG.data(goodchanns,:) = A_hat';
                
                EEG.preprocessed.residualArtifactRemoval.rPCA.check = 'yes';
                EEG.preprocessed.residualArtifactRemoval.rPCA.lambda = lambda;
                EEG.preprocessed.residualArtifactRemoval.rPCA.tol = tol;
                EEG.preprocessed.residualArtifactRemoval.rPCA.maxIter = maxIter;
                
            case 4  % ASR
                disp('Artifact Removal: removing residual artifact using ASR...'); % Artifact Subspace Reconstruction 
                try badchanns = results_QA.badChannels;catch; badchanns = [];end;
                goodchanns = setdiff(selechanns,badchanns);
                % goodchanns = selechanns; % ICA on all selected channels
                
                burst_crit = 5; % Standard deviation cutoff for removal of bursts (via ASR).A quite conservative value is 5.
                burst_crit_refmaxbadchns = 0.075; % This number is the maximum tolerated (0.05-0.3)
                %                             fraction of "bad" channels within a given time window of the recording
                %                             that is considered acceptable for use as calibration data.
                burst_crit_reftolerances = [-3.5,5.5]; %  These are the power tolerances outside of which a channel in a
                %                         given time window is considered "bad", in standard deviations relative to
                %                         a robust EEG power distribution (lower and upper bound). Together with the
                %                         previous parameter this determines how ASR calibration data is be
                %                         extracted from a recording. Can also be specified as 'off' to achieve the
                %                         same effect as in the previous parameter. Default: [-3.5 5.5].
                
                temp_EEG = clean_asr(EEG,burst_crit,[],[],[],burst_crit_refmaxbadchns,burst_crit_reftolerances,[]); 
                EEG.data(goodchanns,:) = temp_EEG.data(goodchanns,:);
                
                EEG.preprocessed.residualArtifactRemoval.ASR.check = 'yes';
                EEG.preprocessed.residualArtifactRemoval.ASR.burst_crit = burst_crit;
                EEG.preprocessed.residualArtifactRemoval.ASR.burst_crit_refmaxbadchns = burst_crit_refmaxbadchns;
                EEG.preprocessed.residualArtifactRemoval.ASR.burst_crit_reftolerances = burst_crit_reftolerances;
        end
    else
        EEG.preprocessed.residualArtifactRemoval.check = 'no';
        disp('Skip residual artifact removal: unused...');
    end
    
    % ---------------------------------------------------------------------
    % [5] bad channel interpolation and re-referencing
    if chanlocsflag == 1 && badChannelInterploateFlag >= 0 && badChannelInterploateFlag < 3
        switch badChannelInterploateFlag           
            case 0  % unused
                disp('Skip interpolating the badChannels: unused');
            case 1 % REST interpolation
                try badchanns = results_QA.badChannels;catch; badchanns = [];end;
                disp(['Detected bad channels are: ',num2str(badchanns)]);
                disp('REST interpolation...');
                EEG = wb_restInterpolate(EEG,badchanns,chanlocsXYZ,selechanns); % interpolate
                EEG.preprocessed.RESITinterpolation.check = 'yes';
                EEG.preprocessed.RESITinterpolation.badchanns = badchanns;
                EEG.preprocessed.RESITinterpolation.comments = 'REST based on 3-concentric spheres head model';
            case 2 % Interpolate the badChannels rows of EEG.data using spherical splines
                try badchanns = results_QA.badChannels;catch; badchanns = [];end;
                disp(['Detected bad channels are: ',num2str(badchanns)]);
                disp('Spherical spline interpolation and re-referencing to REST...');
                if ~isempty(badchanns)
                    temp_EEG = EEG;
                    % unionchanns = union(selechanns,EOGchanns);
                    temp_EEG.data = EEG.data(selechanns,:);
                    temp_EEG.chanlocs = EEG.chanlocs(1,selechanns);
                    temp_EEG.nbchan = length(selechanns);
                    temp_index1 = zeros(length(badchanns),1);
                    for i = 1:length(badchanns)
                        temp_index1(i) = find(selechanns==badchanns(i));
                    end
                    temp_EEG = interpolateChannels(temp_EEG, temp_index1); % interpolate 
                    EEG.data(selechanns,:) = temp_EEG.data;
                end
                
                disp('Re-referencing to REST...');
                EEG = wb_rest_rerefer_concentricspheres(EEG,chanlocsXYZ,selechanns); % re-referencing to REST
                % EEG.data(selechanns,:) = bsxfun ( @minus, EEG.data(selechanns,:), mean(EEG.data(selechanns,:))); % re-referencing to average reference;
                % EEG.ref = 'AVG';
                
                EEG.preprocessed.SphericalSplinesInterpolation.comments = 're-referencing to REST based on 3-concentric spheres head model';
                EEG.preprocessed.SphericalSplinesInterpolation.check = 'yes';
                EEG.preprocessed.SphericalSplinesInterpolation.badchanns = badchanns;
        end
    else
        disp('Skip Interpolating the badChannels: channel location is uncorrect or empty');
        if chanlocsflag == 1
            disp('Re-referencing to REST...');
            EEG = wb_rest_rerefer_concentricspheres(EEG,chanlocsXYZ,selechanns); % re-referencing to REST
            EEG.preprocessed.Interpolation.comments = 're-referencing to REST based on 3-concentric spheres head model';
        end
        EEG.preprocessed.Interpolation.check = 'no';
    end
    % ---------------------------------------------------------------------
    % [6] QA after removing artifacts
    disp('------------------------------------------------------------');
    disp('Quality Assessment of EEG data after artifact removing...');
    % -------------------------------------------------------------------------
    % Quality Assessmenting of EEG data after artifact removing
    [results_QA_clean] = wb_EEG_QA(EEG,WindowSeconds,HighPassband,selechanns,badWindowThreshold,...
        robustDeviationThreshold,PowerFrequency,FrequencyNoiseThreshold,flagNotchFilter,correlationThreshold,...
        ransacCorrelationThreshold,ransacChannelFraction,ransacSampleSize);
    EEG.preprocessed.QA.check = 'yes';
    EEG.preprocessed.QA.comments = 'QA after artifact removal';
    EEG.preprocessed.QA.results = results_QA_clean;
    
    disp(['Overall Data Quality (preprocessed data): ',num2str(results_QA_clean.ODQ)]);
    disp(['Data Quality Rating (preprocessed data): ',results_QA_clean.DataQualityRating]);
    % ---------------------------------------------------------------------
    % [7] Marking bad block with unusually high or low amplitude using zscored STD
    disp('-------------------------------------------------------------');
    disp('Marking residual bad block with unusually high or low amplitude using zscored STD...');
    try
        [EEG,~,data_Z] = wb_mark_EEGblock(EEG,0,0,robustDeviationThreshold,WindowSeconds,selechanns);
        EEG.preprocessed.MarkBadBlock.check = 'yes';
        EEG.preprocessed.MarkBadBlock.comments = 'Marking residual bad block after artifact removal';
        EEG.preprocessed.MarkBadBlock.zscoredGFP = data_Z;
        EEG.preprocessed.MarkBadBlock.STDthreshold = robustDeviationThreshold;
    catch
       EEG.preprocessed.MarkBadBlock.check = 'no';
       disp('Skip marking residual bad block: unused...');
    end
    % ---------------------------------------------------------------------
    if keepUnselectChannsFlag == 0 % do not keep unselected channels if = 0
       EEG.data = EEG.data(selechanns,:);
       try EEG.chanlocs = EEG.chanlocs(1,selechanns);catch;end;
       EEG.nbchan = length(selechanns);
    end
    % ---------------------------------------------------------------------
    % saving results
    if isempty(EEG)
        disp('EEG is empty, subject skipped');
    else
        try
            [~, Origdatfile1, ~] = fileparts(EEG.datfile);
            EEG.datfile = [Origdatfile1,'-prepro'];
        catch
            disp('EEG.datfile may be empty!');
        end
        
        try
            [~, Origdatfile2, ~] = fileparts(EEG.filename);
            EEG.filename = [Origdatfile2,'_prepro'];
        catch
            disp('EEG.filename may be empty!');
        end
        
        try
            [~, Origdatfile3, ~] = fileparts(EEG.setname);
            EEG.setname = [Origdatfile3,'-prepro'];
        catch
            disp('EEG.setname may be empty!');
        end
        
        EEG.comments = 'EEG preprocessed';
        % ---------------------
        EEGprepro = EEG;
    end
else
    try
        disp(['Overall Data Quality of EEG data is ', num2str(results_QA.ODQ),' < threshold of ODQ: ',num2str(thre_ODQ)]);
    catch
        disp('Overall Data Quality of EEG data < threshold of ODQ');
    end;
    warning('Quality of EEG data is too low, skip artifact removing...');
    EEGprepro = [];
end