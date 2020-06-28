function EEG = wb_runICA(EEG,selechanns,ICs,Ntrain,PCs,stop,MaxSteps,sphering)
% run ICA on EEG data based on EEGLAB function 
% runica() - Perform Independent Component Analysis (ICA) decomposition
%            of input data using the logistic infomax ICA algorithm of 
%            Bell & Sejnowski (1995) with the natural gradient feature 
%            of Amari, Cichocki & Yang, or optionally the extended-ICA 
%            algorithm of Lee, Girolami & Sejnowski, with optional PCA 
%            dimension reduction. Annealing based on weight changes is 
%            used to automate the separation process. 
% Input:
%    EEG: EEG structure imported using EEGLAB. EEG.data should be channels
%         X time points OR channels X time points X epoches.
%    selechanns: number with indices of the selected channels
%                   (e.g. [1:4,7:30] or 'all').Default is 'all';
%    ICs: number of ICA components to compute (default -> chans or 'pca' arg) 
%         using rectangular ICA decomposition. This parameter may return 
%         strange results. This is because the weight matrix is rectangular 
%         instead of being square. Do not use except to try to fix the problem. 
%    Ntrain : perform tanh() "extended-ICA" with sign estimation 
%         N training blocks. If N > 0, automatically estimate the 
%         number of sub-Gaussian sources. If N < 0, fix number of 
%         sub-Gaussian comps to -N [faster than N>0] (default=0 -> off)
%    PCs: decompose a principal component (default = 0 -> off)
%          subspace of the data. Value is the number of PCs to retain.
%    stop: stop training when weight-change < this (default is 1e-6
%               if less than 33 channel and 1E-7 otherwise)
%    maxsteps: max number of ICA training steps (default is 512)
%    sphering: ['on'/'off'] flag sphering of data. default is on.
% Output:
%    EEG: The input EEGLAB dataset with new fields icaweights, icasphere 
%             and icachansind (channel indices). 
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC,Lidong@uestc.edu.cn)
% $ 2018.5.11
% -------------------------------------------------------------------------
if nargin < 1
    error ('One input is reqiured at least!!!!!');
elseif nargin == 1
    selechanns = 'all';
    ICs = 'default';
    Ntrain = 0;
    PCs = 'default';
    stop = 1e-6;
    MaxSteps = 512;
    sphering = 'on';
elseif nargin == 2
    ICs = 'default';
    Ntrain = 0;
    PCs = 'default';
    stop = 1e-6;
    MaxSteps = 512;
    sphering = 'on';
elseif nargin == 3
    Ntrain = 0;
    PCs = 'default';
    stop = 1e-6;
    MaxSteps = 512;
    sphering = 'on';
elseif nargin == 4
    PCs = 'default';
    stop = 1e-6;
    MaxSteps = 512;
    sphering = 'on';
elseif nargin == 5
    stop = 1e-6;
    MaxSteps = 512;
    sphering = 'on';
elseif nargin == 6
    MaxSteps = 512;
    sphering = 'on';
elseif nargin == 7
    sphering = 'on';
end

% check inputs:
if isempty(ICs)
    ICs = 'default';
end

if isempty(Ntrain)
    Ntrain = 0;
end

if isempty(PCs)
    PCs = 'default';
end

if isempty(stop)
    stop = 1e-6;
end

if isempty(MaxSteps)
    MaxSteps = 512;
end

if isempty(sphering)
    sphering = 'on';
end

try
    if isempty(EEG.data)
        error('EEG.data is empty!!!!!');
    end
catch
    error('EEG.data is not exist!!!!');
end

% channs
if isequal(selechanns,'all')
    selechanns = 1:size(EEG.data,1);
end
N_channs = length(selechanns); % No. of selected channs
disp(['No. of selected channels: ',num2str(N_channs)]);

% No. of timepoints
Nt = size(EEG.data,2);
disp(['EEG data or epoch legnth:',num2str(Nt)]);

% data

DIM = size(EEG.data);

if length(DIM) == 1
    error('EEG.data is not correct')
elseif length(DIM) == 2
    tempdata = EEG.data(selechanns,:);
elseif length(DIM) == 3
    Ntri = DIM(3);
    tempdata = reshape(EEG.data(selechanns,:,:), N_channs, Nt*Ntri);
    disp('EEG data has been epoched')
    disp(['No. of trials or epoches:',num2str(Ntri)]);
end

tempdata = tempdata - repmat(mean(tempdata,2), [1 size(tempdata,2)]); % zero mean 
tmprank = getrank(double(tempdata(:,1:min(3000, size(tempdata,2)))));
if tmprank ~= size(tempdata,1)
    disp(['It has detected that the rank of your data matrix', ...
          'is lower the number of input data channels. This might',...
          'be because you are including a reference channel or',...
          'because you are running a second ICA decomposition.']);
    disp('reset No. of PCs -> min(PCs, rank of data matrix)');
    disp('reset No. of ICs -> min(ICs, PCs)');
    if isnumeric(PCs)
        PCs = min(tmprank,PCs);
    else
        PCs = tmprank;
    end
    if isnumeric(ICs)
        ICs = min(ICs,PCs);
    end
end
% -------------------------------------------------------------------------
% run ICA
if ~isempty(tempdata)
    if isequal(ICs,'Default') || isequal(ICs,'default')
        if isequal(PCs,'Default') || isequal(PCs,'default')
            [icaweights,icasphere,~,~,~,~,activations] = runica(tempdata,'extended',Ntrain,...
                'stop',stop,'maxsteps',MaxSteps,'sphering',sphering);
        else
            [icaweights,icasphere,~,~,~,~,activations] = runica(tempdata,'pca',PCs,'extended',Ntrain,...
                'stop',stop,'maxsteps',MaxSteps,'sphering',sphering);
        end
    else
        if isequal(PCs,'Default')|| isequal(PCs,'default')
            [icaweights,icasphere,~,~,~,~,activations] = runica(tempdata,'ncomps',ICs,'extended',Ntrain,...
                'stop',stop,'maxsteps',MaxSteps,'sphering',sphering);
        else
            [icaweights,icasphere,~,~,~,~,activations] = runica(tempdata,'ncomps',ICs,'pca',PCs,'extended',Ntrain,...
                'stop',stop,'maxsteps',MaxSteps,'sphering',sphering);
        end
    end
else
    error('data is empty')
end

EEG.icasphere = icasphere;
EEG.icaweights = icaweights;
EEG.icachansind = selechanns;
if length(DIM) == 3
    EEG.activations = activations;
end;
EEG.icawinv = pinv(icaweights * icasphere); % a priori same result as inv

EEG.ICAPara.ICs = ICs;
EEG.ICAPara.Ntrain = Ntrain;
EEG.ICAPara.PCs = PCs;
EEG.ICAPara.stop = stop;
EEG.ICAPara.MaxSteps = MaxSteps;
EEG.ICAPara.sphering = sphering;
end
    
function tmprank2 = getrank(tmpdata)

tmprank = rank(tmpdata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here: alternate computation of the rank by Sven Hoffman
%tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
covarianceMatrix = cov(tmpdata', 1);
[~, D] = eig (covarianceMatrix);
rankTolerance = 1e-7;
tmprank2=sum (diag (D) > rankTolerance);
if tmprank ~= tmprank2
    fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
    tmprank2 = max(tmprank, tmprank2);
end;
end