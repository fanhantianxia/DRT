function MRI = wb_loadMRI(MRIfile,unzip_dir)
% loading fMRI image data. This function will scan all files in folders with its
% subfolders (3 layer max), and load one sample data only (last one).
% Input:
%      MRIfile: Folder OR zip file of one image OR one nii file (only for one subject/sample).
%               If it is a zip file, MRI data will be unzipped to the data path.
%      unzip_dir: temp folder to unzip MRI data. It is required if
%               'MRIfile' is zip file.
% Output:
%      MRI: MRI data
% --------------
% The following MRI file formats are supported
%   CTF - VSM MedTech (*.svl, *.mri version 4 and 5)
%   NIFTi (*.nii) and zipped NIFTi (*.nii.gz)
%   Analyze (*.img, *.hdr)
%   DICOM (*.dcm, *.ima)
%   AFNI (*.head, *.brik)
%   FreeSurfer (*.mgz, *.mgh)
%   MINC (*.mnc)
%   Neuromag - Elekta (*.fif)
%   ANT - Advanced Neuro Technology (*.mri)
%   Yokogawa (*.mrk, incomplete)
%   more details can be seen in ft_read_mri.m from FieldTrip
% -------------------------------------------------------------------------
% Written by Li Dong (Lidong@uestc.edu.cn) $ 2018.1.8
% Revised:
%      can load MATLAB .dat file. $ 20180802 Li Dong
%      can load BIDS-MRI.         $ 20200217 Li Dong
% -------------------------------------------------------------------------
if nargin == 1
    unzip_dir = [];
end

MRI = []; % MRILAB data structure

% scan foulders and load MRI data

if exist(MRIfile,'dir') == 7 % is folder and exist
    % -----------------------------
    listing1 = dir(MRIfile);
    % ----
    for i = 3:length(listing1)
        cur_file = listing1(i).name;
        if exist(fullfile(MRIfile,cur_file),'dir') == 7 % is folder and exist
            listing2 = dir(fullfile(MRIfile,cur_file));
            MRIfile2 = fullfile(MRIfile,cur_file);
            for j = 3:length(listing2)
                cur_file2 = listing2(j).name;
                if exist(fullfile(MRIfile,cur_file,cur_file2),'dir') == 7 % is folder and exist
                    listing3 = dir(fullfile(MRIfile,cur_file,cur_file2));
                    MRIfile3 = fullfile(MRIfile,cur_file,cur_file2);
                    for k = 3:length(listing3)
                        cur_file3 = listing3(k).name;
                        try MRI = readMRI(MRIfile3,cur_file3);catch;end
                    end
                else
                    try MRI = readMRI(MRIfile2,cur_file2);catch;end
                end
            end
        else
            try MRI = readMRI(MRIfile,cur_file);catch;end
        end
    end
    % ------------------------------
else   % is zip file
    if ~isempty(unzip_dir) % if unzip direction is not empty
        [~,name1,ext2] = fileparts(MRIfile); % get filename extension
        if isequal(ext2,'.zip') % is .zip?
            unzip_temp_path = fullfile(unzip_dir,filesep,['unzip_temp_',name1]);
            unzip(MRIfile,unzip_temp_path); % unzip data
            % -----------------------------
            listing1 = dir(unzip_temp_path);
            % ----
            for i = 3:length(listing1)
                cur_file = listing1(i).name;
                if exist(fullfile(unzip_temp_path,cur_file),'dir') == 7 % is folder and exist
                    listing2 = dir(fullfile(unzip_temp_path,cur_file));
                    MRIfile2 = fullfile(unzip_temp_path,cur_file);
                    for j = 3:length(listing2)
                        cur_file2 = listing2(j).name;
                        if exist(fullfile(unzip_temp_path,cur_file,cur_file2),'dir') == 7 % is folder and exist
                            listing3 = dir(fullfile(unzip_temp_path,cur_file,cur_file2));
                            MRIfile3 = fullfile(unzip_temp_path,cur_file,cur_file2);
                            for k = 3:length(listing3)
                                cur_file3 = listing3(k).name;
                                try MRI = readMRI(MRIfile3,cur_file3);catch;end
                            end
                        else
                            try MRI = readMRI(MRIfile2,cur_file2);catch;end
                        end
                    end
                else
                    try MRI = readMRI(unzip_temp_path,cur_file);catch;end
                end
            end
            % ------------------------------
            % remove temp zip folder
            try rmdir(unzip_temp_path,'s'); catch;end;
        end
    else
        error(' Input ''unzip_dir'' is required, if ''MRIfile'' is zip file.');
    end
end

if isempty(MRI)
    warning('failed to load data, please check the MRI file...');
end
% -------------------------------------------------------------------------
% subfunctions

    function MRI = readMRI(MRIfile5,cur_file5)
        % read MRI files using FileTrip function
        % Input:
        %    MRIfile5: input path
        %    cur_file5: MRI file name
        % output:
        %    MRI: MRI data
        
        % disp(cur_file5);  % disp file scaned, for test only
        try
            MRI = ft_read_mri(fullfile(MRIfile5,cur_file5)); % reading in the MRI data;
        catch
        end;
    end
end