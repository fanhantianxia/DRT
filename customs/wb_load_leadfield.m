function G = wb_load_leadfield(leadfield_file)
% load leadfield for EEG data
% Input:
%       leadfield file: A file such as '*\leadfield.dat' or '*\leadfield.mat'. 
%                The leadfield can be a matrix (a matlab *.dat file,
%                output of leadfield.exe,sources X channels) which is 
%                calculated by using the forward theory, based on the 
%                electrode montage, head model and equivalent source model. 
%                It can also be the output of ft_prepare_leadfield.m 
%                (e.g. lf.leadfield,dipoles contain x,y,z-orientations) 
%                based on real head modal (FEM modal) using FieldTrip.
% Output:
%       G: leadfield matrix,sources X channels.
% ------------
% How to calculate leadfield see:
% Li Dong*, Fali Li, Qiang Liu, Xin Wen, Yongxiu Lai, Peng Xu and Dezhong Yao*. 
%              MATLAB Toolboxes for Reference Electrode Standardization Technique (REST) 
%              of Scalp EEG. Frontiers in Neuroscience,  2017:11(601).
% 
% http://www.neuro.uestc.edu.cn/rest/Down.html
% OR 
% Use real head model to calculate leadfield, please see FieldTrip toolbox: 
% http://www.fieldtriptoolbox.org/development/project/example_fem
% http://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_fem?s[]=fem
% 
% -----------
% Written by Li Dong (Lidong@uestc.edu.cn)
% $ 2017.11.27
% -------------------------------------------------------------------------
if ~isempty(leadfield_file)
    leadfield = load(leadfield_file);
else
    error('leadfield is empty.');
end
% get the leafield matrix
if isnumeric(leadfield)
    G = leadfield;
elseif isfield(leadfield,'lf') % ONLY for FieldTrip
    if isstruct(leadfield.lf)
        try
            Npos = size(leadfield.lf.pos,1);
            m = 1;
            for i = 1:Npos
                if ~isempty(leadfield.lf.leadfield{1,i})
                    lf_X(:,m) = leadfield.lf.leadfield{1,i}(:,1); % X orientation of the dipole.
                    lf_Y(:,m) = leadfield.lf.leadfield{1,i}(:,2); % Y orientation of the dipole.
                    lf_Z(:,m) = leadfield.lf.leadfield{1,i}(:,3); % Z orientation of the dipole.
                    m = m + 1;
                end
            end
            G = [lf_X,lf_Y,lf_Z]';
            % the leadfield matrix (sources*3 X chans), which
            % contains the potential or field distributions on all
            % sensors for the x,y,z-orientations of the dipole.
        catch
            error('leadfiled is not calculated by ''ft_prepare_leadfield.m'' (dipoles contain x,y,z-orientations)?');
        end
    end   
elseif isfield(leadfield,'leadfield') % ONLY for FieldTrip
    if iscell(leadfield.leadfield)
        try
            Npos = length(leadfield.leadfield);
            m = 1;
            for i = 1:Npos
                if ~isempty(leadfield.leadfield{1,i})
                    lf_X(:,m) = leadfield.leadfield{1,i}(:,1); % X orientation of the dipole.
                    lf_Y(:,m) = leadfield.leadfield{1,i}(:,2); % Y orientation of the dipole.
                    lf_Z(:,m) = leadfield.leadfield{1,i}(:,3); % Z orientation of the dipole.
                    m = m + 1;
                end
            end
            G = [lf_X,lf_Y,lf_Z]';
            % the leadfield matrix (sources*3 X chans), which
            % contains the potential or field distributions on all
            % sensors for the x,y,z-orientations of the dipole.
        catch
            error('leadfiled is not calculated by ''ft_prepare_leadfield.m'' (dipoles contain x,y,z-orientations)?');
        end
    end
end