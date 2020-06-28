function [lf,headmodel,elec_aligned] = wb_EEG_standardLeadfield_BEM(elecfile,bemfile,elecDirecFlag,gridresolution)

% This funciton is used to generate conduction model of the head based
% on boundary element method (BEM) using standard MRI T1 image, and computes
% the forward model for many dipole locations on a 2D brain mesh or regular
% 3D grid and stores it for efficient inverse modelling using FieldTrip for EEG.
% The coordinates of head model is standard MNI space, and the electrodes
% will be aligned later to the existing standard head model.

% A standard headmodel using 'dipoli' method based on BEM were used in this
% function. The headmodel contains a standard Boundary Element Method volume
% conduction model of the head that can be used for EEG forward and inverse
% computations. The geometry is based on the ¡°colin27¡± template that is
% described further down. The BEM model is expressed in MNI coordinates in mm.
% A very similar BEM volume conduction model (based on the same template data)
% is described and validated by Fuchs et al. in Clin Neurophysiol. 2002 May;113(5):702-12.
% more details see : http://www.fieldtriptoolbox.org/template/headmodel/

% The ¡°colin27¡± anatomical MRI and its relation to the TT and MNI template
% atlas is described in detail on http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach
% The original construction of the averaged MRI is detailed in
% [[http://www.ncbi.nlm.nih.gov/pubmed/9530404| Holmes CJ, Hoge R, Collins L,
% Woods R, Toga AW, Evans AC. Enhancement of MR images using registration for
% signal averaging. J Comput Assist Tomogr. 1998 Mar-Apr;22(2):324-33.]]
% -------------------------------------------------------------------------
% Input:
%    elecfile£º
%       it can be a file of electrode coordinates or a structure array
%       of channel locations (i.e. EEG.chanlocs) from EEG loaded by EEGALB.
%       Files from the following acquisition
%       systems and analysis platforms file formats are supported
%       (more details see readlocs in EEGLAB).
%   bemfile:
%       a mat file contain a standard headmodel using 'dipoli' method based on
%       BEM.
%    elecDirecFlag:
%        0: XYZ coordinates is the electorde array with their Cartesian x
%        (the left ear is defined as -x axis),y (the nasion is the +y axis),
%        z coordinates in three columns.
%        1: XYZ coordinates is the electorde array with their Cartesian x
%        (the nasion is the +x axis),y (the left ear is the +y axis),
%        z coordinates in three columns. Default is 1.
%    gridresolution£º
%        The grid resolution of dipoles (sources) inside the brain.
%        If it is empoty or <=0, the default dipoles are vertices which
%        are little smaller than brain, and the orientations
%        of dipoles are their normal vector directions, i.e.
%        the normals of the brain mesh.If it >0 (unit is mm), the
%        dipoles are distributed on regular 3D grid inside
%        brain mesh.The orientations of dipoles are X, Y and Z
%        oritentations, i.e. there are X, Y and Z oritented dipoles.
%        Default is empty.
% Outputs:
%    lf: leadfield results.
%        lf.leadfieldMatrix: leadfield matrix saved as (channels X sources)
%        lf.label: channel labels.
%        lf.dim: dimension of dipole (source) grid.
%        lf.unit: unit of head model coordinates.
%        lf.pos: positions of dipoles.
%        lf.inside: Boolean value of whether the lf.pos inside the brain.
%        lf.cfg: configuration of leadfield calculation
%        lf.leadfield: leadfield saved as cell.
%    headmodel: headmodel used in the function
%        headmodel.bnd: mesh of scalp, skull and brain.
%        headmodel.cond: conductivity of tissues, order is [scalp,skull and
%               brain].
%        headmodel.type: BEM method used.
%        headmodel.unit:unit of head model coordinates.
%    elec_aligned: aligned eletrode coordinates
%       elec_aligned.elecpos: aligned eletrode postions.
%       elec_aligned.label: channel labels.
%       elec_aligned.cfg: configuration of alignment
% Electrode location file formats:
% The file extension determines its type.
%
%   '.loc' or '.locs' or '.eloc':
%               polar coordinates. Notes: angles in degrees:
%               right ear is 90; left ear -90; head disk radius is 0.5.
%               Fields:   N    angle  radius    label
%               Sample:   1    -18    .511       Fp1
%                         2     18    .511       Fp2
%                         3    -90    .256       C3
%                         4     90    .256       C4
%                           ...
%               Note: In previous releases, channel labels had to contain exactly
%               four characters (spaces replaced by '.'). This format still works,
%               though dots are no longer required.
%   '.sph':
%               Matlab spherical coordinates. Notes: theta is the azimuthal/horizontal angle
%               in deg.: 0 is toward nose, 90 rotated to left ear. Following this, performs
%               the elevation (phi). Angles in degrees.
%               Fields:   N    theta    phi    label
%               Sample:   1      18     -2      Fp1
%                         2     -18     -2      Fp2
%                         3      90     44      C3
%                         4     -90     44      C4
%                           ...
%   '.elc':
%               Cartesian 3-D electrode coordinates scanned using the EETrak software.
%               See readeetraklocs().
%   '.elp':
%               Polhemus-.'elp' Cartesian coordinates. By default, an .elp extension is read
%               as PolhemusX-elp in which 'X' on the Polhemus sensor is pointed toward the
%               subject. Polhemus files are not in columnar format (see readelp()).
%   '.elp':
%               BESA-'.elp' spherical coordinates: Need to specify 'filetype','besa'.
%               The elevation angle (phi) is measured from the vertical axis. Positive
%               rotation is toward right ear. Next, perform azimuthal/horizontal rotation
%               (theta): 0 is toward right ear; 90 is toward nose, -90 toward occiput.
%               Angles are in degrees.  If labels are absent or weights are given in
%               a last column, readlocs() adjusts for this. Default labels are E1, E2, ...
%               Fields:   label      phi  theta
%               Sample:   Fp1        -92   -72
%                         Fp2         92    72
%                         C3         -46    0
%                         C4          46    0
%                           ...
%   '.xyz':
%               Matlab/EEGLAB Cartesian coordinates. Here. x is towards the nose,
%               y is towards the left ear, and z towards the vertex.
%               Fields:   channum   x           y         z     label
%               Sample:   1       .950        .308     -.035     Fp1
%                         2       .950       -.308     -.035     Fp2
%                         3        0           .719      .695    C3
%                         4        0          -.719      .695    C4
%                           ...
%   '.asc', '.dat':
%               Neuroscan-.'asc' or '.dat' Cartesian polar coordinates text file.
%   '.sfp':
%               BESA/EGI-xyz Cartesian coordinates. Notes: For EGI, x is toward right ear,
%               y is toward the nose, z is toward the vertex. EEGLAB converts EGI
%               Cartesian coordinates to Matlab/EEGLAB xyz coordinates.
%               Fields:   label   x           y          z
%               Sample:   Fp1    -.308        .950      -.035
%                         Fp2     .308        .950      -.035
%                         C3     -.719        0          .695
%                         C4      .719        0          .695
%                           ...
%   '.ced':
%               ASCII file saved by pop_chanedit(). Contains multiple MATLAB/EEGLAB formats.
%               Cartesian coordinates are as in the 'xyz' format (above).
%               Fields:   channum  label  theta  radius   x      y      z    sph_theta   sph_phi  ...
%               Sample:   1        Fp1     -18    .511   .950   .308  -.035   18         -2       ...
%                         2        Fp2      18    .511   .950  -.308  -.035  -18         -2       ...
%                         3        C3      -90    .256   0      .719   .695   90         44       ...
%                         4        C4       90    .256   0     -.719   .695  -90         44       ...
%                           ...
%               The last columns of the file may contain any other defined fields (gain,
%               calib, type, custom).
% -------------------------------------------------------------------------
% % usgae of funtion;
% clear all;clc;
% elecfile = 'G:\work_UESTC\Work\WeBrain_Platform\WeBrain_pipelines\EEG-forward-simulation\Curry7_62Channels.loc';
%                     % file of electrode coordinates
%                     % Files from the following acquisition systems and
%                     % analysis platforms file formats are supported (more
%                     % details see readlocs in EEGLAB).
%
% elecDirecFlag = 1;  % 0: XYZ coordinates is the electorde array with their Cartesian x
%                     %    (the left ear is defined as -x axis),y (the nasion is the +y axis),
%                     %    z coordinates in three columns.
%                     % 1: XYZ coordinates is the electorde array with their Cartesian x
%                     %    (the nasion is the +x axis),y (the left ear is the +y axis),
%                     %    z coordinates in three columns. Default is 1.
% gridresolution = 12;   % The grid resolution of dipoles (sources) inside the brain.
%                        % If it is empoty or <=0, the default dipoles are vertices which
%                        % are little smaller than brain, and the orientations
%                        % of dipoles are their normal vector directions, i.e.
%                        % the normals of the brain mesh.If it >0 (unit is mm), the
%                        % dipoles are distributed on regular 3D grid inside
%                        % brain mesh.The orientations of dipoles are X, Y and Z
%                        % oritentations, i.e. there are X, Y and Z oritented dipoles.
%                        % Default is empty.
%
% [lf,headmodel,elec_aligned] = wb_EEG_standardLeadfield_BEM(elecfile,elecDirecFlag,gridresolution);
% -------------------------------------------------------------------------
% Code Summary for working in School of Life Science and Technology,UESTC.
% Author: Li Dong, e-mail: Lidong@uestc.edu.cn
% This template is for non commercial use only.
% It is freeware but not in the public domain.

% Written by Li Dong (Lidong@uestc.edu.cn), UESTC
% $ 2020.2.22
% -------------------------------------------------------------------------
if nargin < 2
    error('At least 2 inputs are requied...');
elseif nargin == 2
    elecDirecFlag = 1;
    gridresolution = [];
elseif nargin == 3
    gridresolution = [];
end


% check inputs
if ~isempty(gridresolution)
    if ~isfinite(gridresolution) || gridresolution <= 0
        warning('The vaule of dipole gridresolution is invalid, use the default settings (set it as empty)');
        gridresolution = [];
    end
end

if isempty(elecDirecFlag) || elecDirecFlag < 0
    elecDirecFlag = 1;
end

%--------------------------------------------------------------------------
% loading the standard head model
disp('----------------------------------');
disp('loading standard head model based on standard MRI...');

load(bemfile); % load headmodel
headmodel = vol;

disp('Tissues: Brain, Skull, Scalp');
disp('Number of vertices: ');
disp([' Brain = ',num2str(size(headmodel.bnd(3).pos,1))]);
disp([' Skull = ',num2str(size(headmodel.bnd(2).pos,1))]);
disp([' Scalp = ',num2str(size(headmodel.bnd(1).pos,1))]);
disp('BEM method: ''diopli''');
disp('Conductivity of tissues: ');
disp([' Brain = ',num2str(headmodel.cond(3))]);
disp([' Skull = ',num2str(headmodel.cond(2))]);
disp([' Scalp = ',num2str(headmodel.cond(1))]);

% h1 = figure;
% title('headmodel(brain,skull,scalp)');
% ft_plot_mesh(headmodel.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.4, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
% hold on;
% ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.3);
% hold on;
% ft_plot_mesh(headmodel.bnd(1),'edgecolor','none','facealpha', 0.2,'facecolor',[0.4 0.6 0.4]);
% camlight;
% hold off;
% -------------------------------------------------------------------------
% aligning the electrodes automatically to an existing head model.

disp('----------------------------------');
disp('aligning the electrodes automatically to an existing head model...');

disp('aligining method: ''project''');
disp('warp: ''globalrescale'' -- apply a rigid-body warp with global rescaling');

% read file of EEG channel locations
if isstruct(elecfile)
    eloc = elecfile;
else
    [eloc] = readlocs(elecfile); % Use EEGLAB function to load channel locations.
end
% elec = ft_read_sens('GSN-HydroCel-129.sfp');
% elec = ft_read_sens('standard_1020.elc'); % 1020 elec
% elec = ft_read_sens('16channs1.xyz');
% [eloc] = readlocs('Curry7-66channels.ced');

% check info of electrode locations (supporting XYZ and -YXZ coordinates).
chanlocsflag = 0;
if ~(isfield(eloc,'X') && isfield(eloc,'Y') && isfield(eloc,'Z') && all([length([eloc.X]),length([eloc.Y]),length([eloc.Z])] > length(eloc)*0.5))
    error('Most of your channels do not have X,Y,Z location measurements, please check your channel location file');
else
    % get the matrix of selected channel locations [3xN]
    [x,y,z] = deal({eloc.X},{eloc.Y},{eloc.Z});
    usable_channels = find(~cellfun('isempty',x) & ~cellfun('isempty',y) & ~cellfun('isempty',z));
    chanlocsXYZ = [cell2mat(x(usable_channels));cell2mat(y(usable_channels));cell2mat(z(usable_channels))];
    chanlocsXYZ = chanlocsXYZ';
    
    % generate elec structrue of Fieldtrip
    elec = [];
    % coordinates
    if elecDirecFlag == 1
        elec.chanpos = [-1*chanlocsXYZ(:,2), chanlocsXYZ(:,1), chanlocsXYZ(:,3)]; % -YXZ
        elec.elecpos = elec.chanpos;
    elseif elecDirecFlag == 0
        elec.chanpos = chanlocsXYZ;  % XYZ
        elec.elecpos = chanlocsXYZ;
    end
    % labels
    try
        for i = 1:length(usable_channels)
            elec.label{i,1} = eloc(usable_channels(i)).labels;
        end
    catch
        disp('Channel labels are invalid, renamed as channel number');
        for i = 1:length(usable_channels)
            elec.label{i,1} = num2str(i);
        end
    end
    % disp info
    chanlocsflag = 1;
    disp(['Number of valid channels: ',num2str(size(chanlocsXYZ,1))]);
    disp('Channel labels: ');
    disp((elec.label)');
end

if chanlocsflag == 1 % if the info of channel location is valid
    % rescale eletrode coordinates
    norm2 = 1;
    if isfield(elec,'elecpos')
        norm2 = sqrt(elec.elecpos(:,1).^2 + elec.elecpos(:,2).^2 + elec.elecpos(:,3).^2);
    elseif isfield(elec,'chanpos')
        norm2 = sqrt(elec.chanpos(:,1).^2 + elec.chanpos(:,2).^2 + elec.chanpos(:,3).^2);
    end
    
    norm1 = sqrt(headmodel.bnd(1).pos(:,1).^2 + headmodel.bnd(1).pos(:,2).^2 + headmodel.bnd(1).pos(:,3).^2);
    max1 = max(norm1);max2 = max(norm2);
    try
        zoom1 = 1 * max1(1)/max2(1);
        elec.elecpos = elec.elecpos * zoom1; % zoom, make sure the electrodes outside the scalp
        elec.chanpos = elec.chanpos * zoom1; % zoom, make sure the electrodes outside the scalp
    catch
        warning('failed to rescale eletrode coordinates');
    end;
    
    %     % plot check
    %     figure;
    %     hold on
    %     ft_plot_vol(headmodel, 'edgecolor', 'none', 'facecolor',[0.5 0.5 0.5],...
    %         'facealpha',0.5,'surfaceonly','false');
    %     camlight
    %     ft_plot_sens(elec,'style', 'b.');
    %     hold off;
    
    % alignment
    cfg               = [];
    cfg.method        ='project';  %  representing the method for aligning
    %  or placing the electrodes;
    % 'interactive'     realign manually using a graphical user interface
    % 'fiducial'        realign using three fiducials (e.g. NAS, LPA and RPA)
    % 'template'        realign the electrodes to match a template set
    % 'headshape'       realign the electrodes to fit the head surface
    % 'project'         projects electrodes onto the head surface
    % 'moveinward'      moves electrodes inward along their normals
    cfg.warp          = 'globalrescale'; % string describing the spatial transformation for the template and headshape methods
    % 'rigidbody'       apply a rigid-body warp (default)
    % 'globalrescale'   apply a rigid-body warp with global rescaling
    % 'traditional'     apply a rigid-body warp with individual axes rescaling
    % 'nonlin1'         apply a 1st order non-linear warp
    % 'nonlin2'         apply a 2nd order non-linear warp
    % 'nonlin3'         apply a 3rd order non-linear warp
    % 'nonlin4'         apply a 4th order non-linear warp
    % 'nonlin5'         apply a 5th order non-linear warp
    % 'dykstra2012'     non-linear wrap only for headshape method, useful for projecting ECoG onto cortex hull
    % 'fsaverage'       surface-based realignment with FreeSurfer fsaverage brain (left->left or right->right)
    % 'fsaverage_sym'   surface-based realignment with FreeSurfer fsaverage_sym left hemisphere (left->left or right->left)
    % 'fsinflated'      surface-based realignment with FreeSurfer individual subject inflated brain (left->left or right->right)
    cfg.headshape     = headmodel.bnd(1); % scalp
    cfg.elec          = elec;        % the electrodes we want to align
    
    elec_aligned = ft_electroderealign(cfg,elec); % realignment
    
    %     % plot check
    %     h2 = figure;
    %     subplot(1,2,1)
    %     hold on
    %     ft_plot_vol(headmodel, 'edgecolor', 'none', 'facecolor',[0.5 0.5 0.5],...
    %         'facealpha',0.5,'surfaceonly','false');
    %     ft_plot_sens(elec,'style', 'b.');
    %     hold off;
    %     title('Before alignement');
    %     subplot(1,2,2)
    %     ft_plot_vol(headmodel, 'edgecolor', 'none', 'facecolor',[0.5 0.5 0.5],...
    %         'facealpha',0.5,'surfaceonly','false');
    %     ft_plot_sens(elec_aligned,'style', 'r.');
    %     title('After alignement');
end
% -------------------------------------------------------------------------
% Calculating leadfield based on the BEM headmodel and aligned eletrodes
disp('----------------------------------');
disp('Calculating leadfield based on the BEM headmodel and aligned eletrodes...');
disp('generate dipole sources...');

if isempty(gridresolution)
    disp('dipoles (sources) are vertices little smaller than brain (~1mm)');
    disp('The orientations of dipoles are their normal vector directions');
    % Compute the normals of a mesh
    xyz_dipOri = wb_eeg_inv_normals(headmodel.bnd(3).pos,headmodel.bnd(3).tri); % Compute the normals of a mesh
    xyz_dipOri = -xyz_dipOri;
    sourcespos = -1*xyz_dipOri + headmodel.bnd(3).pos;  % sources are vertices little smaller than brain (1mm);
    
    cfg                   = [];
    cfg.elec              = elec_aligned;
    cfg.grid.pos          = sourcespos;  % sources
    cfg.grid.mom          = xyz_dipOri'; % normals
    cfg.grid.unit         = 'mm';
    cfg.headmodel         = headmodel;
    cfg.normalize         = 'yes';    %(default = 'no')
    cfg.reducerank        = 3;        % default is 3
    cfg.normalizeparam    = 0.5;      % depth normalization parameter (default = 0.5)
    [lf]        = ft_prepare_leadfield(cfg); % calculate leadfield
else
    disp('dipoles (sources) are vertices within brain ');
    cfg                   = [];
    cfg.elec              = elec_aligned;
    cfg.grid.unit         = 'mm';
    cfg.grid.resolution   = gridresolution; % mm
    cfg.headmodel         = headmodel;
    cfg.normalize         = 'yes';    %(default = 'no')
    cfg.reducerank        = 3;        % default is 3
    cfg.normalizeparam    = 0.5;      % depth normalization parameter (default = 0.5)
    [lf]        = ft_prepare_leadfield(cfg); % calculate leadfield
    
end

if isempty(gridresolution)
    % generate leadfield matrix (channels X sources)
    leadfield = cat(2, lf.leadfield{:});
    lf.leadfieldMatrix = leadfield;
    lf.normals = xyz_dipOri;
else
    % generate leadfield matrix (channels X sources)
    leadfield = [];
    N = size(lf.pos,1);
    m = 1;
    for i = 1:N
        if ~isempty(lf.leadfield{1,i})
            leadfield.X(:,m) = lf.leadfield{1,i}(:,1); % X orientation of the dipole.
            leadfield.Y(:,m) = lf.leadfield{1,i}(:,2); % Y orientation of the dipole.
            leadfield.Z(:,m) = lf.leadfield{1,i}(:,3); % Z orientation of the dipole.
            m = m + 1;
        end
    end
    lf.leadfieldMatrix = leadfield;
end

% figure;
% hold on;
% ft_plot_mesh(headmodel.bnd(3),'vertexcolor',[0,0,0],'vertexsize',1,...
%     'facecolor',[0.2 0.2 0.2],'facealpha',0.2,'edgealpha', 0.2);% camlight;
% quiver3(headmodel.bnd(3).pos(:,1),headmodel.bnd(3).pos(:,2),headmodel.bnd(3).pos(:,3), ...
%     xyz_dipOri(:,1),xyz_dipOri(:,2),xyz_dipOri(:,3),0.5, 'color','r','linewidth',2);
% hold off;
% ---------------------
% plotpos = 200;
% if ~isempty(gridresolution)
%     % plot check (X, Y, Z oritations)
%     figure;
%     hold on;
%     plot3(lf.pos(plotpos,1),lf.pos(plotpos,2),lf.pos(plotpos,3),'.r','MarkerSize',20); % dipole
%     % ft_plot_mesh(lf.cfg.elec.chanpos,'vertexcolor',[1,0,1],'vertexsize',12);camlight; % electrodes
%     ft_plot_mesh(lf.pos(lf.inside,:),'vertexcolor',[0.04,0.14,0.42],'vertexsize',6);camlight;
%     ft_plot_mesh(headmodel.bnd(3),'vertexcolor',[0.04,0.14,0.42],'vertexsize',1,...
%         'facecolor',[0.2 0.2 0.2],'facealpha',0.2,'edgealpha', 0.2);% camlight;
%     hold off;
%
%     figure;
%     for i = 1:3
%         subplot(2,2,i);
%         ft_plot_topo3d(lf.cfg.elec.chanpos,lf.leadfield{plotpos}(:,i));
%         title(num2str(i))
%     end
% else
%
%     % --------------------
%     % plot check (one oritation)
%     figure;
%     subplot(1,2,1);
%     hold on;
%     ft_plot_mesh(lf.cfg.elec.chanpos,'vertexcolor',[1,0,1],'vertexsize',12);camlight;
%     ft_plot_mesh(lf.pos(lf.inside,:),'vertexcolor',[0.04,0.14,0.42],'vertexsize',6);camlight;
%     ft_plot_mesh(headmodel.bnd(3),'vertexcolor',[0.04,0.14,0.42],'vertexsize',1,...
%         'facecolor',[0.2 0.2 0.2],'facealpha',0.2,'edgealpha', 0.2);% camlight;
%     quiver3(lf.pos(plotpos,1),lf.pos(plotpos,2),lf.pos(plotpos,3), ...
%         xyz_dipOri(plotpos,1),xyz_dipOri(plotpos,2),xyz_dipOri(plotpos,3),20, 'color','r','linewidth',2);
%     hold off;
%     subplot(1,2,2);
%
%     hold on;
%     ft_plot_mesh(headmodel.bnd(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.1, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
%     ft_plot_topo3d(lf.cfg.elec.chanpos(1:62,:),lf.leadfield{plotpos});
%     quiver3(lf.pos(plotpos,1),lf.pos(plotpos,2),lf.pos(plotpos,3), ...
%         xyz_dipOri(plotpos,1),xyz_dipOri(plotpos,2),xyz_dipOri(plotpos,3),40, 'color','r','linewidth',1.5);
%     hold off;
% end