function EEG_rest = wb_rest_rerefer_concentricspheres(EEG,chanlocsXYZ,selechanns)
%   Rereferencing to REST (Reference Electrode Standardization Technique)
%   using 3-concentric spheres head model.
%   Input:
%         EEG:  Structure with data and chanlocs fields compatible with
%               an EEGLAB EEG structure (requires .data and .chanlocs fields)
%         chanlocsXYZ: the matrix of selected channel locations [3xN]
%         selechanns: number with indices of the selected channels.[Nx1] or [1xN]. Default is all.
%   Output:
%         EEG_rest: Input EEG structure with the dest_chans rows of
%               EEG.data replaced with their REST referenced values.
%
%   For more see http://www.neuro.uestc.edu.cn/rest/
%   Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
%                       Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305
%  Li Dong*, Fali Li, Qiang Liu, Xin Wen, Yongxiu Lai, Peng Xu and Dezhong Yao*. 
%              MATLAB Toolboxes for Reference Electrode Standardization Technique (REST) 
%              of Scalp EEG. Frontiers in Neuroscience,  2017:11(601).
% -------------------------------------------------------------------------
% Code Summary for working in School of Life Science and Technology,UESTC.
% Author: Li Dong, e-mail: Lidong@uestc.edu.cn
% This template is for non commercial use only.
% It is freeware but not in the public domain.

% Written by Li Dong (Lidong@uestc.edu.cn), UESTC
% $ 2019.10.22
% -------------------------------------------------------------------------
if nargin < 3
    error('3 inputs are required at least!');
end

if length(selechanns) ~= size(chanlocsXYZ,2)
    error('No. of channels with locations not equal to selected channels');
end

if size(EEG.data,1)>size(EEG.data,2)
    warning('No. of channels > No. of time points in data???');
end

% dipoles from 'corti869-3000dipoles.dat'
xyz_dipoles = wb_dipolesXYZ(1);

% Calculate the dipole orientations.
xyz_dipOri = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
xyz_dipOri ( 2601: 3000, 1 ) = 0;
xyz_dipOri ( 2601: 3000, 2 ) = 0;
xyz_dipOri ( 2601: 3000, 3 ) = 1;

% define headmodel
headmodel        = [];
headmodel.type   = 'concentricspheres';
headmodel.o      = [ 0.0000 0.0000 0.0000 ]; % origin coordinates
headmodel.r      = [ 0.8700,0.9200,1]; % radius
headmodel.cond   = [ 1.0000,0.0125,1]; % relative conductivities of Brain,Skull and Scalp
headmodel.tissue = { 'brain' 'skull' 'scalp' };

% calculate leadfield
[G,~] = wb_calc_leadfield3(chanlocsXYZ',xyz_dipoles,xyz_dipOri,headmodel);
%---------------------
Gar = G - repmat(mean(G),size(G,1),1);
temp_data = bsxfun ( @minus, EEG.data(selechanns,:), mean(EEG.data(selechanns,:))); % re-referencing to average reference

data_z = G * pinv(Gar,0.05) * temp_data;  % the value 0.05 is for real data; 
                                     % for simulated data, it may be set as zero.
EEG.data(selechanns,:) = temp_data + repmat(mean(data_z),size(G,1),1); % V = V_avg + AVG(V_0)
EEG.ref = 'REST';

EEG_rest = EEG;

