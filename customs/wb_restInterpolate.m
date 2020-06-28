function [EEG_interp] = wb_restInterpolate(EEG,badchanns,chanlocsXYZ,selechanns)
%   Interpolate the badChannels of EEG.data using Reference Electrode Standardization Technique
%   Input:
%         EEG:  Structure with data and chanlocs fields compatible with
%               an EEGLAB EEG structure (requires .data)
%         badchanns: number with indices of the bad channels.badchanns
%               should belongs to selechanns.
%         chanlocsXYZ: the matrix of selected channel locations [3xN]
%         selechanns: number with indices of the selected channels.[Nx1] or [1xN].Default is all.
%   Output:
%         EEG_interp: Input EEG structure with the dest_chans rows of
%               EEG.data replaced with their interpolated values.
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
elseif nargin == 3
    selechanns = 1:size(EEG.data,1);
end

if length(selechanns) ~= size(chanlocsXYZ,2)
    error('No. of channels with locations not equal to selected channels');
end

% compute leadfield of selected channels
% -------------------
% [X,Y,Z] = sphere(70);              %利用sphere创建矩阵
% xyz1 = [X(:),Y(:),Z(:)];
% xyz1(xyz1(:,3)<0,3)= -0.076;
% % xyz1 = xyz1(xyz1(:,3)>0,:,:);
% figure;plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'.');grid on;

% [theta1,r1,z1] = cart2pol(xyz_dipoles(:,1),xyz_dipoles(:,2),xyz_dipoles(:,3));
% z1 = [-0.076,-0.0757,-0.0182,0.0394,0.0969,0.1539,0.2102,0.2656,0.3199,0.3727,0.424,0.4733,0.5206,0.5655,0.608,0.6478,0.6848,0.7187,0.7495,0.777,0.8011,0.8217,0.8386,0.8519,0.8614,0.8671,0.869];
% n1 = [400,153,154,154,153,151,150,146,144,139,134,129,123,117,110,103,95,86,78,69,60,50,40,31,20,10,1];
% % 400 dipoles on the plane
% n2 = [7,14,22,29,37,43,51,58,66,72,1];
% r1 = 0.0567:0.0868:0.868;
% r1 = [r1,0];
% xyz1 = [];
% for k1 = 1:length(r1)
%     for k2 = 1:length(n2)
%         theta1 = 0:(2*pi/n2(k2)):2*pi;
%         r1 = 0:0.0868:0.868;
%     end
% end
% % 2600 dipoles on the sphere
% theta1 = 0:pi/100:2*pi;
% x = sin(theta1);
% y = cos(theta1);
% figure;plot(x,y,'.');

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
% --------------------
goodchanns = setdiff(selechanns,badchanns);

index1 = zeros(length(goodchanns),1);
for i = 1:length(goodchanns)
    index1(i) = find(selechanns==goodchanns(i));
end

Gar = bsxfun ( @minus, G(index1,:), mean(G(index1,:)));
% Gar = G - repmat(mean(G),size(G,1),1);
temp_data = bsxfun ( @minus, EEG.data(goodchanns,:), mean(EEG.data(goodchanns,:))); % re-referencing to average reference

EEG.data(selechanns,:) = G * pinv(Gar,0.05) * temp_data;  % REST interpolation. the value 0.05 is for real data;for simulated data, it may be set as zero.
%---------------------
EEG.ref = 'REST';
EEG_interp = EEG;
end