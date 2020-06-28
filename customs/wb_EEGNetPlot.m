function wb_EEGNetPlot(M,locFile,shrinkFlag,NodeNameFlag,NodeColor,NodeNameFontSize,NodeMakerSize,titlename,model,ColorbarRange)
% Plot 2D networks across EEG eletrodes.
% For example, EEG connetivity can be calculated by coherence (0-1).
% Input:
%      M:Connection matrix. nodes X nodes
%      locFile: Direction of channel location file. Funtion,'readlocs()', is
%          used to read electrode location coordinates and other infomation.
%          By default the file type is determined using the file extension
%          (More details can be seen in subfuntion 'readlocs()').File type
%          of '*.loc' is suggested.
%          NOTE: In EEGLAB, default setting is that nose alongs +X and
%          left side alongs +Y in the cartesian coordinates read from
%          'locFile'.Therefore, this function wil automatically flip
%          the coordinates. That is,nose alongs +Y, left side alongs -X.
%      shrinkFlag: 1: shrink; 0£¨default£©: do not shrink.
%      NodeFlag: 1 (default): Show eletrode labels; 0: do not show.
%      NodeColor: Node color. Default is [0.5,0.5,0.5].
%      NodeNameFontSize: Font Size of node name. 1 X 1. Default is 8.
%      NodeMakerSize: Default is 4.
%      titlename: Name of the title. It should be a string.Default is
%          empty.
%      model: model = 1(Default): Ploting directional graph (FROM row TO column).
%                Then, red/blue lines represent bidiretional edges, while
%                green lines represent unidirectional edges.If want to plot
%                undiretional graph, M should be symmetric matrix. And the
%                line colors are red (positive) and blue (negative). Default
%                of line windth is from 0 to 2.
%             model = 2: Ploting undirectional graph. The color of edge is
%                corresponding to edge weights (minimal weight -> maximal weight (default),
%                or user defined min -> max). The line width is fixed at 1.
%     ColorbarRange: The colorbar range of edges while setting model = 2.
%                If set model = 2, default is [minimal weight, maximal weight].
%                If input 1 value, the range is input value -> maximal weight.
% Output:
%      h1: handle of the figure. Nose is at top of plot; left is left;
%          right is right.
% Examples:
%       h1 = Dong_EEGNetPlot(M,'~/*.loc');
%       h1 = Dong_EEGNetPlot(M,locFile,...)
%       h1 = Dong_EEGNetPlot(M,locFile,1,1,[],[],NodeMakerSize,titlename)
%       h1 = Dong_EEGNetPlot(M,locFile,1,1,[],[],NodeMakerSize,titlename,1)
%       h1 = Dong_EEGNetPlot(M,locFile,1,1,[],[],NodeMakerSize,titlename,2,[0,1])
% Main function was written by Li Dong (2016/3/17): Li_dong729@163.com
% -------------------------------------------------------------------------
% checking inputs
if nargin < 2
    error('Required at least 2 inputs...');
elseif nargin == 2
    shrinkFlag = 0;
    NodeNameFlag = 1;
    NodeColor = [0.5,0.5,0.5];
    NodeNameFontSize = 8;
    NodeMakerSize = 4;
    titlename = [];
    model = 1;
    ColorbarRange = [];
elseif nargin == 3
    NodeNameFlag = 1;
    NodeColor = [0.5,0.5,0.5];
    NodeNameFontSize = 8;
    NodeMakerSize = 4;
    titlename = [];
    model = 1;
    ColorbarRange = [];
elseif nargin == 4
    NodeColor = [0.5,0.5,0.5];
    NodeNameFontSize = 8;
    NodeMakerSize = 4;
    titlename = [];
    model = 1;
    ColorbarRange = [];
elseif nargin == 5
    NodeNameFontSize = 8;
    NodeMakerSize = 4;
    titlename = [];
    model = 1;
    ColorbarRange = [];
elseif nargin == 6
    NodeMakerSize = 4;
    titlename = [];
    model = 1;
    ColorbarRange = [];
elseif nargin == 7
    titlename = [];
    model = 1;
    ColorbarRange = [];
elseif nargin == 8
    model = 1;
    ColorbarRange = [];
elseif nargin == 9
    ColorbarRange = [];
end

if isempty(shrinkFlag)
    shrinkFlag = 0;
end

if isempty(NodeNameFlag)
    NodeNameFlag = 1;
end

if isempty(NodeColor)
    NodeColor = [0.5,0.5,0.5];
end

if isempty(NodeNameFontSize)
    NodeNameFontSize = 8;
end

if isempty(NodeMakerSize)
    NodeMakerSize = 4;
end
% -------------------------------------------------------------------------
% load channel location
[~,labels,Th,Rd,~] = readlocs(locFile);
Th = pi/180*Th;
[Y,X]   = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates

if NodeNameFlag == 0
    NodeName = [];
else
    NodeName = labels;  % node labels
end
N_node = length(M); % number of nodes
M = M';

% normalize coordinates (change sph_radius to 50);
% 0.85 is the default radius in EEGLAB
ShrinkScale = 80; % shrink fraction %(0 100]
if shrinkFlag == 1
    X = ShrinkScale*X/0.85;
    Y = ShrinkScale*Y/0.85;
else
    X = 100*X/0.85;
    Y = 100*Y/0.85;
end
% -------------------------------------------------------------------------
% Plot nodes
% h1 = figure;
if isempty(NodeName)
    for loc=1:N_node
        plot(X(loc),Y(loc),'o','MarkerEdgeColor','k','MarkerFaceColor',NodeColor,'markersize',NodeMakerSize);
        hold on;
    end
else
    for loc=1:N_node
        plot(X(loc),Y(loc),'o','MarkerEdgeColor','k','MarkerFaceColor',NodeColor,'markersize',NodeMakerSize);
        text(X(loc)+1,Y(loc)-1.5,NodeName{loc},'Fontsize',NodeNameFontSize);
        hold on;
    end
end
% -------------------------------------------------------------------------
% Plot 2D head
theta=0:0.01:2*pi;
m=55*cos(theta);
n=55*sin(theta);
plot(n,m,'k','linewidth',2);
% nose
mm(1)=55*cos(1/2*pi-1/20*pi);mm(2)=55*cos(1/2*pi+1/20*pi);mm(3)=55*cos(1/2*pi);
nn(1)=55*sin(1/2*pi-1/20*pi);nn(2)=55*sin(1/2*pi+1/20*pi);nn(3)=60;%55*sin(1/2*pi)+10;
for i=1:2
    plot([mm(i),mm(3)],[nn(i),nn(3)],'k','linewidth',2);
end
% right ear
alfa=-1/5*pi:0.01:1/5*pi;
mmm1=12.5*cos(alfa)+45.8;
nnn1=12.5*sin(alfa);
plot(mmm1,nnn1,'k','linewidth',2);
% left ear
alfa2=(pi-1/5*pi):0.01:(1/5*pi+pi);
mmm2=12.5*cos(alfa2)-45.8;
nnn2=12.5*sin(alfa2);
plot(mmm2,nnn2,'k','linewidth',2);
% -------------------------------------------------------------------------
% plot edges
if any(M(:))
    switch model
        case 1 % ploting directional graph (FROM row TO column)
            max1 = max(abs(M(:)));
            M = M./max1(1); % normalize to [0,1]
            % Plot bidirectional edge
            sfac = 2;
            temp = M.*M';
            [toInx,fromInx] = find(temp);
            R = [toInx fromInx];
            for i=1:length(toInx),
                xf = X(fromInx(i)); yf = Y(fromInx(i));
                xt = X(toInx(i)); yt = Y(toInx(i));
                h = line([xf xt],[yf yt]);
                set(h,'LineWidth',abs(M(toInx(i),fromInx(i))).*sfac);
                if M(toInx(i),fromInx(i)) > 0
                    set(h,'Color',[1,0.2,0]); % red for reciprocal edge (positive)
                elseif M(toInx(i),fromInx(i)) < 0
                    set(h,'Color',[0,0.4,1]); % blue reciprocal edge (negative)
                end
            end
            
            % Plot unidirectional edge
            [x2,y2] = find(M);
            nonR = [x2 y2];
            nonR = setdiff(nonR,R,'rows');
            toInx = nonR(:,1);
            fromInx = nonR(:,2);
            for i = 1:length(toInx),
                xf = X(fromInx(i)); yf = Y(fromInx(i));
                xt = X(toInx(i)); yt = Y(toInx(i));
                h = line([xf xt],[yf yt]);
                set(h,'LineWidth',abs(M(toInx(i),fromInx(i))).*sfac);
                set(h,'Color',[0.2,0.6,0]); % green for non-recip edges
                cca_arrowh([xf (xf+xt)/2 ],[yf (yf+yt)/2],[0.2,0.6,0],100,60);
            end
        case 2 % Ploting undirectional graph.
            % Plot undirectional edge
            edge_colors = colormap(jet(128));
            weights = M(M~=0);
            if isempty(ColorbarRange);
                minWei = min(weights);
                maxWei = max(weights);
            else
                if length(ColorbarRange) == 1
                    minWei = ColorbarRange(1);
                    maxWei = max(weights);
                else
                    minWei = ColorbarRange(1);
                    maxWei = ColorbarRange(2);
                end
            end
%             colorbar;
%             caxis([minWei(1),maxWei(1)]);
            
            sfac = 1; % linewidth of edges
            % temp = M + M';
            [toInx,fromInx] = find(M);
            for i=1:length(toInx),
                xf = X(fromInx(i)); yf = Y(fromInx(i));
                xt = X(toInx(i)); yt = Y(toInx(i));
                h = line([xf xt],[yf yt]);
                if maxWei(1) == minWei(1)
                    k = 1;
                else
                    k = round(128*(M(toInx(i),fromInx(i))-minWei(1))/(maxWei(1) - minWei(1)));
                end
                set(h,'Color',edge_colors(max(k,1),:),'LineWidth',sfac); % red for reciprocal edge (positive)
            end
            % ------------------------------------------------------------------
    end
end
hold off
title(titlename);
axis equal;
% axis('square');
axis off;
% -------------------------------------------------------------------------
% sub Function
%
%  ARROWH   Draws a solid 2D arrow head in current plot.
%     ARROWH(X,Y,COLOR,SIZE,LOCATION) draws a solid arrow head into
%     the current plot to indicate a direction.  X and Y must contain
%     a pair of x and y coordinates ([x1 x2],[y1 y2]) of two points:
%
%     The first point is only used to tell (in conjunction with the
%     second one) the direction and orientation of the arrow -- it
%     will point from the first towards the second.
%
%     The head of the arrow will be located in the second point.  An
%     example of use is    plot([0 2],[0 4]); ARROWH([0 1],[0 2],'b')
%
%     You may also give two vectors of same length > 2.  The routine
%     will then choose two consecutive points from "about" the middle
%     of each vectors.  Useful if you don't want to worry each time
%     about where to put the arrows on a trajectory.  If x1 and x2
%     are the vectors x1(t) and x2(t), simply put   ARROWH(x1,x2,'r')
%     to have the right direction indicated in your x2 = f(x1) phase
%     plane.
%
%                       (x2,y2)
%                       --o
%                       \ |
%                        \|
%
%
%            o
%        (x1,y1)
%
%     Please note that the following optional arguments need -- if
%     you want to use them -- to be given in that exact order.
%
%     The COLOR argument is exactely the same as for plots, eg. 'r';
%     if not given, blue is default.
%
%     The SIZE argument allows you to tune the size of the arrows.
%
%     The LOCAITON argument only applies, if entire solution vectors
%     have been passed on.  With this argument you can indicate where
%     abouts inside those vectors to take the two points from.
%     Can be a vector, if you want to have more than one arrow drawn.
%
%     Both arguments, SIZE and LOCATION must be given in percent,
%     where 100 means standard size, 50 means half size, respectively
%     100 means end of the vector, 48 means about middle, 0 beginning.
%     Note that those "locations" correspond to the cardinal position
%     "inside" the vector, say "index-wise".
%
%     This routine is mainely intended to be used for indicating
%     "directions" on trajectories -- just give two consecutive times
%     and the corresponding values of a flux and the proper direction
%     of the trajectory will be shown on the plot.  You may also pass
%     on two solution vectors, as described above.
%
%     Note, that the arrow only looks good on the specific axis
%     settings when the routine was actually started.  If you zoom in
%     afterwards, the triangle gets distorted.
%
%     Examples of use:
%     x1 = [0:.2:2]; x2 = [0:.2:2]; plot(x1,x2); hold on;
%     arrowh(x1,x2,'r',100,20);      % passing on entire vectors
%     arrowh([0 1],[0 1],'g',300);   % passing on 2 points

%     Author:       Florian Knorn
%     Email:        florian.knorn@student.uni-magdeburg.de
%     Version:      1.08
%     Filedate:     Oct 28th, 2004
%
%     ToDos:        - More specific shaping-possibilities,
%                   - Keep proportions when zooming or resizing
%     Bugs:         None discovered yet
%
%     If you have suggestions for this program, if it doesn't work for
%     your "situation" or if you change something in it - please send
%     me an email!  This is my very first "public" program and I'd like
%     to improve it where I can -- your help is kindely appreciated!
%     Thank you!

function cca_arrowh(x,y,clr,ArSize,Where);

%-- errors
if exist('x')*exist('y') ~= 1,
    error('Please give enough coordinates!');
end
if ((length(x) < 2) | (length(y) < 2)),
    error('X and Y vectors must each have "length" 2!');
end
if (x(1) == x(2)) & (y(1) == y(2)),
    error('Points superimposed - cannot determine direction!');
end

%-- determine and remember the hold status, toggle if necessary
if ishold == 0,
    WasHold = 0;
    hold on;
else
    WasHold = 1;
end

%-- start for-loop in case several arrows are wanted
for Loop = 1:length(Where)
    %clear ArSize
    
    %-- no errors, move on. if vectors "longer" then 2 are given
    if (length(x) == length(y)) & (length(x) > 2),
        if exist('Where') == 1, %-- and user said, where abouts to put arrow
            j = floor(length(x)*Where(Loop)/100); %-- determine that location
            if j >= length(x), j = length(x) - 1; end
            if j == 0, j = 1; end
        else %-- he didn't tell where, so take the "middle" by default
            j = floor(length(x)/2);
        end %-- now pick the right couple of coordinates
        x1 = x(j); x2 = x(j+1); y1 = y(j); y2 = y(j+1);
        
    else %-- just two points given - must take those
        x1 = x(1); x2 = x(2); y1 = y(1); y2 = y(2);
    end
    
    %-- if no color given, use blue as default
    if exist('clr') ~= 1
        clr = 'b';
    end
    
    %-- determine if size argument was given, if not, set it to default
    if exist('ArSize') ~= 1,
        ArSize = 100 / 10000; %-- 10000 is an arbitrary value...
    else
        ArSize = ArSize / 10000;
    end
    
    %-- get axe ranges and their norm
    OriginalAxis = axis;
    Xextend = abs(OriginalAxis(2)-OriginalAxis(1));
    Yextend = abs(OriginalAxis(4)-OriginalAxis(3));
    
    %-- determine angle for the rotation of the triangle
    if x2 == x1, %-- line vertical, no need to calculate slope
        if y2 > y1
            p = pi/2;
        else
            p= -pi/2;
        end
    else %-- line not vertical, go ahead and calculate slope
        %-- using normed differences (looks better like that)
        m = ( (y2 - y1)/Yextend ) / ( (x2 - x1)/Xextend );
        if x2 > x1, %-- now calculate the resulting angle
            p = atan(m);
        else
            p = atan(m) + pi;
        end
    end
    
    %-- the arrow is made of a transformed "template triangle".
    %-- it will be created, rotated, moved, resized and shifted.
    
    %-- the template triangle (it points "east", centered in (0,0)):
    xt = [1    -sin(pi/6)    -sin(pi/6)];
    yt = [0     cos(pi/6)    -cos(pi/6)];
    
    %-- rotate it by the angle determined above:
    for i=1:3,
        xd(i) = cos(p)*xt(i) - sin(p)*yt(i);
        yd(i) = sin(p)*xt(i) + cos(p)*yt(i);
    end
    
    %-- move the triangle so that its "head" lays in (0,0):
    xd = xd - cos(p);
    yd = yd - sin(p);
    
    %-- stretch/deform the triangle to look good on the current axes:
    xd = xd*Xextend*ArSize;
    yd = yd*Yextend*ArSize;
    
    %-- move the triangle to the location where it's needed
    xd = xd + x2;
    yd = yd + y2;
    
    %-- draw the actual triangle
    patch(xd,yd,clr,'EdgeColor',clr);
    
end % Loops

%-- restore original axe ranges and hold status
axis(OriginalAxis);
if WasHold == 0,
    hold off
end
%-- work done. good bye.
% -------------------------------------------------------------------------
% readlocs() - read electrode location coordinates and other information from a file.
%              Several standard file formats are supported. Users may also specify
%              a custom column format. Defined format examples are given below
%              (see File Formats).
% Usage:
%   >>  eloc = readlocs( filename );
%   >>  EEG.chanlocs = readlocs( filename, 'key', 'val', ... );
%   >>  [eloc, labels, theta, radius, indices] = ...
%                                               readlocs( filename, 'key', 'val', ... );
% Inputs:
%   filename   - Name of the file containing the electrode locations
%                {default: 2-D polar coordinates} (see >> help topoplot )
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'sfp'|'xyz'|'asc'|'polhemus'|'besa'|'chanedit'|'custom']
%                 Type of the file to read. By default the file type is determined
%                 using the file extension (see below under File Formats):
%                  'loc' - an EEGLAB 2-D polar coordinates channel locations file
%                          Coordinates are theta and radius (see definitions below).
%                  'sph' - Matlab spherical coordinates (Note: spherical
%                          coordinates used by Matlab functions are different
%                          from spherical coordinates used by BESA - see below).
%                  'sfp' - EGI Cartesian coordinates (NOT Matlab Cartesian - see below).
%                  'xyz' - Matlab/EEGLAB Cartesian coordinates (NOT EGI Cartesian).
%                          z is toward nose; y is toward left ear; z is toward vertex
%                  'asc' - Neuroscan polar coordinates.
%                  'polhemus' or 'polhemusx' - Polhemus electrode location file recorded
%                          with 'X' on sensor pointing to subject (see below and readelp()).
%                  'polhemusy' - Polhemus electrode location file recorded with
%                          'Y' on sensor pointing to subject (see below and readelp()).
%                  'besa' - BESA-'.elp' spherical coordinates. (Not MATLAB spherical -
%                           see below).
%                  'chanedit' - EEGLAB channel location file created by pop_chanedit().
%                  'custom' - Ascii file with columns in user-defined 'format' (see below).
%   'importmode' - ['eeglab'|'native'] for location files containing 3-D cartesian electrode
%                  coordinates, import either in EEGLAB format (nose pointing toward +X).
%                  This may not always be possible since EEGLAB might not be able to
%                  determine the nose direction for scanned electrode files. 'native' import
%                  original carthesian coordinates (user can then specify the position of
%                  the nose when calling the topoplot() function; in EEGLAB the position
%                  of the nose is stored in the EEG.chaninfo structure). {default 'eeglab'}
%   'format'    - [cell array] Format of a 'custom' channel location file (see above).
%                          {default: if no file type is defined. The cell array contains
%                          labels defining the meaning of each column of the input file:
%                           'channum'   [positive integer] channel number
%                           'labels'    [string] channel name (no spaces)
%                           'theta'     [real degrees] 2-D angle in polar coordinates;
%                                       positive = rotating from nose (0) toward left ear
%                           'radius'    [real] radius for 2-D polar coords; 0.5 is the head
%                                       disk radius and limit for topoplot() plotting)
%                           'X'         [real] Matlab-Cartesian X coordinate (to nose)
%                           'Y'         [real] Matlab-Cartesian Y coordinate (to left ear)
%                           'Z'         [real] Matlab-Cartesian Z coordinate (to vertex)
%                           '-X','-Y','-Z' Matlab-Cartesian coordinates pointing opposite
%                                       to the above.
%                           'sph_theta' [real degrees] Matlab spherical horizontal angle;
%                                       positive = rotating from nose (0) toward left ear.
%                           'sph_phi'   [real degrees] Matlab spherical elevation angle;
%                                       positive = rotating from horizontal (0) upwards.
%                           'sph_radius' [real] distance from head center (unused)
%                           'sph_phi_besa' [real degrees] BESA phi angle from vertical;
%                                       positive = rotating from vertex (0) towards right ear.
%                           'sph_theta_besa' [real degrees] BESA theta horiz/azimuthal angle;
%                                       positive = rotating from right ear (0) toward nose.
%                           'ignore'    ignore column}
%     The input file may also contain other channel information fields
%                           'type'      channel type: 'EEG', 'MEG', 'EMG', 'ECG', others ...
%                           'calib'     [real near 1.0] channel calibration value.
%                           'gain'      [real > 1] channel gain.
%                           'custom1'   custom field #1.
%                           'custom2', 'custom3', 'custom4', etc.    more custom fields
%   'skiplines' - [integer] Number of header lines to skip (in 'custom' file types only).
%                 Note: Characters on a line following '%' will be treated as comments.
%   'readchans' - [integer array] indices of electrodes to read. {default: all}
%   'center'    - [(1,3) real array or 'auto'] center of xyz coordinates for conversion
%                 to spherical or polar, Specify the center of the sphere here, or 'auto'.
%                 This uses the center of the sphere that best fits all the electrode
%                 locations read. {default: [0 0 0]}
% Outputs:
%   eloc        - structure containing the channel names and locations (if present).
%                 It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius'
%                 identical in meaning to the EEGLAB struct 'EEG.chanlocs'.
%   labels      - cell array of strings giving the names of the electrodes. NOTE: Unlike the
%                 three outputs below, includes labels of channels *without* location info.
%   theta       - vector (in degrees) of polar angles of the electrode locations.
%   radius      - vector of polar-coordinate radii (arc_lengths) of the electrode locations
%   indices     - indices, k, of channels with non-empty 'locs(k).theta' coordinate
%
% File formats:
%   If 'filetype' is unspecified, the file extension determines its type.
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
%
% Author: Arnaud Delorme, Salk Institute, 8 Dec 2002 (expanded from the previous EEG/ICA
%         toolbox function)
%
% See also: readelp(), writelocs(), topo2sph(), sph2topo(), sph2cart()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: readlocs.m,v $
% Revision 1.87  2006/06/01 17:40:50  arno
% updating header
%
% Revision 1.86  2006/05/26 15:59:35  scott
% worked on text of instructions for adding a channel type -- NOTE: chaninfo is
% discussed in help message, but NOT implemented!?
%
% Revision 1.85  2006/04/14 21:19:08  arno
% fixing skipping lines
%
% Revision 1.84  2006/03/31 03:11:13  toby
% made '.eloc' equivalent to '.loc' as a filetype
%
% Revision 1.83  2006/02/14 00:01:18  arno
% change xyz format
%
% Revision 1.82  2006/01/20 22:37:08  arno
% default for BESA and polhemus
%
% Revision 1.81  2006/01/12 23:22:39  arno
% fixing indices
%
% Revision 1.80  2006/01/12 22:03:51  arno
% fiducial type
%
% Revision 1.79  2006/01/10 22:56:17  arno
% adding defaultelp option
%
% Revision 1.78  2006/01/10 22:53:49  arno
% [6~[6~changing default besa format
%
% Revision 1.77  2005/11/30 18:31:40  arno
% same
%
% Revision 1.76  2005/11/30 18:29:48  arno
% same
%
% Revision 1.75  2005/11/30 18:28:37  arno
% reformat outputs
%
% Revision 1.74  2005/10/29 03:49:50  scott
% NOTE: there  is no mention of 'chantype' - should at least add a help mention after line  69 -sm
%
% Revision 1.73  2005/09/27 22:08:41  arno
% fixing reading .ced files
%
% Revision 1.72  2005/05/24 17:07:05  arno
% cell2mat - celltomat
%
% Revision 1.71  2005/03/10 17:42:11  arno
% new format for channel location info
%
% Revision 1.70  2005/03/08 23:19:24  arno
% using old function to read asa format
%
% Revision 1.69  2005/03/04 23:17:22  arno
% use fieldtrip readeetrack
%
% Revision 1.65  2004/10/27 01:01:05  arno
% msg format
%
% Revision 1.64  2004/03/23 00:37:56  scott
% clarifying help msg re meaning of 'indices' output
%
% Revision 1.63  2004/03/23 00:22:51  scott
% clarified meaning of output 'indices'
%
% Revision 1.62  2004/02/24 17:17:32  arno
% dbug message
%
% Revision 1.61  2004/01/01 19:12:08  scott
% help message edits
%
% Revision 1.60  2004/01/01 18:57:26  scott
% edit text outputs
%
% Revision 1.59  2004/01/01 01:47:34  scott
% franglais -> anglais
%
% Revision 1.58  2003/12/17 00:55:07  arno
% debug last
%
% Revision 1.57  2003/12/17 00:50:10  arno
% adding index for non-empty electrodes
%
% Revision 1.56  2003/12/05 18:37:56  arno
% debug polhemus x and y fixed
%
% Revision 1.55  2003/12/02 03:21:39  arno
% neuroscan format
%
% Revision 1.54  2003/11/27 00:38:13  arno
% conversion elc
%
% Revision 1.53  2003/11/27 00:31:30  arno
% debuging elc format
%
% Revision 1.52  2003/11/27 00:25:51  arno
% automatically detecting elc files
%
% Revision 1.51  2003/11/05 17:20:23  arno
% first convert spherical instead of carthesian
%
% Revision 1.50  2003/09/18 00:07:05  arno
% further checks for neuroscan
%
% Revision 1.49  2003/07/16 18:52:21  arno
% allowing file type locs
%
% Revision 1.48  2003/06/30 15:00:43  arno
% fixing inputcheck problem
%
% Revision 1.47  2003/05/13 23:31:25  arno
% number of lines to skip in chanedit format
%
% Revision 1.46  2003/05/13 22:09:01  arno
% updating sph format
%
% Revision 1.45  2003/05/13 22:07:07  arno
% removing labels in sfp format
%
% Revision 1.44  2003/05/13 21:14:11  arno
% only write a subset of file format
%
% Revision 1.43  2003/03/10 16:28:12  arno
% removing help for elc
%
% Revision 1.42  2003/03/10 16:26:59  arno
% adding then removing .elc format
%
% Revision 1.41  2003/03/08 17:36:13  arno
% import spherical EGI files correctly
%
% Revision 1.40  2003/03/05 15:38:15  arno
% fixing '.' bug
%
% Revision 1.39  2003/03/04 20:04:44  arno
% adding neuroscan .asc format
%
% Revision 1.38  2003/01/30 16:45:12  arno
% debugging ced format
%
% Revision 1.37  2003/01/10 17:40:11  arno
% removing trailing dots
%
% Revision 1.36  2003/01/03 22:47:00  arno
% typo in warning messages
%
% Revision 1.35  2003/01/03 22:45:48  arno
% adding another warning message
%
% Revision 1.34  2003/01/03 22:41:38  arno
% autodetect format .sfp
%
% Revision 1.33  2003/01/03 22:38:39  arno
% adding warning message
%
% Revision 1.32  2002/12/29 23:04:00  scott
% header
%
% Revision 1.31  2002/12/29 22:37:15  arno
% txt -> ced
%
% Revision 1.30  2002/12/29 22:35:35  arno
% adding coords. info for file format in header, programming .sph, ...
%
% Revision 1.29  2002/12/29 22:00:10  arno
% skipline -> skiplines
%
% Revision 1.28  2002/12/28 23:46:45  scott
% header
%
% Revision 1.27  2002/12/28 02:02:35  scott
% header details
%
% Revision 1.26  2002/12/28 01:32:41  scott
% worked on header information - axis details etcetc. -sm & ad
%
% Revision 1.25  2002/12/27 23:23:35  scott
% edit header msg - NEEDS MORE DETAILS -sm
%
% Revision 1.24  2002/12/27 22:57:23  arno
% debugging polhemus
%
% Revision 1.23  2002/12/27 17:47:32  arno
% compatible with more BESA formats
%
% Revision 1.22  2002/12/26 16:41:23  arno
% new release
%
% Revision 1.21  2002/12/24 02:51:22  arno
% new version of readlocs
%


function [eloc, labels, theta, radius, indices] = readlocs( filename, varargin );

if nargin < 1
    help readlocs;
    return;
end;

% NOTE: To add a new channel format:
% ----------------------------------
% 1) Add a new element to the structure 'chanformat' (see 'ADD NEW FORMATS HERE' below):
% 2)  Enter a format 'type' for the new file format,
% 3)  Enter a (short) 'typestring' description of the format
% 4)  Enter a longer format 'description' (possibly multiline, see ex. (1) below)
% 5)  Enter format file column labels in the 'importformat' field (see ex. (2) below)
% 6)  Enter the number of header lines to skip (if any) in the 'skipline' field
% 7)  Document the new channel format in the help message above.
% 8)  After testing, please send the new version of readloca.m to us
%       at eeglab@sccn.ucsd.edu with a sample locs file.
% The 'chanformat' structure is also used (automatically) by the writelocs()
% and pop_readlocs() functions. You do not need to edit these functions.

chanformat(1).type         = 'polhemus';
chanformat(1).typestring   = 'Polhemus native .elp file';
chanformat(1).description  = [ 'Polhemus native coordinate file containing scanned electrode positions. ' ...
    'User must select the direction ' ...
    'for the nose after importing the data file.' ];
chanformat(1).importformat = 'readelp() function';
% ---------------------------------------------------------------------------------------------------
chanformat(2).type         = 'besa';
chanformat(2).typestring   = 'BESA spherical .elp file';
chanformat(2).description  = [ 'BESA spherical coordinate file. Note that BESA spherical coordinates ' ...
    'are different from Matlab spherical coordinates' ];
chanformat(2).skipline     = -1;
chanformat(2).importformat = { 'type' 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' };
% ---------------------------------------------------------------------------------------------------
chanformat(3).type         = 'xyz';
chanformat(3).typestring   = 'Matlab .xyz file';
chanformat(3).description  = [ 'Standard 3-D cartesian coordinate files with electrode labels in ' ...
    'the first column and X, Y, and Z coordinates in columns 2, 3, and 4' ];
chanformat(3).importformat = { 'channum' '-Y' 'X' 'Z' 'labels'};
% ---------------------------------------------------------------------------------------------------
chanformat(4).type         = 'sfp';
chanformat(4).typestring   = 'BESA or EGI 3-D cartesian .sfp file';
chanformat(4).description  = [ 'Standard BESA 3-D cartesian coordinate files with electrode labels in ' ...
    'the first column and X, Y, and Z coordinates in columns 2, 3, and 4.' ...
    'Coordinates are re-oriented to fit the EEGLAB standard of having the ' ...
    'nose along the +X axis.' ];
chanformat(4).importformat = { 'labels' '-Y' 'X' 'Z' };
chanformat(4).skipline     = -1;
% ---------------------------------------------------------------------------------------------------
chanformat(5).type         = 'loc';
chanformat(5).typestring   = 'EEGLAB polar .loc file';
chanformat(5).description  = [ 'EEGLAB polar .loc file' ];
chanformat(5).importformat = { 'channum' 'theta' 'radius' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(6).type         = 'sph';
chanformat(6).typestring   = 'Matlab .sph spherical file';
chanformat(6).description  = [ 'Standard 3-D spherical coordinate files in Matlab format' ];
chanformat(6).importformat = { 'channum' 'sph_theta' 'sph_phi' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(7).type         = 'asc';
chanformat(7).typestring   = 'Neuroscan polar .asc file';
chanformat(7).description  = [ 'Neuroscan polar .asc file, automatically recentered to fit EEGLAB standard' ...
    'of having ''Cz'' at (0,0).' ];
chanformat(7).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(8).type         = 'dat';
chanformat(8).typestring   = 'Neuroscan 3-D .dat file';
chanformat(8).description  = [ 'Neuroscan 3-D cartesian .dat file. Coordinates are re-oriented to fit ' ...
    'the EEGLAB standard of having the nose along the +X axis.' ];
chanformat(8).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(9).type         = 'elc';
chanformat(9).typestring   = 'ASA .elc 3-D file';
chanformat(9).description  = [ 'ASA .elc 3-D coordinate file containing scanned electrode positions. ' ...
    'User must select the direction ' ...
    'for the nose after importing the data file.' ];
chanformat(9).importformat = 'readeetraklocs';
% ---------------------------------------------------------------------------------------------------
chanformat(10).type         = 'chanedit';
chanformat(10).typestring   = 'EEGLAB complete 3-D file';
chanformat(10).description  = [ 'EEGLAB file containing polar, cartesian 3-D, and spherical 3-D ' ...
    'electrode locations.' ];
chanformat(10).importformat = { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' ...
    'sph_radius' };
chanformat(10).skipline     = 1;
% ---------------------------------------------------------------------------------------------------
chanformat(11).type         = 'custom';
chanformat(11).typestring   = 'Custom file format';
chanformat(11).description  = 'Custom ASCII file format where user can define content for each file columns.';
chanformat(11).importformat = '';
% ---------------------------------------------------------------------------------------------------
% ----- ADD MORE FORMATS HERE -----------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
    'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
    'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };

% ----------------------------------
% special mode for getting the info
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
    eloc = chanformat;
    labels = listcolformat;
    return;
end;

g = finputcheck( varargin, ...
    { 'filetype'	   'string'  {}                 '';
    'importmode'  'string'  { 'eeglab' 'native' } 'eeglab';
    'defaultelp'  'string'  { 'besa'   'polhemus' } 'polhemus';
    'skiplines'   'integer' [0 Inf] 			[];
    'elecind'     'integer' [1 Inf]	    	[];
    'format'	   'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;

if isstr(filename)
    
    % format auto detection
    % --------------------
    if strcmpi(g.filetype, 'autodetect'), g.filetype = ''; end;
    g.filetype = strtok(g.filetype);
    periods = find(filename == '.');
    fileextension = filename(periods(end)+1:end);
    g.filetype = lower(g.filetype);
    if isempty(g.filetype)
        switch lower(fileextension),
            case {'loc' 'locs' }, g.filetype = 'loc';
            case 'xyz', g.filetype = 'xyz';
                fprintf( [ 'WARNING: Matlab Cartesian coord. file extension (".xyz") detected.\n' ...
                    'If importing EGI Cartesian coords, force type "sfp" instead.\n'] );
            case 'sph', g.filetype = 'sph';
            case 'ced', g.filetype = 'chanedit';
            case 'elp', g.filetype = g.defaultelp;
            case 'asc', g.filetype = 'asc';
            case 'dat', g.filetype = 'dat';
            case 'elc', g.filetype = 'elc';
            case 'eps', g.filetype = 'besa';
            case 'sfp', g.filetype = 'sfp';
            otherwise, g.filetype =  '';
        end;
        fprintf('readlocs(): ''%s'' format assumed from file extension\n', g.filetype);
    else
        if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end
        if strcmpi(g.filetype, 'eloc'),  g.filetype = 'loc'; end
    end;
    
    % assign format from filetype
    % ---------------------------
    if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
            & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
        indexformat = strmatch(lower(g.filetype), { chanformat.type }, 'exact');
        g.format = chanformat(indexformat).importformat;
        if isempty(g.skiplines)
            g.skiplines = chanformat(indexformat).skipline;
        end;
        if isempty(g.filetype)
            error( ['readlocs() error: The filetype cannot be detected from the \n' ...
                '                  file extension, and custom format not specified']);
        end;
    end;
    
    % import file
    % -----------
    if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
        eloc = readneurolocs( filename );
        eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
        eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
    elseif strcmp(g.filetype, 'elc')
        eloc = readeetraklocs( filename );
        %eloc = read_asa_elc( filename ); % from fieldtrip
        %eloc = struct('labels', eloc.label, 'X', mattocell(eloc.pnt(:,1)'), 'Y', ...
        %                        mattocell(eloc.pnt(:,2)'), 'Z', mattocell(eloc.pnt(:,3)'));
        eloc = convertlocs(eloc, 'cart2all');
        eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
        eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
    elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
            strcmp(g.filetype, 'polhemus')
        try,
            [eloc labels X Y Z]= readelp( filename );
            if strcmp(g.filetype, 'polhemusy')
                tmp = X; X = Y; Y = tmp;
            end;
            for index = 1:length( eloc )
                eloc(index).X = X(index);
                eloc(index).Y = Y(index);
                eloc(index).Z = Z(index);
            end;
        catch,
            disp('readlocs(): Could not read Polhemus coords. Trying to read BESA .elp file.');
            [eloc, labels, theta, radius, indices] = readlocs( filename, 'defaultelp', 'besa', varargin{:} );
        end;
    else
        % importing file
        % --------------
        if isempty(g.skiplines), g.skiplines = 0; end;
        array = load_file_or_array( filename, max(g.skiplines,0));
        if size(array,2) < length(g.format)
            fprintf(['readlocs() warning: Fewer columns in the input than expected.\n' ...
                '                    See >> help readlocs\n']);
        elseif size(array,2) > length(g.format)
            fprintf(['readlocs() warning: More columns in the input than expected.\n' ...
                '                    See >> help readlocs\n']);
        end;
        
        % removing lines BESA
        % -------------------
        if g.skiplines == -1
            if isempty(array{1,2})
                disp('BESA header detected, skipping three lines...');
                array = load_file_or_array( filename, -2);
            end;
        end;
        
        % removing comments and empty lines
        % ---------------------------------
        indexbeg = 1;
        while isempty(array{indexbeg,1}) | ...
                (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '%' )
            indexbeg = indexbeg+1;
        end;
        array = array(indexbeg:end,:);
        
        % converting file
        % ---------------
        for indexcol = 1:min(size(array,2), length(g.format))
            [str mult] = checkformat(g.format{indexcol});
            for indexrow = 1:size( array, 1)
                if mult ~= 1
                    eval ( [ 'eloc(indexrow).'  str '= -array{indexrow, indexcol};' ]);
                else
                    eval ( [ 'eloc(indexrow).'  str '= array{indexrow, indexcol};' ]);
                end;
            end;
        end;
    end;
    
    % handling BESA coordinates
    % -------------------------
    if isfield(eloc, 'sph_theta_besa')
        if isnumeric(eloc(1).type)
            disp('BESA format detected ( Theta | Phi )');
            for index = 1:length(eloc)
                eloc(index).sph_phi_besa   = eloc(index).labels;
                eloc(index).sph_theta_besa = eloc(index).type;
                eloc(index).labels         = '';
                eloc(index).type           = '';
            end;
            eloc = rmfield(eloc, 'labels');
        elseif isnumeric(eloc(1).labels)
            disp('BESA format detected ( Elec | Theta | Phi )');
            for index = 1:length(eloc)
                eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
                eloc(index).sph_theta_besa = eloc(index).labels;
                eloc(index).labels         = eloc(index).type;
                eloc(index).type           = '';
                eloc(index).radius         = 1;
            end;
        else
            disp('BESA format detected ( Type | Elec | Theta | Phi | Radius )');
        end;
        
        try
            eloc = convertlocs(eloc, 'sphbesa2all');
            eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
        catch, disp('Warning: coordinate conversion failed'); end;
        fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');
        fprintf('          to avoid confusion (these fields can be exported, though)\n');
        eloc = rmfield(eloc, 'sph_phi_besa');
        eloc = rmfield(eloc, 'sph_theta_besa');
        
        % converting XYZ coordinates to polar
        % -----------------------------------
    elseif isfield(eloc, 'sph_theta')
        try
            eloc = convertlocs(eloc, 'sph2all');
        catch, disp('Warning: coordinate conversion failed'); end;
    elseif isfield(eloc, 'X')
        try
            eloc = convertlocs(eloc, 'cart2all');
        catch, disp('Warning: coordinate conversion failed'); end;
    else
        try
            eloc = convertlocs(eloc, 'topo2all');
        catch, disp('Warning: coordinate conversion failed'); end;
    end;
    
    % inserting labels if no labels
    % -----------------------------
    if ~isfield(eloc, 'labels')
        fprintf('readlocs(): Inserting electrode labels automatically.\n');
        for index = 1:length(eloc)
            eloc(index).labels = [ 'E' int2str(index) ];
        end;
    else
        % remove trailing '.'
        for index = 1:length(eloc)
            if isstr(eloc(index).labels)
                tmpdots = find( eloc(index).labels == '.' );
                eloc(index).labels(tmpdots) = [];
            end;
        end;
    end;
    
    % resorting electrodes if number not-sorted
    % -----------------------------------------
    if isfield(eloc, 'channum')
        if ~isnumeric(eloc(1).channum)
            error('Channel numbers must be numeric');
        end;
        allchannum = [ eloc.channum ];
        if any( sort(allchannum) ~= allchannum )
            fprintf('readlocs(): Re-sorting channel numbers based on ''channum'' column indices\n');
            [tmp newindices] = sort(allchannum);
            eloc = eloc(newindices);
        end;
        eloc = rmfield(eloc, 'channum');
    end;
else
    if isstruct(filename)
        eloc = filename;
    else
        disp('readlocs(): input variable must be a string or a structure');
    end;
end;
if ~isempty(g.elecind)
    eloc = eloc(g.elecind);
end;
if nargout > 2
    tmptheta          = { eloc.theta }; % check which channels have (polar) coordinates set
    indices           = find(~cellfun('isempty', tmptheta));
    indbad            = find(cellfun('isempty', tmptheta));
    tmptheta(indbad)  = { NaN };
    theta             = [ tmptheta{:} ];
end;
if nargout > 3
    tmprad            = { eloc.radius };
    tmprad(indbad)    = { NaN };
    radius            = [ tmprad{:} ];
end;
%tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
%disp('Converting channel labels to string');
for index = 1:length(eloc)
    if ~isstr(eloc(index).labels)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
labels = { eloc.labels };
if isfield(eloc, 'ignore')
    eloc = rmfield(eloc, 'ignore');
end;

% process fiducials if any
% ------------------------
fidnames = { 'nz' 'lpa' 'rpa' };
for index = 1:length(fidnames)
    ind = strmatch(fidnames{index}, lower(labels), 'exact');
    if ~isempty(ind), eloc(ind).type = 'FID'; end;
end;

return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skiplines )
if isempty(skiplines),
    skiplines = 0;
end;
if exist(varname) == 2
    array = loadtxt(varname,'verbose','off','skipline',skiplines);
else % variable in the global workspace
    % --------------------------
    try, array = evalin('base', varname);
    catch, error('readlocs(): cannot find the named file or variable, check syntax');
    end;
end;
return;

% check field format
% ------------------
function [str, mult] = checkformat(str)
mult = 1;
if strcmpi(str, 'labels'),         str = lower(str); return; end;
if strcmpi(str, 'channum'),        str = lower(str); return; end;
if strcmpi(str, 'theta'),          str = lower(str); return; end;
if strcmpi(str, 'radius'),         str = lower(str); return; end;
if strcmpi(str, 'ignore'),         str = lower(str); return; end;
if strcmpi(str, 'sph_theta'),      str = lower(str); return; end;
if strcmpi(str, 'sph_phi'),        str = lower(str); return; end;
if strcmpi(str, 'sph_radius'),     str = lower(str); return; end;
if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
if strcmpi(str, 'sph_phi_besa'),   str = lower(str); return; end;
if strcmpi(str, 'gain'),           str = lower(str); return; end;
if strcmpi(str, 'calib'),          str = lower(str); return; end;
if strcmpi(str, 'type') ,          str = lower(str); return; end;
if strcmpi(str, 'X'),              str = upper(str); return; end;
if strcmpi(str, 'Y'),              str = upper(str); return; end;
if strcmpi(str, 'Z'),              str = upper(str); return; end;
if strcmpi(str, '-X'),             str = upper(str(2:end)); mult = -1; return; end;
if strcmpi(str, '-Y'),             str = upper(str(2:end)); mult = -1; return; end;
if strcmpi(str, '-Z'),             str = upper(str(2:end)); mult = -1; return; end;
if strcmpi(str, 'custum1'), return; end;
if strcmpi(str, 'custum2'), return; end;
if strcmpi(str, 'custum3'), return; end;
if strcmpi(str, 'custum4'), return; end;
error(['readlocs(): undefined field ''' str '''']);
% -----------------------------------------------------
% finputcheck() - check Matlab function {'key','value'} input argument pairs
%
% Usage: >> result = finputcheck( varargin, fieldlist );
%        >> [result varargin] = finputcheck( varargin, fieldlist, ...
%                                                         callingfunc, mode );
% Input:
%   varargin  - 'varargin' argument from a function call using 'key', 'value'
%               argument pairs.
%   fieldlist - A 3- to 5-column cell array, one row per 'key'. The first
%               column contains the key string, the second its type,
%               the third the accepted value range, and the fourth the
%               default value.  Allowed types are 'boolean', 'integer',
%               'real', 'string', 'cell' or 'struct'.  For example,
%                       {'key1' 'string' { 'string1' 'string2' } 'defaultval_key1'}
%                       {'key2' 'int' { minint maxint } 'defaultval_key2'}
%  callingfunc - Calling function name for error messages. {default: none}.
%  mode        - ['ignore'|'error'] ignore keywords that are either not specified
%                in the fieldlist cell array or generate an error.
%                {default: 'error'}.
% Outputs:
%   result     - If no error, structure with 'key' as fields and 'value' as
%                content. If error this output contain the string error.
%   varargin   - residual varagin containing unrecognized input arguments.
%                Requires mode 'ignore' above.
%
% Note: In case of error, a string is returned containing the error message
%       instead of a structure.
%
% Example:
%	g = finputcheck(varargin, ...
%               { 'title'         'string'   []       ''; ...
%                 'percent'       'real'     [0 1]    1 ; ...
%                 'elecamp'       'integer'  [1:10]   [] });
% Note:
%   The 'title' argument should be a string. {no default value}
%   The 'percent' argument should be a real number between 0 and 1. {default: 1}
%   The 'elecamp' argument should be an integer between 1 and 10 (inclusive).
%
%   Now 'g.title' will contain the title arg (if any, else the default ''), etc.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: finputcheck.m,v $
% Revision 1.24  2006/03/11 05:37:07  arno
% header
%
% Revision 1.23  2004/11/06 02:54:06  scott
% a few further small edits to the help msg -sm
%
% Revision 1.22  2004/11/05 15:23:37  arno
% ,sg
%
% Revision 1.21  2004/11/05 04:10:44  scott
% help msg. -sm
%
% Revision 1.20  2004/06/09 16:30:42  arno
% adding or if several types
%
% Revision 1.19  2003/10/29 16:35:57  arno
% msg typo
%
% Revision 1.18  2003/07/30 23:53:58  arno
% debug multiple return values
%
% Revision 1.17  2003/07/26 00:21:17  arno
% allowing cell array for values
%
% Revision 1.16  2003/06/30 02:10:10  arno
% strmatch exact
%
% Revision 1.15  2003/01/31 02:35:38  arno
% debugging lowercase/upercase problem
%
% Revision 1.14  2002/11/20 01:05:44  arno
% take into account duplicate parameters
%
% Revision 1.13  2002/11/18 17:15:18  arno
% adding float arg (=real)
%
% Revision 1.12  2002/11/15 02:16:50  arno
% header for web
%
% Revision 1.11  2002/09/30 15:29:23  arno
% autorizing cell arrays for types
%
% Revision 1.10  2002/09/30 00:42:08  arno
% debug input arguments
%
% Revision 1.9  2002/07/29 18:00:53  arno
% debugging for NaN
%
% Revision 1.8  2002/07/29 17:24:22  arno
% header
%
% Revision 1.7  2002/07/20 19:10:41  arno
% debugging output
%
% Revision 1.6  2002/07/19 17:58:11  arno
% returning non-matched 'key' 'val' arguments
%
% Revision 1.5  2002/07/19 17:46:53  arno
% g empty if no varargin
%
% Revision 1.4  2002/07/19 16:27:14  arno
% adding ignore mode
%
% Revision 1.3  2002/07/10 02:18:32  arno
% header info
%
% Revision 1.2  2002/07/10 02:17:27  arno
% debugging error message passing
%
% Revision 1.1  2002/07/10 01:03:19  arno
% Initial revision
%

function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode )

if nargin < 2
    help finputcheck;
    return;
end;
if nargin < 3
    callfunc = '';
else
    callfunc = [callfunc ' ' ];
end;
if nargin < 4
    mode = 'do not ignore';
end;
NAME = 1;
TYPE = 2;
VALS = 3;
DEF  = 4;
SIZE = 5;

varargnew = {};
% create structure
% ----------------
if ~isempty(vararg)
    vararg = removedup(vararg);
    for index=1:length(vararg)
        if iscell(vararg{index})
            vararg{index} = {vararg{index}};
        end;
    end;
    try
        g = struct(vararg{:});
    catch
        g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
    end;
else
    g = [];
end;

for index = 1:size(fieldlist,NAME)
    % check if present
    % ----------------
    if ~isfield(g, fieldlist{index, NAME})
        g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
    end;
    tmpval = getfield( g, {1}, fieldlist{index, NAME});
    
    % check type
    % ----------
    if ~iscell( fieldlist{index, TYPE} )
        res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
            fieldlist{index, VALS}, tmpval, callfunc );
        if isstr(res), g = res; return; end;
    else
        testres = 0;
        tmplist = fieldlist;
        for it = 1:length( fieldlist{index, TYPE} )
            if ~iscell(fieldlist{index, VALS})
                res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                    fieldlist{index, VALS}, tmpval, callfunc );
            else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                    fieldlist{index, VALS}{it}, tmpval, callfunc );
            end;
            if ~isstr(res{it}), testres = 1; end;
        end;
        if testres == 0,
            g = res{1};
            for tmpi = 2:length(res)
                g = [ g 10 'or ' res{tmpi} ];
            end;
            return;
        end;
    end;
end;

% check if fields are defined
% ---------------------------
allfields = fieldnames(g);
for index=1:length(allfields)
    if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
        if ~strcmpi(mode, 'ignore')
            g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
        end;
        varargnew{end+1} = allfields{index};
        varargnew{end+1} = getfield(g, {1}, allfields{index});
    end;
end;


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
NAME = 1;
TYPE = 2;
VALS = 3;
DEF  = 4;
SIZE = 5;
g = [];

switch fieldtype
    case { 'integer' 'real' 'boolean' 'float' },
        if ~isnumeric(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
        end;
        if strcmpi(fieldtype, 'boolean')
            if tmpval ~=0 & tmpval ~= 1
                g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
            end;
        else
            if strcmpi(fieldtype, 'integer')
                if ~isempty(fieldval)
                    if (isnan(tmpval) & ~any(isnan(fieldval))) ...
                            & (~ismember(tmpval, fieldval))
                        g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                    end;
                end;
            else % real or float
                if ~isempty(fieldval)
                    if tmpval < fieldval(1) | tmpval > fieldval(2)
                        g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                    end;
                end;
            end;
        end;
        
        
    case 'string'
        if ~isstr(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
        end;
        if ~isempty(fieldval)
            if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
                g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
            end;
        end;
        
        
    case 'cell'
        if ~iscell(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
        end;
        
        
    case 'struct'
        if ~isstruct(tmpval)
            g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
        end;
        
        
    case '';
    otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
[tmp indices] = unique(cella(1:2:end));
if length(tmp) ~= length(cella)/2
    fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
end;
cella = cella(sort(union(indices*2-1, indices*2)));
% -----------------------------------
% loadtxt() - load ascii text file into numeric or cell arrays
%
% Usage:
%   >> array = loadtxt( filename, 'key', 'val' ...);
%
% Inputs:
%    filename - name of the input file
%
% Optional inputs
%   'skipline' - number of lines to skip {default:0}. If this number is
%                negative the program will only skip non-empty lines
%                (can be usefull for files transmitted from one platform
%                to an other, as CR may be inserted at every lines).
%   'convert'  - 'on' standard text conversion, see note 1
%                'off' no conversion, considers text only
%                'force' force conversion, NaN are returned
%                for non-numeric inputs {default:'on'}
%   'delim'    - ascii character for delimiters. {default:[9 32]
%                i.e space and tab}. It is also possible to enter
%                strings, Ex: [9 ' ' ','].
%   'verbose'  - ['on'|'off'] {default:'on'}
%   'nlines'   - [integer] number of lines to read {default: all file}
%
% Outputs:
%    array - cell array. If the option 'force' is given, the function
%            retrun a numeric array.
%
% Notes: 1) Since it uses cell arrays, the function can handle text input.
%        The function reads each token and then try to convert it to a
%        number. If the conversion is unsucessfull, the string itself
%        is included in the array.
%        2) The function adds empty entries for rows that contains
%        fewer columns than others.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 29 March 2002

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 29 March 2002
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: loadtxt.m,v $
% Revision 1.10  2005/05/24 17:51:16  arno
% remove cell2mat
%
% Revision 1.9  2005/02/14 00:20:57  arno
% fix delim
%
% Revision 1.8  2005/02/04 00:38:50  arno
% upgrading to use finputcheck
%
% Revision 1.7  2004/07/29 21:06:17  arno
% coma separated file debug
%
% Revision 1.6  2004/01/29 01:35:35  arno
% debug numerical conversion
%
% Revision 1.5  2003/11/19 19:28:16  arno
% now reading empty tabs
%
% Revision 1.4  2003/11/07 01:45:32  arno
% adding nline argument
%
% Revision 1.3  2003/01/10 17:28:56  arno
% debug last
%
% Revision 1.2  2003/01/10 17:27:13  arno
% str2num -> str2double
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function array = loadtxt( filename, varargin );

if nargin < 1
    help loadtxt;
    return;
end;
if ~isempty(varargin)
    try, g = struct(varargin{:});
    catch, disp('Wrong syntax in function arguments'); return; end;
else
    g = [];
end;

g = finputcheck( varargin, { 'convert'   'string'   { 'on' 'off' }   'on';
    'skipline'  'integer'  [0 Inf]          0;
    'verbose'   'string'   { 'on' 'off' }   'on';
    'delim'     { 'integer' 'string' } []               [9 32];
    'nlines'    'integer'  []               Inf });
if isstr(g), error(g); end;
g.convert = lower(g.convert);
g.verbose = lower(g.verbose);
g.delim = char(g.delim);

% open the file
% -------------
if exist(filename) ~=2, error( ['file ' filename ' not found'] ); end;
fid=fopen(filename,'r','ieee-le');
if fid<0, error( ['file ' filename ' found but error while opening file'] ); end;

index = 0;
while index < abs(g.skipline)
    tmpline = fgetl(fid);
    if g.skipline > 0 | ~isempty(tmpline)
        index = index + 1;
    end;
end; % skip lines ---------

inputline = fgetl(fid);
linenb = 1;
if strcmp(g.verbose, 'on'), fprintf('Reading file (lines): '); end;
while isempty(inputline) | inputline~=-1
    colnb = 1;
    if ~isempty(inputline)
        switch g.convert
            case 'off',
                while ~isempty(deblank(inputline))
                    % 07/29/04 Petr Janata added following line to
                    % mitigate problem of strtok ignoring leading
                    % delimiters and deblanking residue in the event
                    % of only space existing between delimiters
                    inputline = strrep(inputline,[g.delim g.delim],[g.delim ' ' g.delim]);
                    
                    [array{linenb, colnb} inputline] = strtok(inputline, g.delim);
                    colnb = colnb+1;
                end;
            case 'on',
                while ~isempty(deblank(inputline))
                    [tmp inputline] = mystrtok(inputline, g.delim);
                    if ~isempty(tmp) & tmp(1) > 43 & tmp(1) < 59, tmp2 = str2num(tmp);
                    else tmp2 = []; end;
                    if isempty( tmp2 )  , array{linenb, colnb} = tmp;
                    else                  array{linenb, colnb} = tmp2;
                    end;
                    colnb = colnb+1;
                end;
            case 'force',
                while ~isempty(deblank(inputline))
                    [tmp inputline] = mystrtok(inputline, g.delim);
                    array{linenb, colnb} = str2double( tmp );
                    colnb = colnb+1;
                end;
            otherwise, error('Unrecognised converting option');
        end;
        linenb = linenb +1;
    end;
    inputline = fgetl(fid);
    if linenb > g.nlines
        inputline = -1;
    end;
    if ~mod(linenb,10) & strcmp(g.verbose, 'on'), fprintf('%d ', linenb); end;
end;
if strcmp(g.verbose, 'on'),  fprintf('%d\n', linenb-1); end;
if strcmp(g.convert, 'force'), array = [ array{:} ]; end;
fclose(fid);

% problem strtok do not consider tabulation
% -----------------------------------------
function [str, strout] = mystrtok(strin, delim);
if delim == 9 % tab
    if length(strin) > 1 & strin(1) == 9 & strin(2) == 9
        str = '';
        strout = strin(2:end);
    else
        [str, strout] = strtok(strin, delim);
    end;
else
    [str, strout] = strtok(strin, delim);
end;

% convertlocs() - Convert electrode locations between coordinate systems
%                 using the EEG.chanlocs structure.
%
% Usage: >> newchans = convertlocs( EEG, 'command');
%
% Input:
%   chanlocs  - An EEGLAB EEG dataset OR a EEG.chanlocs channel locations structure
%   'command' - ['cart2topo'|'sph2topo'|'sphbesa2topo'| 'sph2cart'|'topo2cart'|'sphbesa2cart'|
%               'cart2sph'|'sphbesa2sph'|'topo2sph'| 'cart2sphbesa'|'sph2sphbesa'|'topo2sphbesa'|
%               'cart2all'|'sph2all'|'sphbesa2all'|'topo2all']
%                These command modes convert between four coordinate frames: 3-D Cartesian 
%                (cart), Matlab spherical (sph), Besa spherical (sphbesa), and 2-D polar (topo)
%               'auto' -- Here, the function finds the most complex coordinate frame 
%                 and constrains all the others to this one. It searches first for Cartesian 
%                 coordinates, then for spherical and finally for polar. Default is 'auto'.
%
% Optional input
%   'verbose' - ['on'|'off'] default is 'off'.
%
% Outputs:
%   newchans - new EEGLAB channel locations structure
%
% Ex:  CHANSTRUCT = convertlocs( CHANSTRUCT, 'cart2topo');
%      % Convert Cartesian coordinates to 2-D polar (topographic). 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002
%
% See also: readlocs()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function chans = convertlocs(chans, command, varargin);

if nargin < 1
   help convertlocs;
   return;
end;

if nargin < 2
   command = 'auto';
end;
if nargin == 4 && strcmpi(varargin{2}, 'on')
    verbose = 1;
else
    verbose = 0; % off
end;

% test if value exists for default
% --------------------------------
if strcmp(command, 'auto')
    if isfield(chans, 'X') && ~isempty(chans(1).X)
        command = 'cart2all';
        if verbose
            disp('Make all coordinate frames uniform using Cartesian coords');
        end;
    else
        if isfield(chans, 'sph_theta') && ~isempty(chans(1).sph_theta)
            command = 'sph2all';
            if verbose
                disp('Make all coordinate frames uniform using spherical coords');
            end;
        else
            if isfield(chans, 'sph_theta_besa') && ~isempty(chans(1).sph_theta_besa)
                command = 'sphbesa2all';
                if verbose
                    disp('Make all coordinate frames uniform using BESA spherical coords');
                end;
            else
                command = 'topo2all';
                if verbose
                    disp('Make all coordinate frames uniform using polar coords');
                end;
            end;
        end;
    end;
end;

% convert
% -------         
switch command
 case 'topo2sph',
   theta  = {chans.theta};
   radius = {chans.radius};
   indices = find(~cellfun('isempty', theta));
   [sph_phi sph_theta] = topo2sph( [ [ theta{indices} ]' [ radius{indices}]' ] );
   if verbose
       disp('Warning: electrodes forced to lie on a sphere for polar to 3-D conversion');
   end;
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);
   end;
   if isfield(chans, 'sph_radius'),
       meanrad = mean([ chans(indices).sph_radius ]);
       if isempty(meanrad), meanrad = 1; end;
   else
       meanrad = 1;
   end;
   sph_radius(1:length(indices)) = {meanrad};
case 'topo2sphbesa',
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'topo2cart'
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   if verbose
       disp('Warning: spherical coordinates automatically updated');
   end;
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'topo2all',
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sph2cart',
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   if ~isfield(chans, 'sph_radius'), sph_radius(1:length(indices)) = {1};
   else                              sph_radius = {chans.sph_radius};
   end;
   inde = find(cellfun('isempty', sph_radius));
   if ~isempty(inde)
       meanrad = mean( [ sph_radius{:} ]);
       sph_radius(inde) = { meanrad };
   end;
   [x y z] = sph2cart([ sph_theta{indices} ]'/180*pi, [ sph_phi{indices} ]'/180*pi, [ sph_radius{indices} ]');
   for index = 1:length(indices)
      chans(indices(index)).X = x(index);
      chans(indices(index)).Y = y(index);
      chans(indices(index)).Z = z(index);
   end;
case 'sph2topo',
 if verbose
     % disp('Warning: all radii constrained to one for spherical to topo transformation');
 end;
 sph_theta  = {chans.sph_theta};
 sph_phi    = {chans.sph_phi};
 indices = find(~cellfun('isempty', sph_theta));
 [chan_num,angle,radius] = sph2topo([ ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2); % using method 2
 for index = 1:length(indices)
     chans(indices(index)).theta  = angle(index);
     chans(indices(index)).radius = radius(index);
     if ~isfield(chans, 'sph_radius') || isempty(chans(indices(index)).sph_radius)
         chans(indices(index)).sph_radius = 1;
     end;
 end;
case 'sph2sphbesa',
   % using polar coordinates
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2);
   [sph_theta_besa sph_phi_besa] = topo2sph([angle radius], 1, 1);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta_besa  = sph_theta_besa(index);
      chans(indices(index)).sph_phi_besa    = sph_phi_besa(index);
   end;   
case 'sph2all',
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sphbesa2sph',
   % using polar coordinates
   sph_theta_besa  = {chans.sph_theta_besa};
   sph_phi_besa    = {chans.sph_phi_besa};
   indices = find(~cellfun('isempty', sph_theta_besa));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_theta_besa{indices} ]' [ sph_phi_besa{indices} ]' ], 1, 1);
   %for index = 1:length(chans)
   %   chans(indices(index)).theta  = angle(index);
   %   chans(indices(index)).radius = radius(index);
   %   chans(indices(index)).labels = int2str(index);
   %end;   
   %figure; topoplot([],chans, 'style', 'blank', 'electrodes', 'labelpoint');
   
   [sph_phi sph_theta] = topo2sph([angle radius], 2);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);      
   end;
case 'sphbesa2topo',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'sphbesa2cart',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords   
case 'sphbesa2all',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
case 'cart2topo',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'cart2sphbesa',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'cart2sph',
    if verbose
        disp('WARNING: If XYZ center has not been optimized, optimize it using Edit > Channel Locations');
	end;
    X  = {chans.X};
    Y  = {chans.Y};
    Z  = {chans.Z};
    indices = find(~cellfun('isempty', X));
    [th phi radius] = cart2sph( [ X{indices} ], [ Y{indices} ], [ Z{indices} ]);
	for index = 1:length(indices)
		 chans(indices(index)).sph_theta     = th(index)/pi*180;
		 chans(indices(index)).sph_phi       = phi(index)/pi*180;
		 chans(indices(index)).sph_radius    = radius(index);
	end;
case 'cart2all',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
end;

% sph2topo() - Convert from a 3-column headplot file in spherical coordinates
%              to 3-column topoplot() locs file in polar (not cylindrical) coords.
%              Used for topoplot() and other 2-D topographic plotting programs.
%              Assumes a spherical coordinate system in which horizontal angles 
%              have a range [-180,180] deg,  with zero pointing to the right ear. 
%              In the output polar coordinate system, zero points to the nose.
%              See  >> help readlocs
% Usage:
%          >> [chan_num,angle,radius] = sph2topo(input,shrink_factor,method);
%
% Inputs:
%   input         = [channo,az,horiz] = chan_number, azumith (deg), horiz. angle (deg)
%                   When az>0, horiz=0 -> right ear, 90 -> nose 
%                   When az<0, horiz=0 -> left ear, -90 -> nose
%   shrink_factor = arc_length shrinking factor>=1 (deprecated).
%                   1 -> plot edge is 90 deg azimuth {default};
%                   1.5 -> plot edge is +/-135 deg azimuth See 
%                   >> help topoplot(). 
%   method        = [1|2], optional. 1 is for Besa compatibility, 2 is for
%                   compatibility with Matlab function cart2sph(). Default is 2
%
% Outputs:
%   channo  = channel number (as in input)
%   angle   = horizontal angle (0 -> nose; 90 -> right ear; -90 -> left ear)
%   radius  = arc_lengrh from vertex (Note: 90 deg az -> 0.5/shrink_factor);
%             By topoplot() convention, radius=0.5 is the nasion-ear_canal plane.
%             Use topoplot() 'plotrad' to plot chans with abs(az) > 90 deg.
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 6/12/98 
%
% See also: cart2topo(), topo2sph()

% Copyright (C) 6/12/98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% corrected left/right orientation mismatch, Blair Hicks 6/20/98
% changed name sph2pol() -> sph2topo() for compatibility -sm
% 01-25-02 reformated help & license -ad 
% 01-25-02 changed computation so that it works with sph2topo -ad 

function [channo,angle,radius] = sph2topo(input,factor, method)

chans = size(input,1);
angle = zeros(chans,1);
radius = zeros(chans,1);

if nargin < 1
   help sph2topo
   return
end
   
if nargin< 2
  factor = 0;
end
if factor==0
  factor = 1;
end
if factor < 1
  help sph2topo
  return
end

if size(input,2) ~= 3
   help sph2topo
   return
end

channo = input(:,1);
az = input(:,2);
horiz = input(:,3);

if exist('method')== 1 & method == 1
  radius = abs(az/180)/factor;
  i = find(az>=0);
  angle(i) = 90-horiz(i);
  i = find(az<0);
  angle(i) = -90-horiz(i);
else
  angle  = -horiz;
  radius = 0.5 - az/180;
end;

% topo2sph() - convert a topoplot() style 2-D polar-coordinate
%              channel locations file to a 3-D spherical-angle
%              file for use with headplot()
% Usage: 
%   >> [c h] = topo2sph('eloc_file','eloc_outfile', method, unshrink);
%   >> [c h] = topo2sph( topoarray, method, unshrink );
%
% Inputs:
%   'eloc_file'    = filename of polar 2-D electrode locations file used by 
%                    topoplot(). See >> topoplot example or cart2topo()
%   'eloc_outfile' = output file of 3-D electrode locations in spherical angle 
%                    coords. for use in headplot().
%   topoarray      = polar array of 2-D electrode locations, with polar angle
%                    in the first column and radius in the second one.
%   method         = [1|2] 1 is for Besa compatibility, 2 is for
%                    compatibility with Matlab function cart2sph(). {default: 2}
%   unshrink       = [0<real<1] unshrink factor. Enter a shrink factor used
%                    to convert spherical to topo (see sph2topo()). Only 
%                    implemented for 'method' 1 (above). Electrode 'shrink' 
%                    is now deprecated. See >> help topoplot
% Outputs:
%   c = coronal rotation
%   h = horizontal rotation
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 
%
% See also: sph2topo(), cart2topo()

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 3-16-00 changed name to topo2sph() for compatibility with cart2topo() -sm
% 01-25-02 reformated help & license -ad 
% 03-22-02 complete remodeling for returning arguments and taking arrays -ad 

function [c, h] = topo2sph(eloc_locs,eloc_angles, method, unshrink)

MAXCHANS = 1024;

if nargin < 1
    help topo2sph;
    return;
end;
if nargin > 1 && ~isstr(eloc_angles)
	if nargin > 2
		unshrink = method;
	end;
	method = eloc_angles;
else
	method = 2;
end;

if isstr(eloc_locs)
	fid = fopen(eloc_locs);
	if fid<1,
	    fprintf('topo2sph()^G: cannot open eloc_loc file (%s)\n',eloc_locs)
	    return
	end
	E = fscanf(fid,'%d %f %f  %s',[7 MAXCHANS]);
	E = E';
	fclose(fid);
else
    E = eloc_locs;
    E = [ ones(size(E,1),1) E ];
end;
    
if nargin > 1 & isstr(eloc_angles)
	if exist(eloc_angles)==2,
	   fprintf('topo2sph: eloc_angles file (%s) already exists and will be erased.\n',eloc_angles);
	end

	fid = fopen(eloc_angles,'a');
	if fid<1,
	    fprintf('topo2sph()^G: cannot open eloc_angles file (%s)\n',eloc_angles)
	    return
	end
end;

if method == 2
	t = E(:,2); % theta
	r = E(:,3); % radius
	h = -t;  % horizontal rotation
	c = (0.5-r)*180;
else
	for e=1:size(E,1)
		% (t,r) -> (c,h)
		
		t = E(e,2); % theta
		r = E(e,3); % radius
		r = r*unshrink;
		if t>=0
			h(e) = 90-t; % horizontal rotation
		else
			h(e) = -(90+t);
		end
		if t~=0
			c(e) = sign(t)*180*r; % coronal rotation
		else
			c(e) = 180*r;
		end
	end;
	t = t';
	r = r';
end;

for e=1:size(E,1)
   if nargin > 1 & isstr(eloc_angles)
        chan = E(e,4:7);
        fprintf('%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
        fprintf(fid,'%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
   end;     
end

% readeetraklocs() - read 3-D location files saved using the EETrak
%                    digitizing software.
% Usage:
%   >> CHANLOCS = readeetraklocs( filename );
%
% Inputs:
%   filename       - [string] file name
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. 
%                    See help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, Nov 2003
%
% See also: readlocs()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function chanlocs = readeetraklocs( filename )
    
    if nargin < 1
        help readeetraklocs;
        return;
    end;
    
    % read location file
    % ------------------
    locs  = loadtxt( filename );
        
    % get label names
    % ---------------
    indlabels = [];
    indpos    = [];
    for ind = 1:size(locs,1)
        if isstr(locs{ind,1}) 
            if strcmpi(locs{ind,1}, 'Labels')
                indlabels = ind;
            end;
            if strcmpi(locs{ind,1}, 'Positions')
                indpos = ind;
            end;
        end;
    end;
    if isempty(indpos) | isempty(indlabels)
        error('Could not find ''Labels'' or ''Position'' tag in electrode file');
    end;
    
    % get positions
    % -------------
    positions = locs(indpos+1:indlabels-1,1:3);
    labels    = locs(indlabels+1:end,:);
        
    % create structure
    % ----------------
    for index = 1:length(labels)
        chanlocs(index).labels = labels{index};
        chanlocs(index).X      = positions{index,1};
        chanlocs(index).Y      = positions{index,2};
        chanlocs(index).Z      = positions{index,3};
    end;
        
    chanlocs = convertlocs(chanlocs, 'cart2all');

    
    % readneurolocs() - read neuroscan polar location file (.asc)
%
% Usage:
%   >> CHANLOCS = readneurolocs( filename );
%   >> CHANLOCS = readneurolocs( filename, 'key1', val1, ...);
%
% Inputs:
%   filename       - file name or matlab cell array { names x_coord y_coord }
%
% Optional inputs:
%   same as caliblocs()
%   note that if no optional input are provided, re-centering will be
%   performed automatically and re-scaling of coordinates will be
%   performed for '.asc' files (not '.dat' files).
%
% Outputs:
%   CHANLOCS       - EEGLAB channel location data structure. See
%                    help readlocs()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 4 March 2003
%
% See also: readlocs()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function chanlocs = readneurolocs( filename, varargin)

if nargin < 1
    help readneurolocs;
    return;
end;
if nargin < 2
    plottag = 0;
end;

% read location file
% ------------------
if isstr(filename)
    locs  = loadtxt( filename );
end;

if ~isstr(filename) || locs{1,1}(1) == ';' || size(locs,2) < 5
    if ~isstr(filename)
        names = filename{1};
        x     = filename{2};
        y     = filename{3};
    else
        if locs{1,1}(1) == ';'
            % remove trailing control channels
            % --------------------------------
            while isnumeric( locs{end,1} ) & locs{end,1} ~= 0
                locs  = locs(1:end-1,:);
            end;

            % find first numerical index
            % --------------------------
            index = 1;
            while isstr( locs{index,1} )
                index = index + 1;
            end;

            % extract location array
            % ----------------------
            nchans = size( locs, 1 ) - index +1;
            chans  = [locs{end-nchans+1:end, 1:5}];
            chans  = reshape(chans,nchans,5);               %% Added this line in order to get x = chans(:,3)
            names  = locs(end-nchans*2+1: end-nchans, 2);
            for index = 1:length(names)
                if ~isstr(names{index})
                    names{index} = int2str(names{index});
                end;
            end;
            x = chans(:,3);
            y = -chans(:,4);
        else
            [tmp2 tmpind] = sort( [ locs{:,1} ]);
            locs = locs(tmpind,:);
            y      = [ locs{:,end} ];
            x      = [ locs{:,end-1} ];
            x      = x/513.1617*44;
            y      = y/513.1617*44;
            names = locs(:,end-2);
        end;
    end;

    % second solution using angle
    % ---------------------------
    [phi,theta] = cart2pol(x, y);
    phi = phi/pi*180;

    % convert to other types of coordinates
    % -------------------------------------
    labels = names';
    chanlocs = struct('labels', labels, 'sph_theta_besa', mattocell(theta)', 'sph_phi_besa', mattocell(phi)');      %% labels instead of labels(:) 
    chanlocs = convertlocs( chanlocs, 'sphbesa2all');

    for index = 1:length(chanlocs)
        chanlocs(index).labels = num2str(chanlocs(index).labels);
    end;

    % re-calibration
    % --------------
    chanlocs = adjustlocs(chanlocs, 'autoscale', 'on', 'autorotate', 'off', varargin{:});

else % 5 rows, xyz positions
    try
        for index = 1:size(locs,1)
            locs{index,3} = - locs{index,3};
        end;
        chanlocs = struct('labels', locs(:,1), 'type', locs(:,2), 'X', locs(:,4), 'Y', locs(:,3), 'Z', locs(:,5));
        chanlocs = convertlocs( chanlocs, 'cart2all');
    catch
        chanlocs = readlocs(filename, 'filetype', 'custom', 'format', { 'labels' 'ignore' '-Y' 'X' 'Z' });
    end;
end;
    
