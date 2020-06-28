function h1 = wb_NetPlot(M,model,flag1,flag2,nodename,nodecolor,NodenameFontSize,nodeMakerSize)
% Plot 2D networks in an unit circle (raius = 1).
% Input:
%      M:Connection matrix. It should be symmetric matrix or upper
%        triangular matrix.
%      model: model = 1(default): The line width is corresponding to absolute.
%                   edge weights.Line width range from 1 to 3.
%                   This is suitable for less number of nodes.
%             model = 2: The line width is corresponding to absolute edge 
%                   weights,The red line is positive connection, and blue 
%                   line is negative connection.This is suitable for less 
%                   number of nodes (containing positive and negative connection).
%                   Line width range from 1 to 3.
%             model = 3: The color of edge is corresponding to edge
%                   weights (min -> max), and the line width is fixed at
%                   1.This is suitable for large nodes and edges.
%      flag1: flag1 = 0 (default): Undirection graph. M is symmetric
%             matrix or upper triangular matrix.
%             flag1 = 1: direction graph. from i to j.
%      flag2: flag2 = 1 (default): show the node names.
%             flag2 = 0: do not show the node names.
%     nodename: the node names. They are cell strings generated by
%               {'nodename1','nodename2',...'nodenameN'}.
%               Default is the number of nodes.
%     nodecolor: the colors of nodes. row is No. of nodes, colums are RGB
%                values used in the matlab.
%     NodenameFontSize: the Font size of node names. Default is 12.
%     nodeMakerSize:   the node Maker size. Default is 8.
%                      If user inputs various maker sizes, it ranges from 6 to
%                      20, and the No. of nodeMakerSize must be equal to No. of
%                      nodes.
% Output:
%        h1: handle of figure.
% -------------------------------------------------------------------------
% Example: 
% clear all;clc
% M = [0 1 3 0 -2 0;
%      0 0 1 0 0 -1;
%      0 0 0 0 4 2;
%      0 0 0 0 0 0;
%      0 0 0 0 0 1;
%      0 0 0 0 0 0];
% nodename = char('node1','node2','node3','node4','node5','node6');
% nodecolor = [0,0.75,1;
%              0,0.75,1;
%              0,0.75,1;
%              0,0.75,1;
%              1,0.188,0;
%              1,0.188,0;];
% Dong_NetPlot(M);
% Dong_NetPlot(M,1,0,1);
% Dong_NetPlot(M,2,1,1,nodename,nodecolor);
% If you want to save the image, please use below code:
% ---------------------
% frame = getframe(h1);
% im = frame2im(frame);
% imwrite(im,'~\example.tiff','tiff');
% OR
% set(h1,'inverthardcopy','off');
% print(h1,'-dtiff','-r300','~\example.tiff');
% -------------------------------------------------------------------------
% Written by Li Dong,UESTC: Li_dong729@163.com.
% Date: 2015/5/19.
% -------------------------------------------------------------------------
if nargin == 1
    model = 1;
    flag1 = 0;
    flag2 = 1;
    nodename = [];
    nodecolor = [];
    NodenameFontSize = 12;           % the font size of node name.
    nodeMakerSize = 8;               % the maker size of nodes.
elseif nargin == 2
    flag1 = 0;
    flag2 = 1;
    nodename = [];
    nodecolor = [];
    NodenameFontSize = 12;           % the font size of node name.
    nodeMakerSize = 8;               % the maker size of nodes.
elseif nargin == 3
    flag2 = 1;
    nodename = [];
    nodecolor = [];
    NodenameFontSize = 12;           % the font size of node name.
    nodeMakerSize = 8;               % the maker size of nodes.
elseif nargin == 4
    nodename = [];
    nodecolor = [];
    NodenameFontSize = 12;           % the font size of node name.
    nodeMakerSize = 8;               % the maker size of nodes.
elseif nargin == 5
    nodecolor = [];
    NodenameFontSize = 12;           % the font size of node name.
    nodeMakerSize = 8;               % the maker size of nodes.
elseif nargin == 6
    NodenameFontSize = 12;           % the font size of node name.
    nodeMakerSize = 8;               % the maker size of nodes.
elseif nargin == 7
    nodeMakerSize = 8;               % the maker size of nodes.
end

N_node = size(M,1);       % number of nodes.
if any(~isfinite(M))
    M(~isfinite(M)) = 0;  % Toss NaN.
    warning('NaN is contained in your connection matrix!!!');
end
M(eye(N_node)==1)=0; % diagonal is 0.
if flag1 == 0
    M1 = triu(ones(N_node,N_node),1);
    M = M1 .* M;     % extract upper triangular elements.
end
max_weight = max(abs(M(:)));     % maximal values in the M.
min_weight = min(abs(M(:)));     % min values in the M.
Weight_edge = (abs(M) - min_weight(1))/(max_weight(1)-min_weight(1))*2+1; % Line width of edge weights.range from 1 to 3.
if length(nodeMakerSize) == N_node
   maxMakerSize = 20;  % max maker size
   minMakerSize = 6;   % min maker size
   N_Interpolate = 64; % No. of interpolation.
   MakerSizeRange = linspace(minMakerSize,maxMakerSize,N_Interpolate);
   min1 = min(nodeMakerSize);
   max1 = max(nodeMakerSize);
end
r1 = 0.6;                        % distance between origin and middle point of arc.
step1 = 2*pi/N_node;             % step size of degree.
theta = step1:step1:2*pi;        % degrees.
X = cos(theta);                  % x coordinates.
Y = sin(theta);                  % y coordinates.
NodenameColor = 'c';             % The color of node names.
switch model
    case 1
        arrowHeadSize = 0.25;             % the size of arrow head.
    case 2
        arrowHeadSize = 0.25;             % the size of arrow head.
    case 3
        arrowHeadSize = 0.2;             % the size of arrow head.
end
h1 = figure;
xlim([min(X)-0.36,max(X)+0.36]);
ylim([min(Y)-0.36,max(Y)+0.36]);
% -------------------------------------------------------------------------
% plot nodes
hold on;
if isempty(nodecolor) % use default color?
    if length(nodeMakerSize) == 1 % use default nodeMakerSize?
        plot(X,Y,'o','MarkerSize',nodeMakerSize,'MarkerFaceColor',NodenameColor,'MarkerEdgeColor','w');
    else
        for i = 1:N_node
            k1 = round(N_Interpolate*(nodeMakerSize(i)-min1(1))/(max1(1) - min1(1)));
            plot(X(i),Y(i),'o','MarkerSize',MakerSizeRange(max(k1,1)),'MarkerFaceColor',NodenameColor,'MarkerEdgeColor','w');
        end
    end
else
    if length(nodeMakerSize) == 1 % use default nodeMakerSize?
        for i = 1:N_node
            plot(X(i),Y(i),'o','MarkerSize',nodeMakerSize,'MarkerFaceColor',nodecolor(i,:),'MarkerEdgeColor','w');
        end
    else
        for i = 1:N_node
            k1 = round(N_Interpolate*(nodeMakerSize(i)-min1(1))/(max1(1) - min1(1)));
            plot(X(i),Y(i),'o','MarkerSize',MakerSizeRange(max(k1,1)),'MarkerFaceColor',nodecolor(i,:),'MarkerEdgeColor','w');
        end
    end
end

if flag2 == 1 % show the node names?
    if isempty(nodename)
        if isempty(nodecolor)
            for i = 1:N_node
                if theta(i)>pi/2 && theta(i)< 3*pi/2
                    rot = 180*theta(i)/pi - 180;
                    text(sqrt(1.21)*X(i),sqrt(1.21)*Y(i),num2str(i),'Fontsize',NodenameFontSize,...
                        'color',NodenameColor,'HorizontalAlignment','right','rotation',rot);
                else
                    rot = 180*theta(i)/pi;
                    text(sqrt(1.21)*X(i),sqrt(1.21)*Y(i),num2str(i),'Fontsize',NodenameFontSize,...
                        'color',NodenameColor,'HorizontalAlignment','left','rotation',rot);
                end
            end
        else
            for i = 1:N_node
                if theta(i)>pi/2 && theta(i)< 3*pi/2
                    rot = 180*theta(i)/pi - 180;
                    text(sqrt(1.21)*X(i),sqrt(1.21)*Y(i),num2str(i),'Fontsize',NodenameFontSize,...
                        'color',nodecolor(i,:),'HorizontalAlignment','right','rotation',rot);
                else
                    rot = 180*theta(i)/pi;
                    text(sqrt(1.21)*X(i),sqrt(1.21)*Y(i),num2str(i),'Fontsize',NodenameFontSize,...
                        'color',nodecolor(i,:),'HorizontalAlignment','left','rotation',rot);
                end
            end
        end
    else
        if isempty(nodecolor)
            for i = 1:N_node
                 if theta(i)>pi/2 && theta(i)< 3*pi/2
                    rot = 180*theta(i)/pi - 180;
                    text(sqrt(1.21)*X(i),sqrt(1.21)*Y(i),nodename{1,i},'Fontsize',NodenameFontSize,...
                        'color',NodenameColor,'HorizontalAlignment','right','rotation',rot);
                 else
                     rot = 180*theta(i)/pi;
                     text(sqrt(1.21)*X(i),sqrt(1.21)*Y(i),nodename{1,i},'Fontsize',NodenameFontSize,...
                         'color',NodenameColor,'HorizontalAlignment','left','rotation',rot);
                 end
            end
        else
            for i = 1:N_node
                if theta(i)>pi/2 && theta(i)< 3*pi/2
                    rot = 180*theta(i)/pi - 180;
                    text(sqrt(1.44)*X(i),sqrt(1.44)*Y(i),nodename{1,i},'Fontsize',NodenameFontSize,...
                        'color',nodecolor(i,:),'HorizontalAlignment','right','rotation',rot);
                else
                    rot = 180*theta(i)/pi;
                    text(sqrt(1.44)*X(i),sqrt(1.44)*Y(i),nodename{1,i},'Fontsize',NodenameFontSize,...
                        'color',nodecolor(i,:),'HorizontalAlignment','left','rotation',rot);
                end
            end
        end
    end
end
% -------------------------------------------------------------------------
% plot edges
switch model
    case 1 % model 1
        % -----------------------------------------------------------------
        edge_colors = colormap(summer(sum(M(:)~=0)));
        k = 1;                    % counter of edge colors.
        if flag1 == 0
            % undirction graph.
            for i = 1:N_node
                for j = i:N_node
                    if M(i,j) ~= 0
                        x_temp1 = [X(i),sqrt(r1)*(X(i)+X(j))/2,X(j)];
                        y_temp1 = [Y(i),sqrt(r1)*(Y(i)+Y(j))/2,Y(j)];
                        values = spcrv([[x_temp1(1),x_temp1,x_temp1(end)];[y_temp1(1),y_temp1,y_temp1(end)]],3);
                        plot(values(1,:),values(2,:),'-','linewidth',abs(Weight_edge(i,j)),'Color',edge_colors(k,:));
                        k = k+1;
                    end
                end
            end
        else
            % direction gragh
            for i = 1:N_node
                for j = 1:N_node
                    if M(i,j) ~= 0 && i ~= j
                        % lines
                        x_temp1 = [X(i),sqrt(r1)*(X(i)+X(j))/2,X(j)];
                        y_temp1 = [Y(i),sqrt(r1)*(Y(i)+Y(j))/2,Y(j)];
                        values = spcrv([[x_temp1(1),x_temp1,x_temp1(end)];[y_temp1(1),y_temp1,y_temp1(end)]],3);
                        plot(values(1,:),values(2,:),'-','linewidth',Weight_edge(i,j),'Color',edge_colors(k,:));
                        d1 = sqrt((X(i)-X(j)).^2 + (Y(i)-Y(j)).^2);
                        [coor1,coor2,coor3] = dong_draw_arrow([x_temp1(2),y_temp1(2)],...
                            [X(j),Y(j)],arrowHeadSize + 0.6 - 0.3*d1);
                        fill([coor1(1),coor2(1),coor3(1)],[coor1(2),coor2(2),coor3(2)],edge_colors(k,:));% this fills the arrowhead
                        k = k+1;
                    end
                end
            end
        end
    case 2 % model 2
        % -----------------------------------------------------------------
        edge_colors = [0,0.75,1;1,0.2,0];
        if flag1 == 0
            % undirction graph.
            for i = 1:N_node
                for j = i:N_node
                    if M(i,j) ~= 0
                        x_temp1 = [X(i),sqrt(r1)*(X(i)+X(j))/2,X(j)];
                        y_temp1 = [Y(i),sqrt(r1)*(Y(i)+Y(j))/2,Y(j)];
                        values = spcrv([[x_temp1(1),x_temp1,x_temp1(end)];[y_temp1(1),y_temp1,y_temp1(end)]],3);
                        if M(i,j)<0
                            plot(values(1,:),values(2,:),'-','linewidth',Weight_edge(i,j),'Color',edge_colors(1,:));
                        else
                            plot(values(1,:),values(2,:),'-','linewidth',Weight_edge(i,j),'Color',edge_colors(2,:));
                        end
                    end
                end
            end
        else
            % direction gragh
            for i = 1:N_node
                for j = 1:N_node
                    if M(i,j) ~= 0 && i ~= j
                        % lines
                        x_temp1 = [X(i),sqrt(r1)*(X(i)+X(j))/2,X(j)];
                        y_temp1 = [Y(i),sqrt(r1)*(Y(i)+Y(j))/2,Y(j)];
                        values = spcrv([[x_temp1(1),x_temp1,x_temp1(end)];[y_temp1(1),y_temp1,y_temp1(end)]],3);
                        d1 = sqrt((X(i)-X(j)).^2 + (Y(i)-Y(j)).^2);
                        [coor1,coor2,coor3] = dong_draw_arrow([x_temp1(2),y_temp1(2)],...
                            [X(j),Y(j)],arrowHeadSize + 0.6 - 0.3*d1);
                        if M(i,j) < 0
                            plot(values(1,:),values(2,:),'-','linewidth',Weight_edge(i,j),'Color',edge_colors(1,:));
                            fill([coor1(1),coor2(1),coor3(1)],[coor1(2),coor2(2),coor3(2)],edge_colors(1,:));% this fills the arrowhead
                        else
                            plot(values(1,:),values(2,:),'-','linewidth',Weight_edge(i,j),'Color',edge_colors(2,:));
                            fill([coor1(1),coor2(1),coor3(1)],[coor1(2),coor2(2),coor3(2)],edge_colors(2,:));% this fills the arrowhead
                        end
                    end
                end
            end
        end
    case 3 % model 3
        % -----------------------------------------------------------------
        edge_colors = colormap(summer(128));
        weights = M(M~=0);
        minWei = min(weights);
        maxWei = max(weights);
        colorbar;
        caxis([minWei(1),maxWei(1)]);
        if flag1 == 0
            % undirction graph.
            for i = 1:N_node
                for j = i:N_node
                    if M(i,j) ~= 0
                        x_temp1 = [X(i),sqrt(r1)*(X(i)+X(j))/2,X(j)];
                        y_temp1 = [Y(i),sqrt(r1)*(Y(i)+Y(j))/2,Y(j)];
                        values = spcrv([[x_temp1(1),x_temp1,x_temp1(end)];[y_temp1(1),y_temp1,y_temp1(end)]],3);
                        k = round(128*(M(i,j)-minWei(1))/(maxWei(1) - minWei(1)));
                        plot(values(1,:),values(2,:),'-','linewidth',1,'Color',edge_colors(max(k,1),:));
                    end
                end
            end
        else
            % direction gragh
            for i = 1:N_node
                for j = 1:N_node
                    if M(i,j) ~= 0 && i ~= j
                        % lines
                        x_temp1 = [X(i),sqrt(r1)*(X(i)+X(j))/2,X(j)];
                        y_temp1 = [Y(i),sqrt(r1)*(Y(i)+Y(j))/2,Y(j)];
                        values = spcrv([[x_temp1(1),x_temp1,x_temp1(end)];[y_temp1(1),y_temp1,y_temp1(end)]],3);
                        k = round(128*(M(i,j)-minWei(1))/(maxWei(1) - minWei(1)));
                        plot(values(1,:),values(2,:),'-','linewidth',1,'Color',edge_colors(max(k,1),:));
                        
                        [coor1,coor2,coor3] = dong_draw_arrow([x_temp1(2),y_temp1(2)],...
                            [X(j),Y(j)],arrowHeadSize);
                        fill([coor1(1),coor2(1),coor3(1)],[coor1(2),coor2(2),coor3(2)],edge_colors(max(k,1),:));% this fills the arrowhead
                    end
                end
            end
        end
end
axis equal
hold off
axis off
set(h1,'color',[0,0,0]);
% -------------------------------------------------------------------------
% subFunction
function [coor1,coor2,coor3] = dong_draw_arrow(startpoint,endpoint,headsize)
% Calcualte head of 2D arrow from point [x1,y1] to [x2,y2].
% accepts two [x y] coords and one double headsize.
% Input:
%      startpoint: coordinate of start point.
%      endpoint:   coordinate of end point.
%      headsize:   head size of arrow.
% Output:
%     coor1: coordinate of end point.
%     coor2: coordinate of vector 1
%     coor3: coordinate of vector 2
% -------------------------------------------------------------------------
% Written by Li Dong.
% Date: 2015/5/14.
% -------------------------------------------------------------------------
v1 = headsize * (startpoint - endpoint)/2.5;

theta1 = 20 * pi/180;
theta2 = -20 * pi/180;
rotMatrix1 = [cos(theta1),-sin(theta1);sin(theta1),cos(theta1)]; % rotation matrix
rotMatrix2 = [cos(theta2),-sin(theta2);sin(theta2),cos(theta2)]; % rotation matrix

v2 = v1 * rotMatrix1;
v3 = v1 * rotMatrix2;
coor1 = endpoint;
coor2 = coor1 + v2;
coor3 = coor1 + v3;