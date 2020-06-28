function NetMeasure = wb_calculate_BasicNetIndices(M,type)
% Calculate network measures based on  graph theory using BCT toolbox 
% (BCT toolbox: https://sites.google.com/site/bctnet/Home)
% Input:
%    M: binary/weighted/directed/undirected connection matrix.For weighted 
%          network, value of edge may belong to [0,1].
%          directed: from cloumn -> to row (j -> i).
%    type: the type of compatible associated network.
%          'BU': binary undirected network;
%          'BD': binary directed network;
%          'WU': weighted undirected network;
%          'WD': weighted directed network.
% Output:
%    NetMeasure: a structure array which contains network measures.

% NetMeasure.degree: node degree;
% NetMeasure.indegree: node indegree;
% NetMeasure.outdegree: node outdegree;
% NetMeasure.Kn: mean degree of network;
% NetMeasure.Kcost: cost of network;
% NetMeasure.Kn_in: mean indegree of network;
% NetMeasure.Kn_out: mean outdegree of network;

% NetMeasure.strength: node strength;
% NetMeasure.strength_m: mean node strength of network;
% NetMeasure.instrength: node instrength;
% NetMeasure.outstrength: node outstrength;
% NetMeasure.instrength_m: mean node instrength of network;
% NetMeasure.outstrength_m: mean node outstrength of network;

% NetMeasure.Cn: node clustering coefficient;
% NetMeasure.Cn_m: mean clustering coefficient of network;
% NetMeasure.Ln: characteristic path length of network;
% NetMeasure.Eglob: global efficiency of network;
% NetMeasure.Eloc: node local efficiency;
% NetMeasure.Eloc_m: mean local efficiency of network;

% NetMeasure.BC: node betweenness centrality;
% NetMeasure.BC_m: mean node betweenness centrality of network;

% NetMeasure.assort_coef: assortativity coefficient (undirected graph: strength/strength correlation);
% NetMeasure.assort_coef1: assortativity coefficient (directed graph: out-strength/in-strength correlation);
% NetMeasure.assort_coef2: assortativity coefficient (directed graph: in-strength/out-strength correlation);
% NetMeasure.assort_coef3: assortativity coefficient (directed graph: out-strength/out-strength correlation);
% NetMeasure.assort_coef4: assortativity coefficient (directed graph: in-strength/in-strength correlation);

% NetMeasure.rich_club.rich_coef: rich-club coefficients at level(degree)k;1 X levels 
% NetMeasure.rich_club.proportion: proportions (0<p<1) of the strongest weights; 1 X proportions
% NetMeasure.rich_club.rich_coef_proportion: rich-club coefficients of each proportion; proportions X levels
% NetMeasure.rich_club.Nk: number of nodes with degree > k;
% NetMeasure.rich_club.Ek: number of edges remaining in subgraph with degree > k.

% reference:
% Rubinov, M. and O. Sporns (2010). "Complex network measures of brain connectivity
%     : uses and interpretations." Neuroimage 52(3): 1059-1069.
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Li_dong729@163.com)
% $ 2018.1.29
% -------------------------------------------------------------------------
if nargin < 2
    error ('two inputs are reqiured');
end
M(~isfinite(M)) = 0;  % Toss NaN's
N = size(M,1); % number of nodes;
NetMeasure = [];
NetMeasure.nettype = type;

switch type
    case 'BU' % binary undirected
        % degree
        degree = degrees_und(M); % node degree
        Kn = mean(degree); % mean degree of network
        Kcost = sum(degree)/(N*(N-1)); % cost of network
        
        NetMeasure.degree = degree;
        NetMeasure.Kn = Kn;
        NetMeasure.Kcost = Kcost;

        % clustering coefficient
        Cn = clustering_coef_bu(M);
        Cn_m = mean(Cn);
        NetMeasure.Cn = Cn;
        NetMeasure.Cn_m = Cn_m;
        
        % characteristic path length
        D1 = distance_bin(M);
        Ln = charpath(D1);
        NetMeasure.Ln = Ln;
        
        % global efficiency
        Eglob = efficiency_bin(M);
        NetMeasure.Eglob = Eglob;
        
        % local efficiency
        Eloc = efficiency_bin(M,1);
        Eloc_m = mean(Eloc);
        
        NetMeasure.Eloc = Eloc;
        NetMeasure.Eloc_m = Eloc_m;
        
        % node betweenness centrality
        BC = betweenness_bin(M);
        BC = BC./((N-1)*(N-2));% normalize
        BC_m = mean(BC);
        
        NetMeasure.BC = BC;
        NetMeasure.BC_m = BC_m;
        
        % assortativity coefficient
        assort_coef = assortativity_bin(M,0);
        NetMeasure.assort_coef = assort_coef;
        
        % rich club
        [rich_coef,Nk,Ek] = rich_club_bu(M,length(M)-1);
        rich_coef(~isfinite(rich_coef)) = 0;
        NetMeasure.rich_club.rich_coef = rich_coef;
        NetMeasure.rich_club.Nk = Nk(isfinite(rich_coef));
        NetMeasure.rich_club.Ek = Ek(isfinite(rich_coef));
        % % topological overlap
        % GTOM = gtom(M,numSteps); % generalized topological overlap
        % %               measure (GTOM) matrix for number of steps, numSteps.
        % % Neighborhood overlap
        % [EC,ec,~] = edge_nei_overlap_bu(M); % edge neighborhood overlap matrix
        % % Matching index
        % M0 = matching_ind_und(M); % matching index matrix
        
    case 'BD' % binary directed
        % degree
        [indegree,outdegree,degree] = degrees_dir(M); % In/Out/node degree
        Kn_in = mean(indegree); % mean indegree of network
        Kn_out = mean(outdegree);% mean outdegree of network
        Kn = mean(degree); % mean degree of network
        Kcost = sum(degree)/(N*(N-1)); % cost of network
        
        NetMeasure.degree = degree;
        NetMeasure.indegree = indegree;
        NetMeasure.outdegree = outdegree;
        NetMeasure.Kn_in = Kn_in;
        NetMeasure.Kn_out = Kn_out;
        NetMeasure.Kn = Kn;
        NetMeasure.Kcost = Kcost;
        
        % clustering coefficient
        Cn = clustering_coef_bd(M);
        Cn_m = mean(Cn);
        
        NetMeasure.Cn = Cn;
        NetMeasure.Cn_m = Cn_m;
        
        % characteristic path length
        D1 = distance_bin(M);
        Ln = charpath(D1);
        NetMeasure.Ln = Ln;
        
        % global efficiency
        Eglob = efficiency_bin(M);
        NetMeasure.Eglob = Eglob;
        
        % local efficiency
        Eloc = efficiency_bin(M,1);
        Eloc_m = mean(Eloc);
        
        NetMeasure.Eloc = Eloc;
        NetMeasure.Eloc_m = Eloc_m;
        
        % node betweenness centrality
        BC = betweenness_bin(M);
        BC = BC./((N-1)*(N-2));
        BC_m = mean(BC);
        
        NetMeasure.BC = BC;
        NetMeasure.BC_m = BC_m;
        
        % assortativity coefficient
        assort_coef1 = assortativity_bin(M,1);
        assort_coef2 = assortativity_bin(M,2);
        assort_coef3 = assortativity_bin(M,3);
        assort_coef4 = assortativity_bin(M,4);
        NetMeasure.assort_coef1 = assort_coef1;
        NetMeasure.assort_coef2 = assort_coef2;
        NetMeasure.assort_coef3 = assort_coef3;
        NetMeasure.assort_coef4 = assort_coef4;
        
        % rich club
        [rich_coef,Nk,Ek] = rich_club_bd(M,2*(length(M)-1));
        rich_coef(~isfinite(rich_coef)) = 0;
        NetMeasure.rich_club.rich_coef = rich_coef;
        NetMeasure.rich_club.Nk = Nk(isfinite(rich_coef));
        NetMeasure.rich_club.Ek = Ek(isfinite(rich_coef));
        
    case 'WU' % weighted undirected
        A = M > 0; % adjacency matrix
        L = M;
        L(A) = 1 ./ L(A); % mapping from weight to distance (usually weight inversion)
        
        % strength
        strength = strengths_und(M); % node strength
        strength_m = mean(strength); % mean strength of network
        
        NetMeasure.strength = strength;
        NetMeasure.strength_m = strength_m;
        
        % clustering coefficient
        Cn = clustering_coef_wu(M);
        Cn_m = mean(Cn); % mean clustering coefficient
        
        NetMeasure.Cn = Cn;
        NetMeasure.Cn_m = Cn_m;
        
        % characteristic path length
        D1 = distance_wei(L);
        Ln = charpath(D1);
        NetMeasure.Ln = Ln;
        
        % Global efficiency
        Eglob = efficiency_wei(M);
        NetMeasure.Eglob = Eglob;
        
        % Local efficiency
        Eloc = efficiency_wei(M,2);
        Eloc_m = mean(Eloc);
        NetMeasure.Eloc = Eloc;
        NetMeasure.Eloc_m = Eloc_m;
        
        % node betweenness centrality
        BC = betweenness_wei(L);
        BC = BC./((N-1)*(N-2));
        BC_m = mean(BC);
        
        NetMeasure.BC = BC;
        NetMeasure.BC_m = BC_m;
        
        % Assortativity coefficient
        assort_coef = assortativity_wei(M,0);
        NetMeasure.assort_coef = assort_coef;
        
        % rich club
        Nnode = length(M);
        temp1 = (M~=0);
        if sum(temp1(:)) == Nnode*(Nnode-1) % full connection matrix
            thre = 0.1:0.1:0.9;
            rich_coef = [];
            for i = 1:length(thre)
                Mp = threshold_proportional(M,thre(i));
                ri = rich_club_wu(Mp,length(Mp)-1);
                % ri(~isfinite(ri)) = 0;
                rich_coef(i,:) = ri;
            end
            NetMeasure.rich_club.proportion = thre;
            NetMeasure.rich_club.rich_coef_proportion = rich_coef;
            NetMeasure.rich_club.rich_coef = nansum(rich_coef)./size(rich_coef,1); % mean across proportions
        else
            rich_coef = rich_club_wu(M);
            rich_coef(~isfinite(rich_coef)) = 0;
            NetMeasure.rich_club.rich_coef = rich_coef;
        end
        
    case 'WD' % weighted directed
        A = M > 0; % adjacency matrix
        L = M;
        L(A) = 1 ./ L(A); % mapping from weight to distance (usually weight inversion)
        
        % strength
        [instrength,outstrength,strength] = strengths_dir(M); % node strength
        instrength_m = mean(instrength); % mean instrength of network
        outstrength_m = mean(outstrength); % mean outstrength of network
        strength_m = mean(strength); % mean strength of network
        
        NetMeasure.strength = strength;
        NetMeasure.instrength = instrength;
        NetMeasure.outstrength = outstrength;
        NetMeasure.strength_m = strength_m;
        NetMeasure.instrength_m = instrength_m;
        NetMeasure.outstrength_m = outstrength_m;
        
        % clustering coefficient
        Cn = clustering_coef_wd(M);
        Cn_m = mean(Cn);
        
        NetMeasure.Cn = Cn;
        NetMeasure.Cn_m = Cn_m;
        
        % characteristic path length
        D1 = distance_wei(L);
        Ln = charpath(D1);
        NetMeasure.Ln = Ln;
        
        % Global efficiency
        Eglob = efficiency_wei(M);
        NetMeasure.Eglob = Eglob;
        
        % Local efficiency
        Eloc = efficiency_wei(M,2);
        Eloc_m = mean(Eloc);
        
        NetMeasure.Eloc = Eloc;
        NetMeasure.Eloc_m = Eloc_m;
        
        % node betweenness centrality
        BC = betweenness_wei(L);
        BC = BC./((N-1)*(N-2));
        BC_m = mean(BC);
        
        NetMeasure.BC = BC;
        NetMeasure.BC_m = BC_m;
        
        % rich club
        Nnode = length(M);
        temp1 = (M~=0);
        if sum(temp1(:)) == Nnode*(Nnode-1) % full connection matrix
            thre = 0.1:0.1:0.9;
            rich_coef = [];
            for i = 1:length(thre)
                Mp = threshold_proportional(M,thre(i));
                ri = rich_club_wd(Mp,2*(length(Mp)-1));
                % ri(~isfinite(ri)) = 0;
                rich_coef(i,:) = ri;
            end
            NetMeasure.rich_club.proportion = thre;
            NetMeasure.rich_club.rich_coef_proportion = rich_coef;
            NetMeasure.rich_club.rich_coef = nansum(rich_coef)./size(rich_coef,1); % mean across proportions
        else
            rich_coef = rich_club_wd(M);
            rich_coef(~isfinite(rich_coef)) = 0;
            NetMeasure.rich_club.rich_coef = rich_coef;
        end
       
end
        

