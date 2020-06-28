function [PSI,PLV] = wb_calculate_PLV(data)
% Calculate Phase Synchronization Index
% Input:
%      data: The EEG potentials,channels X time points. e.g. 
%            62 channels X 5000 time points. EEG data should be filtered
%            first (passband filtering) to calculate PSI/PLV in specific band.
% Output:
%      PSI: connection matrix of PSI
%      PLV: connection matrix of PLI
% references:
% Bob, P., et al. (2008). "EEG phase synchronization in patients with paranoid schizophrenia." Neurosci Lett 447(1): 73-77.
% Edagawa, K. and M. Kawasaki (2017). "Beta phase synchronization in the frontal-temporal-cerebellar network during auditory-to-motor rhythm learning." Scientific Reports 7.
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC, Li_dong729@163.com)
% $ 2017.12.26
% -------------------------------------------------------------------------

if nargin < 1
    error ('One input is reqiured at least.');
end
% -----------------------
[Nchanns,~] = size(data);

PSI = zeros(Nchanns,Nchanns);% connection matrix 
PLV = zeros(Nchanns,Nchanns);% connection matrix 

for i1 = 1:Nchanns-1
     for j1 = (i1+1):Nchanns
        [PSI(i1,j1),PLV(i1,j1)] = calculate_PLV(data(i1,:),data(j1,:));
     end 
end

function [psi,plv] = calculate_PLV(xr,yr)
% Phase synchronization
% input:  
%      xr,yr are two different sequence, xr and yr should be filtered by bandpass 
%      filtering
% output:
%        PLV: phase synchronization index or phase locking
% -----------------------
% Discrete-time analytic signal using Hilbert transform
x = hilbert(xr);        
y = hilbert(yr);

%获取信号的瞬时相位 the instantaneous phases of two signal components
x_theta = angle(x);     
y_theta = angle(y);

%瞬时相位差 phase differences
xy_theta = x_theta-y_theta;

%phase locking/phase synchronization index
plv = abs(sum(exp(i*xy_theta)))/length(xr); % phase locking value
psi = sqrt((mean(cos(xy_theta)))^2+(mean(sin(xy_theta)))^2); % phase synchronization index
end
end





