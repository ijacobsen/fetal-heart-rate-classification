% Author: Ian Jacobsen

%% ----- Input Arguments -----
% data = data
% regularizationMethod = 0 for kezi, 1 for stdnrml
% kezi: (x - xmin)/(xmax - xmin)
% stdnrml: (x-xmean)/(stdx)

%% ----- Output Values -----
% regularizedData = data matrix with features on the same order of magnitude
% paramNum = numerator parameter to pass back (dMin for kezi, dMean for stdnrml)
% paramDen = denominator parameter to pass back (dMax-dMin for kezi, dstd for stdnrml)


%% ----- Unit Test -----
% run this as a test script
% result should be:
% mean(mean(data - check)) = 0 for both checks

% close all; clear all; clc;
% 
% % create test cases
% sl(:, :, 1) = [2 4; 1 3];
% sl(:, :, 2) = [4 7; 9 8];
% sl(:, :, 3) = [8 2; 9 1];
% sl(:, :, 4) = [7 3; 7 1];
% 
% % hand calculate statistics
% col1 = [2 1 4 9 8 9 7 7]';
% muCol1 = mean(col1);
% stdCol1 = std(col1);
% minCol1 = min(col1);
% maxCol1 = max(col1);
% 
% % hand calculate statistics
% col2 = [4 3 7 8 2 1 3 1]';
% muCol2 = mean(col2);
% stdCol2 = std(col2);
% minCol2 = min(col2);
% maxCol2 = max(col2);
% 
% % make matrices
% mu = repmat([muCol1, muCol2], 2, 1, 4);
% std = repmat([stdCol1, stdCol2], 2, 1, 4);
% maxim = repmat([maxCol1, maxCol2], 2, 1, 4);
% minim = repmat([minCol1, minCol2], 2, 1, 4);
% 
% % check stdnrml
% method = 1;
% check = (sl - mu) ./ (std);
% data = doRegularizeData(sl, method);
% mean(mean(data - check), 3)
% 
% % check kezi
% method = 0;
% check = (sl - minim) ./ (maxim - minim);
% data = doRegularizeData(sl, method);
% mean(mean(data - check), 3)

%% Function Beginning
function [regularizedData, paramNum, paramDen] = doRegularizeData(data, regularizationMethod)

% constants
M = length(data(:, 1, 1)); % number of training examples
D = length(data(1, :, 1)); % number of features
N = length(data(1, 1, :)); % number of observations

% reshape matrix to find min and max
dataTemp(1:M, 1:D) = data(:, :, 1);
for iter = 2:N
    dataTemp = [dataTemp; data(:, :, iter)];
end

dMax = max(dataTemp); % find max
dMin = min(dataTemp); % find min
dMean = mean(dataTemp); % find mean
dstd = std(dataTemp); % find std

dMaxRep = repmat(dMax, M, 1, N); % resize for operations
dMinRep = repmat(dMin, M, 1, N); % resize for operations
dMeanRep = repmat(dMean, M, 1, N); % resize for operations
dstdRep = repmat(dstd, M, 1, N); % resize for operations

if (regularizationMethod == 0)
    regularizedData = (data - dMinRep)./(dMaxRep-dMinRep); % kezi's method
    paramNum = dMin;
    paramDen = (dMax - dMin);
else
    regularizedData = (data - dMeanRep)./(dstdRep); % standard normal data
    paramNum = dMean;
    paramDen = dstd;
end

end