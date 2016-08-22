% Author: Ian Jacobsen
%% ----- Input Arguments -----
% data = matrix of data in format (M X D X N)
% numFolds = number of folds (5 sounds good)
% targets = targets for data ... 1 means healthy, 0 means unhealthy

%% ----- Output Values -----
% foldedData = (floor(M/K) x D x N x K) where K is the number of folds
% foldedTargets = (floor(M/K) x 1 x K), where the first (M/K)*posPercent are 1, and the remaining are 0

%% ----- Unit Test -----
% randHealthyExamp(1, :, :) - healthyExamp(healthyPerm(1), :, :) 
% should return 0 matrix

%% ----- Begin Function -----
function [foldedData, foldedTargets] = doCreateFolds(data, numFolds, targets)

%find size of data
m = length(data(:, 1, 1));

% find length of feature vector
d = length(data(1, :, 1));

% find number of observations
N = length(data(1, 1, :));

% find number of healthy and unhealthy
healthySize = sum(targets(:, 1));
unhealthySize = m - healthySize;

% determine split
healthyPercent = healthySize/m;
unhealthyPercent = 1 - healthyPercent;

% determine size for each fold
numHealthyPerFold = floor(healthyPercent*m/numFolds);
numUnhealthyPerFold = floor(unhealthyPercent*m/numFolds);

% find indices for healthy and unhealthy
healthyExampIndex = find(targets(:, 1) == 1);
unhealthyExampIndex = find(targets(:, 1) == 0);

% form healthy and unhealthy matrices
healthyExamp = data(healthyExampIndex, :, :);
unhealthyExamp = data(unhealthyExampIndex, :, :);

% form random permutation
healthyPerm = randperm(healthySize);
unhealthyPerm = randperm(unhealthySize);

% shuffle healthy and unhealthy examples
randHealthyExamp = healthyExamp(healthyPerm, :, :);
randUnhealthyExamp = unhealthyExamp(unhealthyPerm, :, :);

% create folds and associated targets
foldedData = zeros(numHealthyPerFold+numUnhealthyPerFold, d, N, numFolds);
foldedTargets = zeros(numHealthyPerFold+numUnhealthyPerFold, 1, numFolds);
for iter = 1:numFolds
    hIndexMin = numHealthyPerFold*(iter-1) + 1;
    hIndexMax = numHealthyPerFold*iter;
    uIndexMin = numUnhealthyPerFold*(iter-1) + 1;
    uIndexMax = numUnhealthyPerFold*iter;
    foldedData(:, :, :, iter) = [randHealthyExamp(hIndexMin:hIndexMax, :, :); randUnhealthyExamp(uIndexMin:uIndexMax, :, :)];
    foldedTargets(:, 1, iter) = [ones(numHealthyPerFold, 1); zeros(numUnhealthyPerFold, 1)];
end


end