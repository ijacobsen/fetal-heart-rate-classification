% Author: Ian Jacobsen

%% ----- Input Arguments -----
% data = data
% keepVar = number between 0 and 1 that serves as a percentage of how much
% informatino we want to retain in the reduced matrix

%% ----- Output Values -----
% reducedData = reduced dimension data matrix
% projectionMatrix = matrix that can be used to rotate future fetus

%% Function Beginning
function [reducedData, projectionMatrix] = doPCA(data, keepVar)

% constants
M = length(data(:, 1, 1)); % number of training examples
D = length(data(1, :, 1)); % number of features
N = length(data(1, 1, :)); % number of observations

% reshape matrix to find singular values
dataTemp(1:M, 1:D) = data(:, :, 1);
for iter = 2:N
    dataTemp = [dataTemp; data(:, :, iter)];
end

% PCA
covMatrix = 1/M * (dataTemp' * dataTemp); % form covariance matrix
[U, S, V] = svd(covMatrix); % perform singular value decomposition

% retain keepVar% variance
for iter = 1 : length(S)
    varRet(iter) = sum(sum(S(:, 1:iter)), 2)/sum(diag(S));
    if varRet(iter) >= keepVar
        k = iter;
        break;
    end
end

% project and rotate
projectionMatrix = U(:, 1:k);
tempReducedData = dataTemp*projectionMatrix; % transformed features stored in data matrix

% reshape data matrix into (M x D x N)
for iter = 1:N
    reducedData(:, :, iter) = tempReducedData((iter-1)*M+1:iter*M, :);
end


end