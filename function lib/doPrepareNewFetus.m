% Author: Ian Jacobsen

%% ----- Input Arguments -----
% newFet = matrix of new fetus' to classify
% numParam = subtract from data
% denomParam = divide above difference by denomParam ie (d - numParam)/(denomParam)
% projectionMatrix = reduction matrix from PCA 

%% ----- Output Values -----
% data = prepared data matrix (regularized and reduced)

%% Function Beginning
function data = doPrepareNewFetus(newFet, numParam, denomParam, projectionMatrix)

% constants
M = length(newFet(:, 1, 1));
D = length(newFet(1, :, 1));
N = length(newFet(1, 1, :));

numeratorParam = repmat(numParam, M, 1, N);
denominatorParam = repmat(denomParam, M, 1, N);

% regularize data
newFet = (newFet - numeratorParam)./(denominatorParam);

% reshape matrix to prepare for reduction  (M*N x D)
dataTemp(1:M, 1:D) = newFet(:, :, 1);
for iter = 2:N
    dataTemp = [dataTemp; newFet(:, :, iter)];
end

% rotate and reduce data
newFetRed = dataTemp * projectionMatrix;

% reshape matrix for return (M x D x N)
for iter = 1:N
    data(:, :, iter) = newFetRed((iter-1)*M+1:iter*M, :);
end

end