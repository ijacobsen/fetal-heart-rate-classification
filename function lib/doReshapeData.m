% Author: Ian Jacobsen
% Puts in form (D x N x M) for use with Kevin Murphy's HMM Library

%% ----- Input Arguments -----
% data in format MxDxN

%% ----- Output Values -----
% data in format DxNxM

%% Function Beginning
function shapedData = doReshapeData(data)

numExamples = length(data(:, 1, 1));
numObservations = length(data(1, 1, :));
numFeatures = length(data(1, :, 1));

shapedData = zeros(numFeatures, numObservations, numExamples);

for iter = 1:numExamples
    shapedData(:, : , iter) = data(iter, :, :);
end

end