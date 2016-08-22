% Author: Ian Jacobsen
% Puts in form (D x N x M) for use with Kevin Murphy's HMM Library

%% ----- Input Arguments -----
% data in format MxDxNxK

%% ----- Output Values -----
% data in format DxNxM

%% Function Beginning
function shapedData = doReshape_4DtoFHR(data)

numFolds = length(data(1, 1, 1, :));
tempData = data(:, :, :, 1);

for iter = 2:numFolds
    tempData = cat(1, tempData, data(:, :, :, iter));
end

shapedData = tempData;

end