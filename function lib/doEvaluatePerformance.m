% Author: Ian Jacobsen

%% ----- Input Arguments -----
% results in format Mx1
% targets in format Mx1

%% ----- Output Values -----
% performance in size Mx3 with format:
% performance = [percentCorrect, percentTruePositive, percentTrueNegative]

%% Function Beginning
function performance = doEvaluatePerformance(results, targets)

    testSize = length(results(:, 1));

    % evaluate performance
    error = abs(targets(1:testSize, 1) - results);
    sumError = sum(error);
    ratio = sumError/testSize;
    percentCorrect = 100*(1-ratio);
    
    
    % true positive
    healthyIndex = find(targets(1:testSize, 1) == 1);
    error = sum(abs(results(healthyIndex) - targets(healthyIndex, 1)));
    truePositive = length(healthyIndex) - error;
    percentTruePositive = 100*truePositive/length(healthyIndex);
    
    
    % true negative
    unhealthyIndex = find(targets(1:testSize, 1) == 0);
    error = sum(abs(~results(unhealthyIndex) - ~targets(unhealthyIndex, 1)));
    trueNegative = length(unhealthyIndex) - error;
    percentTrueNegative = 100*trueNegative/length(unhealthyIndex);
    
    
    performance = [percentCorrect, percentTruePositive, percentTrueNegative];
    
end