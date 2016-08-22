% Author: Ian Jacobsen

%% ----- Input Arguments -----
% fileNameData is the file to load data from
% fileNameTargets is the file to load targets from

% *** IMPORTANT! ***
% the name of the variable that's loaded from 'fileNameData' must be 'computedFeatures'
% likewise with 'fileNameTargets' being "outcome_metrics"

%% ----- Output Values -----
% data = data ....
% targets = target values

%% ----- Function Beginning -----
function [data, targets] = doLoadData(fileNameData, fileNameTargets)

% load data and targets
load(fileNameData)
load(fileNameTargets)

data = computedFeatures;
targets = outcome_metrics{1};

end