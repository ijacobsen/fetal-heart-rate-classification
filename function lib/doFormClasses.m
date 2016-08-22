% Author: Ian Jacobsen
% WARNING this is hardcoded for only 2 classes

%% ----- Input Arguments -----
% target = target values of pH
% threshold = threshold level which defines two classes

%% ----- Output Values -----
% targets = (M x 2) matrix of the form [healthy, nothealthy] 

%% Function Beginning
function targets = doFormClasses(target, threshold)

numClasses = 2;
M = length(target);

% initialize targets to 0
targets = zeros(M, numClasses);

% discriminator of classes
targets((target>threshold), 1) = 1;
targets(~(targets(:, 1)), 2) = 1;

end