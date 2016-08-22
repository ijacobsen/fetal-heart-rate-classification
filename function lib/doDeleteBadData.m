% Author: Ian Jacobsen
%
% This function has a dependency on the .mat file "deleteTheseIndices.mat".
%
% This file contains indices from the original data matrix of fetuses which
% we believe may contain too much artificial data.

%% ----- Input Arguments -----
% data in format MxDxN
% targets in format Mx1

%% ----- Output Values -----
% data in format MxDxN
% targets in format Mx1

%% Function Beginning

function [data, targets] = doDeleteBadData(data, targets)

load('deleteTheseIndices.mat');

% delete bad data
data(del', :, :) = [];

% delete bad targets
targets(del', :, :) = [];

end