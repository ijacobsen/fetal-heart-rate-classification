
% data matrix format
% [med, medDev, stv, sti, msti, poinstd1, poinstd2, poinccm higuch, psd1, psd2, psd3, psd4, psd5, ltv, delta, sampEn, fuzzyEn, mFHR, sdFHR, LTI, STV,  II]
% [ 1     2      3    4     5      6          7         8      9     10    11    12    13    14    15    16     17      18      19     20    21,  22,  23]  

clear all; close all;

%% load data and outputs

% file names
dataFile = 'ian_24mins_noLast5_withOverlap.mat';
targetFile = 'outcome_metrics.mat';

% parameters (we can vary these)
thresholdpH = 7.1501; % below is negative class (unhealthy), above is positive class (healthy)

% get data from files
[data, target] = doLoadData(dataFile, targetFile);

% delete bad indices
[data, target] = doDeleteBadData(data, target);

%data = data(:, [6 7 18 23 20], :);

% form classes (only 2 classes!)
targets = doFormClasses(target, thresholdpH);

% clean up workspace
clear target;
clear dataFile;
clear targetFile;
clear thresholdpH;
clear healthyIndices;
clear unhealthyIndices;
clear healthyTrainingData;
clear unhealthyTrainingData;
clear healthyTrainingSize;
clear unhealthyTrainingSize;
clear testSize;

%% create k-folds

numFolds = 10; % this can be changed

[kFoldData, kFoldTargs] = doCreateFolds(data, numFolds, targets(:, 1));

%% perform cross validation
PCA = 0; % 1 for PCA

performance = doCrossValidationNB(kFoldData, kFoldTargs, PCA);


