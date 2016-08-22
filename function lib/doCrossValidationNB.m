% Author: Ian Jacobsen
%
% This script contains the full procedure for performing cross-validation
% for the Naive Bayes model,
%
%% ----- Input Arguments -----
% data in format MxDxNxK
% targets in format Mx1xK
% PCA as a boolean flag (1 to perform PCA, 0 to not perform PCA)

%% ----- Output Values -----
% performance in the form [percentCorrect, percentTruePositive, percentTrueNegative]

%% Function Beginning
function performance = doCrossValidationNB(data, targets, PCA)

numFolds = length(data(1, 1, 1, :));

% train on k-1, test on 1
numCombs = nchoosek(numFolds, numFolds-1);

% find all combinations
combSeq = combnk(1:numCombs, numCombs-1);

for iter = 1:numCombs    
    fprintf('fold number: %d \n', iter)
    % ******** label training and testing data ********
    trainingData = doReshape_4DtoFHR(data(:, :, :, [combSeq(iter, :)]));
    trainingTargets = repmat(targets(:, :, iter), numCombs-1, 1); % CHECK THIS
    testingData = data(:, :, :, numCombs+1-iter); % because combination is in order
    testingTargets = targets(:, :, iter);
    

    if (PCA == 1)
        % ******** regularize data ********
        
        % parameter (we can change this)
        regMethod = 1; % either 0 for Kezi regularization, or 1 for standard normal
        
        % get data on same order of magnitude
        [trainingData, numeratorParam, denominatorParam] = doRegularizeData(trainingData, regMethod);
        
        % ******** principal components analysis ********
        
        % parameters
        keepVar = .98; % keep 98% of variance
        
        % get reduced feature matrix and the projection matrix
        [trainingData, projectionMatrix] = doPCA(trainingData, keepVar);
        
        % clean up workspace
        clear keepVar;
    end
    
    % ******** train model ******** 
    
    % estimate parameters for naive bayes model
    trainingTargets = [trainingTargets, ~trainingTargets];
    [prior, muHatHealthy, sigmaHatHealthy, muHatUnHealthy, sigmaHatUnHealthy] = doTrainNaiveBayes(trainingData, trainingTargets);
    % ******** test model ******** 

    if (PCA == 1)
        % regularize data
        testingData = doPrepareNewFetus(testingData, numeratorParam, denominatorParam, projectionMatrix);
    end
    
    % calculate log-likelihood and classify new fetus
    [results, likelihood] = doClassifyNaiveBayes(testingData, prior, muHatHealthy, sigmaHatHealthy, muHatUnHealthy, sigmaHatUnHealthy);

    
    performance(iter, :) = doEvaluatePerformance(results, testingTargets);
    
    
end

end
