% Author: Ian Jacobsen
%
% This script contains the full procedure for performing cross-validation
% for the HMM model,
%
%% ----- Input Arguments -----
% data in format MxDxNxK
% targets in format Mx1xK
% M as number of mixtures
% Q as number of hidden states
% PCA as a boolean flag (1 to perform PCA, 0 to not perform PCA)

%% ----- Output Values -----
% performance in the form [percentCorrect, percentTruePositive, percentTrueNegative]

%% Function Beginning
function performance = doCrossValidationHMM(data, targets, M, Q, PCA)

numFolds = length(data(1, 1, 1, :));

% train on k-1, test on 1
numCombs = nchoosek(numFolds, numFolds-1);

% find all combinations
combSeq = combnk(1:numCombs, numCombs-1);

% initialize size
performance = zeros(numCombs, 3);

for iter = 1:numCombs    
    fprintf('fold number: %d \n', iter)
    % ******** label training and testing data ********
    trainingData = doReshape_4DtoFHR(data(:, :, :, [combSeq(iter, :)]));
    trainingTargets = repmat(targets(:, :, iter), numCombs-1, 1); 
    testingData = data(:, :, :, numCombs+1-iter); % because combination is in order
    testingTargets = targets(:, :, iter);
    

    if (PCA == 1)
        % ******** regularize data ********
        
        % parameter (we can change this)
        regMethod = 0; % either 0 for Kezi regularization, or 1 for standard normal
        
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
    
    % reshape data for hmm library
    shapedData = doReshapeData(trainingData);
    shapedTestData = doReshapeData(testingData);
    
    % learn parameters using EM algorithm
    
    % healthy model
    posExamples = find(trainingTargets(:, 1) == 1);
    posData = shapedData(:, :, posExamples);
    [h_LL, h_prior1, h_transmat1, h_mu1, h_Sigma1, h_mixmat1] = doTrainHMM(posData, M, Q);
    
    % unhealthy model
    negExamples = find(trainingTargets(:, 1) == 0);
    negData = shapedData(:, :, negExamples);
    [u_LL, u_prior1, u_transmat1, u_mu1, u_Sigma1, u_mixmat1] = doTrainHMM(negData, M, Q);
    
    % ******** test model ******** 

    if (PCA == 1)
        % regularize data
        testingData = doPrepareNewFetus(testingData, numeratorParam, denominatorParam, projectionMatrix);
    end
    
    % reshape data
    testingData = doReshapeData(testingData);
    
    % estimate likelihood
    numTest = length(testingData(1, 1, :));
    logLikHealthy = zeros(numTest, 1);
    logLikUnhealthy = zeros(numTest, 1);
    for jter = 1:numTest
        
        logLikHealthy(jter) = mhmm_logprob(testingData(:, :, jter), h_prior1, h_transmat1, h_mu1, h_Sigma1, h_mixmat1);
        logLikUnhealthy(jter) = mhmm_logprob(testingData(:, :, jter), u_prior1, u_transmat1, u_mu1, u_Sigma1, u_mixmat1);
        
    end
    
    results = (logLikHealthy > logLikUnhealthy);
    
    performance(iter, :) = doEvaluatePerformance(results, testingTargets);
    
    
end

end
