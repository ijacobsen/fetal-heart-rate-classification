% Author: Ian Jacobsen

%% ----- Input Arguments -----
% data = data
% targets = classes (one column for each class)

%% ----- Output Values -----
% q = class priors
% mu = sample means
% sigma = sample standard deviations

%% ----- Unit Test -----
% % run this as a test script
% % outputs should be:
% % mean(mu(Un)Healthy,3) = approx [1 2 3]
% % diag(mean(sigma(Un)Healthy, 3))' = approx [0.5 0.5 0.5]
% clear all; close all; clc;
% % create 3x3x5 data matrix
% for iter = 1:5
%     testMat(:, :, iter) = normrnd([1 2 3;1 2 3;1 2 3; 1 2 3],.5,4,3);
% end
% % create targets
% outcome = [[1 0 0 1]', ~[1 0 0 1]'];
% % function call
% [q, muHealthy, sigmaHealthy, muUnHealthy, sigmaUnHealthy] = doTrainNaiveBayes(testMat, outcome);
% % test output
% mean(muHealthy, 3)
% diag(mean(sigmaHealthy, 3))'
% mean(muUnHealthy, 3)
% diag(mean(sigmaUnHealthy, 3))'

%% ----- Function Beginning -----
function [q, muHealthy, sigmaHealthy, muUnHealthy, sigmaUnHealthy] = doTrainNaiveBayes(data, targets)

% constant declarations
M = length(data(:, 1, 1)); % find number of training examples
D = length(data(1, :, 1)); % find length of feature vector
N = length(data(1, 1, :)); % find number of observations
numClass = length(targets(1, :)); % find number of classes

% find healthy and unhealthy index
healthyIndex = (targets(:, 1) == 1);
unhealthyIndex = (targets(:, 2) == 1);

% probability initializations
q = zeros(1, numClass); % initialize priors
muHealthy = zeros(1, D, N); % initialize mu
sigmaHealthy = zeros(D, D, N); % initialize sigma
muUnHealthy = zeros(1, D, N); % initialize mu
sigmaUnHealthy = zeros(D, D, N); % initialize sigma

% find priors
q = sum(targets)/M;

% estimate parameters
muHealthy = mean(data(healthyIndex, :, :));                     % (1 x D x N)
for iter = 1:N                                                  % (D x D x N) Covariance Matrix
    sigmaHealthy(:, :, iter) = cov(data(healthyIndex, :, iter));
end 

muUnHealthy = mean(data(unhealthyIndex, :, :));                 % (1 x D x N)
for iter = 1:N                                                  % (D x D x N) Covariance Matrix
    sigmaUnHealthy(:, :, iter) = cov(data(unhealthyIndex, :, iter));
end 

end