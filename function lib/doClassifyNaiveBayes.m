% Author: Ian Jacobsen
% NOTE: This script has a dependence on the function "logmvnpdf" that was
% written by Benjamin Dichter. This function can be found here:
% http://www.mathworks.com/matlabcentral/fileexchange/34064-log-multivariate-normal-distribution-function
%% ----- Input Arguments -----
% newFet = matrix of new fetuses to classify ... in the form (M x D x N)
% prior
% muHat(Un)Healthy = estimate of mu for multivariate gaussian distr for the (Un)Healthy class (1 x D x N)
% sigmaHat(Un)Healthy = estimate of sigma for multivariate gaussian distr for the (Un)Healthy class (D x D x N)

%% ----- Output Values -----
% estimatedClass = 1 for healthy, 0 for unhealthy
% loglikelihood = loglikelihood for estimated class

%% Function Beginning
function [estimatedClass, loglikelihood] = doClassifyNaiveBayes(newFet, prior, muHatHealthy, sigmaHatHealthy, muHatUnHealthy, sigmaHatUnHealthy)

% constants
M = length(newFet(:, 1, 1));
D = length(newFet(1, :, 1));
N = length(newFet(1, 1, :));

% priors
% healthyPrior = prior(1);
% unhealthyPrior = prior(2);
% *** better performance with uniform priors
healthyPrior = 0.5;
unhealthyPrior = 0.5; 

% likelihood functions p(x|y)
healthyLikelihood = zeros(M, 1);
unhealthyLikelihood = zeros(M, 1);


for fetusNum = 1:M
    for observNum = 1:N
        sumTermHealthy = logmvnpdf(newFet(fetusNum, :, observNum), muHatHealthy(1, :, N), sigmaHatHealthy(:, :, N));
        sumTermUnHealthy = logmvnpdf(newFet(fetusNum, :, observNum), muHatUnHealthy(1, :, N), sigmaHatUnHealthy(:, :, N));
        
        healthyLikelihood(fetusNum) = healthyLikelihood(fetusNum) + sumTermHealthy;
        unhealthyLikelihood(fetusNum) = unhealthyLikelihood(fetusNum) + sumTermUnHealthy;
    end
end



% form probabilities
healthyEstimate = healthyPrior * healthyLikelihood;
unhealthyEstimate = unhealthyPrior * unhealthyLikelihood;

% classify
estimatedClass = healthyEstimate>unhealthyEstimate;
loglikelihood(estimatedClass) = healthyEstimate(estimatedClass);
loglikelihood(~estimatedClass) = unhealthyEstimate(~estimatedClass);

% reshape
loglikelihood = loglikelihood';

end