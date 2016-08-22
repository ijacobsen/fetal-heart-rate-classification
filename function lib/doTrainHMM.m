% Author: Ian Jacobsen

% This function has dependencies on the HMM Toolbox written by Kevin Murphy.
% This Toolbox can be downloaded here:
% https://www.cs.ubc.ca/~murphyk/Software/HMM/hmm.html

% O is the number of features
% T is the number of observations
% nex is the number of examples

%% ----- Input Arguments -----
% data is data ...
% M is the number of mixtures
% Q is the number of states
%% ----- Output Values -----
% LL is a vector of likelihoods for each iteration
% prior1 is the estimated prior
% transmat1 is the estimated transition matrix
% mu1 is the estimated mean
% sigma1 is the estimated covariance
% mixmat1 is the estimated weights for mixtures

%% Function Beginning

function [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = doTrainHMM(data, M, Q)

O = length(data(:, 1, 1));                                                  % vector dimension
T = length(data(1, :, 1));                                                  % length of one time series
nex = length(data(1, 1, :));                                                % number of time series
cov_type = 'full';                                                          

prior0 = normalise(rand(Q,1));                                              % some random prior
transmat0 = mk_stochastic(rand(Q,Q));                                       % more random initial values 

[mu0, Sigma0] = mixgauss_init(Q*M, reshape(data, [O T*nex]), cov_type);     % initial estimates for parameters mu and sigma for each state
mu0 = reshape(mu0, [O Q M]);                                                % " reshape
Sigma0 = reshape(Sigma0, [O O Q M]);                                        % " reshape
mixmat0 = mk_stochastic(rand(Q,M));                                         % "

[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = mhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 25); % EM algorithm



end