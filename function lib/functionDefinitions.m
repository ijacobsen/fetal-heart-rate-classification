%% Power Spectrum Density
function varagout = EstPsd(fhr, f1, f2, f3, fs)
    N = size(fhr,2);
    m = size(fhr, 1);
    p0 = zeros(m, 1);
    p1 = zeros(m, 1);
    p2 = zeros(m, 1);
    p3 = zeros(m, 1);
    ratio = zeros(m, 1);
    for iter = 1:m

        p0(iter) = bandpower(fhr(iter, :),fs,[0; f1])';
        p1(iter) = bandpower(fhr(iter, :),fs,[f1; f2])';
        p2(iter) = bandpower(fhr(iter, :),fs,[f2; f3])';
        p3(iter) = bandpower(fhr(iter, :),fs,[f3; fs/2])';
        ratio(iter) = p1(iter)/(p2(iter)+p3(iter));

    end
    
    varagout = [p0, p1, p2, p3, ratio]
end
% input:
%   fhr: FHR series in matrix form: row is FHR index, column is time index
%   f: vector of cut-off frequencies
%       f[0] = f_vlf
%       f[1] = f_lf
%       f[2] = f_mf
%       f[3] = f_hf
%   fs: Sampling Frequency of the data
% output:
%   p0: power in VLF band
%   p1: power in LF band
%   p2: power in MF band
%   p3: power in HF band
%   ratio: p1/(p2+p3)

%% Delta
function feature_delta = delta(fhr, fs)
    feature_delta = zeros(size(fhr,1),1);
    M = size(fhr,2)/(60*fs);
    for j = 1:M
        temp = fhr(:,60*fs*j-(60*fs)+1:60*fs*(j));
        feature_delta = feature_delta + (max(temp,[],2)-min(temp,[],2));
    end
    feature_delta = feature_delta./M;
end
% input:
%   fhr: FHR series in matrix form: row is FHR index, column is time index
% output: 
%   feature_delta: the Delta value defined in 4(a) in column vector form: row is FHR index

%% Long Term Variability
function feature_lti = ltv(fhr)
    squared_fhr = fhr.^2;
    squared_sum = squared_fhr(:,1:end-1)+squared_fhr(:,2:end);
    squared_sum = squared_sum.^.5;
    feature_lti = iqr(squared_sum,2);
end
% input:
%   fhr: FHR series in matrix form: row is FHR index, column is time index
% output: 
%   feature_delta: long-term irregularity defined in 4(b) in column vector form: row is FHR index

%% Short Term Variability
function feature_stv = stv(fhr, fs)

n = fs*60; % number of samples per minute
data = fhr;

diffMatrix = [(abs(data(:, 2:end) - data(:, 1:end-1))), zeros(length(data(:, 1)), 1)];

for iter=0:length(data(1, :))/n-1
    
    vs(:, iter+1) = (1/n).*sum(diffMatrix(:, iter*n+1:(iter+1)*n), 2);
    
end

feature_stv = mean(vs, 2);

end

%% Fuzzy Entropy
function feature_fuzzen = FuzzyEn(fhr, m, n)

mOrig = m;
N = length(fhr(1, :)); % time length
fetusSize = length(fhr(:, 1));
rMult = 0.2;
feature_fuzzen = zeros(fetusSize, 1);

for fetusNum = 1:fetusSize

    data = fhr(fetusNum, :);
    r = rMult*std(data);
    
    % set U vector
    m = mOrig;
    u = zeros(N-m+1, m);
    for iter = 1:N-m+1
        u(iter, :) = data(iter:iter+m-1) - repmat(mean(data(iter:iter+m-1)), 1, m);
    end
    
    % set C vector
    fuzzDist = 0;
    
    for iter = 1:N-m+1
        uLagged = circshift(u, iter);
        maxDiffMat = max(abs((u - uLagged)), [], 2);
        fuzzDist = fuzzDist + sum(exp((-(maxDiffMat).^n)/r));
    end
    B = (1/(N-m))*(1/(N-m+1))*fuzzDist;


     % set U+1 vector
     m = m+1;
    u = zeros(N-m+1, m);
    for iter = 1:N-m+1
        u(iter, :) = data(iter:iter+m-1) - repmat(mean(data(iter:iter+m-1)), 1, m);
    end
    
    % set C vector
    fuzzDist = 0;
    
    for iter = 1:N-m+1
        uLagged = circshift(u, iter);
        maxDiffMat = max(abs((u - uLagged)), [], 2);
        fuzzDist = fuzzDist + sum(exp((-(maxDiffMat).^n)/r));
    end
    A = (1/(N-m))*(1/(N-m+1))*fuzzDist;
    
    
    feature_fuzzen(fetusNum, 1) = log(B) - log(A);
end



end

%% Modified Short Term Irregularity
function feature_msti = msti(fhr, fs)

n = fs*60;
data = fhr;

differenceMatrix = [data(:, 2:end) - data(:, 1:end-1), zeros(length(data(:, 1)), 1)];

for iter=0:length(data(1, :))/n-1
    iqr_arctan_mat(:, iter+1) = iqr(atand(differenceMatrix(:, iter*n+1:(iter+1)*n)), 2);
end

feature_msti = mean(iqr_arctan_mat, 2);

end

%% Median Deviation
function feature_deviation = MedianDeviation(fhr)

N = length(fhr(1, :));

s_bar = median(fhr, 2);

feature_deviation = median(fhr - repmat(s_bar, 1, N), 2);

%% Sample Entropy
function feature_sampen = SampEn(fhr, m)

mOrig = m;
N = length(fhr(1, :)); % time length
fetusSize = length(fhr(:, 1));
rMult = 0.2;
feature_sampen = zeros(fetusSize, 1);

for fetusNum = 1:fetusSize

    data = fhr(fetusNum, :);
    r = rMult*std(data);
    
    % set U vector
    m = mOrig;
    u = zeros(N-m+1, m);
    for iter = 1:N-m+1
        u(iter, :) = data(iter:iter+m-1);
    end
    
    % set C vector
    count = 0;
    
    for iter = 1:N-m+1
        uLagged = circshift(u, iter);
        normMat = sqrt(sum((u - uLagged).^2, 2));
        count = sum(normMat <= r) + count;
    end
    B = (1/(N-m+1))*count;


     % set U+1 vector
     m = m+1;
    u = zeros(N-m+1, m);
    for iter = 1:N-m+1
        u(iter, :) = data(iter:iter+m-1);
    end
    
    % set C+1 vector
    count = 0;
    
    for iter = 1:N-m+1
        uLagged = circshift(u, iter);
        normMat = sqrt(sum((u - uLagged).^2, 2));
        count = sum(normMat <= r) + count;
    end
    A = (1/(N-m+1))*count;
    
    
    feature_sampen(fetusNum, 1) = -log(A/B);
    
end



end

%% Poincare
function varagout = poincare(fhr, m)

    dataOrig = fhr;
    
 % create RR interval and lag by one
    N = length(dataOrig(1, :));
    x_n = 60./dataOrig(:, 1:end-1);
    x_nPlusOne = 60./dataOrig(:, 2:end);
    data = 60 ./ dataOrig;
    %scatter(x_n, x_nPlusOne)

    SDSD = std(x_n - x_nPlusOne, 0, 2); % standard deviation of successive difference
    SDRR = std(x_n, 0, 2); % standard deviation of RR interval

    SD1 = sqrt(0.5) * SDSD; % as taken from Complex Correlation Measure (Karmakar)
    SD2 = sqrt(2*SDRR.^2 - 0.5*SDSD.^2); % ^^

    Cn = pi*(SD1.*SD2);

    consts = data(:, N-m).*data(:, N-1) + data(:, 2).*data(:, m+1) - data(:, N-1-m).*data(:, N) - data(:, 1).*data(:, m+2);

    sum1 = sum(data(:, 3:N-m) .* data(:, 1+m:N-2), 2);
    sum2 = sum(data(:, 2:N-m) .* data(:, 1+m:N-1), 2);
    sum3 = sum(data(:, 1:N-1-m) .* data(:, 2+m:N), 2);
    sum4 = sum(data(:, 1:N-2-m) .* data(:, 3+m:N), 2);

    CCM = (1./(2*Cn*(N-2))) .* (consts + sum1 + sum2 + sum3 + sum4);

    returnMat = [SD1, SD2, CCM];
    
    feature_sd1 = returnMat(:, 1);
    feature_sd2 = returnMat(:, 2);
    feature_ccm = returnMat(:, 3);
    
    varagout = [feature_sd1, feature_sd2, feature_ccm]
end

%% Higuchi
function [feature_hfd] = higuchi(fhr)
    kmax = floor(log2(size(fhr,2)))-1;
    feature_hfd = zeros(size(fhr,1),1);
    N = size(fhr,2);
    L=zeros(size(fhr,1),fix(log2(kmax)));
for i = 1:(log2(kmax))
    k = 2^(i);
    L(:,i) = zeros(size(fhr,1),1);
for s=1:k
    temp = fhr(:,s:k:N);
    L(:,i) = L(:,i) + (sum(abs(diff(temp,1,2)),2)*((N-1)/floor((N-s)/k)))/k;
end
    L(:,i) = L(:,i)/k;
end
for j=1:size(feature_hfd)
    pf = polyfit(1:log2(kmax),log2(L(j,:)),1);
    feature_hfd(j) = -pf(1);
end
end
% input:
%   fhr: FHR series in matrix form: row is FHR index, column is time index
% output: 
%   feature_hfd: Higuchi dimension defined in 5 in column vector form: row is FHR index
%   kmax: column vector with maximum of k value used for each FHR

%% Short Term Irregularity
function feature_sti = sti(fhr, fs)

n = fs*60; % number of samples per minute
data = fhr;

quotientMatrix = [data(:, 2:end)./data(:, 1:end-1), ones(length(data(:, 1)), 1)];

for iter=0:length(data(1, :))/n-1
    iqr_arctan_mat(:, iter+1) = iqr(atand(quotientMatrix(:, iter*n+1:(iter+1)*n)), 2);
end

feature_sti = mean(iqr_arctan_mat, 2);

end

