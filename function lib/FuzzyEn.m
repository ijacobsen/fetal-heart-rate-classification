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