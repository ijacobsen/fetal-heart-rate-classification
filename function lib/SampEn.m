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