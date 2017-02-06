function [C,index] = kmeanspp(X,k)
S = RandStream.getGlobalStream;
n = size(X,2);
% Select the first seed by sampling uniformly at random
index = zeros(1,k);
[C(:,1), index(1)] = datasample(S,X,1,2);
minDist = inf(n,1);

% Select the rest of the seeds by a probabilistic model
for ii = 2:k
    minDist = min(minDist,distfun(X,C(:,ii-1)));
    denominator = sum(minDist);
    if denominator==0 || isinf(denominator) || isnan(denominator)
        C(:,ii:k) = datasample(S,X,k-ii+1,2,'Replace',false);
        break;
    end
    sampleProbability = minDist/denominator;
    [C(:,ii), index(ii)] = datasample(S,X,1,2,'Replace',false,...
        'Weights',sampleProbability);
end

end


function d = distfun(X,C)
d = zeros(size(X,2),1);
for i = 1:size(X,2)
    d(i) = sum((X(:,i)-C).^2);
end

end