function [NearestCentroid, idx] = getNearestCluster(y, Centroids)
% if size(y,1)~= 1 && size(y,2) ~=1
%     y = y';
% end

% Parameters
nCluster = size(Centroids,1);

if size(y,1) == 1 % just one time instant
    % Euclidean distance
    dist    = sqrt(sum((Centroids-(repmat(y,[nCluster, 1]))).^2,2));
    [~,idx] = min(dist);
    NearestCentroid = Centroids(idx,:);
else
    %y = y';
    M = size(y,1);
    dist = zeros(M,nCluster);
    for iCluster = 1:nCluster
       dist(:,iCluster) = sqrt(sum((  repmat(Centroids(iCluster,:),[M 1])-y  ).^2,2));  
    end
    [~,idx] = min(dist,[],2);
    NearestCentroid = Centroids(idx,:);
end
end