function [w, Sw,Sb, D, V] = LDA(Data, centroids, labels, Data_mean) 
%% Linear Discriminant Analysis
% Determine discrimant vectors of centroids
disp('COMPUTING: LDA ...')

if size(Data,1) == size(centroids,2)
    Data = Data';
end

if nargin < 5
    disp('WARNING: Compute mean for LDA based on input data ...')
    Data_mean = mean(Data);
end


%% Parameters
Nclusters = size(centroids,1);
Npod      = size(centroids,2);

if size(Data_mean,1) == size(centroids,2)
    Data_mean = Data_mean';
end

%% LDA
% 1) Means: are given by centroids

% 2) Scatter matrices
% within-cluster scatter matrix
Sw = zeros(Npod,Npod);
for iCluster = 1:Nclusters
    idx                   = find(labels==iCluster);
    Nk                    = length(idx);
    fluct_cluster         = Data(idx,:)-repmat(centroids(iCluster,:),[Nk 1]);
    C                     = (fluct_cluster'*fluct_cluster); %(1/(Nk-1))*
    Sw                    = Sw + C; %(Nk-1)*
end
    
% between-cluster scatter matrix
Sb = zeros(Npod,Npod);
for iCluster = 1:Nclusters
    idx                   = find(labels==iCluster);
    Nk                    = length(idx);
    C                     = Nk*((centroids(iCluster,:) - Data_mean)'*(centroids(iCluster,:) - Data_mean));
    Sb                    = Sb + C;
end

% 3) Solving generalized eigenvalue problem for Sw^-1*Sb
%[V,D] = eig(Sw\Sb);
%[V,D] = eig(Sw,Sb);
% [V,D]  = eigs(pinv(Sw)*Sb, Nclusters-1);
% [D,IX] = sort(diag(D),'descend');
% V      = V(:,IX);

try
    A = pinv(Sw)*Sb;
catch
    A = 0;
    Sw
end
if size(A,2) < Nclusters-1
    [V,D]  = eig(A);
    p = size(A,2);
else
    [V,D]  = eig(A);
    p = Nclusters - 1;
end
[D,IX] = sort(diag(D),'descend');
V      = V(:,IX(1:p));

% Check
eig_err = zeros(size(V,2),1);
for iEV = 1:size(V,2)
    v               = V(:,iEV);
    eig_err(iEV)    = max(abs(Sb*v - D(iEV).*Sw*v));
end

% 4) Selection of linear discriminants (based on D ranking)
eigv_sum  = sum(D);
eigv_norm = real(D./eigv_sum);


% Normalize w vectors, 
%w = V(:,1:p);
for i = 1:p 
    w(:,i) = V(:,i) / sqrt(V(:,i)'*V(:,i));
end;


end