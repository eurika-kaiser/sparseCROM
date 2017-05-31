function kld = KLD(P, Q)
% Kullback-Leibler divergence
% Check equal size
if all( (size(P)-size(Q)) == 0) == 0
    disp('ERROR: Wrong size')
    return
end

% Check vec or mat
if any(size(P)==1) == 1
    type = 'vec';
    if size(P,1) == 1
        P = P';
    end
    if size(Q,1) == 1
        Q = Q';
    end
    
elseif all(size(P)>1)
    type = 'mat';
end



% Parameters
nCluster = size(P,1);

% Compute divergence
switch type
    
    case 'mat'
        KLDtmp = zeros(size(P,1),size(P,2));
        for iCluster = 1:nCluster
            for jCluster = 1:nCluster
                if P(iCluster,jCluster)==0
                    KLDtmp(iCluster,jCluster) = 0;
                else
                    KLDtmp(iCluster,jCluster) = P(iCluster,jCluster)*log(P(iCluster,jCluster)/Q(iCluster,jCluster));
                end
            end
        end
        kld = sum(sum(KLDtmp));
        
    case 'vec'
        KLDtmp = zeros(size(P,1),1);
        for iCluster = 1:nCluster
                if P(iCluster,1)==0
                    KLDtmp(iCluster,1) = 0;
                else
                    KLDtmp(iCluster,1) = P(iCluster,1)*log(P(iCluster,1)/Q(iCluster,1));
                end
        end
        kld = sum(KLDtmp);
        
end

end
