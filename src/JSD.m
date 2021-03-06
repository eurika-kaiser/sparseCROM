function jsd = JSD(P, Q)
% Jensen-Shannon divergence
% Check equal size
if all( (size(P)-size(Q)) == 0) == 0
    disp('ERROR: Wrong size')
    return
end

% Parameters
M = 1/2*(P+Q);

% Compute divergence
jsd = 1/2*KLD(P,M) + 1/2*KLD(Q,M);

end

