function [select_sensors,threshold] = applyThresholdToSensors(s,Nfeatures,Nclusters)
% threshold               = norm(s,'fro')/options.Npod/2;%max(max(abs(s)))/10^6 %
% select_sensors          = find(sum(abs(s),2) > threshold); % find locations above threshold

threshold               = norm(s,'fro')/(2*Nfeatures*Nclusters);
select_sensors          = [];
for i = 1:size(s,1)
    for j = 1:size(s,2)
        if abs(s(i,j)) > threshold
            select_sensors = [select_sensors; i];
        end
    end
end
select_sensors = unique(select_sensors);