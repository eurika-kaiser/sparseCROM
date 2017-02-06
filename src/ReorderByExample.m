function [LabelsUpdated,index] = ReorderByExample(LabelsRef,Labels)
N = max(unique(Labels));
M = length(Labels);

index = zeros(N,1);
for i = 1:N
    idx = find(LabelsRef==i);
    is_set = 0;
    for j = 1:length(idx)
        if index(Labels(idx(j)))==0
            index(Labels(idx(j))) = i;
            is_set = 1;
            break;
        end
    end
    if is_set == 0
        idx = find(index==0);
        index(idx(1)) = i;
    end
end
LabelsUpdated = zeros(M,1);
for i = 1:N
   idx = find(Labels==i);
   LabelsUpdated(idx)=index(i);
end