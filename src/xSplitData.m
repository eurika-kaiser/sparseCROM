function Results = xSplitData(M,traindata_perc,testdata_perc)
% Assumes Data.ts R^(NxM), N spatial/modes, M time

%% Splitting into test and validation data
disp(['---> Split data into training and test data sets ... ', num2str(traindata_perc),'/',num2str(testdata_perc)])
if traindata_perc < 1
    Results.M_traindata                                     = floor(traindata_perc*M);
    Results.M_testdata                                      = M - Results.M_traindata;
    Results.selection_traindata                             = sort(randperm(M,Results.M_traindata)','ascend');
    Results.selection_testdata                              = [1:M]';
    Results.selection_testdata(Results.selection_traindata) = [];

else
    Results.M_traindata                                     = M;
    Results.M_testdata                                      = M;
    Results.selection_traindata                             = [1:M];
    Results.selection_testdata                              = [1:M]';

end

end