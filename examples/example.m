%% sparseCROM of the periodic double gyre flow
clear all, close all, clc

addpath('../src/');

options.path2save = 'output/';
mkdir(options.path2save)
load('DoubleGyre.mat');

%% Parameters
Nclusters   = 4;           % number of clusters
Nkmeans     = 100;         % repetitions of k-means
rmax        = Data.N3-1;   % number of features to keep
cmap        = bone(100);
Data.C      = randn(100,Data.N1*Data.N2); % Gaussian ramdom measurement matrix // Change to random rows of eye for random point measurements

%% Plot vorticity
cmin = -2; cmax = 2; cvals = [cmin:0.5:cmax]; m = 10;
filename = ['Vorticity_m',sprintf('%03g',m)];
fhandle = figure;
contourf(Data.x,Data.y,Data.Omega(:,:,m),cvals);caxis([cmin cmax]),hold on
quiver(Data.x,Data.y,Data.U(1:Data.N1,:,m),Data.U(Data.N1+1:end,:,m),'k')
colormap(cmap),
xlabel('x'), ylabel('y')
box on
daspect([1 1 1])
set(gca,'FontSize',14,'LineWidth',1)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', fullfile(options.path2save,filename));
close(fhandle);

%% POD via SVD
X           = reshape(Data.Omega,[Data.N1*Data.N2 Data.N3]);
Xp          = X - repmat(mean(X,2),[1 Data.N3]);
Xp          = [Xp];
[U,S,V]     = svd(Xp,'econ'); U = U(:,1:rmax);
ai          = (U'*Xp)';
Ntimes      = size(ai,1);

%% (1) --- Compute CROM from full-state features as benchmark
rng('default')      % set seed for
Data2cluster        = ai;
c0_Labels           = cell(Nkmeans,1);
c0_Labels_cpp       = cell(Nkmeans,1);
c0_Centroids        = cell(Nkmeans,1);
c0_Centroids_index  = cell(Nkmeans,1);
c0_Centroids_cpp    = cell(Nkmeans,1);
sumD                = cell(Nkmeans,1);
D                   = cell(Nkmeans,1);
kmeans_eval         = zeros(Nkmeans,1);
for i = 1:Nkmeans
    [c0_Centroids_cpp{i},c0_Centroids_index{i}] = kmeanspp(Data2cluster',Nclusters);
    c0_Centroids_cpp{i} = c0_Centroids_cpp{i}';
    Cpp = c0_Centroids_cpp{i};
    [c0_Labels{i}, c0_Centroids{i}, sumD{i}, D{i}] = kmeans(Data2cluster, Nclusters,'Start',Cpp,'Maxiter',5000);
    kmeans_eval(i)  = sum(sumD{i});
    disp(['i = ',num2str(i)])
end
[~,idx] = min(kmeans_eval);

ClusteringResults.c0_Labels         = c0_Labels{idx};
ClusteringResults.c0_Centroids      = c0_Centroids{idx};
ClusteringResults.c0_Centroids_cpp   = c0_Centroids_cpp{idx};
ClusteringResults.sumD              = sumD{idx};
ClusteringResults.D                 = D{idx};
ClusteringResults.c0_Centroids_index = c0_Centroids_index{idx};

figure,
plot(Data2cluster(:,1),Data2cluster(:,2),'ok')
hold on
plot(ClusteringResults.c0_Centroids_cpp(:,1),ClusteringResults.c0_Centroids_cpp(:,2),'xr')

% Find ICs // Check
D = zeros(Ntimes,Nclusters);
for i = 1:Ntimes
    for j = 1:Nclusters
        D(i,j) = sum((Data2cluster(i,:)-ClusteringResults.c0_Centroids_cpp(j,:)).^2);
    end
end
ClusteringResults.ICidx = zeros(Nclusters,1);
for i = 1:Nclusters
    idx = find(D(:,i)==0);
    ClusteringResults.ICidx(i) = idx;
end
all((ClusteringResults.c0_Centroids_index'-ClusteringResults.ICidx)==0)
save(fullfile(options.path2save,'ClusteringResults.mat'),'ClusteringResults');
clear c0_Labels c0_Labels_cpp c0_Centroids c0_Centroids_index c0_Centroids_cpp sumD kmeans_eval

% Transition matrix
[CROM.P, CROM.c1_Labels, CROM.c1_Centroids] = determineClusterTransitionMat(ClusteringResults.c0_Labels, ...
    ClusteringResults.c0_Centroids, Ntimes,'basic', []);

save(fullfile(options.path2save,'CROMResults.mat'),'CROM');

%% (2) --- Take measurements
Data.Omega = reshape(Data.Omega,[Data.N1*Data.N2 Data.N3]); % Reshape vorticity field into vectors
Data.Y = Data.C*Data.Omega; Data.Y = Data.Y';

%% (3A) --- Recompute Labels from measurements for existing clustering
name_post       = 'keep_clustering/';
path2save       = [options.path2save,name_post];
mkdir(path2save)

% Params
selection2plot  = [1,5,10,50,100];
select_Nmeas    = [1:200]; NmeasMax = length(select_Nmeas);
M               = size(Data.Y,1); % number of snapshots

% Init
Labels          = CROM.c1_Labels;
classerr        = zeros(NmeasMax,1);

for i = 1:NmeasMax
    Nmeas       = select_Nmeas(i); % number of measurements for this run
    Nfeatures   = size(Data.Y,2);  % Reset number of features
    if Nmeas<Nfeatures             % Check if more measurements are taken than available
        Nfeatures = Nmeas;
    end
    
    % Recompute centroids from measurements
    CentroidsY = zeros(Nclusters,Nfeatures);
    for iC = 1:Nclusters
        CentroidsY(iC,:) = mean(Data.Y(Labels == iC,1:Nfeatures));
    end
    
    % Determine labels from measurements and recomputed centroids using
    % Nearest-Centroid method
    [~,LabelsY] = getNearestCluster(Data.Y(:,1:Nfeatures),CentroidsY);
    
    % Classification error
    classerr(i) = (M - length(find(Labels==LabelsY)))./M;
    
    % Plot
    if any(selection2plot==select_Nmeas(i))
        
        f1 = figure;
        plot(Labels,'.-k','LineWidth',2), hold on
        plot(LabelsY,'.--r','LineWidth',2)
        xlabel('Snapshot number')
        ylabel('Cluster')
        axis tight
        ylim([0 Nclusters+1])
        set(gca,'yTick',[1,Nclusters]);
        set(gca,'FontSize',14)
        set(gcf,'Position',[50 50 600 200])
        set(gcf,'PaperPositionMode','auto')
        print('-depsc2', '-loose', [path2save,['Classification_Meas_',num2str(select_Nmeas(i)),'.eps']]);
        close(f1);
        
    end
end

% Plot classification error
f1 = figure;
plot(select_Nmeas,classerr,'-k','LineWidth',2)
xlabel('Number of measurements')
ylabel('Classification error')
axis([0 M 0 1])
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2save,['ClassificationError.eps']]);
close(f1);

%% (3B) --- sparseCROM from measurement labels
cmap = (gray(12));

% Determine and plot transition probabilisty matrix from full-state fetaures
DataSet         = TransitionMatrixIdentification(Labels,1);
Q               = DataSet.P;
tmp = Q; % Prepare for log10 plots
tmp(tmp(:)<0.001) = -4; idx = find(tmp(:)>=0.001);
tmp(idx) = log10(tmp(idx));

f1 = figure;
p1 = imagesc(tmp); axis equal, axis tight
colormap(cmap)
cb = colorbar;
ax1 = get(gca); Pos1 = ax1.Position;
cb.Ticks = [-4,-3.5:1:-0.5,0];
cb.TickLabels = [0,0.0005,0.005,0.05,0.5,1];
cb.Position = [0.8 cb.Position(2) 0.025 cb.Position(4)];
set(gca,'Position',[ax1.Position(1)-0.01 ax1.Position(2:4)])
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 250 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2save,['Q.eps']]);
close(f1);

% Init
jsd         = zeros(NmeasMax,1);

for i = 1:NmeasMax
    Nmeas = select_Nmeas(i);
    Nfeatures = size(Data.Y,2);
    if Nmeas<Nfeatures
        Nfeatures = Nmeas;
    end
    M         = size(Data.Y,1);
    
    % Centroids
    CentroidsY = zeros(Nclusters,Nfeatures);
    for iC = 1:Nclusters
        idx = find(Labels == iC);
        CentroidsY(iC,:) = mean(Data.Y(idx,1:Nfeatures));
    end
    
    % Labels
    [~,LabelsY] = getNearestCluster(Data.Y(:,1:Nfeatures),CentroidsY);
    
    % CROM
    DataSet = TransitionMatrixIdentification(LabelsY,1);
    Py = DataSet.P;
    
    % Jensen-Shannon divergence
    jsd(i) = JSD(Py,Q);
    
    % Plot
    if any(selection2plot==select_Nmeas(i))
        cmap = (gray(12)); %flipud
        tmp = Py;
        tmp(tmp(:)<0.001) = -4;
        idx = find(tmp(:)>=0.001);
        tmp(idx) = log10(tmp(idx));
        
        f1 = figure;
        p1 = imagesc(tmp); axis equal, axis tight
        colormap(cmap)
        cb = colorbar;
        ax1 = get(gca); Pos1 = ax1.Position;
        cb.Ticks = [-4,-3.5:1:-0.5,0];
        cb.TickLabels = [0,0.0005,0.005,0.05,0.5,1];
        cb.Position = [0.8 cb.Position(2) 0.025 cb.Position(4)];
        set(gca,'Position',[ax1.Position(1)-0.01 ax1.Position(2:4)])
        set(gca,'FontSize',14)
        set(gcf,'Position',[0 0 250 200])
        set(gcf,'PaperPositionMode','auto')
        
        print('-depsc2', '-loose', [path2save,['P_',num2str(select_Nmeas(i)),'.eps']]);
        close(f1);
    end
end

f1 = figure;
plot(select_Nmeas,jsd,'-k','LineWidth',2)
xlabel('Number of measurements')
ylabel('JSD')
axis([0 M 0 max(jsd)])
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2save,['JSD.eps']]);
close(f1);


%% (4) --- sparseCROM from clustered measurement data
name_post = 'new_clustering/';
path2save = [options.path2save,name_post];
mkdir(path2save)

% Init
classerr        = zeros(NmeasMax,1);
Labels_new      = zeros(Ntimes,NmeasMax); 
Centroids_new   = cell(NmeasMax,1);
IndexSorting    = zeros(Nclusters,NmeasMax);
jsd             = zeros(NmeasMax,1);

for i = 1:NmeasMax
    Nmeas = select_Nmeas(i);
    Nfeatures = size(Data.Y,2);
    if Nmeas<Nfeatures
        Nfeatures = Nmeas;
    end
    M         = size(Data.Y,1);
    
    % Clustering
    Data2cluster = Data.Y(:,1:Nfeatures);
    Cpp          = Data2cluster(ClusteringResults.c0_Centroids_index,:);
    [Labels_new(:,i), Centroids_new{i}, ~, ~] = kmeans(Data2cluster, Nclusters,'Start',Cpp,'Maxiter',5000);
    
    % Reorder labels + update centroids
    [Labels_new(:,i),IndexSorting(:,i)] = ReorderByExample(Labels,Labels_new(:,i));
    Centroids_new{i} = Centroids_new{i}(IndexSorting(:,i),:);
    
    % Transition probability matrix
    DataSet = TransitionMatrixIdentification(Labels_new(:,i),1);
    Py = DataSet.P;
    
    % Classification error
    classerr(i) = (M - length(find(Labels==Labels_new(:,i))))./M;
    
    % Plot
    if any(selection2plot==select_Nmeas(i))
        
        f1 = figure;
        plot(Labels,'.-k','LineWidth',2), hold on
        plot(Labels_new(:,i),'.--r','LineWidth',2)
        xlabel('Snapshot number')
        ylabel('Cluster')
        ylim([0 Nclusters+1])
        set(gca,'yTick',[1,Nclusters]);
        set(gca,'FontSize',14)
        set(gcf,'Position',[50 50 600 200])
        set(gcf,'PaperPositionMode','auto')
        print('-depsc2', '-loose', [path2save,['Classification_Meas_',num2str(select_Nmeas(i)),'.eps']]);
        close(f1);
    end

    % Jensen-Shannon divergence
    jsd(i) = JSD(Py,Q);
    
    % Plot
    if any(selection2plot==select_Nmeas(i))
        cmap = (gray(12)); %flipud
        tmp = Py;
        tmp(tmp(:)<0.001) = -4;
        idx = find(tmp(:)>=0.001);
        tmp(idx) = log10(tmp(idx));
        
        f1 = figure;
        p1 = imagesc(tmp); axis equal, axis tight
        colormap(cmap)
        cb = colorbar;
        ax1 = get(gca); Pos1 = ax1.Position;
        cb.Ticks = [-4,-3.5:1:-0.5,0];
        cb.TickLabels = [0,0.0005,0.005,0.05,0.5,1];
        cb.Position = [0.8 cb.Position(2) 0.025 cb.Position(4)];
        set(gca,'Position',[ax1.Position(1)-0.01 ax1.Position(2:4)])
        set(gca,'FontSize',14)
        set(gcf,'Position',[0 0 250 200])
        set(gcf,'PaperPositionMode','auto')
        print('-depsc2', '-loose', [path2save,['P_',num2str(select_Nmeas(i)),'.eps']]);
        close(f1);
    end
end

f1 = figure;
plot(select_Nmeas,classerr,'-k','LineWidth',2)
xlabel('Number of measurements')
ylabel('Classification error')
axis([0 M 0 1])
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2save,['ClassificationError.eps']]);
close(f1);


f1 = figure;
plot(select_Nmeas,jsd,'-k','LineWidth',2)
xlabel('Number of measurements')
ylabel('JSD')
axis([0 M 0 max(jsd)])
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2save,['JSD.eps']]);
close(f1);

%% (5) --- Find optimized point measurements tailored to full-stated CROM
name_post = 'optimized_locations/';
path2save = [options.path2save,name_post];
mkdir(path2save)

% Parameters
recompute_classifier    = 1;
recompute_svd           = 1;
traindata_perc          = 1;%0.8; % set <1 for cross-validation
testdata_perc           = 1-traindata_perc;
Nfeatures               = 3;  
Nrepetitions            = 100;
M                       = Data.N3;
lambda                  = 10;
epsilon                 = 10^(-10);

% Split data into training and test set for crossvalidation (if traindata_perc < 1)
Results_splitting = xSplitData(M,traindata_perc,testdata_perc);

% POD via SVD
[U,S,~] = svd(bsxfun(@minus, Data.Omega(:,Results_splitting.selection_traindata), mean(Data.Omega(:,Results_splitting.selection_traindata),2)),'econ');
S = diag(S); sigma = S;
evals = sigma.^2/Results_splitting.M_traindata; 
Psi = U;
ai  = Psi'*bsxfun(@minus, Data.Omega,mean(Data.Omega(:,Results_splitting.selection_traindata),2)); ai = ai';

% Use full-state clustering results (could also use Labels from (1), but here recomputed from the test set for cross-validation)
Data2cluster = ai(Results_splitting.selection_traindata,:);
Cpp          = Data2cluster(ClusteringResults.c0_Centroids_index,:);
[Labels_new, Centroids_new, ~, ~] = kmeans(Data2cluster, Nclusters,'Start',Cpp,'Maxiter',5000);

% Reorder labels + update centroids
[Labels_new,IndexSorting] = ReorderByExample(Labels,Labels_new);
Centroids_new = Centroids_new(IndexSorting,:);
    
% Linear discriminant analysis
[w, ~, ~, ~, ~] = LDA(ai(Results_splitting.selection_traindata,1:Nfeatures), Centroids_new(:,1:Nfeatures), Labels_new) ;

% SSPOC
s = SSPOC(Psi(:,1:Nfeatures),w,'lambda',lambda, 'epsilon',epsilon);
[sensors_locations,threshold] = applyThresholdToSensors(s,Nfeatures,Nclusters);
Nsensors = length(sensors_locations);

% Take measurements
Data.Y = Data.Omega(sensors_locations,:)';

% Recompute classifier
CentroidsY = zeros(Nclusters,Nsensors);
for iC = 1:Nclusters
    idx = find(Labels_new == iC);
    CentroidsY(iC,:) = mean(Data.Y(Results_splitting.selection_traindata(idx),1:Nsensors));
end
[wY, ~, ~, ~, ~] = LDA(Data.Y(Results_splitting.selection_traindata,:), CentroidsY, Labels_new) ;

% Classify test set
[~,LabelsYtest_centroids] = getNearestCluster(Data.Y(Results_splitting.selection_testdata,:),CentroidsY);

CentroidsYdecision = CentroidsY*wY;
DataYdecision      = Data.Y*wY;
[~,LabelsYtest_lda] = getNearestCluster(DataYdecision(Results_splitting.selection_testdata,:),CentroidsYdecision);

% Plot classification
f1=figure;
plot(Labels,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
plot(Results_splitting.selection_testdata,Labels(Results_splitting.selection_testdata),'o','Color',0.8.*ones(1,3),'MarkerFaceColor',0.8.*ones(1,3),'MarkerSize',10)
plot(Results_splitting.selection_testdata,LabelsYtest_centroids,'^b','MarkerSize',8)
plot(Results_splitting.selection_testdata,LabelsYtest_lda,'xr','MarkerSize',7)
legend('full state', 'full state (test)', 'centroids', 'lda')
axis([0 M+100 0 Nclusters+1])
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 600 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2save,['ClassificationTimeSeries.eps']]);
close(f1);

TF = LabelsYtest_centroids == Labels(Results_splitting.selection_testdata);
ClassAccuracy_centroids = length(TF(TF==1))/length(Results_splitting.selection_testdata)

TF = LabelsYtest_lda == Labels(Results_splitting.selection_testdata);
ClassAccuracy_lda = length(TF(TF==1))/length(Results_splitting.selection_testdata)

% Plot vorticity with sensors
Omega = reshape(Data.Omega,[Data.N1 Data.N2 Data.N3]); Omega = Omega(:,:,m);
cmin = -2; cmax = 2; cvals = [cmin:0.5:cmax]; m = 10;
filename = ['Vorticity_m',sprintf('%03g',m),'_SSPOCsensors'];
fhandle = figure;
contourf(Data.x,Data.y,Omega,cvals);caxis([cmin cmax]),hold on
quiver(Data.x,Data.y,Data.U(1:Data.N1,:,m),Data.U(Data.N1+1:end,:,m),'k')
[X,Y] = meshgrid(Data.x,Data.y); 
X = reshape(X,[Data.N1*Data.N2, 1]); Y = reshape(Y,[Data.N1*Data.N2, 1]);
scatter(X(sensors_locations), Y(sensors_locations),5,'r')
colormap(cmap),
xlabel('x'), ylabel('y')
box on
daspect([1 1 1])
set(gca,'FontSize',14,'LineWidth',1)
set(gcf,'Position',[100 100 300 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', fullfile(path2save,filename));
close(fhandle);