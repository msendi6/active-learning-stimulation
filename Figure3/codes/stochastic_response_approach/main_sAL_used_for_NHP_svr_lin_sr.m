%% Implementation of the AL algorithms in:
%%
%% D. Wu, "Pool-based sequential active learning for regression," IEEE Trans. on Neural 
%% Networks and Learning Systems, 30(5), pp. 1348-1359, 2019.
%%
%% Compare 10 methods:
%% 1. BL
%% 2. EMCM
%% 3. RD + EMCM
%% 4. QBC
%% 5. RD + QBC
%% 6. GS
%% 7. RD + GS
%% 8. RD only
%% 9. EQBC
%% 10. EEMCM

%% Mo Sendi msendi@mclean.harvard.edu adapted from a code develpoed by Dongrui WU, drwu@hust.edu.cn
method_name={'RS', 'EMCM', 'RD+EMCM','QBC','RD+QBC','GS','RD+GS','RD','EQBC','EEMCM'};
clc; 
clearvars; 
close all;   rng('default');
datasets={'NHP-theta'};
nRepeat=1000; % number of repeats to get statistically significant results;  100 was used in the paper
rr=.1; % RR parameter
nBoots=4;
numAlgs=10;
MSEs=cell(1,length(datasets)); CCs=MSEs;

load('~/Figure3/data/stochastic_response_data.mat')
pool_train_n=80;
test_n=20;
warning off

%%

for s=1:length(datasets)   


    
    numY0=pool_train_n;
    minN=3; % mininum number of training samples
    maxN=20; % maximum number of training samples
    MSEs{s}=nan(numAlgs,maxN,nRepeat); CCs{s}=MSEs{s};
    
    %% Iterative approaches; random sampling
    for r=1:nRepeat
        [s r]

        y=[];
        for i=1:length(data)
        y{i,1} = datasample(data(i,4:end),1);
        end
        y=cell2mat(y);
      
        dataU=horzcat(data(:,1:3),y);
        X_all=dataU(:,1:3); y_all=dataU(:,4);
        % Get the total number of samples
        totalSamples = size(X_all, 1);
        % Specify the number of samples you want to select
        numSamplesToSelect = test_n;

        % Generate random indices for selecting samples
        idsTest = randperm(totalSamples, numSamplesToSelect);

        % Use the random indices to select the samples from X and y
        XTest = X_all(idsTest, :);
        YTest = y_all(idsTest, :);

        % Create a logical index for the samples that were not selected
        logicalIndexNotSelected = true(totalSamples, 1);
        logicalIndexNotSelected(idsTest) = false;

        % Use the logical index to select the samples not selected
        X_pool_train = X_all(logicalIndexNotSelected, :);
        Y_pool_train = y_all(logicalIndexNotSelected, :);

        X0=X_pool_train;
        Y0=Y_pool_train;
         %% Pre-compute distances for GS
        distX0=zeros(numY0);
       for i=1:numY0
        M=X0-repmat(X0(i,:),numY0,1);
        distX0(:,i)=sqrt(sum(M.^2,2));
        end
        %% random effect; 80% samples
        ids=datasample(1:numY0,round(numY0*0.8),'Replace',false);
        X=X0(ids,:); Y=Y0(ids); numY=length(Y);
        
        numYTest=length(YTest);
        idsTrain=repmat(datasample(1:numY,maxN,'Replace',false),numAlgs,1);
        distX=distX0(ids,ids);
        
        %% For EQBC and EEMCM
        idxOutlierGroup=1; idsGood9=1:numY; numEle=zeros(1,minN);
        while ~isempty(idxOutlierGroup)
            ids=kmeans(X(idsGood9,:),minN);
            for i=1:minN
                numEle(i)=sum(ids==i);
            end
            idxOutlierGroup=find(numEle<.02*length(idsGood9));
            idsOutliers=[];
            for i=1:length(idxOutlierGroup)
                idsOutliers=cat(2,idsOutliers,find(ids==idxOutlierGroup(i)));
            end
            idsGood9(idsOutliers)=[];
        end
        idsGood10=idsGood9;
        
        for n=minN:maxN
            %% 1. BL: No AL
            svrMdl1=fitrsvm(X(idsTrain(1,1:n),:),Y(idsTrain(1,1:n)),'Standardize',true,'KernelFunction','linear');
            Y1=predict(svrMdl1,XTest);
            MSEs{s}(1,n,r)=sqrt(mean((Y1-YTest).^2));
            CCs{s}(1,n,r)=corr(Y1,YTest);
            
            %% 2. EMCM; the first minN samples were obtained randomly

            svrMdl2=fitrsvm(X(idsTrain(2,1:n),:),Y(idsTrain(2,1:n)),'Standardize',true,'KernelFunction','linear');
            Y2=predict(svrMdl2,XTest);
            MSEs{s}(2,n,r)=sqrt(mean((Y2-YTest).^2));
            CCs{s}(2,n,r)=corr(Y2,YTest);
            %% Select new samples by EMCM
            idsUnlabeled=1:numY; idsUnlabeled(idsTrain(2,1:n))=[];
            C=max(1,ceil(n*rand(nBoots,n)));
            Ys=repmat(Y,1,nBoots);
            for i=1:nBoots
            svrMdl=fitrsvm(X(idsTrain(2,C(i,:)),:),Y(idsTrain(2,C(i,:))),'Standardize',true,'KernelFunction','linear');
            Ys(idsUnlabeled,i)=predict(svrMdl,X(idsUnlabeled,:));

            end
            EMCM=zeros(1,length(idsUnlabeled));
            Y2=Y; Y2(idsUnlabeled)=predict(svrMdl2,X(idsUnlabeled,:));
            
            for i=1:length(idsUnlabeled)
                for j=1:nBoots
                    EMCM(i)=EMCM(i)+norm((Ys(idsUnlabeled(i),j)-Y2(idsUnlabeled(i)))*X(idsUnlabeled(i),:));
                end
            end
            [~,idx]=max(EMCM);
            idsTrain(2,n+1)=idsUnlabeled(idx);
            
            %% 3. RD + EMCM
            
            svrMdl3=fitrsvm(X(idsTrain(3,1:n),:),Y(idsTrain(3,1:n)),'Standardize',true,'KernelFunction','linear');
            Y3=predict(svrMdl3,XTest);
            MSEs{s}(3,n,r)=sqrt(mean((Y3-YTest).^2));
            CCs{s}(3,n,r)=corr(Y3,YTest);
            %% Select the new sample by RD + EMCM
            idsUnlabeled=1:numY; idsUnlabeled(idsTrain(3,1:n))=[];
            Y3=Y; Y3(idsUnlabeled)=predict(svrMdl3,X(idsUnlabeled,:));
            [ids,~,~,D]=kmeans(X,n+1,'Replicates',10); %%% Study the effect of n+1 and n+2
            idsP=cell(1,n+1);
            numP=zeros(1,n+1);
            for i=1:n+1
                idsP{i}=find(ids==i);
                numP(i)=length(idsP{i});
            end
            [numP,ids2]=sort(numP,'descend');
            idsP=idsP(ids2); D=D(:,ids2);
            for k=1:n+1
                if sum(ismember(idsTrain(3,1:n),idsP{k}))==0
                    Ys=repmat(Y(idsP{k}),1,nBoots);
                    for i=1:nBoots

                        svrMdl=fitrsvm(X(idsTrain(3,C(i,:)),:),Y(idsTrain(3,C(i,:))),'Standardize',true,'KernelFunction','linear');
                        Ys(:,i)=predict(svrMdl,X(idsP{k},:));

                    end
                    EMCM=zeros(1,numP(k));
                    for i=1:numP(k)
                        for j=1:nBoots
                            EMCM(i)=EMCM(i)+norm((Ys(i,j)-Y3(idsP{k}(i)))*X(idsP{k}(i),:));
                        end
                    end
                    [~,idx]=max(EMCM);
                    idsTrain(3,n+1)=idsP{k}(idx);
                    break;
                end
            end
            
            %% 4. QBC; the first minN samples were obtained randomly
            svrMdl4=fitrsvm(X(idsTrain(4,1:n),:),Y(idsTrain(4,1:n)),'Standardize',true,'KernelFunction','linear');
            Y4=predict(svrMdl4,XTest);
            MSEs{s}(4,n,r)=sqrt(mean((Y4-YTest).^2));
            CCs{s}(4,n,r)=corr(Y4,YTest);
            %% Select new samples by QBC
            Ys=repmat(Y,1,nBoots);
            idsUnlabeled=1:numY; idsUnlabeled(idsTrain(4,1:n))=[];
            for i=1:nBoots            
            svrMdl=fitrsvm(X(idsTrain(4,C(i,:)),:),Y(idsTrain(4,C(i,:))),'Standardize',true,'KernelFunction','linear');
            Ys(idsUnlabeled,i)=predict(svrMdl,X(idsUnlabeled,:));
            end
            QBC=zeros(1,length(idsUnlabeled));
            for i=1:length(idsUnlabeled)
                QBC(i)=var(Ys(idsUnlabeled(i),:));
            end
            [~,idx]=max(QBC);
            idsTrain(4,n+1)=idsUnlabeled(idx);
            
            
            %% 5. RD + QBC

            svrMdl5=fitrsvm(X(idsTrain(5,1:n),:),Y(idsTrain(5,1:n)),'Standardize',true,'KernelFunction','linear');
            Y5=predict(svrMdl5,XTest);
            MSEs{s}(5,n,r)=sqrt(mean((Y5-YTest).^2));
            CCs{s}(5,n,r)=corr(Y5,YTest);
            %% Select the new sample by RD + QBC
            for k=1:n+1
                if sum(ismember(idsTrain(5,1:n),idsP{k}))==0
                    Ys=repmat(Y(idsP{k}),1,nBoots);
                    for i=1:nBoots
                        svrMdl=fitrsvm(X(idsTrain(5,C(i,:)),:),Y(idsTrain(5,C(i,:))),'Standardize',true,'KernelFunction','linear');
                        Ys(:,i)=predict(svrMdl,X(idsP{k},:));    
 
                    end
                    QBC=zeros(1,numP(k));
                    for i=1:numP(k)
                        QBC(i)=var(Ys(i,:));
                    end
                    [~,idx]=max(QBC);
                    idsTrain(5,n+1)=idsP{k}(idx);
                    break;
                end
            end
            
            
            %% 6. GS
            svrMdl6=fitrsvm(X(idsTrain(6,1:n),:),Y(idsTrain(6,1:n)),'Standardize',true,'KernelFunction','linear');
            Y6=predict(svrMdl6,XTest); 
            MSEs{s}(6,n,r)=sqrt(mean((Y6-YTest).^2));
            CCs{s}(6,n,r)=corr(Y6,YTest);
            %% Select new samples by GS
            idsUnlabeled=1:numY; idsUnlabeled(idsTrain(6,1:n))=[];
            dist=min(distX(idsUnlabeled,idsTrain(6,1:n)),[],2);
            [~,idx]=max(dist);
            idsTrain(6,n+1)=idsUnlabeled(idx);
            
            %% 7. GS with new initialization

            svrMdl7=fitrsvm(X(idsTrain(7,1:n),:),Y(idsTrain(7,1:n)),'Standardize',true,'KernelFunction','linear');
            Y7=predict(svrMdl7,XTest);
            MSEs{s}(7,n,r)=sqrt(mean((Y7-YTest).^2));
            CCs{s}(7,n,r)=corr(Y7,YTest);
            %% Select the new sample by RD + GS
            for k=1:n+1
                if sum(ismember(idsTrain(7,1:n),idsP{k}))==0
                    dists=zeros(numP(k),n);
                    for i=1:n
                        M=X(idsP{k},:)-repmat(X(idsTrain(7,i),:),numP(k),1);
                        dists(:,i)=sqrt(sum(M.^2,2));
                    end
                    dist=min(dists,[],2);
                    [~,idx]=max(dist);
                    idsTrain(7,n+1)=idsP{k}(idx);
                    break;
                end
            end
            
            %% 8. RD only
           
            svrMdl8=fitrsvm(X(idsTrain(8,1:n),:),Y(idsTrain(8,1:n)),'Standardize',true,'KernelFunction','linear');
            Y8=predict(svrMdl8,XTest);
            MSEs{s}(8,n,r)=sqrt(mean((Y8-YTest).^2));
            CCs{s}(8,n,r)=corr(Y8,YTest);
            %% Select the new sample by RD
            for k=1:n+1
                if sum(ismember(idsTrain(8,1:n),idsP{k}))==0
                    dist=D(idsP{k},k);
                    [~,idx]=min(dist);
                    idsTrain(8,n+1)=idsP{k}(idx);
                    break;
                end
            end
            
            %% 9. EQBC
            svrMdl9=fitrsvm(X(idsTrain(9,1:n),:),Y(idsTrain(9,1:n)),'Standardize',true,'KernelFunction','linear');
            Y9=predict(svrMdl9,XTest);
            MSEs{s}(9,n,r)=sqrt(mean((Y9-YTest).^2));
            CCs{s}(9,n,r)=corr(Y9,YTest);
            %% Select a new sample by EQBC
            Ys=repmat(Y,1,nBoots);
            for i=1:nBoots
               
             svrMdl=fitrsvm(X(idsTrain(9,C(i,:)),:),Y(idsTrain(9,C(i,:))),'Standardize',true,'KernelFunction','linear');
             Ys(idsGood9,i)=predict(svrMdl,X(idsGood9,:));

            end
            QBC=zeros(1,length(idsGood9));
            for i=1:length(idsGood9)
                QBC(i)=var(Ys(idsGood9(i),:));
            end
            [~,idx]=max(QBC);
            idsTrain(9,n+1)=idsGood9(idx);
            idsGood9(idx)=[];
            
            %% 10. EEMCM
%             if n==minN
%                 idsTrain(10,1:n)=idsTrain(9,1:n); idsGood10=idsGood9;
%             end
            svrMdl10=fitrsvm(X(idsTrain(10,1:n),:),Y(idsTrain(10,1:n)),'Standardize',true,'KernelFunction','linear');
            Y10=predict(svrMdl10,XTest);
            MSEs{s}(10,n,r)=sqrt(mean((Y10-YTest).^2));
            CCs{s}(10,n,r)=corr(Y10,YTest);
            %% Select a new sample by EEMCM
            idsUnlabeled=1:numY; idsUnlabeled(idsTrain(10,1:n))=[];
            Y10=Y; Y10(idsUnlabeled)=predict(svrMdl10,X(idsUnlabeled,:)); 
            C=max(1,ceil(n*rand(nBoots,n)));
            Ys=repmat(Y,1,nBoots);
            for i=1:nBoots
             svrMdl=fitrsvm(X(idsTrain(10,C(i,:)),:),Y(idsTrain(10,C(i,:))),'Standardize',true,'KernelFunction','linear');
             Ys(idsGood10,i)=predict(svrMdl,X(idsGood10,:));    
 

            end
            EMCM=zeros(1,length(idsGood10));
            for i=1:length(idsGood10)
                for j=1:nBoots
                    EMCM(i)=EMCM(i)+norm((Ys(idsGood10(i),j)-Y10(idsGood10(i)))*X(idsGood10(i),:));
                end
            end
            [~,idx]=max(EMCM);
            idsTrain(10,n+1)=idsGood10(idx);
            idsGood10(idx)=[];
            
        end
    end
    
    %% Plot results
 close all
    AUCMSE=zeros(numAlgs,length(datasets)); AUCCC=AUCMSE;
    linestyle={'k-','r-','r--','g-','g--','b-','b--','k--','m-','m--'};
    ids=[1 2 3 4 5 6 7 8 9 10];
    color{1}=[0,0,255]/255;
    color{2}=[191,191,0]/255;
    color{3}=[237,177,32]/255;
    color{4}=[0,127,0]/255;
    color{5}=[217,83,25]/255;
    color{6}=[191,0,191]/255;
    color{7}=[255,0,0]/255;
    color{8}=[0,204,204]/255;
    color{9}=[153,51,0]/255;
    color{10}=[255,0,255]/255;

    numY0=length(Y0);
    minN=size(X0,2); % mininum number of training samples
    maxN=min(60,max(20,ceil(.1*numY0))); % maximum number of training samples
    mMSEs=squeeze(nanmean(MSEs{s},3)); AUCMSE(:,s)=nansum(mMSEs,2);
    mCCs=squeeze(nanmean(CCs{s},3)); AUCCC(:,s)=nansum(mCCs,2);
    
    figure;
    set(gcf,'DefaulttextFontName','times new roman','DefaultaxesFontName','times new roman');
    subplot(121); 
    hold on;
    for i=ids
        plot(minN:maxN,mMSEs(i,minN:maxN),'Color',color{i},'linewidth',3);
    end
    h=legend('RS','EMCM','RD+EMCM','QBC','RD+QBC','GS','RD+GS','RD','EQBC','EEMCM','location','northeast');
    set(h,'fontsize',15);
    axis tight; box off;
    xlabel('N','fontsize',12);     ylabel('RMSE','fontsize',12);
    
    subplot(122); hold on;
    for i=ids
        plot(minN:maxN,mCCs(i,minN:maxN),'Color',color{i},'linewidth',3);
    end
    set(h,'fontsize',11);
    axis tight; box off ;
    xlabel('N','fontsize',12);     
    ylabel('CC','fontsize',12);

end

save('~/Figure3/results/NHP_SVR_LIN_stochastic_response','MSEs','CCs') %please change the name to avoid overwritting the shared results

