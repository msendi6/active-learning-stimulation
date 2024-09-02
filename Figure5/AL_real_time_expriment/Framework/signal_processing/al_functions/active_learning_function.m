function[next_stim_param]=active_learning_function(params_space,old_stim_params,brain_res, reg_type,query_type)

%% inputs
% params_space:the entire parameter space we are going to serach in
%old_stim_params: the latest  stimulation parameters we adjusted
% brain_res: the post-stimulation signal recored from the brain
% reg_type: 1-ridge, 2-lasso, 3-Decision tree regression 
%query_type:1-EMCM, 2-QBC, 3- GS, 4-RD, 5-RD+EMCM, 6-RD+QBC, 7-RD+GS

%%outputs
%next_stim_param: the next stimulation parameters that AL suggested
%new_params_space: is the parameters space for the next search 
 
X_obs=old_stim_params; 
y_obs=brain_res;

if query_type==1
%% 1-Select new samples by EMCM
nBoots=4; % number of repeats to get statistically significant results
[next_stim_param]=EMCM(X_obs,y_obs,params_space,reg_type,nBoots);
elseif query_type==2
%% 2-Select new samples by QBC
nBoots=4; % number of repeats to get statistically significant results
[next_stim_param]=QBC(X_obs,y_obs,params_space,reg_type,nBoots);
elseif query_type==3
%% 3-Select new samples by Greedy_sampling
[next_stim_param]=Greedy_sampling(X_obs,params_space);
elseif query_type==4
%% 4-Select new samples by RD
[next_stim_param]=RD_naive(X_obs,params_space);
elseif query_type==5
%% 5-Select new sample by RD+EMCM
nBoots=4;
[next_stim_param]=RD_EMCM(X_obs,y_obs,params_space,reg_type,nBoots);
elseif query_type==6
%% 6-Select new sample by RD+QBC
nBoots=4;
[next_stim_param]=RD_QBC(X_obs,y_obs,params_space,reg_type,nBoots);
elseif query_type==7
%% 7-Select new sample by RD+GS
[next_stim_param]=RD_greedy_sampling(X_obs,params_space);
end

%% Query Strategy
%EMCM
function[next_stim]=EMCM(X,y,future_X,regression_type,n_of_Boots)
if regression_type==1 % ridge regression
k=0.1;
[coef]=ridge_reg(y,X,k);
y_pred=future_X*coef(2:end) + coef(1);
 
elseif regression_type==2 % SVR
svrMdl=svr_reg(y,X);
y_pred=predict(svrMdl,future_X);
 
elseif regression_type==3 % gpr
[gprMdl]=gpr_reg(y,X);
y_pred=predict(gprMdl,future_X); 

elseif regression_type==4 % Decision tree regression
 [treeMdl]=Dtree_reg(y,X);
 y_pred=predict(treeMdl,future_X);
end


all_obs_data=horzcat(X,y);
for i=1:n_of_Boots
 y_bootstraped= datasample(all_obs_data,size(all_obs_data,1),'Replace',true); 
 
 if regression_type==1 % ridge regression
 k=0.1;
 [coef]=ridge_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1),k);
 y_pred_B(:,i)=future_X*coef(2:end) + coef(1);

 elseif regression_type==2 % SVR
 svrMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(svrMdl,future_X);

 elseif regression_type==3 % gpr
 gprMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(gprMdl,future_X);

 elseif regression_type==4 % Decision tree regression
 [treeMdl]=Dtree_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(treeMdl,future_X);
 end
 
end
MEMC=zeros(1,length(future_X));

for i=1:length(future_X)
                for j=1:n_of_Boots
                    MEMC(i)=MEMC(i)+norm((y_pred_B(i,j)-y_pred(i))*future_X(i,:));
                end
end
[~,ids]=sort(MEMC,'descend');
next_stim=future_X(ids(1),:);
% future_X(ids(1),:)=[];
% new_future_X=future_X;
end
%QBC
function[next_stim]=QBC(X,y,future_X,regression_type,n_of_Boots)
all_obs_data=horzcat(X,y);
for i=1:n_of_Boots
 y_bootstraped= datasample(all_obs_data,size(all_obs_data,1),'Replace',true); 
 
 if regression_type==1 % ridge regression
 k=0.1;
 [coef]=ridge_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1),k);
 y_pred_B(:,i)=future_X*coef(2:end) + coef(1);

 elseif regression_type==2 % SVR
 svrMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(svrMdl,future_X);

 elseif regression_type==3 % gpr
 gprMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(gprMdl,future_X);

 elseif regression_type==4 % Decision tree regression
 [treeMdl]=Dtree_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(treeMdl,future_X);
 end
 
end
QBC_var=zeros(1,length(future_X));

for i=1:length(future_X)
QBC_var(i)=var(y_pred_B(i,:));
end
[~,ids]=sort(QBC_var,'descend');
next_stim=future_X(ids(1),:);
% future_X(ids(1),:)=[];
% new_future_X=future_X;

end
%GS
function[next_stim]=Greedy_sampling(X,future_X)

GS=zeros(length(future_X),length(X));

for i=1:length(future_X)
for   j=1:length(X)
    GS(i,j)=sqrt(sum((future_X(i,:) -X(j,:)) .^ 2));
end
end

dist=min(GS,[],2);
[~,ids]=sort(dist,'descend');
next_stim=future_X(ids(1),:);

end
%RD
function[next_stim]=RD_naive(X,future_X)

[target_space]=RD_target_space(X,future_X);

[~,~,~,D]=kmeans(target_space,1,'Replicates',10);
[~,ids]=sort(D,'descend');
next_stim=target_space(ids(1),:);


end
%RD+EMCM
function [next_stim]=RD_EMCM(X,y,future_X,regression_type,n_of_Boots)
[target_space]=RD_target_space(X,future_X);

if regression_type==1 % ridge regression
k=0.1;
[coef]=ridge_reg(y,X,k);
y_pred=target_space*coef(2:end) + coef(1);
 
elseif regression_type==2 % SVR
svrMdl=svr_reg(y,X);
y_pred=predict(svrMdl,target_space);
 
elseif regression_type==3 % gpr
[gprMdl]=gpr_reg(y,X);
y_pred=predict(gprMdl,target_space); 

elseif regression_type==4 % Decision tree regression
 [treeMdl]=Dtree_reg(y,X);
 y_pred=predict(treeMdl,target_space);
end

all_obs_data=horzcat(X,y);
for i=1:n_of_Boots
 y_bootstraped= datasample(all_obs_data,size(all_obs_data,1),'Replace',true); 
 
 if regression_type==1 % ridge regression
 k=0.1;
 [coef]=ridge_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1),k);
 y_pred_B(:,i)=target_space*coef(2:end) + coef(1);

 elseif regression_type==2 % SVR
 svrMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(svrMdl,target_space);

 elseif regression_type==3 % gpr
 gprMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(gprMdl,target_space);

 elseif regression_type==4 % Decision tree regression
 [treeMdl]=Dtree_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(treeMdl,target_space);
 end
 
end


MEMC=zeros(1,length(target_space));

for i=1:length(target_space)
                for j=1:n_of_Boots
                    MEMC(i)=MEMC(i)+norm((y_pred_B(i,j)-y_pred(i))*target_space(i,:));
                end
end
[~,ids]=sort(MEMC,'descend');
next_stim=target_space(ids(1),:);


end
%RD+QBC
function [next_stim]=RD_QBC(X,y,future_X,regression_type,n_of_Boots)
[target_space]=RD_target_space(X,future_X);

all_obs_data=horzcat(X,y);
for i=1:n_of_Boots
 y_bootstraped= datasample(all_obs_data,size(all_obs_data,1),'Replace',true); 
 
 if regression_type==1 % ridge regression
 k=0.1;
 [coef]=ridge_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1),k);
 y_pred_B(:,i)=target_space*coef(2:end) + coef(1);

 elseif regression_type==2 % SVR
 svrMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(svrMdl,target_space);

 elseif regression_type==3 % gpr
 gprMdl=svr_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(gprMdl,target_space);

 elseif regression_type==4 % Decision tree regression
 [treeMdl]=Dtree_reg(y_bootstraped(:,end),y_bootstraped(:,1:end-1));
 y_pred_B(:,i)=predict(treeMdl,target_space);
 end
 
end
QBC_var=zeros(1,length(target_space));

for i=1:length(target_space)
QBC_var(i)=var(y_pred_B(i,:));
end
[~,ids]=sort(QBC_var,'descend');
next_stim=target_space(ids(1),:);

end
%RD+GS
function [next_stim]=RD_greedy_sampling(X,future_X)
[target_space]=RD_target_space(X,future_X);

GS=zeros(length(target_space),length(future_X));
if length(target_space)==1
next_stim=target_space
else 
for i=1:length(target_space)
for   j=1:length(future_X)
    GS(i,j)=sqrt(sum((target_space(i,:) - future_X(j,:)) .^ 2));
end
end
dist=min(GS,[],2);
[~,ids]=sort(dist,'descend');
next_stim=target_space(ids(1),:);
end
end

%% RD to find the target space
function [search_area]=RD_target_space(X,future_X)
d=size(X,1);
all_param=[X;future_X];
[ids,~,~,D]=kmeans(all_param,d+1,'Distance','cityblock'); %%% 
idsP=cell(1,d+1);
numP=zeros(1,d+1);
for i=1:d+1
idsP{i}=find(ids==i);
end
search_area=[];
for i=1:d+1
    A= all_param(idsP{1,i},:);
    S=0;
    for j=1:d
    I = ismember(A,X(j,:),'rows');
    S=S+sum(I);   
    
    end
    if S==0
     search_area=vertcat(search_area,A);
    end
    
end

end
%% Regression models
%ridge, lasso, elastic net
function[coef0,coef]=lasso_reg(y,X,alpha,cross_validation)
[B,FitInfo] = lasso(X,y,'Alpha',alpha,'CV',cross_validation);
idxLambda1SE = FitInfo.Index1SE;
coef = B(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);
end
% decision tree
function[treeMdl]=Dtree_reg(y,X)
treeMdl = fitrtree(X,y);
end

% ridge regression
function[b]=ridge_reg(y,X,k)
b = ridge(y,X,k,0);
end

% Support vector regression
function[svrMdl]=svr_reg(y,X)
svrMdl=fitrsvm(X,y,'Standardize',true,'KernelFunction','gaussian');
end

% Guassian process regression
function[gprMdl]=gpr_reg(y,X)

gprMdl = fitrgp(X,y,'Basis','linear','FitMethod','exact','PredictMethod','exact');

end
end
