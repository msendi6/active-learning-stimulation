%Run each section separately 
clear
%Figure 3C
close all
clc
load('~/Figure2/results/NHP_RR_S1.mat')
all_reg_cc=[]
all_reg_rmse=[]

all_pval_rmse=[]
all_pval_cc=[]

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end


pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');



all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

load('~/Figure3/results/NHP_RR_brain_state_dynamic.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));


load('~/Figure3/results/NHP_RR_stochastic_response.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

%Fig.2a 
figure
data=all_reg_rmse-all_reg_rmse(:,1);
data(:,1)=[]

data(all_pval_rmse>0.05/9)=0;
imagesc(data)
hold on

plot([11 0],[1.5 1.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[3.5 3.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([1.5 1.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([2.5 2.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([3.5 3.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([4.5 4.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([5.5 5.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([6.5 6.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([7.5 7.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([8.5 8.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([9.5 9.5],[11 0],'color',[1 1 1], 'linewidth',7)
xticks([])
yticks([])
width=800;
height=250;
set(gcf,'position',[0,0,width,height])
caxis([-2 ,2])
colormap(redblue)
set(gca,'visible','off')
set(gca,'XColor', 'none','YColor','none')



%%
clear
%Figure 3D
close all
clc
load('~/Figure2/results/NHP_SVR_LIN_S1.mat')
all_reg_cc=[]
all_reg_rmse=[]

all_pval_rmse=[]
all_pval_cc=[]

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end


pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');



all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

load('~/Figure3/results/NHP_SVR_LIN_brain_state_dynamic.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));


load('~/Figure3/results/NHP_SVR_LIN_stochastic_response.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

%Fig.2a 
figure
data=all_reg_rmse-all_reg_rmse(:,1);
data(:,1)=[]

data(all_pval_rmse>0.05/9)=0;
imagesc(data)
hold on

plot([11 0],[1.5 1.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[3.5 3.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([1.5 1.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([2.5 2.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([3.5 3.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([4.5 4.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([5.5 5.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([6.5 6.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([7.5 7.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([8.5 8.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([9.5 9.5],[11 0],'color',[1 1 1], 'linewidth',7)
xticks([])
yticks([])
width=800;
height=250;
set(gcf,'position',[0,0,width,height])
caxis([-2 ,2])
colormap(redblue)
set(gca,'visible','off')
set(gca,'XColor', 'none','YColor','none')



%%
%Figure 3E
clear
close all
clc
load('~/Figure2/results/NHP_SVR_RBF_S1.mat')
all_reg_cc=[]
all_reg_rmse=[]

all_pval_rmse=[]
all_pval_cc=[]

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end


pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');



all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

load('~/Figure3/results/NHP_SVR_RBF_brain_state_dynamic.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));


load('~/Figure3/results/NHP_SVR_RBF_stochastic_response.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

%Fig.2a 
figure
data=all_reg_rmse-all_reg_rmse(:,1);
data(:,1)=[]

data(all_pval_rmse>0.05/9)=0;
imagesc(data)
hold on

plot([11 0],[1.5 1.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[3.5 3.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([1.5 1.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([2.5 2.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([3.5 3.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([4.5 4.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([5.5 5.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([6.5 6.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([7.5 7.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([8.5 8.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([9.5 9.5],[11 0],'color',[1 1 1], 'linewidth',7)
xticks([])
yticks([])
width=800;
height=250;
set(gcf,'position',[0,0,width,height])
caxis([-2 ,2])
colormap(redblue)
set(gca,'visible','off')
set(gca,'XColor', 'none','YColor','none')


%%
%Figure 3F
clear
close all
clc
load('~/Figure2/results/NHP_GPR_S1.mat')
all_reg_cc=[]
all_reg_rmse=[]

all_pval_rmse=[]
all_pval_cc=[]

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end


pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');



all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

load('~/Figure3/results/NHP_GPR_brain_state_dynamic.mat')

auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end
%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

load('~/Figure3/results/NHP_GPR_stochastic_response.mat')
load('/Users/msendi/Library/CloudStorage/OneDrive-MassGeneralBrigham/TReNDS work/active learning paperv2/For_manuscript/Figure3/results/NHP_GPR_dyn2.mat')
auc_rmse_all=[]
auc_cc_all=[]
for i=1:10
mse=MSEs{1,1};

corcof=CCs{1,1};

auc_rmse_all(:,i)=squeeze(nansum(mse(i,1:end,:),2));

auc_cc_all(:,i)=squeeze(nansum(corcof(i,1:end,:),2));
end

auc_rmse_all(175:end,:)=[]
auc_cc_all(175:end,:)=[]

%removing outliear from rmse
for i=1:10

    %%removing outlier
data=auc_rmse_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_rmse_all(:,i)=data;
end
%removing outliear from cc
for i=1:10

    %%removing outlier
data=auc_cc_all(:,i);
% Calculate the first quartile (Q1) and third quartile (Q3)
Q1 = quantile(data, 0.25);
Q3 = quantile(data, 0.75);
% Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1;
% Define a multiplier for identifying outliers (e.g., 1.5 times IQR)
multiplier = 1.5;
% Find the indices of outliers
outlier_indices = (data < Q1 - multiplier * IQR) | (data > Q3 + multiplier * IQR);
% Remove outliers from the data
data(outlier_indices)=NaN
auc_cc_all(:,i)=data;
end

pval=[]
for j=2:10
    [h,p]=ttest2(auc_rmse_all(:,1),auc_rmse_all(:,j));
    pval(j-1)=p
end
pval=pval'

all_pval_rmse=vertcat(all_pval_rmse,pval');

pval=[]
for j=2:10
    [h,p]=ttest2(auc_cc_all(:,1),auc_cc_all(:,j));
    pval(j-1)=p
end
pval=pval'
all_pval_cc=vertcat(all_pval_cc,pval');

all_reg_cc=vertcat(all_reg_cc,nanmean(auc_cc_all));
all_reg_rmse=vertcat(all_reg_rmse,nanmean(auc_rmse_all));

%Fig.2a 
figure
data=all_reg_rmse-all_reg_rmse(:,1);
data(:,1)=[]

data(all_pval_rmse>0.05/9)=0;
imagesc(data)
hold on

plot([11 0],[1.5 1.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[3.5 3.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([1.5 1.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([2.5 2.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([3.5 3.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([4.5 4.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([5.5 5.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([6.5 6.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([7.5 7.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([8.5 8.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([9.5 9.5],[11 0],'color',[1 1 1], 'linewidth',7)
xticks([])
yticks([])
width=800;
height=250;
set(gcf,'position',[0,0,width,height])
caxis([-2 ,2])
colormap(redblue)
set(gca,'visible','off')
set(gca,'XColor', 'none','YColor','none')

figure
data=all_reg_cc-all_reg_cc(:,1);
data(:,1)=[]
data(all_pval_cc>0.05/9)=0;
imagesc(data)
hold on

plot([11 0],[1.5 1.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[3.5 3.5],'color',[1 1 1], 'linewidth',7)
plot([11 0],[2.5 2.5],'color',[1 1 1], 'linewidth',7)
plot([1.5 1.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([2.5 2.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([3.5 3.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([4.5 4.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([5.5 5.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([6.5 6.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([7.5 7.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([8.5 8.5],[11 0],'color',[1 1 1], 'linewidth',7)
plot([9.5 9.5],[11 0],'color',[1 1 1], 'linewidth',7)
xticks([])
yticks([])
width=800;
height=250;
set(gcf,'position',[0,0,width,height])
caxis([-0.5 ,.5])
%colormap(othercolor('BuOrR_14'))
colormap(redblue)
set(gca,'visible','off')
set(gca,'XColor', 'none','YColor','none')
