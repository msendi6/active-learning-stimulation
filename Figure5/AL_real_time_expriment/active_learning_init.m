function[init_stim_param,sample_idx]=active_learning_init(params_space)
% a code find the best parameters to initialize the sampling
%. This initialization ensures representativeness
% this function introduce n (is the dimension of parameter space) initializing sample. 
%Inputs:
% Param_space: is the entire stim param space we are going to sample from!
%Outputs:
%init_stim_param: initial parameters we will start with. 

n=size(params_space,2)+1;

[ids,~,~,D]=kmeans(params_space,n,'Distance','cityblock','Replicates',10);

 idsP=cell(1,n);
 numP=zeros(1,n);
 for i=1:n
 idsP{i}=find(ids==i);
 numP(i)=length(idsP{i});
 end
 [~,ids2]=sort(numP,'descend');
 idsP=idsP(ids2);
 sample_idx=zeros(1,n);
 for k=1:n
 dist=D(ids==ids2(k),ids2(k));
 [~,idx]=min(dist);
 sample_idx(k)=idsP{k}(idx);
 end
 
 init_stim_param=params_space(sample_idx,:);
end
