function [ net,rsq,sig,prsq ] = regression_network(inputs,targets,nhidden,outlayer,est)
%REGRESSION_NETWROK generate and train feedforward regression NN 
numLayers=2;
if(nargin<4)
    outlayer='tansig';
end
if(nargin<5)
    est=0;
end

%% remove nans from inputs and targets
% ind=find(~isnan(sum(inputs,1)) & ~isnan(targets));
% inputs=inputs(:,ind);
% targets=targets(ind);
%% generate network
net = feedforwardnet(nhidden);
% net = patternnet(nhidden);
% net.layers{2}.transferFcn = outlayer;
%% train network
net = train(net,inputs,targets,'useParallel','yes');%,'useGPU','yes');
outputs=net(inputs);
rsq=nancorr(outputs(:),targets(:))^2;
sig=nanstd(outputs(:)-targets(:));

%% estimate degrees of freedom
if(est)
%     enet = feedforwardnet(nhidden);
%     enet.layers{2}.transferFcn = outlayer;
    N=size(inputs,1);    
    for i=1:N %go over each input
        I=inputs;
        I(i,:)=nanmedian(I(i,:));
%         I(i,:)=I(i,randperm(size(I,2))); %shuffle
        outputs=net(I);
        prsq(i)=(rsq-nancorr(outputs(:),targets(:))^2)/rsq;
%         ind=setdiff(1:N,i);
%         enet = train(enet,inputs(ind,:),targets,'useParallel','yes');%,'useGPU','yes');
%         eoutputs=enet(inputs(ind,:));
%         ersq=nancorr(eoutputs(:),targets(:))^2;
%         prsq(i)=(rsq-ersq)/rsq;
    end    
else
    prsq=[];    
end

end