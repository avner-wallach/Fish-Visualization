function [ net,rsq,sig,prsq ] = ni_network(inputs,targets,nhidden,outlayer,est)
%NI_NETWROK generate and train negative image regression. hidden layer is
%random and static, output layer is performing gradient decent
if(nargin<4)
    outlayer='tansig';
end
if(nargin<5)
    est=0;
end
%% generate granular stage
% x=[1;1]*[1:nhidden];
% x=x(:);
y=[]; x=[];
% if(size(inputs,1)>2)
    for i=1:nhidden %each row = one GC
        J=randperm(size(inputs,1),randi(min(size(inputs,1),3)));
        y=[y;J'];
        x=[x;i*ones(size(J'))];
    end
    w=randn(size(x));
    W=sparse(x,y,w);    
% else    
%     W=randn(nhidden,size(inputs,1));
% end
B=randn(nhidden,1);
mf=W*inputs+B; %'synaptic input' into GC
if(nhidden<=2e3)
    gc=tansig(mf);  %GC output activity
else
    gc=zeros(size(mf));
    for i=1:size(mf,2)
        gc(:,i)=tansig(mf(:,i));
    end
end
%% generate network[
net1 = fitnet(1,'trainbr');
% net.layers{1}.transferFcn = outlayer;
net1.inputs{1}.processFcns={};
net1.outputs{2}.processFcns={'mapminmax'};
net1=configure(net1,gc,targets);
net1.biases{2}.learn=0;
net1.b{2}=0;
net1.layerWeights{2,1}.learn=0;
net1.LW{2,1}=1;
%% train network
net1 = train(net1,gc,targets,'useParallel','yes');%,'useGPU','yes');

%% generate output network
net = fitnet(nhidden);
net.inputs{1}.processFcns={};
net.outputs{2}.processFcns={'mapminmax'};
net.layers{2}.transferFcn = outlayer;
net = configure(net,inputs,targets);
net.iw{1,1}=W;
net.b{1}=B;
net.b{2}=net1.b{1};
net.LW{2,1}=net1.iw{1,1};

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