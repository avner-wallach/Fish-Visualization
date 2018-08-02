% function [strct]=tuning1d(x_vec,z_vec,x_edges,minsamp,q)
function [ F ] = plot_tuning1d( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% params
params.headfix='on';
params.x_lim=[0 55];   %default x bin limits for headfixed view
params.x_nbins=24;  %default x bin number
params.y_nbins=24;  %default y bin number
params.minsamp=1;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
params.objidx=1;    %object 
params.freqmode='on';
params.t_col=0;     %targe column
params.t_mean=0;  %target mean 
params.quant=[0.25 0.75];
% params.clim='auto';

%% get varins
n=numel(varargin);
if(n<3 | n>3&~mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};
params.ind=1:size(data_struct.data,1);
params.x_col=find(cellfun(@(x) (strcmp(x,'iei')),data_struct.fnames));
if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end

%%  produce bins
if(params.x_col~=find(cellfun(@(x) (strcmp(x,'iei')),data_struct.fnames)))
    x=data_struct.data(:,params.x_col(1));    
    Ix=find(~isoutlier(x));
    params.x_lim=[nanmin(x(Ix)) nanmax(x(Ix))];
end
params.x_edges=linspace(params.x_lim(1),params.x_lim(2),params.x_nbins+1);
params.x_bins=edge2bin(params.x_edges);

%% get data vectors
x=data_struct.data(params.ind,params.x_col(1));
z=data_struct.data(params.ind,z_col);
if(params.t_col~=0)
    t=data_struct.data(params.ind,params.t_col);
    if(strcmp(data_struct.fnames{params.t_col},'iei') & strcmp(params.freqmode,'on'))
        dt=nanmean(diff(sort(t)));
        t=1./(t+dt*randn(size(t)));
    end
end
if(strcmp(data_struct.fnames{z_col},'iei') & strcmp(params.freqmode,'on'))
    z=1./z;
end
if(strcmp(data_struct.fnames{params.x_col},'iei') & strcmp(params.freqmode,'on'))
    x=1./x;
end
% %% headfix
% if(strcmp(params.headfix,'on'))
%     a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct.fnames));
%     a=data_struct.data(params.ind,a_col);
%     objx=file_struct.objects(params.objidx).x;
%     objy=file_struct.objects(params.objidx).y;
%     [x,y]=allo2ego(objx,objy,a,x,y); %obj coordinates rel. to LED
% end    
%%
[N0,Mz,Sz,Qz,I]=tuning1d();
AL=min(N0,quantile(N0(:),params.maxa));
F=gcf;
A=gca;
hold on;
H=plot(x(I),z(I),'.');
set(H,'MarkerSize',6,'Color',0.75*[1 1 1]);
H=plot(params.x_bins,Mz);
% H=area(params.x_bins,Qz);
H=plot(params.x_bins,Qz(:,1),':k');
H=plot(params.x_bins,Qz(:,2),':k');
% set(H(2),'FaceAlpha',0);
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
% set(A,'ALim',[0 quantile(N0(:),params.maxa)]);
% view(0,90);
set(gca,'FontSize',params.fontsize)
hold off;

[c,p]=get_corr_pval(x(I),z(I))

function [N0,Mz,Sz,Qz,I]=tuning1d()
Mx=numel(params.x_edges)-1;
N0=zeros(params.x_nbins,1);
Mz=nan(params.x_nbins,1);
Sz=nan(params.x_nbins,1);
Qz=nan(params.x_nbins,numel(params.quant));
I=[];
for i=1:params.x_nbins
    idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1));
    if(numel(idx)>0 & params.t_col>0) %there is a target variable
        inds=get_subpop_indices(t(idx),params.t_mean);
        idx=idx(inds);
    end
    N0(i)=numel(idx);
    if(numel(idx)>=params.minsamp)
        Mz(i)=nanmean(z(idx));
        Sz(i)=nanstd(z(idx));
        Qz(i,:)=quantile(z(idx),params.quant);
    end
    I=[I;idx];
end
end

%% get correlation and p-value
    function [c,p]=get_corr_pval(xx,zz)
        ind=find(~isnan(xx) & ~isnan(zz) & ~isinf(xx) & ~isinf(zz));
        xx=xx(ind);
        zz=zz(ind);
        c=corr(xx,zz);
        M=1000;
        for m=1:M
            c_bs(m)=corr(xx,zz(randperm(numel(zz))));
        end
        p=sum(abs(c_bs)>=abs(c))/M;
    end
end


