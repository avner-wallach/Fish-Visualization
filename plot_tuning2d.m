function [ F ] = plot_tuning2d( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% params
params.headfix='on';
params.bg='on';
params.x_lim=[-100 100];   %default x bin limits for headfixed view
params.y_lim=[-400 200];   %default x bin limits for headfixed view
params.x_nbins=24;  %default x bin number
params.y_nbins=24;  %default y bin number
params.minsamp=1;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
params.objidx=1;    %object 
params.freqmode='on';
params.t_col=0;     %targe column
params.t_mean=0;  %target mean 
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
params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct.fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct.fnames))];%default is spatial mapping
if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
%% produce bins
if(strcmp(params.headfix,'off'))
    x=data_struct.data(:,params.xy_cols(1));
    y=data_struct.data(:,params.xy_cols(2));
    Ix=find(~isoutlier(x));
    Iy=find(~isoutlier(y));
    params.x_lim=[nanmin(x(Ix)) nanmax(x(Ix))];
    params.y_lim=[nanmin(y(Iy)) nanmax(y(Iy))];
end
params.x_edges=linspace(params.x_lim(1),params.x_lim(2),params.x_nbins+1);
params.y_edges=linspace(params.y_lim(1),params.y_lim(2),params.y_nbins+1);
params.x_bins=edge2bin(params.x_edges);
params.y_bins=edge2bin(params.y_edges);

%% get data vectors
x=data_struct.data(params.ind,params.xy_cols(1));
y=data_struct.data(params.ind,params.xy_cols(2));
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
%% headfix
if(strcmp(params.headfix,'on'))
    a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct.fnames));
    a=data_struct.data(params.ind,a_col);
    objx=file_struct.objects(params.objidx).x;
    objy=file_struct.objects(params.objidx).y;
    [x,y]=allo2ego(objx,objy,a,x,y); %obj coordinates rel. to LED
end    
%%
[N0,Mz,Sz]=tuning2d();
AL=min(N0,quantile(N0(:),params.maxa));
F=gcf;
A=gca;
if(strcmp(params.headfix,'off') & strcmp(params.bg,'on'))
    imshow(file_struct.BG);
end
hold on;
S=surf(params.x_bins,params.y_bins,ones(size(Mz')),...
    'CData',Mz',...
    'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
    ,'AlphaData',AL','AlphaDataMapping','scaled');
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
set(A,'ALim',[0 quantile(N0(:),params.maxa)]);
view(0,90);
set(gca,'FontSize',params.fontsize)
hold off;
colorbar;

function [N0,Mz,Sz]=tuning2d()
Mx=numel(params.x_edges)-1;
My=numel(params.y_edges)-1;
N0=zeros(params.x_nbins,params.y_nbins);
Mz=nan(params.x_nbins,params.y_nbins);
Sz=nan(params.x_nbins,params.y_nbins);
for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                 y>params.y_edges(j) & y<=params.y_edges(j+1));
%         if()
        if(numel(idx)>0 & params.t_col) %there is a target variable
            inds=get_subpop_indices(t(idx),params.t_mean);
            idx=idx(inds);
        end
        N0(i,j)=numel(idx);
        if(numel(idx)>=params.minsamp)
            Mz(i,j)=nanmean(z(idx));
            Sz(i,j)=nanstd(z(idx));
        end
    end
end
end

end

