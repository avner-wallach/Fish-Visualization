function [ F,Mz ] = plot_tuning_dist( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.headfix='on';
params.bg='on';
% params.x_lim=[-150 150];   %default x bin limits for headfixed view
% params.y_lim=[-400 200];   %default x bin limits for headfixed view
params.r_lim=[1 150];   %distance from RF center
% params.x_nbins=25;  %default x bin number
% params.y_nbins=25;  %default y bin number
params.r_nbins=5;
params.rf=[-15 0]; %RF center coordinates
params.minsamp=1;
params.maxa=0.3; %maximal alpha percentile
params.fontsize=16;
params.objidx=1;    %object 
params.freqmode='on';
params.t_col=0;     %targe column
params.t_edges=0;  %target mean / boundaries
% params.clim='auto';
%% get varins
n=numel(varargin);
if(n<3 | n>3&~mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};
% params.ind=1:size(data_struct(1).data,1);
params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
params.K=numel(data_struct);

%% get data vectors
x=[]; y=[]; z=[]; t=[]; r=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1))*pxl2mm;
    y1=data_struct(k).data(:,params.xy_cols(2))*pxl2mm;
    % headfix
    if(strcmp(params.headfix,'on'))
        a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(k).fnames));
        a=data_struct(k).data(:,a_col);
        objx=file_struct(k).objects(params.objidx).x*pxl2mm;
        objy=file_struct(k).objects(params.objidx).y*pxl2mm;
        [x1,y1]=allo2ego(objx,objy,a,x1,y1); %obj coordinates rel. to LED
        r1=((x1-params.rf(1)).^2 + (y1-params.rf(2)).^2).^(0.5);
    end    
    
    x=[x;x1];
    y=[y;y1];
    r=[r;r1];
    z=[z;data_struct(k).data(:,z_col)];
    if(params.t_col~=0)
        t=[t;data_struct(k).data(:,params.t_col)];
    end    
end

if(isfield(params,'func') & params.t_col~=0)
    if(numel(params.t_col)>1)
        t=params.func(t,2);
    else
        t=params.func(t);
    end
end

if(params.t_col~=0 & numel(params.t_col)==1 & ...
    strcmp(data_struct(1).fnames{params.t_col},'iei') & strcmp(params.freqmode,'on'))
    dt=nanmean(diff(sort(t)));
    t=1./(t+dt*randn(size(t)));
end
if(strcmp(data_struct(1).fnames{z_col},'iei') & strcmp(params.freqmode,'on'))
    z=1./z;
end
    
%% produce bins
params.r_edges=linspace(params.r_lim(1),params.r_lim(2),params.r_nbins+1);
params.r_bins=edge2bin(params.r_edges);
params.t_nbins=numel(params.t_edges)-1;
params.t_bins=edge2bin(params.t_edges);
%% tuning    
[N0,Mz,Sz]=tuning2d();
AL=min(N0,quantile(N0(:),params.maxa));
F=gcf;
A=gca;

hold on;
S=surf(params.r_bins,params.t_bins,ones(size(Mz')),...
    'CData',Mz',...
    'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
    ,'AlphaData',AL','AlphaDataMapping','scaled');
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
% set(A,'ALim',[0 max(quantile(N0(:),params.maxa),10)]);
view(0,90);
set(gca,'FontSize',params.fontsize)
set(gca,'Xlim',params.r_bins([1 end]),'Ylim',params.t_bins([1 end]));
hold off;
colorbar;

function [N0,Mz,Sz]=tuning2d()
N0=zeros(params.r_nbins,params.t_nbins);
Mz=nan(params.r_nbins,params.t_nbins);
Sz=nan(params.r_nbins,params.t_nbins);
for i=1:params.r_nbins
    for j=1:params.t_nbins
        idx=find(r>params.r_edges(i) & r<=params.r_edges(i+1) &...
                 t>=params.t_edges(j) & t<params.t_edges(j+1) & ...
                    y>=0 );   %only ipsi-front quadrant
        N0(i,j)=numel(idx);
        if(numel(idx)>=params.minsamp)
            Mz(i,j)=nanmean(z(idx));
            Sz(i,j)=nanstd(z(idx));
        end
    end
end
end

end

