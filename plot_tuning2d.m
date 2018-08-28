function [ F ] = plot_tuning2d( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.headfix='on';
params.bg='on';
params.x_lim=[-150 150];   %default x bin limits for headfixed view
params.y_lim=[-400 200];   %default x bin limits for headfixed view
params.x_nbins=25;  %default x bin number
params.y_nbins=25;  %default y bin number
params.minsamp=1;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
params.objidx=1;    %object 
params.freqmode='on';
params.t_col=0;     %targe column
params.t_val=0;  %target mean / boundaries
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
x=[]; y=[]; z=[]; t=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    % headfix
    if(strcmp(params.headfix,'on'))
        a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(k).fnames));
        a=data_struct(k).data(:,a_col);
        objx=file_struct(k).objects(params.objidx).x;
        objy=file_struct(k).objects(params.objidx).y;
        [x1,y1]=allo2ego(objx,objy,a,x1,y1); %obj coordinates rel. to LED
    end    
    
    x=[x;x1];
    y=[y;y1];
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
if(strcmp(params.headfix,'off'))
    Ix=find(~isoutlier(x));
    Iy=find(~isoutlier(y));
    if(params.xy_cols==[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
             find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))])
         params.x_lim=[1 size(file_struct(1).BG,1)];
         params.y_lim=[1 size(file_struct(1).BG,2)];
    else
        params.x_lim=[nanmin(x(Ix)) nanmax(x(Ix))];
        params.y_lim=[nanmin(y(Iy)) nanmax(y(Iy))];
    end
end
params.x_edges=linspace(params.x_lim(1),params.x_lim(2),params.x_nbins+1);
params.y_edges=linspace(params.y_lim(1),params.y_lim(2),params.y_nbins+1);
params.x_bins=edge2bin(params.x_edges);
params.y_bins=edge2bin(params.y_edges);


% x=data_struct.data(params.ind,params.xy_cols(1));
% y=data_struct.data(params.ind,params.xy_cols(2));
% z=data_struct.data(params.ind,z_col);
%%
[N0,Mz,Sz]=tuning2d();
AL=min(N0,quantile(N0(:),params.maxa));
F=gcf;
A=gca;
if(strcmp(params.headfix,'off') & (strcmp(params.bg,'on') | strcmp(params.bg,'image')))
    sX=size(file_struct(1).BG,2);
    sY=size(file_struct(1).BG,1);
    ix=([1:sX])*pxl2mm;
    iy=([1:sY])*pxl2mm;    
    if(strcmp(params.bg,'on'))
        imagesc(ix,iy,file_struct(1).BG);
    else
        imagesc(ix,iy,params.image);
    end
end
if(strcmp(params.headfix,'on') & isfield(params,'image'))
    sX=size(params.image,2);
    sY=size(params.image,1);
    ix=([1:sX]-sX/2)*pxl2mm;
    iy=(-[1:sY]+sY/3)*pxl2mm;
    imagesc(ix,iy,params.image);
end

hold on;
S=surf(params.x_bins*pxl2mm,params.y_bins*pxl2mm,ones(size(Mz')),...
    'CData',Mz',...
    'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
    ,'AlphaData',AL','AlphaDataMapping','scaled');
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
set(A,'ALim',[0 max(quantile(N0(:),params.maxa),10)]);
view(0,90);
set(gca,'FontSize',params.fontsize)
set(gca,'Xlim',params.x_bins([1 end])*pxl2mm,'Ylim',params.y_bins([1 end])*pxl2mm);
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
        if(numel(idx)>0 & params.t_col) %there is a target variable
            if(numel(params.t_val)==1)  %target mean                                    
                inds=get_subpop_indices(t(idx),params.t_val);
            else %boundaries
                inds=find(t(idx)>=params.t_val(1) & t(idx)<params.t_val(2));
            end
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

