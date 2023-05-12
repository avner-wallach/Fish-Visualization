function [ Izt,Izt_xy,Izxy,Izxy_t,izt_xy ] = compute_info( varargin )
%COMPUTE_INFO compute mutual info I(z,t|x,y)
%   
pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.headfix='on';
params.bg='on';
params.bgcol='w';
params.x_lim=[-550 550];   %default x bin limits for headfixed view
params.y_lim=[-550 450];   %default x bin limits for headfixed view
params.x_nbins=15;  %default x bin number
params.y_nbins=15;  %default y bin number
params.info_nbins=5; %max number of bins for info calculations
params.minsamp=35;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
% params.objidx=1;    %object 
params.freqmode='on';
params.t_col=0;     %targe column
params.t_val=0;  %target mean / boundaries
params.ind_lim=[];
params.mfunc='mean'; %function to plot: mean-mean response in each bin; slope-sensitivity to t; offset- intersect at t=0
params.maxz=inf;
params.plotting=0;
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

a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));

params.K=numel(data_struct);

%% get data vectors
x=[]; y=[]; a=[]; z=[]; t=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    a1=data_struct(k).data(:,a_col);
    
    x=[x;x1];
    y=[y;y1];
    a=[a;a1];    
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

if(isfield(params,'zfunc'))
    z=params.zfunc(z,2);
elseif(strcmp(data_struct(1).fnames{z_col},'iei') & strcmp(params.freqmode,'on'))
    z=1./z;
end

if(params.t_col~=0 & numel(params.t_col)==1 & ...
    strcmp(data_struct(1).fnames{params.t_col},'iei') & strcmp(params.freqmode,'on'))
    dt=nanmean(diff(sort(t)));
    t=1./(t+dt*randn(size(t)));
end

z(abs(z)>params.maxz)=nan;
%% default object idx
if(~isfield(params,'objidx'))
    params.objidx=numel(file_struct(1).objects);
end

%% reduce indices
if(numel(params.ind_lim))
    if(numel(params.ind_lim)==2)
        inds=inrange(1:numel(z),params.ind_lim);
    else
        inds=params.ind_lim;
    end
    z=z(inds);
    if(numel(t))
        t=t(inds,:);
    end
    x=x(inds);
    y=y(inds);
    a=a(inds);
    
%     if(params.roinum>0 & strfind(data_struct(1).fnames{z_col},'sc'))
%         raster=raster(inrange(raster(:,2),params.ind_lim),:);
%         raster(:,2)=raster(:,2)-params.ind_lim(1);
%     end
end
%% coodintate transforms
if(strcmp(params.headfix,'on'))
    if(params.objidx==0)
        convert_to_wall;
    else
        convert_to_object;
    end
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

%% compute information
%remove nans;
idx=find(~isnan(x+y+z+t));
x=x(idx);
y=y(idx);
z=z(idx);
t=t(idx);

% discretize
% if(numel(unique(z))>params.info_nbins) %not discrete
%     [z,ez]=discretize(z,params.info_nbins);
%     z_bins=edge2bin(ez);
% else    
%     z_bins=unique(z);
%     ez=bin2edge(z_bins);
% end
% 
% if(numel(unique(t))>params.info_nbins) %not discrete
%     [t,et]=discretize(t,params.info_nbins);
%     t_bins=edge2bin(et);
% else
%     t_bins=unique(t);
%     et=bin2edge(t_bins);
% end

[x,ex]=discretize(x,params.x_edges);
[y,ey]=discretize(y,params.y_edges);
[tmp1,tmp2,xy]=unique(1000*x+y); %convert to single discrete vector
xy(isnan(x+y))=nan;

%reafference info
Izt=mutInfo(z,t)/entropy(z);
Idx=[];
for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x==i & y==j);
%          idx=find(x>=params.x_edges(i) & x<params.x_edges(i+1) &...
%              y>=params.y_edges(j) & y<params.y_edges(j+1));        
        N0(i,j)=numel(idx);
        if(numel(idx)>=params.minsamp)
            izt_xy(i,j)=mutInfo(z(idx),t(idx),params.info_nbins);
            c(i,j)=nancorr(z(idx),t(idx));
            ez_xy(i,j)=entropy(z(idx),params.info_nbins);
            Idx=[Idx;idx];
        else            
            izt_xy(i,j)=nan;
            c(i,j)=nan;
            ez_xy(i,j)=nan;            
        end
    end
end
izt_xy=izt_xy./ez_xy;
% izt_xy=izt_xy/entropy(z);
Izt_xy=nansum(izt_xy(:).*N0(:))/sum(N0(:));
% Izt_xy=nanmean(izt_xy(:));
        
%exafference info
idx=find(~isnan(xy));
z=z(idx);
t=t(idx);
xy=xy(idx);

Izxy=mutInfo(z,xy)/entropy(z,params.info_nbins);
% discretize
% if(numel(unique(z))>params.info_nbins) %not discrete
%     [z,ez]=discretize(z,params.info_nbins);
%     z_bins=edge2bin(ez);
% else    
%     z_bins=unique(z);
%     ez=bin2edge(z_bins);
% end
% 
if(numel(unique(t))>params.info_nbins) %not discrete
    [t,et]=discretize(t,params.info_nbins);
    t_bins=edge2bin(et);
else
    t_bins=unique(t);
    et=bin2edge(t_bins);
end

for i=1:numel(t_bins)
    idx=find(t==i);
    M0(i)=numel(idx);
    if(numel(idx)>=params.minsamp)
        izxy_t(i)=mutInfo(z(idx),xy(idx))/entropy(z(idx),params.info_nbins);
    else
        izxy_t(i)=nan;
    end
end

Izxy_t=nansum(izxy_t(:).*M0(:))/sum(M0(:));

%% plotting
if(params.plotting)
    F=figure;
    A=gca;
    if(isfield(params,'image'))
        sX=size(params.image,2);
        sY=size(params.image,1);
        ix=([1:sX]-sX/2)*pxl2mm;
        iy=(-[1:sY]+sY/3)*pxl2mm;    
        if(params.bgcol=='w')
            imagesc(ix,iy,params.image);
        else
            imagesc(ix,iy,255-params.image);
        end
    end
    hold on;
    S=surf(params.x_bins*pxl2mm,params.y_bins*pxl2mm,zeros(size(izt_xy')),...
        'CData',izt_xy',...
        'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
        ,'AlphaData',N0','AlphaDataMapping','scaled');
    if(isfield(params,'clim'))
        set(A,'CLim',params.clim);
    end
    a=max(quantile(N0(:),params.maxa),10);
    set(A,'ALim',[0 a]);
    view(0,90);
    set(gca,'FontSize',params.fontsize)
    set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',params.x_bins([1 end])*pxl2mm,'Ylim',params.y_bins([1 end])*pxl2mm);
    hold off;
    colorbar;
    if(params.bgcol=='k')
        set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);
        set(gcf,'Color',[0 0 0]);
        colorbar(gca,'Color',[1 1 1]);
    end
end

    function convert_to_wall()
%         file_struct(1).circle = file_struct(1).circle*pxl2mm;
        rx=file_struct(1).circle(3)/2; ry=file_struct(1).circle(4)/2;    %ellipse radii
        x0=file_struct(1).circle(1) + rx; y0=file_struct(1).circle(2) + ry; %center point
        phi=atan2((y-y0),(x-x0));    %azimuth in tank
        R1=hypot((x-x0),(y-y0)); %distance of fish from center
        R2=(rx*ry)./sqrt((ry*cos(phi)).^2 + (rx*sin(phi)).^2); %distance of nearest point from cetner
        ra=R2-R1;    %distance of fish to nearest wall
        ra(ra<0)=nan;
        th=phi-a;   %egocentric angle of closest wall
        x=ra.*sin(th);
        y=ra.*cos(th);
%         th(th>pi)=th(th>pi)-2*pi;
%         th(th<-pi)=th(th<-pi)+2*pi;    
    end

    function convert_to_object()
          objx=file_struct(k).objects(params.objidx).x;%*pxl2mm;
          objy=file_struct(k).objects(params.objidx).y;%*pxl2mm;
          [x,y]=allo2ego(objx,objy,a,x,y); %obj coordinates rel. to LED
%           th=atan2(x,y);    %azimuth of object
%           r=hypot(x,y); %distance of object from fish's head      
    end

end