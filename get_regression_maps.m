function [ F ] = get_regression_maps( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.headfix='on';
params.bg='on';
params.x_lim=[-250 250];   %default x bin limits for headfixed view
params.y_lim=[-400 200];   %default x bin limits for headfixed view
params.x_nbins=15;  %default x bin number
params.y_nbins=15;  %default y bin number
params.minsamp=1;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
params.objidx=1;    %object 
params.freqmode='on';
t_cols=0;     %targe column
params.t_val=0;  %target mean / boundaries
params.pth=0.05; %pvalue threshold for significance
% params.clim='auto';
%% get varins
n=numel(varargin);
if(n<4 | n>4&mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};
t_cols=varargin{4};
% params.ind=1:size(data_struct(1).data,1);
params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
if(n>4)
    i=5;
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
    if(strcmp(params.headfix,'home'))
        a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(k).fnames));
        a=data_struct(k).data(:,a_col);
    end    
    
    x=[x;x1];
    y=[y;y1];
    z=[z;data_struct(k).data(:,z_col)];
    t=[t;data_struct(k).data(:,t_cols)];
        
end

if(isfield(params,'func'))
    if(numel(t_cols)>1)
        t=params.func(t,2);
    else
        t=params.func(t);
    end
end

iei_ind=find(cellfun(@(x) strcmp(x,'iei'),data_struct(1).fnames(t_cols)));
if(numel(iei_ind) & strcmp(params.freqmode,'on'))
    dt=nanmean(diff(sort(t(:,iei_ind))));
    t(:,iei_ind)=1./(t(:,iei_ind)+dt*randn(size(t,iei_ind),1));
end
if(strcmp(data_struct(1).fnames{z_col},'iei') & strcmp(params.freqmode,'on'))
    z=1./z;
end
%% home analysis
if(strcmp(params.headfix,'home'))
    [x0,y0,a0]=get_home_angle(file_struct.home);
    
    %convert to home-cenetered coordinates
    x=x-x0;
    y=y0-y;
    %rotate
    x_=(cos(a0).*x-sin(a0).*y);
    y_=(sin(a0).*x+cos(a0).*y);
    a_=a0-a;
    ind=find(sin(a_)<0); %fish facing 'south'
    x_(ind)=-x_(ind);    y_(ind)=-y_(ind);  %mirror fish
    x=x_;
    y=y_;
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

%% get regression maps
[pval,c,M,B,P,T,N]=tuning2d();
%% plot base response map
figure;
AL=min(N,quantile(N(:),params.maxa));
plot_map(B,AL);

%% plot correlation maps
for m=1:numel(t_cols)
    figure;
    AL=-log10(100*pval(:,:,m));
    plot_map(c(:,:,m),AL);
    colormap(modified_jet);
    set(gca,'CLim',[-1 1]);
end
%% get regression
ind=find(x>params.x_edges(1) & x<=params.x_edges(end) &...
         y>params.y_edges(1) & y<=params.y_edges(end));      %all smaples in maps
b_est=interp2(params.x_bins,params.y_bins,B,x(ind),y(ind));
for m=1:numel(t_cols)
    m_est(:,m)=interp2(params.x_bins,params.y_bins,M(:,:,m),x(ind),y(ind));
end
z_est=b_est+sum(m_est.*(t(ind,:)-repmat(T,numel(ind),1)),2);
figure;
plot(z(ind),z_est,'.');


function plot_map(Mz,A0)
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
if(strcmp(params.headfix,'home') & strcmp(params.bg,'on') )
    frame=file_struct(1).BG;
    sX=size(frame,2);
    sY=size(frame,1);
    frame=imtranslate(frame,[sX/2-x0,sY/2-y0],'FillValues',128); %center image around home
    frame=imrotate(frame,rad2deg(a0),'bilinear','crop');
    frame=imcrop(frame,[sX/2+params.x_lim(1),sY/2+params.y_lim(1),diff(params.x_lim),diff(params.y_lim)]);        
    
    ix=([params.x_lim(1):params.x_lim(2)])*pxl2mm;
    iy=([params.y_lim(1):params.y_lim(2)])*pxl2mm;    
    imagesc(ix,iy,frame);
end

hold on;
S=surf(params.x_bins*pxl2mm,params.y_bins*pxl2mm,ones(size(Mz')),...
    'CData',Mz',...
    'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
    ,'AlphaData',A0','AlphaDataMapping','scaled');
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
set(A,'ALim',[0 mean(A0(:))]);
view(0,90);
set(gca,'FontSize',params.fontsize)
set(gca,'Xlim',params.x_bins([1 end])*pxl2mm,'Ylim',params.y_bins([1 end])*pxl2mm);
set(gca,'YDir','normal');
axis('image');
hold off;
colorbar;
end

function [pval,c,M,B,P,T,N]=tuning2d()
B=nan(params.x_nbins,params.y_nbins); %response mean map 
N=zeros(params.x_nbins,params.y_nbins); %sample map 
c=zeros(params.x_nbins,params.y_nbins,size(t,2));
pval=ones(params.x_nbins,params.y_nbins,size(t,2));
M=zeros(params.x_nbins,params.y_nbins,size(t,2));
P=nan(params.x_nbins,params.y_nbins,size(t,2));
T=nanmean(t); %overall mean posture

for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                 y>params.y_edges(j) & y<=params.y_edges(j+1));             
        if(numel(idx)>0) 
            Ind=idx(~isnan(sum(t(idx,:),2)) & ~isnan(z(idx))); %remove nans            
            N(i,j)=numel(Ind);
            B(i,j)=mean(z(Ind));
            for m=1:size(t,2) %go over all targets
                P(i,j,m)=mean(t(Ind,m)); %posture map
                [c(i,j,m),pval(i,j,m)]=corr(t(Ind,m),z(Ind));
                [r,M(i,j,m),b]=regression(t(Ind,m)-P(i,j,m),z(Ind),'one');
            end
        end
    end
end
% remove regression where correlations are insignificant
M(pval>params.pth)=0;

%correct base response map according to posture maps
TT=repmat(reshape(T,1,1,[]),params.x_nbins,params.y_nbins,1);
B=B+sum(M.*(TT-P),3);

end

    function [x0,y0,a0]=get_home_angle(home)
        x0=mean(home(:,1));
        y0=mean(home(:,2));
        D=diff([home;home(1,:)]).*[1;1;-1;-1];
        V=hypot(D(:,1),D(:,2));
        a1=atan2(D(:,2),D(:,1));
        a0=mean(a1(V>median(V)))+pi/2; %direction of main axes
        params.x_lim=[-2*min(V) 2*min(V)];
        params.y_lim=[-max(V) max(V)];
    end    



end

