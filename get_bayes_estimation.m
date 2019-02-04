function [ pxy,pz,pz_xy,bins,pm,ps ] = get_bayes_estimation( varargin )
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
params.z_nbins=15;  %default y bin number
params.minsamp=1;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
params.objidx=1;    %object 
params.freqmode='on';
params.t_col=0;     %targe columns
params.visualize=0; %plot poterior maps
% params.pth=0.05; %pvalue threshold for significance
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
    if(strcmp(params.headfix,'home'))
        a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(k).fnames));
        a=data_struct(k).data(:,a_col);
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
if(params.t_col~=0)
    iei_ind=find(cellfun(@(x) strcmp(x,'iei'),data_struct(1).fnames(params.t_col)));
    if(numel(iei_ind) & strcmp(params.freqmode,'on'))
        dt=nanmean(diff(sort(t(:,iei_ind))));
        t(:,iei_ind)=1./(t(:,iei_ind)+dt*randn(size(t,iei_ind),1));
    end
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
Iz=find(~isoutlier(z));
params.z_lim=[nanmin(z(Iz)) nanmax(z(Iz))];
params.x_edges=linspace(params.x_lim(1),params.x_lim(2),params.x_nbins+1);
params.y_edges=linspace(params.y_lim(1),params.y_lim(2),params.y_nbins+1);
params.z_edges=linspace(params.z_lim(1),params.z_lim(2),params.z_nbins+1);
params.x_bins=edge2bin(params.x_edges);
params.y_bins=edge2bin(params.y_edges);
params.z_bins=edge2bin(params.z_edges);
%% get regression maps
[beta,rsq]=get_regressions();

%% get prior and likelyhood - naive
[pxy,pz_xy,pz,pm,ps]=get_dsitributions();
%compute posterior (normalized)
% pz_xy=pz_xy./repmat(permute(pz,[3,1,2]),params.x_nbins,params.y_nbins,1);

%% get posterior and likelyhood - corrected
% [pz_xy_c,pm_c,ps_c]=get_corrected_dsitributions();


if(params.visualize==1)
    pxy_z=pz_xy.*repmat(pxy,1,1,params.z_nbins);
    figure;
    for i=1:params.z_nbins
        plot_map(pxy_z(:,:,i));
    %     set(gca,'Alim',[0 max(pxy_z(:))]);
        fig(i)=getframe;
    end
%     make_video_from_imgs(fig,'posteriors','',3);
end
%% up-smaple maps

    
bins={params.x_bins,params.y_bins,params.z_bins};

function plot_map(Mz)
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
A0=Mz;
hold on;
S=surf(params.x_bins*pxl2mm,params.y_bins*pxl2mm,ones(size(Mz')),...
    'CData',Mz',...
    'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
    ,'AlphaData',A0','AlphaDataMapping','scaled');
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
set(A,'ALim',[0 quantile(A0(:),0.99)]);
view(0,90);
set(gca,'FontSize',params.fontsize)
set(gca,'Xlim',params.x_bins([1 end])*pxl2mm,'Ylim',params.y_bins([1 end])*pxl2mm);
set(gca,'YDir','normal');
axis('image');
hold off;
% colorbar;
end

function [pxy,pz_xy,pz,pm,ps]=get_dsitributions();
pxy=zeros(params.x_nbins,params.y_nbins); %prior distribution
pz_xy=zeros(params.x_nbins,params.y_nbins,params.z_nbins); %likelyhood func.
% Mz=NaN(params.x_nbins,params.y_nbins,params.z_nbins); %likelyhood func.

for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                 y>params.y_edges(j) & y<=params.y_edges(j+1));             
        if(numel(idx)>0) 
            Ind=idx(~isnan(z(idx))); %remove nans            
            N(i,j)=numel(Ind);
            pz_xy(i,j,:)=histcounts(z(Ind),params.z_edges,'Normalization','probability');
            pm(i,j)=mean(z(Ind));
            ps(i,j)=std(z(Ind));
        end
    end
end
pxy=N/sum(N(:));

%normalization (measurement distribution)
idx=find(x>params.x_edges(1) & x<=params.x_edges(end) &...
         y>params.y_edges(1) & y<=params.y_edges(end));      %all smaples in maps
pz=histcounts(z(idx),params.z_edges,'Normalization','probability');
end

function [pz_xy_c,pm_c,ps_c]=get_corrected_dsitributions()

pz_xy_c=zeros(params.x_nbins,params.y_nbins,params.z_nbins); %likelyhood func.

X=[t ones(size(z))]; %add units

%get prediction coefficients
for i=1:size(beta,3)
    b(i,:)=interp2(params.x_bins,params.y_bins,beta(:,:,i),x,y);
end
%get deviations from predicted response
z_=dot(X',b)';

for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                 y>params.y_edges(j) & y<=params.y_edges(j+1));             
        if(numel(idx)>0) 
            Ind=idx(~isnan(z_(idx))); %remove nans            
            N(i,j)=numel(Ind);
            pz_xy_c(i,j,:)=histcounts(z_(Ind),params.z_edges,'Normalization','probability');
            pm_c(i,j)=mean(z_(Ind));
            ps_c(i,j)=std(z_(Ind));
        end
    end
end

end


function [beta,rsq]=get_regressions()
    X=[t ones(size(z))]; %add units
    
    for i=1:params.x_nbins
        for j=1:params.y_nbins
            idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                     y>params.y_edges(j) & y<=params.y_edges(j+1));             
            if(numel(idx)>0) 
                Ind=idx(~isnan(z(idx))); %remove nans            
                N(i,j)=numel(Ind);
                [b,bint,r,rint,stats]=regress(z(Ind),X(Ind,:));
                beta(i,j,:)=b;
                rsq(i,j)=stats(1);
            end
        end
    end      
end
end

