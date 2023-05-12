function [ Mz,N0,out_arg ] = plot_tuning2d( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.headfix='on';
params.bg='on';
params.bgcol=[1 1 1];
params.x_lim=[-550 550];   %default x bin limits for headfixed view
params.y_lim=[-550 450];   %default x bin limits for headfixed view
params.r_max=0;
params.x_nbins=25;  %default x bin number
params.y_nbins=25;  %default y bin number
params.minsamp=10;
params.zscore=1;
params.maxa=0.4; %maximal alpha percentile
params.fontsize=16;
params.freqmode='on';
params.t_col=0;     %targe column
params.ind_lim=[];
params.mfunc='mean'; %function to plot: mean-mean response in each bin; slope-sensitivity to t; offset- intersect at t=0
params.maxz=inf;
params.maxr=inf;
params.roinum=0;
params.psth_lines=5;
params.fit_mode='slope';
params.binsize=1;
params.upsamp=1;
params.roc=0;
params.plotting=1;

out_arg=[];
%% get varins
n=numel(varargin);
if(n<3 | n>3&~mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};
% if pca was not performed yet
if(numel(data_struct(1).data)==0)
    for k=1:numel(data_struct)
        data_struct(k).data=data_struct(k).data0;
        data_struct(k).fnames=data_struct(k).fnames0;
    end
end
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

if(isfield(params,'rastind'))
    raster=data_struct(1).raster{params.rastind};
end

if(params.t_col~=0 & numel(params.t_col)==1 & ...
    strcmp(data_struct(1).fnames{params.t_col},'iei') & strcmp(params.freqmode,'on'))
    dt=nanmean(diff(sort(t)));
    t=1./(t+dt*randn(size(t)));
end

z(abs(z)>params.maxz)=nan;
%% perform spike count
if(isfield(params,'rastind') & isfield(params,'scwindow'))
        indt=find(inrange(raster(:,1),params.scwindow));                        
        z=histcounts(raster(indt,2),bin2edge([1:numel(z)]));
end
%% zscore/dezscore data 
if(isfield(params,'zscale'))
    z=z*params.zscale(2) + params.zscale(1);
else
    if(params.zscore)
        z=(z-nanmean(z))/nanstd(z);
        t=(t-nanmean(t))/nanstd(t);
    end
end
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
    if(exist('a'))
        a=a(inds);
    end
    if(params.roinum>0 & strfind(data_struct(1).fnames{z_col},'sc'))
        raster=raster(inrange(raster(:,2),params.ind_lim),:);
        raster(:,2)=raster(:,2)-params.ind_lim(1);
    end
end

%% remove outliers
if(params.r_max>0)
    rx=file_struct(1).circle(3)/2; ry=file_struct(1).circle(4)/2;    %ellipse radii
    x0=file_struct(1).circle(1) + rx; y0=file_struct(1).circle(2) + ry; %center point
    R1=hypot((x-x0),(y-y0)); %distance of fish from center
    N0=numel(z);
    inds=find(R1<params.r_max);
    z=z(inds);
    if(numel(t))
        t=t(inds,:);
    end
    x=x(inds);
    y=y(inds);
    if(exist('a'))
        a=a(inds);
    end
    if(params.roinum>0 & strfind(data_struct(1).fnames{z_col},'sc'))
        ind=find(ismember(raster(:,2),inds));
        N=numel(inds);
        q=zeros(size(N0,1));
        q(r)=[1:N];
        drast=raster(ind,:); 
        drast(:,2)=q(drast(:,2));
        raser=drast;
    end
end    
    
%% 
if(strcmp(params.headfix,'on'))
    if(params.objidx==0)
        convert_to_wall;
    else
        convert_to_object;
    end
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

%%
if(strcmp(params.mfunc,'mean') | strcmp(params.mfunc,'std'))
    [N0,Mz,Sz]=tuning2d();
    if(strcmp(params.mfunc,'std'))
        Mz=Sz;
    end
    AL=N0;
    AL=min(AL,10);
%     clim=quantile(Mz(:),[0.1 0.9]);
else
    [N0,P1,P2,Pv1]=fitting2d();
    if(strcmp(params.mfunc,'slope'))
        Mz=P1;
        AL=Pv1;
        AL(AL>2)=2;
    else
        Mz=P2;
        AL=N0;
        AL=min(AL,10);
    end
end
clim=[-1 1];

if(params.upsamp>1)
    [Y,X]=meshgrid(params.y_bins,params.x_bins);
    F=scatteredInterpolant(X(~isnan(Mz)),Y(~isnan(Mz)),Mz(~isnan(Mz)),'linear','none');
%     F=scatteredInterpolant(X(:),Y(:),Mz(:));
    params.x_bins=linspace(params.x_bins(1),params.x_bins(end),params.x_nbins*params.upsamp);
    params.y_bins=linspace(params.y_bins(1),params.y_bins(end),params.y_nbins*params.upsamp);
    [Yq,Xq]=meshgrid(params.y_bins,params.x_bins);
    Mz=F(Xq,Yq);    
    N0(isnan(N0))=0;
    AL(isnan(AL))=0;
    N0=interp2(Y,X,N0,Yq,Xq);
    AL=interp2(Y,X,AL,Yq,Xq);
end

if(params.plotting)
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
        if(params.bgcol(1)==0)
            imagesc(ix,iy,255-params.image);  
        else
            imagesc(ix,iy,params.image);  
        end
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
    S=surf(params.x_bins*pxl2mm,params.y_bins*pxl2mm,zeros(size(Mz')),...
        'CData',Mz',...
        'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
        ,'AlphaData',AL','AlphaDataMapping','scaled');
    % a=max(quantile(N0(:),params.maxa),10);
    set(A,'ALim',[0 15]);
    view(0,90);
    set(gca,'FontSize',params.fontsize)
    set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',params.x_bins([1 end])*pxl2mm,'Ylim',params.y_bins([1 end])*pxl2mm);
    hold off;
    cb=colorbar;
    set(gca,'Color','none','XColor','none','YColor','none');
    set(gcf,'Color',params.bgcol);
    set(cb,'Color',1-params.bgcol,'Location','east','AxisLocation','out');
%     cb.Position = cb.Position + [0.015 0 -.01 0];
%     cb.Ticks = [-2 -1 0 1 2];
    if(strcmp(params.mfunc,'slope'))
        if(params.bgcol==[1 1 1])
%             colormap(gca,modified_jet);
            colormap(gca,flipud(brewermap(64,'RdBU')));
        else
            colormap(gca,invert_map(flipud(brewermap(64,'RdBU'))));
        end
        set(gca,'alim',[1 2.1]);    
    else
        colormap(gca,'parula');
    end
    if(isfield(params,'clim'))
        set(A,'CLim',params.clim);
    else
        set(A,'Clim',clim);
    end
    
    if(isfield(params,'poly'))
        params.roinum=numel(params.poly);
    end
    
    if(params.roinum>0)
        if(~isfield(params,'cmaps'))
            C=brewermap(10,'Set1');
        else
            C=params.cmaps;
        end
        r=1;
        d=1e2;
        if(isfield(params,'poly'))
            p=params.poly;
        else
            for i=1:params.roinum %get roi PSTH
                h=imellipse(A);
                p{i}=h.getVertices;                
                h.setColor(C(i,:));
            end        
        end
        for i=1:numel(p)
            in{i} = find(inpolygon(x,y,p{i}(:,1)/pxl2mm,p{i}(:,2)/pxl2mm));
            %rr(i)=min(ceil(numel(in{i})/500),10);
            rr(i)=params.psth_lines;
            d=min(d,numel(in{i}));
            if(params.t_col>0)
                cmap{i}=colormap_graded(C(i,:));
            else
                cmap{i}=repmat(C(i,:),64,1);
                t=randn(size(z));
            end
        end
        
        if(~isfield(params,'aux_axes'))
            figure;
            params.aux_axes=axes;
        end
        
        if(isfield(params,'rastind'))
%             d=min(d,1000);
            [out_arg.a_out,out_arg.b_out]=plot_sorted_raster(raster,t,params.ops,in,d,rr,params.binsize,cmap,params.bgcol,params.aux_axes,params.fit_mode);
        else %lfp
            trnum=params.ops.seg(1).LFPgroups{str2num(data_struct(1).fnames{z_col}(end))};
            out_arg.a_out=plot_sorted_lfp(z,t,in,data_struct(1).avtraces(:,:,trnum,1),params.ops,params.bgcol,params.aux_axes);
        end
        
        out_arg.p=p;
        out_arg.inds=in;
        
    end

    if(params.roc==1)
        C=[1 .4921 0;0 0.4921 .754;0.4921 0.746 0.4];    
        if(isfield(params,'poly'))
            p=params.poly;
        else
            h=imellipse(A);
            p=h.getVertices;
            h.setColor(C(1,:));
        end
        in = find(inpolygon(x,y,p(:,1)/pxl2mm,p(:,2)/pxl2mm));    
        out = find(~inpolygon(x,y,p(:,1)/pxl2mm,p(:,2)/pxl2mm));    
        targets=zeros(size(z));
        targets(in)=1;
        figure, cdfplot(z(in));
        hold on;
        cdfplot(z(out));
        figure, plotroc(targets',rescale(z'))
    end
end

function [N0,Mz,Sz]=tuning2d()
Mx=numel(params.x_edges)-1;
My=numel(params.y_edges)-1;
N0=zeros(params.x_nbins,params.y_nbins);
Mz=nan(params.x_nbins,params.y_nbins);
Sz=nan(params.x_nbins,params.y_nbins);
all_idx=find(x>params.x_edges(1) & x<=params.x_edges(end) &...
         y>params.y_edges(1) & y<=params.y_edges(end));

for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                 y>params.y_edges(j) & y<=params.y_edges(j+1));
        if(numel(idx)>0)
            if(isfield(params,'t_col') & isfield(params,'t_val')) %there is a target variable
                if(numel(params.t_val)==1)  %target mean                                    
                    inds=get_subpop_indices(t(idx),params.t_val);
                elseif(numel(params.t_val)==2) %boundaries
                    inds=find(t(idx)>=params.t_val(1) & t(idx)<params.t_val(2));
                else
                    inds=1:numel(idx);
                end            
                idx=idx(inds);
            end
        end
        N0(i,j)=numel(idx);
%         pval(i,j)=boottest_subset(z(idx),z(all_idx),1e3,@mean,'single');
        if(numel(idx)>=params.minsamp)
            Mz(i,j)=nanmean(z(idx));
            Sz(i,j)=nanstd(z(idx));
        end
    end
end
end

function [N0,P1,P2,Pv1]=fitting2d()
N0=zeros(params.x_nbins,params.y_nbins);
P1=nan(params.x_nbins,params.y_nbins);
P2=nan(params.x_nbins,params.y_nbins);
Pv1=nan(params.x_nbins,params.y_nbins);
% Pv2=nan(params.x_nbins,params.y_nbins);
% zz=(z-nanmean(z))/nanstd(z);
% tt=(t-nanmean(t))/nanstd(t);
zz=z;
tt=t;
for i=1:params.x_nbins
    for j=1:params.y_nbins
        idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                 y>params.y_edges(j) & y<=params.y_edges(j+1));
        idx=idx(~isnan(tt(idx)+zz(idx)));             
        N0(i,j)=numel(idx);
        if(numel(idx)>=max(params.minsamp,3))
            [f,g]=fit(tt(idx),zz(idx),'poly1');
            P1(i,j)=f.p1;
            P2(i,j)=f.p2;
            [rp,pval1]=corr(tt(idx),zz(idx));
%             [h0,pval2]=ttest2(zz(idx),zz);
            Pv1(i,j)=-log10(pval1);
%             Pv2(i,j)=-log10(pval2);
        end
    end
end

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

    function convert_to_wall()
%         file_struct(1).circle = file_struct(1).circle*pxl2mm;
        rx=file_struct(1).circle(3)/2; ry=file_struct(1).circle(4)/2;    %ellipse radii
        x0=file_struct(1).circle(1) + rx; y0=file_struct(1).circle(2) + ry; %center point
        phi=atan2((y-y0),(x-x0));    %azimuth in tank
        R1=hypot((x-x0),(y-y0)); %distance of fish from center
        R2=(rx*ry)./sqrt((ry*cos(phi)).^2 + (rx*sin(phi)).^2); %distance of nearest point from cetner
        ra=R2-R1;    %distance of fish to nearest wall
        ra(ra<0 | ra>params.maxr)=nan;        
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

