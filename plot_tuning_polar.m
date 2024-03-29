function [ Mz,N0,out_arg ] = plot_tuning_polar( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% params
params.bg='on';
params.r_lim=[0 180]; %radius limits
params.th_lim=[-pi pi]; %azimouth limits
params.r_nbins=15;  %default r bin number
params.th_nbins=15;  %default th bin number
params.minsamp=25;
params.fontsize=16;
params.freqmode='on';
params.t_col=0;     %targe column
% params.t_val=0;  %target mean / boundaries
params.objidx=0;    %object: 0=no object, wall analysis; >0= object centered
params.tank='circ'; %rectangular or circular tank
params.mfunc='mean'; %function to plot: mean-mean response in each bin; std-std response; slope-sensitivity to t; offset- intersect at t=0
params.bgcol=[1 1 1]; %background color
params.ind_lim=[];
params.zscore=1;
params.roinum=0;
params.binsize=.5;
params.upsamp=.25;
r=[]; th=[];
out_arg=[];
% params.clim='auto';
%% get varins
n=numel(varargin);
if(n<3 | n>3&~mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};

%convert ungrouped data
if(numel(data_struct.data)==0)
    data_struct.fnames=data_struct.fnames0;
    data_struct.data=data_struct.data0;
end
% params.ind=1:size(data_struct(1).data,1);
params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
params.a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));

if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
params.K=numel(data_struct);
%% get pixel size
if(isfield(params,'ops'))
    pxl2mm=params.ops.pxlsize;
else
    pxl2mm=str2num(getenv('PXLSIZE'));
end
%% get data vectors
x=[]; y=[]; a=[]; z=[]; t=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    a1=data_struct(k).data(:,params.a_col);
    
    x=[x;x1];
    y=[y;y1];
    a=[a;a1];
    z=[z;data_struct(k).data(:,z_col)];
    if(params.t_col~=0)
        t=[t;data_struct(k).data(:,params.t_col)];
    end    
end
x=x*pxl2mm;
y=y*pxl2mm;

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
%% reduce indices
if(numel(params.ind_lim))
    inds=inrange(1:numel(z),params.ind_lim);
    z=z(inds);
    if(numel(t))
        t=t(inds,:);
    end
    x=x(inds);
    y=y(inds);
    a=a(inds);
    if(params.roinum>0 & strfind(data_struct(1).fnames{z_col},'sc'))
        raster=raster(inrange(raster(:,2),params.ind_lim),:);
        raster(:,2)=raster(:,2)-params.ind_lim(1);
    end
end
%% convert to wall coordinates
R=params.r_lim(2);
if(params.objidx==0)
    convert_to_wall;
else
    convert_to_object;
end

%% produce bins
params.r_edges=linspace(params.r_lim(1),params.r_lim(2),params.r_nbins+1);
params.th_edges=linspace(params.th_lim(1),params.th_lim(2),params.th_nbins+1);
params.r_bins=edge2bin(params.r_edges);
params.th_bins=edge2bin(params.th_edges);
%% zscore/dezscore data 
if(isfield(params,'zscale'))
    z=z*params.zscale(2) + params.zscale(1);
else
    if(params.zscore)
        z=(z-nanmean(z))/nanstd(z);
        t=(t-nanmean(t))/nanstd(t);
    end
end
%%
if(strcmp(params.mfunc,'mean') | strcmp(params.mfunc,'std'))
    [N0,Mz,Sz]=tuning2d();
    if(strcmp(params.mfunc,'std'))
        Mz=Sz;
    end
    AL=min(N0,15);
else
    [N0,P1,P2,Pv]=fitting2d();
    if(strcmp(params.mfunc,'slope'))
        Mz=P1;
        AL=min(Pv,-log10(0.01));
    else
        Mz=P2;
        AL=min(N0,15);
    end
end

Mzz=[Mz Mz(:,1)];
AL=[AL AL(:,1)];
F=gcf;
A=gca;
if(isfield(params,'image'))
    sX=size(params.image,2);
    sY=size(params.image,1);
    ix=([1:sX]-sX/2)*pxl2mm;
    iy=(-[1:sY]+sY/3)*pxl2mm;
    if(params.bgcol==[1 1 1])
        imagesc(ix,iy,params.image);
    else
        imagesc(ix,iy,255-params.image);
    end
    hold on;
end

S=polarplot3d(Mzz,AL,'AngularRange',[-pi pi],'plottype','surfa','RadialRange',params.r_lim,'AxisLocation','off','MeshScale',params.upsamp*[1 1]);
if(strcmp(params.mfunc,'slope'))
    set(A,'clim',[-1 1]);
    set(A,'ALim',[-log10(0.1) -log10(0.01)+1]);
else
    set(A,'ALim',[0 20]);
end
if(isfield(params,'alim'))
    set(A,'ALim',params.alim);
end
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
view(0,90);
set(gca,'FontSize',params.fontsize)
set(gca,'XLim',params.r_lim(2)*[-1 1],'YLim',params.r_lim(2)*[-1 1]);
set(gca,'Ydir','normal');
if(strcmp(params.mfunc,'slope'))
    if(params.bgcol==[0 0 0])
        colormap(gca,invert_map(flipud(brewermap(64,'RdBu'))))
    else
%         colormap(gca,modified_jet);
        colormap(gca,flipud(brewermap(64,'RdBU')));
    end
else
    colormap(gca,'parula');
end

axis('image');
hold off;
set(gca,'Color','none','XColor','none','YColor','none','XGrid','off','YGrid','off');
colorbar(gca,'Color',1-params.bgcol,'FontSize',params.fontsize,'FontWeight','normal','LineWidth',1);
if(isfield(params,'poly'))
    params.roinum=numel(params.poly);
end

if(params.roinum>0)
    if(~isfield(params,'cmaps'))
        C=brewermap(10,'Set1');
    else
        C=params.cmaps;
    end
    rr=1;
    d=1e2;
    for i=1:params.roinum %get roi PSTH
        if(isfield(params,'poly'))
            p{i}=params.poly{i};
        else
            h=imellipse(A);
            p{i}=h.getVertices;
    %         h=impoly(A);
    %         p=h.getPosition;
            h.setColor(C(i,:));
        end
        x=r.*sin(th); y=r.*cos(th);
        in{i} = find(inpolygon(x,y,p{i}(:,1),p{i}(:,2)));
        rr(i)=3;%min(ceil(numel(in{i})/1e3),10);
        d=min(d,numel(in{i}));
       cmap{i}=colormap_graded(C(i,:));        
    end
    if(~isfield(params,'aux_axes'))
        figure;
        params.aux_axes=axes;
    end            
    if(strfind(data_struct(1).fnames{z_col},'sc')) %spikes
        a_out=plot_sorted_raster(raster,t,params.ops,in,d,rr,params.binsize,cmap,params.bgcol,params.aux_axes);
    else %lfp
        trnum=params.ops.seg(1).LFPgroups{str2num(data_struct(1).fnames{z_col}(end))};
        a_out=plot_sorted_lfp(z,t,in,data_struct(1).avtraces(:,:,trnum,1),params.ops,params.bgcol,params.aux_axes);
    end
%     if(isfield(params,'poly'))
        out_arg.inds=in;
%     else
        out_arg.p=p;
        out_arg.a_out=a_out;
%     end
end


function [N0,Mz,Sz]=tuning2d()
Mx=numel(params.r_edges)-1;
My=numel(params.th_edges)-1;
N0=zeros(params.r_nbins,params.th_nbins);
Mz=nan(params.r_nbins,params.th_nbins);
Sz=nan(params.r_nbins,params.th_nbins);
% zz=(z-nanmean(z))/nanstd(z);
zz=z;
for i=1:params.r_nbins
    for j=1:params.th_nbins
        idx=find(r>params.r_edges(i) & r<=params.r_edges(i+1) &...
                 th>params.th_edges(j) & th<=params.th_edges(j+1));
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
        if(numel(idx)>=params.minsamp)
            Mz(i,j)=nanmean(zz(idx));
            Sz(i,j)=nanstd(zz(idx));
        end
    end
end
end

function [N0,P1,P2,Pv]=fitting2d()
N0=zeros(params.r_nbins,params.th_nbins);
Pv=zeros(params.r_nbins,params.th_nbins);
P1=nan(params.r_nbins,params.th_nbins);
P2=nan(params.r_nbins,params.th_nbins);
zz=z;%(z-nanmean(z))/nanstd(z);
tt=t;%(t-nanmean(t))/nanstd(t);
ff=figure;
for i=1:params.r_nbins
    for j=1:params.th_nbins
        idx=find(r>params.r_edges(i) & r<=params.r_edges(i+1) &...
            th>params.th_edges(j) & th<=params.th_edges(j+1));
        idx=idx(~isnan(tt(idx)+zz(idx)));             
        N0(i,j)=numel(idx);
        if(numel(idx)>=max(params.minsamp,3))
            f=fit(tt(idx),zz(idx),'poly1');
            P1(i,j)=f.p1;
            P2(i,j)=f.p2;
            [rp,pval]=corr(tt(idx),zz(idx));
            Pv(i,j)=-log10(pval);
        end
    end
end
close(ff);
end

function convert_to_wall()
    if(strcmp(params.tank,'rect'))
        % corner locations
        corners=sort(file_struct(1).corners);
        x_cor=[mean(corners(1:2,1)) mean(corners(3:4,1))]*pxl2mm;
        y_cor=[mean(corners(1:2,2)) mean(corners(3:4,2))]*pxl2mm;

        ind1=find(x<(x_cor(1)+R) & y>(y_cor(1)+R) & y<(y_cor(2)-R)); %west
        ind2=find(x>(x_cor(2)-R) & y>(y_cor(1)+R) & y<(y_cor(2)-R)); %east
        ind3=find(y<(y_cor(1)+R) & x>(x_cor(1)+R) & x<(x_cor(2)-R)); %north
        ind4=find(y>(y_cor(2)-R) & x>(x_cor(1)+R) & x<(x_cor(2)-R)); %south

        r=[x(ind1)-x_cor(1);...
            x_cor(2)-x(ind2);...
            y(ind3)-y_cor(1);...
            y_cor(2)-y(ind4)];
        th=[-pi-a(ind1);...
            0-a(ind2);...
            -pi/2-a(ind3);...
            pi/2-a(ind4)];
        th(th>pi)=th(th>pi)-2*pi;
        th(th<-pi)=th(th<-pi)+2*pi;
        z=[z(ind1);z(ind2);z(ind3);z(ind4)];
        if(params.t_col~=0)
            t=[t(ind1);t(ind2);t(ind3);t(ind4)];
        end
    else %circular tank
        file_struct(1).circle = file_struct(1).circle*pxl2mm;
        rx=file_struct(1).circle(3)/2; ry=file_struct(1).circle(4)/2;    %ellipse radii
        x0=file_struct(1).circle(1) + rx; y0=file_struct(1).circle(2) + ry; %center point
        phi=atan2((y-y0),(x-x0));    %azimuth in tank
        R1=hypot((x-x0),(y-y0)); %distance of fish from center
        R2=(rx*ry)./sqrt((ry*cos(phi)).^2 + (rx*sin(phi)).^2); %distance of nearest point from cetner
        r=R2-R1;    %distance of fish to nearest wall
%         r(r<=10)=nan;
        th=phi-a;   %egocentric angle of closest wall
        th(th>pi)=th(th>pi)-2*pi;
        th(th<-pi)=th(th<-pi)+2*pi;    
    end   
end

function convert_to_object()
      objx=file_struct(k).objects(params.objidx).x*pxl2mm;
      objy=file_struct(k).objects(params.objidx).y*pxl2mm;
      [x,y]=allo2ego(objx,objy,a,x,y); %obj coordinates rel. to LED
      th=atan2(x,y);    %azimuth of object
      r=hypot(x,y); %distance of object from fish's head      
end
    
end

