function [ Mz,N0 ] = plot_tuning_polar( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.bg='on';
params.r_lim=[0 180]; %radius limits
params.th_lim=[-pi pi]; %azimouth limits
params.r_nbins=15;  %default r bin number
params.th_nbins=15;  %default th bin number
params.minsamp=5;
params.mina=0.15;%minimal alpha percentile
params.maxa=0.3; %maximal alpha percentile
params.fontsize=16;
params.freqmode='on';
params.t_col=0;     %targe column
params.t_val=0;  %target mean / boundaries
params.tank='circ'; %rectangular or circular tank
params.mfunc='mean'; %function to plot: mean-mean response in each bin; slope-sensitivity to t; offset- intersect at t=0
params.bgcol='w'; %background color, w=white, b=black
params.ind_lim=[];
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
params.a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));

if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
params.K=numel(data_struct);

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
end
%% convert to wall coordinates
R=params.r_lim(2);
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
    th=phi-a;   %egocentric angle of closest wall
    th(th>pi)=th(th>pi)-2*pi;
    th(th<-pi)=th(th<-pi)+2*pi;    
end
    

%% produce bins
params.r_edges=linspace(params.r_lim(1),params.r_lim(2),params.r_nbins+1);
params.th_edges=linspace(params.th_lim(1),params.th_lim(2),params.th_nbins+1);
params.r_bins=edge2bin(params.r_edges);
params.th_bins=edge2bin(params.th_edges);

%%
if(strcmp(params.mfunc,'mean'))
    [N0,Mz,Sz]=tuning2d();
    AL=N0;
else
    [N0,P1,P2,Pv]=fitting2d();
    if(strcmp(params.mfunc,'slope'))
        Mz=P1;
        AL=Pv;
    else
        Mz=P2;
        AL=N0;
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
    if(params.bgcol=='w')
        imagesc(ix,iy,params.image);
    else
        imagesc(ix,iy,255-params.image);
    end
    hold on;
end

S=polarplot3d(Mzz,AL,'AngularRange',[-pi pi],'plottype','surfa','RadialRange',[0 150]);
if(isfield(params,'clim'))
    set(A,'CLim',params.clim);
end
if(strcmp(params.mfunc,'slope'))
    set(A,'ALim',[-log10(0.05) -log10(0.01)]);
else
    set(A,'ALim',[quantile(N0(:),params.mina) quantile(N0(:),params.maxa)]);
end
view(0,90);
set(gca,'FontSize',params.fontsize)
set(gca,'XLim',params.r_lim(2)*[-1 1],'YLim',params.r_lim(2)*[-1 1]);
set(gca,'Ydir','normal');
axis('image');
hold off;
colorbar;
if(params.bgcol=='k')
    set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);
    set(gcf,'Color',[0 0 0]);
    colorbar(gca,'Color',[1 1 1]);
end

function [N0,Mz,Sz]=tuning2d()
Mx=numel(params.r_edges)-1;
My=numel(params.th_edges)-1;
N0=zeros(params.r_nbins,params.th_nbins);
Mz=nan(params.r_nbins,params.th_nbins);
Sz=nan(params.r_nbins,params.th_nbins);
for i=1:params.r_nbins
    for j=1:params.th_nbins
        idx=find(r>params.r_edges(i) & r<=params.r_edges(i+1) &...
                 th>params.th_edges(j) & th<=params.th_edges(j+1));
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

function [N0,P1,P2,Pv]=fitting2d()
N0=zeros(params.r_nbins,params.th_nbins);
P1=nan(params.r_nbins,params.th_nbins);
P2=nan(params.r_nbins,params.th_nbins);
zz=(z-nanmean(z))/nanstd(z);
tt=(t-nanmean(t))/nanstd(t);
ff=figure;
for i=1:params.r_nbins
    for j=1:params.th_nbins
        idx=find(r>params.r_edges(i) & r<=params.r_edges(i+1) &...
            th>params.th_edges(j) & th<=params.th_edges(j+1));
        idx=idx(~isnan(tt(idx)+zz(idx)));             
        N0(i,j)=numel(idx);
        if(numel(idx)>=params.minsamp)
            f=fit(tt(idx),zz(idx),'poly1');
%             if(i==8 & j==2)
%                 figure(ff);
%                 H=plot(tt(idx),zz(idx),'.');
%                 set(H,'MarkerSize',18);
%                 hold on;
%                 H=plot(f);
%                 set(H,'Color',[0 0 0],'LineWidth',3);
%                 legend('off');
%                 hold off;
%                 title(['i=',num2str(i),'; j=',num2str(j)]);
%                 set(gcf,'Color',[1 1 1]);
%                 drawnow;
%             end
            P1(i,j)=f.p1;
            P2(i,j)=f.p2;
            [rp,pval]=corr(tt(idx),zz(idx));
%             t_slope=rp*sqrt((numel(idx)-2)/(1-rp^2));
%              N0(i,j)=-log10(pval);            
            Pv(i,j)=-log10(pval);
        end
    end
end
close(ff);
end

end

