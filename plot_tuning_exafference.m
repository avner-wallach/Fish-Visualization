function [ rsq ] = plot_tuning_exafference( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% params
params.bg='on';
params.r_lim=[0 180]; %radius limits
params.th_lim=[-pi pi]; %azimouth limits
params.r_nbins=25;  %default r bin number
params.th_nbins=25;  %default th bin number
params.minsamp=5;
params.mina=0.15;%minimal alpha percentile
params.maxa=0.3; %maximal alpha percentile
params.fontsize=16;
params.freqmode='on';
params.t_col=[1 5:14];     %targe column
% params.t_val=0;  %target mean / boundaries
params.objidx=0;    %object: 0=no object, wall analysis; >0= object centered
params.tank='circ'; %rectangular or circular tank
params.mfunc='mean'; %function to plot: mean-mean response in each bin; std-std response; slope-sensitivity to t; offset- intersect at t=0
params.bgcol=[1 1 1]; %background color,
params.ind_lim=[];
params.nhidden=10;
params.fb_cols=[];

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
x=[]; y=[]; a=[]; z=[]; t=[]; fb=[];
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
    if(numel(params.fb_cols) & params.fb_cols~=0)
        fb=[fb;data_struct(k).data(:,params.fb_cols)];
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
    if(numel(params.ind_lim)==2)
        inds=inrange(1:numel(z),params.ind_lim);
    elseif(numel(params.ind_lim==1))
        inds=randperm(numel(z),params.ind_lim);
    end
    z=z(inds);
    if(numel(t))
        t=t(inds,:);
    end
    x=x(inds);
    y=y(inds);
    a=a(inds);
%     if(params.roinum>0)
%         raster=raster(inrange(raster(:,2),params.ind_lim),:);
%         raster(:,2)=raster(:,2)-params.ind_lim(1);
%     end
end
%% convert to wall coordinates
R=params.r_lim(2);
r=[];
th=[];
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
%% filter
z=medfilt1(z,3);
% fb=medfilt1(fb,3);
% 
loc=[r th];
%% get regression models
% [ net_ex,rsq_ex,sig_ex] = ni_network(loc',z');
% [ net_re,rsq_re,sig_re] = ni_network(t',z');
% % z_ex=net_ex(loc'); %ex-afference signal
% % z_re=z'-z_ex;   %reafference component
% % 
% % % get all model
% [ net_all,rsq_all,sig_all] = ni_network([loc t]',z');
% %scramble loc
% L=loc(randperm(size(loc,1)),:);
% [ net_scr,rsq_ex_scr,sig_ex_scr]= ni_network([L t]',z');
% %scramble mot
% T=t(randperm(size(t,1)),:); %scramble all motor
% [ net_scr,rsq_mot_scr,sig_mot_scr]= ni_network([loc T]',z');
% [ net_scr,rsq_all_scr,sig_all_scr]= ni_network([L T]',z');
% rsq.ex=rsq_ex;
% rsq.re=rsq_re;
% rsq.all=rsq_all;
% rsq.ex_ctl=rsq_ex_scr;
% rsq.re_ctl=rsq_mot_scr;
% rsq.ctl=rsq_all_scr;
% itemize and scramble motor info
% for i=1:size(t,2)
%     itemize
%     [net_i,rsq_i(i),sig_i(i)] = ni_network(t(:,i)',z');
%     T=t;
%     T(:,i)=T(randperm(size(T,1)),i); %scramble
%     [ net_scr,rsq_scr(i),sig_scr(i)]= ni_network([loc T]',z');
% end
% % NI models
% [ net_ni,rsq_ni,sig_ni ] = ni_network([t]',z_re); %motor only
% [ net_lfb,rsq_lfb,sig_lfb] = ni_network(z',z_re); %local only
% [ net_fb,rsq_fb,sig_fb] = ni_network([fb]',z_re); %global only
% [ net_lfbni,rsq_lfbni,sig_lfbni] = ni_network([t z]',z_re); %local + motor
% [ net_fbni,rsq_fbni,sig_fbni] = ni_network([t fb]',z_re); %local + motor
% prior pdf p(x,y)
% [X,Y]=meshgrid(params.r_bins,params.th_bins);
% [f,xi] = ksdensity(loc,[X(:) Y(:)]);
% F=log(f);
% F(isinf(F))=min(F(~isinf(F)));
% [ net_prior,rsq_prior,sig_prior] = ni_network([X(:) Y(:)]',F');
%% get exafference map
zz=(z-nanmean(z))/nanstd(z);
[N0,Mz,Sz]=tuning2d();
AL=N0;

%interpolate exafference 
z_ex=interp2(params.r_bins,params.th_bins,Mz',r,th);

%plot data vs. exafference
S=scatter(zz,z_ex,25,0*[1 1 1],'filled');
set(S,'MarkerFaceAlpha',0.1,'MarkerFaceColor',1-params.bgcol,'MarkerEdgeAlpha',0);
hold on;
f1=fit(zz(~isnan(zz.*z_ex)),z_ex(~isnan(zz.*z_ex)),'poly1');
v=[nanmin(z_ex) nanmax(z_ex)];
% H=plot(v,v,':k');
% set(H,'LineWidth',3);
H=plot(f1);
COL=brewermap(9,'Set1');
set(H,'LineWidth',3,'Color',COL(2,:),'LineStyle',':');
rsq=nancorr(zz,z_ex).^2;
set(gca,'Color',params.bgcol,'XColor',1-params.bgcol,'YColor',1-params.bgcol,'FontSize',16,'XLabel',[],'YLabel',[]);
legend('off');

% ind=find(~isnan(z_ex) & ~isnan(zz));
% [f1,g1]=fit(z_ex(ind),zz(ind),'poly1');
% plot(f1);
% rsq=g1.rsquare;
% set(gca,'FontSize',params.fontsize)
% set(gca,'XLim',params.r_lim(2)*[-1 1],'YLim',params.r_lim(2)*[-1 1]);
% set(gca,'Ydir','normal');
% if(params.bgcol=='k')
%     set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);
%     set(gcf,'Color',[0 0 0]);
%     colorbar(gca,'Color',[1 1 1]);
% end


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

