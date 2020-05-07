function [ strct ] = get_polar_model( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.mode='cart'; %cart=cartesian ; polar=polar
params.nhidden=10;
params.cornersize=180;
params.freqmem=0;
params.freqmode='on';
params.t_col=[1 5:14];     %input columns
params.tank='circ'; %rectangular or circular tank
%% get varins
n=numel(varargin);
if(n<3 | n>3&~mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};  %output variable

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
fnames=data_struct(1).fnames(params.t_col);
iiei=find(cellfun(@(x) numel(strfind(x,'iei')),fnames));

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
        if(params.freqmem>0)
            iei1=t(:,iiei);
            for i=1:params.freqmem
                iei1=[nan;iei1(1:end-1)];
                t=[iei1 t];
                fnames={['iei',num2str(i)],fnames{:}};
            end
            iiei=find(cellfun(@(x) numel(strfind(x,'iei')),fnames));
        end
        
    end    
end
x=x*pxl2mm;
y=y*pxl2mm;
%z=detrend(z);

if(isfield(params,'func') & params.t_col~=0)
    if(numel(params.t_col)>1)
        t=params.func(t,2);
    else
        t=params.func(t);
    end
end

if(params.t_col~=0  & ...
    numel(iiei) & strcmp(params.freqmode,'on'))
    dt=nanmean(diff(sort(t(:,iiei(1)))));
    t(:,iiei)=1./(t(:,iiei)+dt*randn(size(t,1),params.freqmem+1));    
end
if(strcmp(data_struct(1).fnames{z_col},'iei') & strcmp(params.freqmode,'on'))
    z=1./z;
end

%% convert to wall coordinates
R=params.cornersize;
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
        t=[t(ind1,:);t(ind2,:);t(ind3,:);t(ind4,:)];
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
%% convert to cartesian if in cartesian mode
M=100;

if(strcmp(params.mode,'cart'))
    x=r.*cos(th);
    y=r.*sin(th);
    loc=[x y];
    x_edges=linspace(min(x),max(x),M);
    y_edges=linspace(min(y),max(y),M+1);
    x_bins=edge2bin(x_edges);
    y_bins=edge2bin(y_edges);
    loc_bins={x_bins,y_bins};
else
    loc=[r th];
    r_edges=linspace(min(r),max(r),M);
    th_edges=linspace(-pi,pi,M+1);
    r_bins=edge2bin(r_edges);
    th_bins=edge2bin(th_edges);
    loc_bins={r_bins,th_bins};
end

%% get regression model
% t0=[1 nanmean(t,1)]; %mean posture
% thx=cos(th);
% thy=sin(th);
z(z<-3.3)=nan;
% [ net_rth,rsq_rth,sig_rth ] = regression_network(loc',z',params.nhidden);
% z_ex=net_rth(loc'); %ex-afference signal
% z_re=z'-z_ex;   %reafference component
% [ net_ni,rsq_ni,sig_ni ] = regression_network([t]',z_re,params.nhidden);
[ net_rthp,rsq_rthp,sig_rthp ] = regression_network([loc t]',z',params.nhidden);

%prior pdf p(x,y)
[X,Y]=meshgrid(loc_bins{1},loc_bins{2});
[f,xi] = ksdensity(loc,[X(:) Y(:)]);
F=log(f);
F(isinf(F))=min(F(~isinf(F)));
[ net_prior,rsq_prior,sig_prior] = regression_network([X(:) Y(:)]',F',params.nhidden);

%normalization pdf p(z)
xi=linspace(nanmin(z),nanmax(z),1e4);
[f,xi] = ksdensity(z,xi(:));
[ net_norm,rsq_norm,sig_norm ] =regression_network(xi(:)',f(:)',params.nhidden);

%cond. normalization pdf p(z|z_-1)
xi=linspace(nanmin(z),nanmax(z),300);
xj=linspace(nanmin(z),nanmax(z),301);
[X,Y]=meshgrid(xi,xj);
[fj,xj] = ksdensity(z,xj(:));
[fij,xij]= ksdensity([z(1:end-1) z(2:end)],[X(:) Y(:)]); %joint prob.
[ net_joint,rsq_joint,sig_joint] =regression_network(xij(:)',fij(:)',50);
fc=reshape(fij,301,300)./repmat(fj,1,300);

% N0=reshape(f,params.th_nbins,params.r_nbins)';
% p_rth=reshape(net_rth([R(:) TH(:)]'),params.th_nbins,params.r_nbins)';
%% output
strct.model_loc.net=net_rth;
strct.model_loc.rsq=rsq_rth;
strct.model_loc.sig=sig_rth;

strct.model_ni.net=net_ni;
strct.model_ni.rsq=rsq_ni;
strct.model_ni.sig=sig_ni;

strct.model_loc_mot.net=net_rthp;
strct.model_loc_mot.rsq=rsq_rthp;
strct.model_loc_mot.sig=sig_rthp;

strct.norm.net=net_norm;
strct.norm.rsq=rsq_norm;
strct.norm.sig=sig_norm;

strct.prior.net=net_prior;
strct.prior.rsq=rsq_prior;
strct.prior.sig=sig_prior;

strct.loc.max=max(loc,[],1);
strct.loc.min=min(loc,[],1);
strct.loc.mean=[nanmean(loc(:,1)) nanmean(loc(:,2))];

strct.mot.max=max(t,[],1);
strct.mot.min=min(t,[],1);
strct.mot.mean=nanmean(t,1);

end

