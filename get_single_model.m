function [ strct ] = get_single_model( varargin )
%GET_SINGLE_MODEL create non-linear regression NN for eod response

pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.mode='cart'; %cart=cartesian ; polar=polar
params.hidmode='BP'; %BP=backprop; GC=granule-cell (random)
params.nhidden=100;
params.cornersize=180;
params.x_lim=[-150 150];   %default x bin limits for headfixed view
params.y_lim=[-150 150];   %default x bin limits for headfixed view
params.freqmem=0;
params.freqmode='on';
params.t_col=[1 5:14];     %input columns
params.ref='circ'; %circ- circular tank; rect- rectangular tank; pole- pole object
params.rmax=inf; %max distance from object to use for modeling
params.objidx=1;
params.sampnum=inf;
params.priorbins=200;
% params.fb_cols=[];
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
params.lfb_cols=z_col;

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
    if(strcmp(params.ref,'pole'))       
        objx=file_struct(k).objects(params.objidx).x;
        objy=file_struct(k).objects(params.objidx).y;
        [y1,x1]=allo2ego(objx,objy,a1,x1,y1); %obj coordinates rel. to LED
    end    
    
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
if(strcmp(params.ref,'rect'))
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
elseif(strcmp(params.ref,'circ')) %circular tank
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
M=params.priorbins;
p=.75;
if(~strcmp(params.ref,'pole')) %wall
    if(strcmp(params.mode,'cart'))
        x=r.*cos(th);
        y=r.*sin(th);
        loc=[x y];
        x_edges=linspace(min(x)*p,max(x)*p,M);
        y_edges=linspace(min(y),max(y),M+1);
        x_bins=edge2bin(x_edges);
        y_bins=edge2bin(y_edges);
        loc_bins={x_bins,y_bins};
    else
        loc=[r th];
        r_edges=linspace(0,max(r)*p,M);
        th_edges=linspace(-pi,pi,M+1);
        r_bins=edge2bin(r_edges);
        th_bins=edge2bin(th_edges);
        loc_bins={r_bins,th_bins};
    end
else %pole object
    loc=[x y];
    r=sqrt(x.^2+y.^2);
    x_edges=linspace(params.x_lim(1),params.x_lim(2),M);
    y_edges=linspace(params.y_lim(1),params.y_lim(2),M+1);
    x_bins=edge2bin(x_edges);
    y_bins=edge2bin(y_edges);
    loc_bins={x_bins,y_bins};
end    

%% get regression model
z=medfilt1(z,3);
% fb=medfilt1(fb,3);
%% downsample
if(~isinf(params.sampnum) & numel(z)>params.sampnum)
    ind=randperm(numel(z),params.sampnum);
    z=z(ind);
    loc=loc(ind,:);
    t=t(ind,:);
    r=r(ind);
end
%remove far samples
ind=find(r<=params.rmax);
z=z(ind);
loc=loc(ind,:);
if(numel(t))
    t=t(ind,:);
end

%% model
[ net_locp,rsq_locp,sig_locp,prsq] = regression_network([loc t]',z',params.nhidden,'tansig',1);

%prior pdf p(x,y)
[X,Y]=meshgrid(loc_bins{1},loc_bins{2});
[f,xi] = ksdensity(loc,[X(:) Y(:)]);
F=log(f);
F(isinf(F))=min(F(~isinf(F)));
[ net_prior,rsq_prior,sig_prior] = regression_network([X(:) Y(:)]',F',params.nhidden);

%% output
strct.z_col=z_col;
strct.t_col=params.t_col;

strct.model_single.net=net_locp;
strct.model_single.rsq=rsq_locp;
strct.model_single.sig=sig_locp;
strct.model_single.cols=params.t_col;
strct.model_single.prsq=prsq;

strct.prior.net=net_prior;
strct.prior.rsq=rsq_prior;
strct.prior.sig=sig_prior;

strct.loc.max=[max(loc_bins{1}),max(loc_bins{2})];
strct.loc.min=[min(loc_bins{1}),min(loc_bins{2})];
strct.loc.mean=[nanmean(loc(:,1)) nanmean(loc(:,2))];

strct.mot.max=max(t,[],1);
strct.mot.min=min(t,[],1);
strct.mot.mean=nanmean(t,1);

end

