function [ strct ] = get_encoding_model( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.mode='cart'; %cart=cartesian ; polar=polar
params.nhidden=20;
params.cornersize=180;
params.freqmem=0;
params.freqmode='on';
params.t_col=[1 5:14];     %input columns
params.ref='circ'; %circ- circular tank; rect- rectangular tank; pole- pole object
params.objidx=1;
params.fb_cols=[];
%% get varins
n=numel(varargin);
if(n<3 | n>4&mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};  %output variable
lfp_nets=varargin{4};  %net struct for lfp processing

params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
params.a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));
params.lfb_cols=z_col;

if(n>4)
    i=5;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
params.K=numel(data_struct);
fnames=data_struct(1).fnames(params.t_col);
iiei=find(cellfun(@(x) numel(strfind(x,'iei')),fnames));

%% get data vectors
x=[]; y=[]; a=[]; z=[]; t=[]; fb=[]; lfb=[];
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
    if(numel(params.fb_cols) & params.fb_cols~=0)
        fb=[fb;data_struct(k).data(:,params.fb_cols)];
    end
    
    if(numel(params.lfb_cols) & params.lfb_cols~=0)
        lfb=[lfb;data_struct(k).data(:,params.lfb_cols)];
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
    
    %alternative cart calculation
    %get wall allocentric position
%     xw=x0+R2.*cos(phi);  yw=y0+R2.*sin(phi);
%     [y1,x1]=allo2ego(xw,yw,a,x,y); %obj coordinates rel. to LED    

end
%% convert to cartesian if in cartesian mode
M=300;
p=1.5;
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
    x_edges=linspace(min(x)*p,max(x)*p,M);
    y_edges=linspace(min(y)*p,max(y)*p,M+1);
    x_bins=edge2bin(x_edges);
    y_bins=edge2bin(y_edges);
    loc_bins={x_bins,y_bins};
end    

%% get regression model
% z=medfilt1(z,3);
fb=medfilt1(fb,3);
lfb=medfilt1(lfb,3);

% lfp->spikes
[ net_lfb,rsq_lfb,sig_lfb] = regression_network([lfb]',z',params.nhidden);

% all fb->spikes
[ net_fb,rsq_fb,sig_fb] = regression_network([fb]',z',params.nhidden);

% motor->spikes
[ net_ni,rsq_ni,sig_ni] = regression_network([t]',z',params.nhidden);

%feeadback + motor-> spikes
[ net_fbni,rsq_fbni,sig_fbni,prsq] = regression_network([t fb]',z',params.nhidden,'tansig',1);

%process lfp using network
lfb_ex=lfb-lfp_nets.model_fbni.net([t fb]')';
fb_ex=fb;
fb_ex(:,params.fb_cols==params.lfb_cols)=lfb_ex;

% lfp->spikes
[ net_lfb_f,rsq_lfb_f,sig_lfb_f] = regression_network([lfb_ex]',z',params.nhidden);

% all fb->spikes
[ net_fb_f,rsq_fb_f,sig_fb_f] = regression_network([fb_ex]',z',params.nhidden);

%feeadback + motor-> spikes
[ net_fbni_f,rsq_fbni_f,sig_fbni_f,prsq_f] = regression_network([t fb_ex]',z',params.nhidden,'tansig',1);

% ni filtering -> spikes
%% output
strct.z_col=z_col;
strct.fb_col=params.fb_cols;

strct.model_lfb.net=net_lfb;
strct.model_lfb.rsq=rsq_lfb;
strct.model_lfb.sig=sig_lfb;
strct.model_lfb.cols=params.lfb_cols;

strct.model_fb.net=net_fb;
strct.model_fb.rsq=rsq_fb;
strct.model_fb.sig=sig_fb;
strct.model_fb.cols=params.fb_cols;

strct.model_fbni.net=net_fbni;
strct.model_fbni.rsq=rsq_fbni;
strct.model_fbni.sig=sig_fbni;
strct.model_fbni.cols=[params.t_col params.fb_cols];

end

