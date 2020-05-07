function [  ] = evaluate_decoding_model( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=str2num(getenv('PXLSIZE'));
figure;
COL=colormap('lines');
% COL(4,:)=0.5*[1 1 1];
% fcolor=[0 0 0];
fcolor=[1 1 1];
close(gcf);
%% params
params.mode='cart'; %cart=cartesian ; polar=polar
params.nhidden=100;
params.cornersize=180;
params.freqmem=0;
params.freqmode='on';
params.fontsize=14;
params.msize=8;
params.t_col=[1 5:14];     %input columns
params.ref='circ'; %rect=rectangular; circ=circular tank; pole=pole object
params.r_max=0;
%% get varins
n=numel(varargin);
if(n<4 | n>4&mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};  %output variable
net_struct=varargin{4};
fb_col=net_struct.fb_col;

params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
params.a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));

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
x=[]; y=[]; a=[]; z=[]; d=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    a1=data_struct(k).data(:,params.a_col);
    
    x=[x;x1];
    y=[y;y1];
    a=[a;a1];
    z=[z;data_struct(k).data(:,z_col)];
    d=[d;data_struct(k).data];
%     if(params.t_col~=0)
%         t=[t;data_struct(k).data(:,params.t_col)];
%         if(params.freqmem>0)
%             iei1=t(:,iiei);
%             for i=1:params.freqmem
%                 iei1=[nan;iei1(1:end-1)];
%                 t=[iei1 t];
%                 fnames={['iei',num2str(i)],fnames{:}};
%             end
%             iiei=find(cellfun(@(x) numel(strfind(x,'iei')),fnames));
%         end
%         
%     end    
%     fb=[fb;data_struct(k).data(:,fb_col)];
    
end
x=x*pxl2mm;
y=y*pxl2mm;
%z=detrend(z);
z=medfilt1(z,3);
% fb=medfilt1(fb,3);
if(strcmp(params.freqmode,'on'))
    d(:,1)=1./d(:,1);
end

% if(isfield(params,'func') & params.t_col~=0)
%     if(numel(params.t_col)>1)
%         t=params.func(t,2);
%     else
%         t=params.func(t);
%     end
% end

% if(params.t_col~=0  & ...
%     numel(iiei) & strcmp(params.freqmode,'on'))
%     dt=nanmean(diff(sort(t(:,iiei(1)))));
%     t(:,iiei)=1./(t(:,iiei)+dt*randn(size(t,1),params.freqmem+1));    
% end
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
else %pole
    %object position
    xp = file_struct(1).objects.x*pxl2mm;
    yp = file_struct(1).objects.y*pxl2mm;
    %egocentric object position
    [y,x]=allo2ego(xp,yp,a,x,y); 

end
%% convert to cartesian if in cartesian mode
M=30;
if(~strcmp(params.ref,'pole'))
    if(strcmp(params.mode,'cart'))
        x=r.*cos(th);
        y=r.*sin(th);
        loc=[x y];
        x_edges=linspace(min(x),max(x),M);
        y_edges=linspace(min(y),max(y),M+1);
        x_bins=edge2bin(x_edges);
        y_bins=edge2bin(y_edges);
        loc_edges={x_edges,y_edges};    
        loc_bins={x_bins,y_bins};
    else
        loc=[r th];
        r_edges=linspace(min(r),max(r),M);
        th_edges=linspace(-pi,pi,M+1);
        r_bins=edge2bin(r_edges);
        th_bins=edge2bin(th_edges);
        loc_edges={x_edges,y_edges};        
        loc_bins={r_bins,th_bins};
    end
else %pole
    loc=[x y];
    r=sqrt(x.^2+y.^2);
    x_edges=linspace(min(x),max(x),M);
    y_edges=linspace(min(y),max(y),M+1);
    x_bins=edge2bin(x_edges);
    y_bins=edge2bin(y_edges);
    loc_edges={x_edges,y_edges};    
    loc_bins={x_bins,y_bins};    
end
%% limit range
if(params.r_max>0)
    ind=find(r<=params.r_max);
    loc=loc(ind,:);
    z=z(ind);
    d=d(ind,:);
end

%% evaluate decoding using regression model
z_ctl=z(randperm(numel(z)));
F1=figure;
F2=figure;
F3=figure;
acclag=50;
pq=[0.25 0.5 0.75];

prior=exp(net_struct.prior.net(loc')');
normz=net_struct.norm.net(z')';
normz_ctl=net_struct.norm.net(z_ctl')';
zjoint=[nan;net_struct.joint.net([z(1:end-1) z(2:end)]')']; %joint prob. z_n,z_n-1
zcond=zjoint./normz;

% motor agnostic model
z_reg=net_struct.model_loc.net(loc')';
sig_reg=net_struct.model_loc.sig;
l_loc=(1/sqrt(2*pi*sig_reg^2))*exp(-(z_reg-z).^2/(2*sig_reg^2));
l_loc_ctl=(1/sqrt(2*pi*sig_reg^2))*exp(-(z_reg-z_ctl).^2/(2*sig_reg^2));
p_loc=l_loc.*prior./normz;
% p_loc_ctl=l_loc_ctl.*prior./normz_ctl;
p_loc_ctl=prior;
[pval_loc,N0,pmean_loc,pmean_loc_ctl]=pval2d(p_loc,p_loc_ctl);

lratio=l_loc_ctl./zcond;
pshift=[nan;p_loc_ctl(1:end-1)];
for i=1:acclag
    p_acc=pshift.*movprod(lratio,i)./p_loc_ctl;
    pshift=[nan;pshift(1:end-1)];
    acc_ratio(i,:)=quantile(p_acc,pq);
    pconv(i)=sum(p_acc>1)/numel(p_acc);
end
plotconv(F3,pconv,0);

lratio=l_loc./zcond;
pshift=[nan;p_loc(1:end-1)];
for i=1:acclag
    p_acc=pshift.*movprod(lratio,i)./p_loc;
    pshift=[nan;pshift(1:end-1)];
    acc_ratio(i,:)=quantile(p_acc,pq);
    pconv(i)=sum(p_acc>1)/numel(p_acc);
end
plotconv(F3,pconv,1);
    
Mz=log(pmean_loc);
pval_loc(pval_loc<0.01)=0.01;
AL=-log(pval_loc);
% AL(pval_loc==0)=2*nanmax(-log(pval_loc(pval_loc>0)));
figure, plot2d;
if(nanmin(pval_loc(:))<0.05)
    set(gca,'Alim',[-log(0.05) nanmax(AL(:))]);
else
    set(gca,'Alim',[nanmax(AL(:)) -log(0.05)]);
end
ploterror(z_ctl-z_reg(randperm(numel(z_reg))),F1,F2,0);
ploterror(z-z_reg,F1,F2,1);

% nimodels={'ni','lfb','fb','fbni'};
nimodels={'ni','fbni'};
for m=1:numel(nimodels)
    net=net_struct.(['model_',nimodels{m}]).net;
    cols=net_struct.(['model_',nimodels{m}]).cols;
    e_reg=net(d(:,cols)')';
    z_loc=net_struct.(['model_',nimodels{m},'loc']).net(loc')';
    sig_reg=net_struct.(['model_',nimodels{m},'loc']).sig;
    z_reg=z-e_reg;
    z_reg_ctl=z_ctl-e_reg;
    l_ni=(1/sqrt(2*pi*sig_reg^2))*exp(-(z_reg-z_loc).^2/(2*sig_reg^2));
    l_ni_ctl=(1/sqrt(2*pi*sig_reg^2))*exp(-(z_reg_ctl-z_loc).^2/(2*sig_reg^2));
    p_ni=l_ni.*prior./normz;
    p_ni_ctl=prior;
    [pval_ni,N0,pmean_ni,pmean_ni_ctl]=pval2d(p_ni,p_ni_ctl);

    lratio=l_ni./zcond;
    pshift=[nan;p_ni(1:end-1)];
    for l=1:acclag
        p_acc=pshift.*movprod(lratio,l)./p_ni;
        pshift=[nan;pshift(1:end-1)];
        acc_ratio(l,:)=quantile(p_acc,pq);
        pconv(l)=sum(p_acc>1)/numel(p_acc);
    end
    plotconv(F3,pconv,m+1);

    Mz=log(pmean_ni);
    pval_ni(pval_ni<0.01)=0.01;
    AL=-log(pval_ni);
    figure, plot2d;
    if(nanmin(pval_ni(:))<0.05)
        set(gca,'Alim',[-log(0.05) nanmax(AL(:))]);
    else
        set(gca,'Alim',[nanmax(AL(:)) -log(0.05)]);
    end

    ploterror(z_reg-z_loc,F1,F2,m+1);
end

% % SMC model
% z_reg=net_struct.model_loc_mot.net([loc t]')';
% sig_reg=net_struct.model_loc_mot.sig;
% l_loc_mot=(1/sqrt(2*pi*sig_reg^2))*exp(-(z_reg-z).^2/(2*sig_reg^2));
% l_loc_mot_ctl=(1/sqrt(2*pi*sig_reg^2))*exp(-(z_reg-z_ctl).^2/(2*sig_reg^2));
% p_loc_mot=l_loc_mot.*prior./normz;
% % p_loc_mot_ctl=l_loc_mot_ctl.*prior./normz_ctl;
% p_loc_mot_ctl=prior;
% [pval_loc_mot,N0,pmean_loc_mot,pmean_loc_mot_ctl]=pval2d(p_loc_mot,p_loc_mot_ctl);
% 
% lratio=l_loc_mot./zcond;
% pshift=[nan;p_loc_mot(1:end-1)];
% for i=1:acclag
%     p_acc=pshift.*movprod(lratio,i)./p_loc_mot;
%     pshift=[nan;pshift(1:end-1)];
%     acc_ratio(i,:)=quantile(p_acc,pq);
%     pconv(i)=sum(p_acc>1)/numel(p_acc);
% end
% plotconv(F3,pconv,3);
% 
% Mz=log(pmean_loc_mot);
% pval_loc_mot(pval_loc_mot<0.01)=0.01;
% AL=-log(pval_loc_mot);
% % AL(pval_loc_mot==0)=2*nanmax(-log(pval_loc_mot(pval_loc_mot>0)));
% figure, plot2d;
% set(gca,'Alim',[-log(0.05) nanmax(AL(:))]);
% ploterror(z-z_reg,F1,F2,3);
% 

    function [pval,N0,pmean,pmean_ctl]=pval2d(p,p_ctl)
        Mx=numel(loc_edges{1})-1;
        My=numel(loc_edges{2})-1;
        N=1e3;
        N0=zeros(Mx,My);
        pval=nan(Mx,My);
        pmean=zeros(Mx,My);

        for i=1:Mx
            for j=1:My
                idx=find(inrange(loc(:,1),[loc_edges{1}(i) loc_edges{1}(i+1)]) &...
                         inrange(loc(:,2),[loc_edges{2}(j) loc_edges{2}(j+1)]));
                if(numel(idx)>0 & sum(~isnan(p(idx))) & sum(~isnan(p_ctl(idx))))
%                     pval(i,j)=boottest(p(idx),p_ctl(idx),N);
                    pval(i,j)=signrank(p(idx),p_ctl(idx),'tail','right');
                    pmean(i,j)=nanmean(p(idx)./p_ctl(idx));
                    pmean_ctl(i,j)=nanmean(p_ctl(idx));
                end
                N0(i,j)=numel(idx);
            end
        end        
    end

    function plot2d()
        if(isfield(params,'image'))
            sX=size(params.image,2);
            sY=size(params.image,1);
            ix=([1:sX]-sX/2)*pxl2mm;
            iy=(-[1:sY]+sY/3)*pxl2mm;
            imagesc(ix,iy,params.image);
        end

        hold on;
        S=surf(loc_bins{2},loc_bins{1},ones(size(Mz)),...
            'CData',Mz,...
            'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
            ,'AlphaData',AL,'AlphaDataMapping','scaled');
        if(isfield(params,'clim'))
            set(gca,'CLim',params.clim);
        end
        view(0,90);
        set(gca,'FontSize',params.fontsize)
        set(gca,'Xlim',loc_bins{1}([1 end]),'Ylim',loc_bins{2}([1 end]));
        set(gca,'YDir','normal');
        axis('image');
        hold off;
        set(gcf,'Color',fcolor);
        set(gca,'Color',fcolor,'XColor',fcolor,'YColor',fcolor);
%         colorbar;
    end

    function ploterror(e,F1,F2,k)
        if(k==0)
            C=0.5*[1 1 1];
        else
            C=COL(k,:);
        end
        e_edges=linspace(-2,2,1e2);
        e_bins=edge2bin(e_edges);
        h=histcounts(e,e_edges,'Normalization','probability');
        figure(F1);
        H=plot(e_bins,h);
        set(H,'Color',C,'LineWidth',3);
        hold on;
        set(gca,'FontSize',params.fontsize,'YLim',[0 0.07])        
        set(gcf,'Color',fcolor);
        set(gca,'Color',fcolor,'XColor',1-fcolor,'YColor',1-fcolor);
        
        evar=nanvar(e);
        
        lags=[0:50];
        for i=1:numel(lags)
            xc(i)=nancorr(e(1:end-lags(i)),e(1+lags(i):end));
        end
        figure(F2);
        H=plot(lags,xc*evar);
        set(H,'Color',C,'LineWidth',3,'MarkerSize',params.msize);        
        hold on;
        set(gca,'FontSize',params.fontsize,'Ylim',[0 0.5])        
        set(gcf,'Color',fcolor);
        set(gca,'Color',fcolor,'XColor',1-fcolor,'YColor',1-fcolor);
    end

    function plotconv(F3,ratio,k)
        if(k==0)
            C=0.5*[1 1 1];
        else
            C=COL(k,:);
        end        
        figure(F3);
        H=plot(ratio,'-');
        set(H,'Color',C,'LineWidth',3);        
        hold on;
        set(gca,'FontSize',params.fontsize)        
        set(gcf,'Color',fcolor);
        set(gca,'Color',fcolor,'XColor',1-fcolor,'YColor',1-fcolor);        
    end

end

