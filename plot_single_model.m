function [ Mz ] = plot_single_model( varargin )
%GET_SINGLE_MODEL plot tuning of eod 

pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.bg='on';
params.fcolor=[0 0 0];
% params.r_lim=[0 200]; %radius limits
% params.th_lim=[-pi pi]; %azimouth limits
params.rx_nbins=200;  %default r/x bin number
params.thy_nbins=201;  %default th/y bin number
params.fontsize=16;
params.t_col=[];            %condition columns
params.t_val=[];
params.model='loc_mot';
params.mode='cart';
params.maxa=0.25;
% params.tank='circ'; %rectangular or circular tank
%% get varins
n=numel(varargin);
if(n<1 | n>1&~mod(n,2))
    error('too few input arguments!');
end
net_struct=varargin{1};
if(n>1)
    i=2;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end

%% get regression model
net=net_struct.(params.model).net;

if(strcmp(params.mode,'polar'))
    r_edges=linspace(net_struct.loc.min(1),net_struct.loc.max(1),params.rx_nbins+1);
    th_edges=linspace(-pi,pi,params.thy_nbins+1);
    r_bins=edge2bin(r_edges);
    th_bins=edge2bin(th_edges);
    [R,TH]=meshgrid(r_bins,th_bins);
    loc=[R(:) TH(:)];
else %cartesian
    x_edges=linspace(net_struct.loc.min(1),net_struct.loc.max(1),params.rx_nbins+1);
    y_edges=linspace(net_struct.loc.min(2),net_struct.loc.max(2),params.thy_nbins+1);    
    x_bins=edge2bin(x_edges);
    y_bins=edge2bin(y_edges);    
    [X,Y]=meshgrid(x_bins,y_bins);
    loc=[X(:) Y(:)];
end

prior=(reshape(net_struct.prior.net(loc'),params.thy_nbins,params.rx_nbins)');
prior=prior/abs(nanmax(prior(:)));
% prior(prior<0)=0;
P=repmat(net_struct.mot.mean,size(loc,1),1);
P(:,params.t_col)=params.t_val;
z=net([loc P]');    
Mz=reshape(z,params.thy_nbins,params.rx_nbins)';
% Mz=prior;
if(strcmp(params.mode,'polar'))
    plotsurf_polar();
else
    plotsurf_cart();
end
%% output

%%
    function plotsurf_polar()
        AL=min(prior,quantile(prior(:),params.maxa));
        Mz=[Mz Mz(:,1)];
        AL=[AL AL(:,1)];
        F=gcf;
        A=gca;
        if(isfield(params,'image'))
            sX=size(params.image,2);
            sY=size(params.image,1);
            ix=([1:sX]-sX/2)*pxl2mm;
            iy=(-[1:sY]+sY/3)*pxl2mm;
            imagesc(ix,iy,params.image);
            hold on;
        end

        S=polarplot3d(Mz,AL,'AngularRange',[-pi pi],'plottype','surfa','RadialRange',r_bins);        
        if(isfield(params,'clim'))
            set(A,'CLim',params.clim);
        end
        % set(A,'ALim',[quantile(N0(:),params.mina) quantile(N0(:),params.maxa)]);
        view(0,90);
        set(gca,'FontSize',params.fontsize)
        set(gca,'XLim',net_struct.loc.max(1)*[-1 1],'YLim',net_struct.loc.max(1)*[-1 1]);
        set(gca,'Ydir','normal');
        axis('image');
        hold off;
        set(gcf,'Color',params.fcolor);
        set(gca,'Color',params.fcolor,'XColor',[0 0 0],'YColor',[0 0 0]);
%         colorbar;
    end

    function plotsurf_cart()
        F=gcf;
        A=gca;
        if(isfield(params,'image'))
            sX=size(params.image,2);
            sY=size(params.image,1);
            ix=([1:sX]-sX/2)*pxl2mm;
            iy=(-[1:sY]+sY/3)*pxl2mm;
            imagesc(ix,iy,params.image);
            hold on;
        end

%         AL=min(prior,quantile(prior(:),params.maxa));
        AL=prior;
        hold on;
        S=surf(y_bins,x_bins,ones(size(Mz)),...
            'CData',Mz,...
            'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
            ,'AlphaData',AL,'AlphaDataMapping','scaled');
        if(isfield(params,'clim'))
            set(A,'CLim',params.clim);
        end
        view(0,90);
        set(gca,'FontSize',params.fontsize)
        set(gca,'XLim',net_struct.loc.max(1)*[-1 1],'YLim',net_struct.loc.max(1)*[-1 1]);
        set(gca,'Ydir','normal');
        axis('image');
        hold off;
        set(gcf,'Color',params.fcolor);
        set(gca,'Color',params.fcolor,'XColor',[0 0 0],'YColor',[0 0 0]);

    end
end

