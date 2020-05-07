function [ Mz ] = plot_polar_model( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=str2num(getenv('PXLSIZE'));
%% params
params.bg='on';
params.fcolor=[0 0 0];
% params.cornersize=180;
params.r_lim=[0 200]; %radius limits
params.th_lim=[-pi pi]; %azimouth limits
params.r_nbins=40;  %default r bin number
params.th_nbins=41;  %default th bin number
% params.minsamp=50;
params.fontsize=16;
% params.freqmode='on';
% params.p_col=[1 5:14];     %posture columns
params.t_col=[];            %condition columns
params.t_val=[];
params.model='loc_mot';
params.mode='cart';
% params.tank='circ'; %rectangular or circular tank
%% get varins
n=numel(varargin);
if(n<1 | n>1&~mod(n,2))
    error('too few input arguments!');
end
net_struct=varargin{1};
% file_struct=varargin{2};
% z_col=varargin{3};  %output variable
if(n>1)
    i=2;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end

%% get regression model
r_edges=linspace(params.r_lim(1),params.r_lim(2),params.r_nbins+1);
th_edges=linspace(-pi,pi,params.th_nbins+1);
r_bins=edge2bin(r_edges);
th_bins=edge2bin(th_edges);
% r_bins=linspace(net_struct.loc.min(1),net_struct.loc.max(1),params.r_nbins);
% th_bins=linspace(net_struct.loc.min(2),net_struct.loc.max(2),params.th_nbins);
[R,TH]=meshgrid(r_bins,th_bins);
if(strcmp(params.mode,'polar'))
    loc=[R(:) TH(:)];
else
    loc=[R(:).*cos(TH(:)) R(:).*sin(TH(:))];
end
net=net_struct.(['model_',params.model]).net;
prior=exp(reshape(net_struct.prior.net(loc'),params.th_nbins,params.r_nbins)');
% prior(prior<0)=0;
if(~strcmp(params.model,'loc_mot'))
    z=net(loc');    
else
    P=repmat(net_struct.mot.mean,numel(R(:)),1);
    P(:,params.t_col)=params.t_val;
    z=net([loc P]');    
end
Mz=reshape(z,params.th_nbins,params.r_nbins)';
% Mz=prior;
plotsurf();
%% output

%%
    function plotsurf()
        AL=min(prior,quantile(prior(:),0.25));
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

end

