% script for visualizing bayesian decoding of wall position
fpath='Z:\mormyrid_data';
samplerate=30e3;
setenv('DATAPATH',fpath);
setenv('SAMPLERATE','30000');
setenv('FRAMERATE','50');
pxl2mm=0.423;
setenv('PXLSIZE',num2str(pxl2mm)); %pxl to mm conversion
clim=[-2 2];
alim=[];
hfix='off';
netmode='fbni'; %loc,ni,loc_mot,fb,none
accmode=0;
accnum=5;
dispmode='tuning'; %tuning,likelihood,posterior
%hpf filter
lp=100*2/samplerate;
[NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
[Bhp,Ahp] = butter(NN,Wn,'high');
%% general parameters 
date=20190131;
segs=[6:42 44 49:53];
K=1;
% date=20190201;
% segs=[6:10 16:19 25:33];
% K=3;
fname=[fpath,'\',num2str(date),'\'];
fnum=1; %file to take
setenv('SESSDATE',num2str(date));
varname='group3';
datafile=[fname,'data_',num2str(segs(fnum))];
% frameind=23180+[1:300];
frameind=200+[1:300];

%% get video and data files
loading=0;
if(loading)
    load('Z:\mormyrid_data\analyzed data\20190131_alldata.mat');
    
    vid = VideoReader([fname,'video_',num2str(segs(fnum)),'.avi']);

%     load([fname,'networks_pole_17.mat']);
    load([fname,'networks_wall_19.mat']);
    
    load(datafile);        %for traces
    % get images
end
gframes=0;
if(gframes)
    F=figure;
    set(F,'Color',[0 0 0],'Position',[1 41 1920 964]);
    [frames,ax_grid]=get_fish_image(segs(fnum),frameind,'headfix',hfix,'vid',vid,'data',data,'track','on');
end
%% select network
switch netmode
    case 'loc_mot'
        net=net_struct.model_loc_mot.net;
        sig=net_struct.model_loc_mot.sig;
    case 'loc'
        net=net_struct.model_loc.net;
        sig=net_struct.model_loc.sig;
    case 'ni'
        net=net_struct.model_niloc.net;
        sig=net_struct.model_niloc.sig;
        net_ni=net_struct.model_ni.net;
        cols=net_struct.model_ni.cols;
    case 'lfb'
        net=net_struct.model_lfbloc.net;
        sig=net_struct.model_lfbloc.sig;
        net_ni=net_struct.model_lfb.net; 
        cols=net_struct.model_lfb.cols;
    case 'fb'
        net=net_struct.model_fbloc.net;
        sig=net_struct.model_fbloc.sig;
        net_ni=net_struct.model_fb.net; 
        cols=net_struct.model_fb.cols;
    case 'fbni'
        net=net_struct.model_fbniloc.net;
        sig=net_struct.model_fbniloc.sig;
        net_ni=net_struct.model_fbni.net; 
        cols=net_struct.model_fbni.cols;
end

pnet=net_struct.prior.net;
rmin=0;
rmax=max(net_struct.loc.max);
%% get data stats
% lvar=find(cellfun(@(x) strcmp(x,varname),peod(K).fnames)); %variable in processed data to use
lvar=net_struct.z_col;
cvar=13;    %number of channel to plot

f_ind=[frame(K).ind0(fnum):frame(K).ind0(fnum+1)-1]; %frame indices for file
f_ind=f_ind(frameind);
fdata=pframe(K).data(f_ind,:);    %frame data for file
ft=frame(K).t(f_ind,:);  %frame time
eind=[peod(K).ind0(fnum):peod(K).ind0(fnum+1)-1];   %eod indices for file
% if(fnum>1)
%     eind=eind+1;    %remove after this is fixed in collect_data
% end
edata=peod(K).data(eind,:);   %eod data for file
et=peod(K).t(eind,:);    %eod time
%% get variables
fx=fdata(:,1)*pxl2mm;  fy=fdata(:,2)*pxl2mm;  fa=fdata(:,3);

%object position
xp = file(K).objects.x*pxl2mm;
yp = file(K).objects.y*pxl2mm;

circle=file(K).circle*pxl2mm;
rx=circle(3)/2; ry=circle(4)/2;    %ellipse radii
x0=circle(1) + rx; y0=circle(2) + ry; %center point

%egocentric object position
[fy1,fx1]=allo2ego(xp,yp,fa,fx,fy); 

% find frames with eod
clear frind Tr idx;
iii=find(et>=ft(1) & et<=ft(end));
e=1; frind=[]; idx=[];
while(e<=numel(iii))
    [m,i]=min(abs(ft-et(iii(e)))); %find frame closest to eod
    frind=[frind i];
%     ii=find(data.EOD.t==et(iii(e)));
%     idx=[idx ii];
    Tr(e,:)=data.EOD.traces(iii(e),1:180,cvar);
    e=e+1;
end
% frind(frind>frameind(end))=[];
% idx(end)=[];
idx=iii;
edata=edata(idx,:);
et=et(idx);
ex=edata(:,2)*pxl2mm;  ey=edata(:,3)*pxl2mm;  ea=edata(:,4);
epose=edata(:,[1 5:14]);
epose(:,1)=1./epose(:,1);
ez=edata(:,lvar);
efb=edata(:,cols);
if(ismember(1,cols))
    efb(:,cols==1)=1./efb(:,cols==1);
end

fisi=interp1(et,epose(:,1),ft,'spline');
fpose=[fisi fdata(:,4:13)];
Tr=filtfilt(Bhp,Ahp,Tr')';
% Tr=my_detrend(Tr,420);
Trx=[1:size(Tr,2)]/samplerate*1e3;

%% grid 
%allocentric
Rcirc=300;
[X,Y]=meshgrid(ax_grid{1}*pxl2mm,ax_grid{2}*pxl2mm);
mask=(((X-x0).^2/rx^2+(Y-y0).^2/ry^2)<1);
immask=repmat(uint8(((X-x0).^2+(Y-y0).^2)<Rcirc^2)*255,1,1,3);

%% play decoder

F1=figure;
COL=colormap('lines');
COL(4,:)=[236 157 118]/255;
COL(5,:)=[57 181 74]/255;
COL(6,:)=[86 124 141]/255;
colormap('parula');
set(F1,'Color',[0 0 0],'Position',[1 41 1920 964]);
A1=axes('XTick',[],'YTick',[]);
set(A1,'Position',[0.05 0.11 0.775 0.815]);
A2=axes('Position',[ 0.65    0.7029    0.1461    0.2019],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0]);
hold on;

A3=axes('Position',[ 0.65    0.1029    0.1461    0.2019],'XTick',[],'YTick',[-2 -1 1 2],'Color',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0]);
xlabel('EOD');
hold on;
N=10; %number of LFP amps to display
H=[];
z_net=[];
z_aff=[];
z_ni=[]; ni0=[];
zt=[];
if(~isinf(accnum))
    Zl=ones([size(X),accnum]);
else
    Zl=ones([size(X),1]);
end
Zp=ones(size(X));
for j=1:numel(frames)
    %rotate map
%     if(accmode)
%         for k=1:size(Zl,3)
% %         Zl(:,:,k)=rotateAround(Zl(:,:,k),x0/pxl2mm,y0/pxl2mm,-rad2deg(dphi(j)*R1(j)/R2(j)));
%         end
%     end
    
    axes(A1);
    cla(A1);
    hold on;
    imagesc(ax_grid{1}([1 end])*pxl2mm,ax_grid{2}([1 end])*pxl2mm,(frames(j).cdata));
    imagesc(ax_grid{1}([1 end])*pxl2mm,ax_grid{2}([1 end])*pxl2mm,immask,'AlphaData',~immask(:,:,1))
    view(0,90);    
    Hobj=plot(xp,yp,'*');
    set(Hobj,'MarkerSize',12,'Color',COL(2,:));

%     set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',ax_grid{1}([1 end])*pxl2mm,'Ylim',ax_grid{2}([1 end])*pxl2mm);
    

    if((ismember(j,frind) & ~strcmp(dispmode,'none'))) %|...
%             strcmp(dispmode,'tuning') & ~strcmp(netmode,'fb'))   %EOD emitted in this frame
    
        ei=find(frind==j);  %eod index
        
        %polar grid- shifted and tilted
        R=((X-fx(j)).^2+(Y-fy(j)).^2).^(.5);
        TH=atan2((Y-fy(j)),(X-fx(j)))-fa(j);    
        TH(TH>pi)=TH(TH>pi)-2*pi;
        TH(TH<-pi)=TH(TH<-pi)+2*pi; 
        X_=R.*cos(TH);   Y_=R.*sin(TH);
        loc=[X_(:) Y_(:)];

        if(strcmp(netmode,'loc_mot'))
            z=net([loc repmat(fpose(j,:),size(loc,1),1)]');
            z0=net([fx1 fy1 fpose(j,:)]');
        else
            z=net(loc');
            z0=net([fx1(j) fy1(j)]');
            if(~strcmp(netmode,'loc'))
%                 ni=net_ni(fpose(j,:)'); %negative image correction
%                 z=z+ni;
%                 ni0=ni;
%             end
%             if(strcmp(netmode,'fb'))
%                 ni=net_ni([fpose(j,:) efb(ei,:) ez(ei)]'); %negative image correction
                ni=net_ni([efb(ei,:)]'); %negative image correction
                z=z+ni;
                ni0=ni;
            end
        end
        Z=reshape(z,size(X,1),size(X,2));                            
        
        prior=pnet(loc');
        P=reshape(prior,size(X,1),size(X,2));
        P(R>=rmax | R<=rmin)=min(P(:));
%         P(~mask)=min(P(:));
        if(~strcmp(dispmode,'tuning'))
            L=1/(sqrt(2*pi*sig^2))*exp(-(ez(ei)-Z).^2/(2*sig^2));
            if(accmode)
                if(isinf(accnum))
                    Zl=Zl.*L;
                    Zl=Zl/nansum(Zl(:));
                    L=Zl;
                else
                    Zl(:,:,1:end-1)=Zl(:,:,2:end);%shift
                    Zl(:,:,end)=L;
                    L=prod(Zl,3);
                end
            end
            if(strcmp(dispmode,'posterior'))
                Zp=L.*exp(P);
            else
                Zp=L;
            end
            Zp=Zp./nansum(Zp(:)); %normalize            
            
            [mle,ii]=nanmax(Zp(:));
            mlx=X(ii); mly=Y(ii);
%             Hml=plot(mlx,mly,'+m');
            AL=Zp;
            clim=[0 1e-4];
            alim=[0 max(AL(:))/10];
        else
            Zp=Z;
            AL=P;
            clim=[-1 1.5];
            alim=quantile(P(P>min(P(:))),[0.2 1]);            

        end
        
        S=surf(ax_grid{1}*pxl2mm,ax_grid{2}*pxl2mm,ones(size(Z)),...
            'CData',Zp,'AlphaData',AL,...
            'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
            ,'AlphaDataMapping','scaled');
            set(gca,'Alim',alim,'Clim',clim);
            hold off;
    end
    set(gca,'YDir','reverse');
        
    axes(A2);
    TRM=max(Tr(:));
    TRm=min(Tr(:));
    if(ismember(j,frind) & ~strcmp(dispmode,'tuning') )
        cla(A2);
        trind=find(frind==j)+[-3:0];
        trind(trind<1)=[];
        H=plot(Trx,Tr(trind,:));
        set(H(1:end-1),'LineWidth',2,'Color',0.5*[1 1 1]);
        set(H(end),'LineWidth',3,'Color',COL(4,:));
        HH=plot([4 5],TRm*0.8*[1 1]);
        set(HH,'LineWidth',2,'Color',[1 1 1]);
        HH=plot([4 4],TRm*0.4*[1 2]);
        set(HH,'LineWidth',2,'Color',[1 1 1]);
        T=text(4.1,TRm*0.9,'1ms');
        set(T,'Color',[1 1 1],'FontSize',14);
        T=text(3.7,TRm*0.8,'.5mV');
        set(T,'Color',[1 1 1],'FontSize',14,'Rotation',90);
    end
    set(A2,'Xlim',[min(Trx) max(Trx)],'Ylim',[min(Tr(:)) max(Tr(:))]);
    
    axes(A3);
    if(ismember(j,frind) & ~strcmp(dispmode,'tuning') & ~strcmp(dispmode,'none'))
        cla(A3);
        z_net=[z_net z0];
        z_aff=[z_aff ez(ei)];
        z_ni=[z_ni ni0];
        zt=[zt ei];
        if(numel(z_net)>N)
            zt(1)=[];
            z_net(1)=[];
            z_aff(1)=[];
            if(numel(z_ni))
                z_ni(1)=[];
            end
        end
        if(numel(z_ni))
%             H1=plot(zt,-z_ni,':+');
%             set(H1,'LineWidth',2,'Color',COL(6,:));
%             hold on;
            H1=plot(zt,z_aff-z_ni,'-*');
            set(H1,'LineWidth',2,'Color',COL(4,:));
            hold on;
        else
            H1=plot(zt,z_aff,'-*');
            set(H1,'LineWidth',2,'Color',COL(4,:));
            hold on;
        end

        H1=plot(zt,z_net,':o');
        set(H1,'LineWidth',2,'Color',[1 1 1]);
        set(A3,'XLim',[zt(1) zt(1)+N],'YLim',[-3 3],...
            'XColor',[1 1 1],'YColor',[1 1 1],'XTick',zt,'XTickLabel','','TickLength',[.05 .05]);    
    end

    frames1(j)=getframe(F1);
end
