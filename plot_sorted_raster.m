function [ha,hb] = plot_sorted_raster( raster,svar,ops,winds,lines,Ns,binsize,cmap,bgc,axes_in,fit_mode)
if(nargin<10)
    figure;
    axes_in=axes;
end
if(nargin<11)
    fit_mode='slope';
end
%bgc=[0 0 0];
% bgc=[1 1 1];
fontsize=7;
ksmooth=10;
% msize=1;
msize=1.67;
nline=.5;
wline=2;
if(nargin<7)
    binsize=10;
end
zscored=0;
slope_only=1;
swin=ops.scwin*1e3; %spike window in ms
vbin=0.05;
% swin=[0 0.04]*1e3;
% swin=[0 30];
if(isfield(ops,'eodblankpre'))
    bpre=ops.eodblankpre;
    bpost=ops.eodblankpost;
else
    bpre=ops.blankpre;
    bpost=ops.blankpost;
end

if(nargin<8 | numel(cmap)==0)
    col=brewermap(10,'Set1');
    for i=1:10
        cmap{i}=colormap_graded(col(i,:));
        if(bgc==[0 0 0])
            cmap(i)=flipud(cmap(i));
        end
    end    
end
    
DR={};
h={};
C=[];
M=0;
X=cell(numel(Ns),1);
Y=cell(numel(Ns),1);
G=[]; g=1;
% V=[];
Sc=cell(numel(Ns),max(Ns));
Vc=cell(numel(Ns),max(Ns));
S=cell(numel(Ns),max(Ns));
V=cell(numel(Ns),max(Ns));

edges=[];
bins=[];
%% get SC vector
edges=[-ops.rasterpre*1e3:binsize:(-bpre*1e3)  (bpost*1e3):binsize:ops.rasterpost*1e3];
bins=edge2bin(edges); 
indt=find(inrange(raster(:,1),swin/1e3));                        
sc_all=histcounts(raster(indt,2),bin2edge([1:numel(svar)]));
sc_m=nanmean(sc_all);
sc_s=nanstd(sc_all);

%%
if(iscell(winds))
   for j=1:numel(winds)
       N=Ns(j);
       wind=winds{j};
       COL=cmap{j};
       get_rasters_psth;
       Clines{j}=COL;
%        cmap{j}=COL;
       C=[C;COL(1:N,:)];
   end
else
    j=1;
    N=Ns;
    wind=winds;
    COL=colormap(gca,cmap{1});
    get_rasters_psth;
    Clines{1}=COL;
    
    C=[C;COL(1:N,:)];
end
indblank=find(inrange(bins,[-bpre bpost]*1e3));
%% plot rasters and psth
d=0.025;
xlim=[-20 50];
[a1,p1]=tight_subplot(1, 3, [.0 .04],[.0 .00],[.00 .0],[],[.4 .3 .3],axes_in); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
a1(2).Position(1)=p1{2}(1)-0.025;
[a11,p11]=tight_subplot(numel(Ns), 1, [.02 .0],[.0 .00],[.00 .0],[],[],a1(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[a12,p12]=tight_subplot(numel(Ns), 1, [.02 .0],[.0 .00],[.00 .0],[],[],a1(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
ha=[a11;a12;a1(3)];

% [ha,pa]=tight_subplot(numel(Ns),2,[d d],[0 0],[0 0],[],[],axes_in); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
hm=45;

for i=1:numel(Ns)
    [hb(i,:),pb]=tight_subplot(2,1,[0 0],[0 0],[0 0],[.35 .65],[1],a11(i));
    axes(hb(i,1));
    hb(i,1).Tag='Raster';
    if(numel(DR{i})>0)
        scatter(DR{i}(:,1),DR{i}(:,2),msize,Vc{i}(DR{i}(:,2)),'filled');
    end
    M=lines;
    hold on;
    A=area([-bpre bpost]*1e3,[M+1 M+1],'EdgeColor','none','FaceColor',.8*[1 1 1],'FaceAlpha',1,'ShowBaseLine','off');
    set(gca,'Xlim',xlim,'Ylim',[0 M+1]);   
    set(gca,'Color',bgc,'XColor','none','YColor',1-bgc,'FontSize',fontsize,'XTick',[],'YTick',[]);    
    set(gca,'Colormap',cmap{i});    
    
    axes(hb(i,2));
    hb(i,2).Tag='PSTH';
    for j=1:size(h{i},1)
        hh=smooth(h{i}(j,:),ksmooth);
        hh(indblank)=nan;
        H=plot(bins,hh);
        hm=max(hm,nanmax(hh));
        set(H,'Color',Clines{i}(j,:),'LineWidth',nline);
        hold on;        
    end
    if(i==numel(Ns))
        x1=30;
        y1=25;
        H=plot(x1+[0 20],y1+[0 0],'k',x1+[0 0],y1+[0 20],'k');
        set(H,'LineWidth',nline);        
    end

    A=area([-bpre bpost]*1e3,[1e3 1e3],'EdgeColor','none','FaceColor',.8*[1 1 1]);        
    set(gca,'Color','none','XColor',1-bgc,'YColor',1-bgc,'FontSize',fontsize,'box','off');
%     set(gca,'Xlim',xlim,'Ylim',[0 hm*1.1],'YTick',[]);    
end
set(hb(:,2),'Xlim',xlim,'Ylim',[0 hm*1.2],'YTick',[]);
set(hb(1:end,2),'XTick',[]);
%% plot input-output

vedges=[-2:vbin:2];
vbins=edge2bin(vedges);
for j=1:numel(Ns)
    axes(a12(j));
    a12(j).Tag='Regression';
    for k=1:numel(vedges)-1
        scs(k)=nanmean(S{j}(inrange(V{j},vedges([k:k+1]))));
    end
    H=scatter(vbins,scs,msize,vbins,'filled');
%     H=scatter(Vc{j},smooth(Vc{j},Sc{j},10),msize,Vc{j},'filled');
    set(H,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);
    hold on;
    f=fit(V{j},S{j},'poly1');
    if(slope_only)
        
        p = predint(f,vedges,0.95,'functional','on');
        if(strcmp(fit_mode,'slope'))
            clim=[-.5 .5];
            if(bgc==[0 0 0])
                C1=colormap(gca,invert_map(flipud(brewermap(64,'RdBu'))))
            else
                C1=colormap(gca,flipud(brewermap(64,'RdBU')));
            end            
            cval=min(max((f.p1-clim(1))/diff(clim),0),1);    
            MCOL=interp1(linspace(0,1,size(C1,1)),C1,cval);            
            A1=patch([vedges fliplr(vedges)],[p(:,2)' fliplr(p(:,1)')],MCOL,'LineStyle','none','FaceColor',MCOL,'FaceAlpha',0.5);        
            H=plot(vedges,p(:,1)); set(H,'Color',0.5*[1 1 1],'LineWidth',nline); legend('off');
            H=plot(vedges,p(:,2)); set(H,'Color',0.5*[1 1 1],'LineWidth',nline); legend('off');
            H=plot(vedges,f(vedges)); set(H,'Color',MCOL,'LineWidth',wline); legend('off');        
            set(gca,'Clim',[-2.5 2.5],'XAxisLocation','origin','YAxisLocation','origin',...
                'Xlim',[-2.5 2.5],'Ylim',[-1.75 1.75],'YTick',[],'XTick',[],'FontSize',fontsize);    
            
        else
            clim=[0 120];
            C1=colormap(gca,'parula');
            cval=min(max((f.p2-clim(1))/diff(clim),0),1);    
            MCOL=interp1(linspace(0,1,size(C1,1)),C1,cval);            
            H=plot(vedges,f(vedges)); set(H,'Color',[0 0 0],'LineWidth',wline); legend('off');        
            H=plot(0,f.p2,'o');                        
            set(H,'MarkerFaceColor',MCOL,'MarkerSize',msize*2,'MarkerEdgeColor',1-bgc);
            set(gca,'Clim',[-2.5 2.5],'XAxisLocation','origin','YAxisLocation','origin',...
                'Xlim',[-2.5 2.5],'Ylim',clim*1.25,'YTick',[],'XTick',[],'FontSize',fontsize);    
            
        end
        
        ci(1:2,:,j)=confint(f,0.95);
        ci(3,1,j)=f.p1; ci(3,2,j)=f.p2;
        mc(j,:)=MCOL;        

    else
        C1=colormap(gca,'parula');        
        Hf=plot(V{j},f(V{j}));
        set(Hf,'Color',0.5*[1 1 1],'LineWidth',wline);
        H=plot(0,f(0),'o');
        cval=min(max((f(0)-clim(1))/diff(clim),0),1);    
        MCOL=interp1(linspace(0,1,size(C1,1)),C1,cval);
        set(H,'MarkerFaceColor',MCOL,'MarkerSize',msize*2,'MarkerEdgeColor',1-bgc);
        set(gca,'Clim',[-2.5 2.5],'XAxisLocation','origin','YAxisLocation','origin',...
            'Xlim',[-2.5 2.5],'Ylim',clim*1.25,'YTick',[],'XTick',[],'FontSize',fontsize);    
        
    end
    if(j==numel(Ns))
        x1=-2;
        y1=1.5;
        H=plot(x1+[0 1],y1+[0 0],'k',x1+[0 0],y1+[0 1],'k');
        set(H,'LineWidth',nline);        
%         set(gca,'XTick',[-2 -1 0 1 2]);
    end
    set(gca,'Color',bgc,'YColor',1-bgc,'XColor',1-bgc,'Colormap',cmap{j});    
end

axes(ha(end));
bw=.2;
if(strcmp(fit_mode,'slope'))
    j=1;
else
    j=2;
end
for i=1:size(mc,1)
    H=rectangle('Position',[i-bw/2 ci(1,j,i) bw diff(ci(1:2,j,i))]);
    set(H,'FaceColor',[mc(i,:),0.5],'EdgeColor',0.5*[1 1 1]);
    hold on;
    H=plot(i,ci(3,j,i),'o');
    set(H,'MarkerEdgeColor',1-bgc,'MarkerFaceColor',mc(i,:),'MarkerSize',5);        
end
set(gca,'Color','none','YColor',1-bgc,'XColor',1-bgc,'Xlim',[.5 i+.5],'Ylim',clim,'YTickLabelMode','auto','XTick',[],'XAxisLocation','origin','Box','off','FontSize',fontsize);    

    function get_rasters_psth
        if(numel(wind)==2)
            ind=find(inrange(raster(:,2),wind));
            rast=raster(ind,:); 
            rast(:,2)=rast(:,2)-wind(1);
            if(wind(1)<1)
                wind(1)=1;
            end
            if(isinf(wind(2)) | wind(2)>numel(svar))
                wind(2)=numel(svar);
            end
            sv=svar(inrange([1:numel(svar)],wind));
%             sc1=sc_all(inrange([1:numel(svar)],wind));
        elseif(numel(wind)>2)
            ind=find(ismember(raster(:,2),wind));
            q=zeros(numel(svar),1);
            q(wind)=[1:numel(wind)];
            rast=raster(ind,:); 
            rast(:,2)=q(rast(:,2));        
            sv=svar(wind);    
%             sc1=sc_all(wind);
        else
            rast=raster;
            sv=svar;
%             sc1=sc_all;
        end
        [rast,v]=sort_raster(rast,sv);

        %convert to ms
        rast(:,1)=rast(:,1)*1e3;
        
        %% PSTH
        redges=round(linspace(1,numel(v)+1,N+1));
%         edges=[-ops.rasterpre*1e3:binsize:ops.rasterpost*1e3];
        edges=[-ops.rasterpre*1e3:binsize:(-bpre*1e3)  (bpost*1e3):binsize:ops.rasterpost*1e3];
        bins=edge2bin(edges);        
        bind=find(inrange(bins,swin));
        L=floor(lines/N);
        K=0;
        Sc{j}=[];
        Vc{j}=[];
        DR{j}=[];
        S{j}=[];
        V{j}=[];
        h{j}=[];
        for i=1:N
            ind=find(rast(:,2)>=redges(i) & rast(:,2)<redges(i+1));
            h{j}=[h{j};histcounts(rast(ind,1),edges)/diff(redges(i:i+1))/binsize*1e3];
            indt=find(inrange(rast(:,1),swin));                        
            vind=[redges(i):redges(i+1)-1]';
            dv=v(vind);
%             dsc=sc1(vind);
%             sc=histcounts(rast(indt,2),bin2edge([redges(i):redges(i+1)]));             
%             Sc{j}=[Sc{j};sc'];

            
            vvind=unique(round(linspace(vind(1),vind(end),L)));
            [drast]=dilute_rast(rast(ind,:),vvind);
            drast(:,2)=drast(:,2)+K;
            K=K+numel(vvind);
            DR{j}=[DR{j};drast];    
            Vc{j}=[Vc{j};v(vvind)];                        

            iind=ind(inrange(rast(ind,1),swin));
            sc=histcounts(rast(iind,2),bin2edge(vind')); 
            if(zscored)
                sc=(sc-sc_m)/sc_s;
            else
                sc=sc*40;
            end
            
            S{j}=[S{j};sc'];
            V{j}=[V{j};v(vind)];
        end
        
        if(N==1)
            COL(1,:)=COL(32,:);
        else
            COL(1:N,:)=interp1([1:size(COL,1)],COL,linspace(1,size(COL,1),N));
        end
%         M=M+numel(r);
    end

    function [drast]=dilute_rast(rrast,indin)
%         r=unique(round(linspace(n1,n2,lines)));
        ind=find(ismember(rrast(:,2),indin));
        q=zeros(indin(end)-indin(1)+1,1);
        q(indin)=[1:numel(indin)];
        drast=rrast(ind,:); 
        drast(:,2)=q(drast(:,2));
%         M=M+lines;
    end

end

