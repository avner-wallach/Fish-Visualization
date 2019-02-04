% script for visualizing bayesian decoding of object position
fpath='Z:\mormyrid_data';
samplerate=30e3;
setenv('DATAPATH',fpath);
setenv('SAMPLERATE','30000');
setenv('FRAMERATE','50');
setenv('PXLSIZE','0.423'); %pxl to mm conversion
x_lim=[-250 250];   %default x bin limits for headfixed view
y_lim=[-400 200];   %default x bin limits for headfixed view

% load('Z:\mormyrid_data\analyzed data\20180829_alldata.mat');
dates={20180829,20180829,20180830,[20180830 20180831] 20180831 20180901 20180902};
segs={[1 17],[18 50],[0 18],[19 25;0 10],[14 36],[0 44],[0 20]};
K_obj=[2 4];
I=16; J=2;
% get probability maps
[ pxy,pz,pz_xy,bins,pm,ps]  = get_bayes_estimation(peod(K_obj),file(K_obj),I,'headfix','on');
% pxy_z=pz_xy.*repmat(pxy,1,1,size(pz_xy,3));
L=pz_xy./repmat(permute(pz,[3,1,2]),size(pz_xy,1),size(pz_xy,2),1);% likelihood
%% get video file and data
k=1;
D=dates{K_obj(k)};
seg=segs{K_obj(k)};
d=1;
fname=[fpath,'\',num2str(D(d)),'\'];
filenum=[18,23:50];%[seg(d,1):seg(d,2)];
f=2;
setenv('SESSDATE',num2str(D(d)));

% vid = VideoReader([fname,'video_',num2str(filenum(f)),'.avi']);
datafile=[fname,'data_',num2str(filenum(f))];
% load(datafile);        
f_ind=[pframe(K_obj(k)).ind0(f):pframe(K_obj(k)).ind0(f+1)-1];
fdata=pframe(K_obj(k)).data(f_ind,:);
ft=pframe(K_obj(k)).t(f_ind,:);
eind=[peod(K_obj(k)).ind0(f):peod(K_obj(k)).ind0(f+1)-1];
edata=peod(K_obj(k)).data(eind,:);
et=peod(K_obj(k)).t(eind,:);
%% get segments with objects
[fx,fy]=allo2ego(data.FILE.objects.x,data.FILE.objects.y,fdata(:,3),fdata(:,1),fdata(:,2)); %obj coordinates 
[ex,ey]=allo2ego(file(K_obj(k)).objects.x,file(K_obj(k)).objects.y,edata(:,4),edata(:,2),edata(:,3)); %obj coordinates- eod 
% ind=find(inrange(fx,x_lim) & inrange(fy,y_lim));
 ind=find(sqrt(fx.^2 + fy.^2)<200);
seg_edges = find_discont(ind); %segments of object passing
s=63;
II=seg_edges(s,1):seg_edges(s,2);
% [frames,ax_grid]=get_fish_image(filenum(f),II,'objects','on','headfix','on','vid',vid,'data',data);

%compute changes of fish position between frames
dx=[0;diff(fdata(II,1))];
dy=[0;diff(fdata(II,2))];
a=[0;fdata(II(1:end-1),3)];
[dx,dy]=allo2ego(dx,dy,a,0,0);
da=[0;diff(fdata(II,3))];

% find frames with eod
eind=find(inrange(et,ft(seg_edges(s,:))));
for e=1:numel(eind)
    [m,frind(e)]=min(abs(ft(II)-et(eind(e)))); %find frame closest to eod
    ii=find(data.EOD.t==et(eind(e)));
    Tr(e,:)=data.EOD.traces(ii,:,J);
end
Tr=my_detrend(Tr,420);
Trx=[1:size(Tr,2)]/samplerate*1e3;

Pz_xy=my_interp3(bins{2},bins{1},bins{3},L,edata(eind,I));
% Pxy=interp2(bins{2},bins{1},pxy,ax_grid{2},ax_grid{1}','spine');
Pxy=pxy;
%% play decoder
posterior=1;
recursive=1;

F1=figure;
set(F1,'Color',[0 0 0]);
A1=axes('XTick',[],'YTick',[]);
A2=axes('Position',[ 0.7946    0.7029    0.1461    0.2019],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0]);
hold on;
H=[];
Mz=ones(numel(bins{1}),numel(bins{2}));
flag=1;
for j=1:numel(frames)
    axes(A1);
    cla(A1);
    hold on;
    imagesc(ax_grid{1},ax_grid{2},flipud(frames(j).cdata));
    view(0,90);    
    Mz=rotate_map(Mz,-dx(j),-dy(j),-da(j),bins,ax_grid);    
    jj=find(frind==j);
    if(numel(jj) & flag)
        for q=1:numel(jj)
            if(recursive)
                Mz=Pz_xy(:,:,jj(q)).*Mz;
            else
                Mz=Pz_xy(:,:,jj(q));
            end
        end
%         flag=0;
    elseif(~recursive)
        Mz(:)=0;       
    end
    A0=(Mz);    
    if(posterior)
        A0=A0.*Pxy;
    end
    A0=A0/sum(A0(:)); %normalize
    S=surf(bins{1},bins{2},ones(size(A0')),...
        'CData',A0',...
        'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
        ,'AlphaData',A0','AlphaDataMapping','scaled');
%     set(gca,'Clim',[0 0.
    if(max(A0(:))>0)
        set(gca,'ALim',[0 quantile(A0(:),0.93)]);
    else
        set(gca,'ALim',[0 1]);
    end

    set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',ax_grid{1}([1 end]),'Ylim',ax_grid{2}([1 end]));
    hold off;
    axes(A2);
    if(numel(H))
        set(H,'LineWidth',1,'Color',0.5*[1 1 1]);
    end
    if(ismember(j,frind))
        H=plot(Trx,Tr(find(frind==j),:));
        set(H,'LineWidth',3,'Color',[1 1 1]);
    end
    set(A2,'Xlim',[min(Trx) max(Trx)],'Ylim',[min(Tr(:)) max(Tr(:))]);
    set(A1,'XTick',[],'YTick',[]);
    frames1(j)=getframe(F1);
end

%% run with surrogate data
% generate serogate data
posterior=1;
recursive=1;
N=1;

Pm=interp2(bins{1},bins{2},pm,fx(II),fy(II),'spline');
Ps=interp2(bins{1},bins{2},ps,fx(II),fy(II),'spline');
for i=1:numel(Pm)
    sdata(i,1:N)=randn(1,N)*Ps(i)+Pm(i);
end
% Sdata=mean(sdata,2);
for n=1:N
    Pz_xy_s(:,:,:,n)=interp3(bins{2},bins{1},bins{3},L,bins{2},bins{1},sdata(:,n),'spline');
end
Pz_xy_m(:,:,:)=interp3(bins{2},bins{1},bins{3},L,bins{2},bins{1},Pm,'spline');
frind_s=1:numel(sdata);

F1=figure;
set(F1,'Color',[0 0 0]);
A1=axes('XTick',[],'YTick',[]);
A2=axes('Position',[ 0.7946    0.7029    0.1461    0.2019],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0]);
hold on;
H=[];
flag=1;
Mz=ones(numel(bins{1}),numel(bins{2}));
for j=1:numel(frames)
    axes(A1);
    cla(A1);
    hold on;
    imagesc(ax_grid{1},ax_grid{2},flipud(frames(j).cdata));
    view(0,90);    
    Mz=rotate_map(Mz,-dx(j),-dy(j),-da(j),bins(1:2),ax_grid(1:2));    
    jj=find(frind_s==j);
    
    if(numel(jj) & flag)
        if(~recursive)
            Mz=ones(numel(bins{1}),numel(bins{2}));
        end
        for q=1:5
            Mz=Pz_xy_s(:,:,jj,q).*Mz;
%         plot(ex(eind(find(frind==j))),ey(eind(find(frind==j))),'*w');
        end
        
%         Mz=Pz_xy_m(:,:,jj).*Mz;
%         flag=1;
%     else
    elseif(~recursive)
        Mz(:)=0;
    end
    if(posterior)
        A0=Mz.*pxy;   
    else
        A0=Mz;
    end
    A0=A0/sum(A0(:)); %normalize    
%     A0=real((A0)^(1));
    S=surf(bins{1},bins{2},ones(size(A0')),...
        'CData',A0',...
        'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
        ,'AlphaData',A0','AlphaDataMapping','scaled');
%     set(gca,'Clim',[0 0.
    a=[0.25 0.99];
    a1=quantile(A0(:),a)
    if(max(a1)>0)
        set(gca,'ALim',[0 a1(2)]);
    else
        set(gca,'ALim',[0 1]);
    end

    set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',ax_grid{1}([1 end]),'Ylim',ax_grid{2}([1 end]));
    hold off;
%     axes(A2);
%     if(numel(H))
%         set(H,'LineWidth',1,'Color',0.5*[1 1 1]);
%     end
%     if(ismember(j,frind))
%         H=plot(Trx,Tr(find(frind==j),:));
%         set(H,'LineWidth',3,'Color',[1 1 1]);
%     end
%     set(A2,'Xlim',[min(Trx) max(Trx)],'Ylim',[min(Tr(:)) max(Tr(:))]);
    set(A1,'XTick',[],'YTick',[]);
    frames2(j)=getframe(F1);
end
%% show likelyhood and posterior maps
posterior=0;
N=100
sdata=linspace(min(edata(eind,I)),max(edata(eind,I)),N);
for n=1:N
    Pz_xy_s(:,:,n)=interp3(bins{2},bins{1},bins{3},L,bins{2},bins{1},sdata(n),'spline');
end
frind_s=1:numel(sdata);

F1=figure;
set(F1,'Color',[0 0 0]);
A1=axes('XTick',[],'YTick',[]);
A2=axes('Position',[ 0.7946    0.7029    0.1461    0.2019],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0]);
hold on;
H=[];
flag=1;
Mz=ones(numel(bins{1}),numel(bins{2}));
for j=1:N
    axes(A1);
    cla(A1);
    hold on;
    imagesc(ax_grid{1},ax_grid{2},flipud(frames(1).cdata));
    view(0,90);    
    Mz=Pz_xy_s(:,:,j);
    if(posterior)
        A0=Mz.*pxy;   
    else
        A0=Mz;
    end
    A0=A0/sum(A0(:)); %normalize    
%     A0=real((A0)^(1));
    S=surf(bins{1},bins{2},ones(size(A0')),...
        'CData',A0',...
        'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
        ,'AlphaData',A0','AlphaDataMapping','scaled');
    a=[0.25 0.99];
    a1=quantile(A0(:),a);
    if(max(a1)>0)
        set(gca,'ALim',[0 a1(2)]);
    else
        set(gca,'ALim',[0 1]);
    end

    set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',ax_grid{1}([1 end]),'Ylim',ax_grid{2}([1 end]));
    hold off;
    set(A1,'XTick',[],'YTick',[]);
    frames3(j)=getframe(F1);
end
