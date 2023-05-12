% adir='20190131'; seg=2; fnum=15; ss=300; spkgroup=1; chid=4; rastind=[2:4];
% adir='20200123'; seg=3; fnum=15; ss=1500; spkgroup=1; chid=2; rastind=[2:4];
adir='20190624'; seg=6; fnum=9; ss=13760; spkgroup=1; chid=1; rastind=[1:2];
% script for visualizing freely moving recordings
apath='Z:\analysis_data';
opsload=1;
if(opsload)
    load([apath,filesep,adir,filesep,'output'],'ops','eod','file');
    load([apath,filesep,adir,filesep,adir,'_archieve'],'frame');
end
setenv('PXLSIZE',num2str(ops.pxlsize)); %pxl to mm conversion
pxl2mm=ops.pxlsize;
% netmode='none'; %loc,ni,loc_mot,fb,none
% accmode=0;
% accnum=5;
dispmode='none'; %tuning,likelihood,posterior,none
%hpf filter
lp=100*2/ops.samplerate;
[NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
[Bhp,Ahp] = butter(NN,Wn,'high');
%% general parameters 
fdate=num2str(ops.seg(seg).dates(1));
dname=[ops.datapath,'\',fdate,'\'];
fname=ops.seg(seg).files(fnum); %file to take
cvar=ops.seg(seg).LFPgroups{ops.seg(seg).Spkgroups(spkgroup)}(chid);
% setenv('SESSDATE',num2str(date));
% varname='group4';
datafile=[dname,'data_',num2str(fname)];
framenum=800; %number of frames
frameind=ss+[0:framenum];

%% get video and data files
loading=1;
if(loading)
    load(datafile);        %for traces
end

gamp=1;
if(gamp)
    chns=find(ops.outchans);
    ch=chns(cvar);
    sesspath=[ops.datapath,'\',fdate,'\'];
    [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(fname)],ops.samplerate,ops.chan_num,ops.blocksize,ch,'adc');
    t=t-data.FILE.offset;
end

gframes=1;
if(gframes)
    vidname=[dname,'video_',num2str(fname)];
    trackfile=[dname,'video_',num2str(fname),ops.trackname];
    vid = VideoReader([vidname,'.avi']);
    F=figure;
    set(F,'Color',[0 0 0],'Position',[1 41 1920 964]);
    setenv('DATAPATH',ops.datapath);
    setenv('SESSDATE',fdate);
    setenv('FRAMESCALE',num2str(ops.framescale));
    [frames,ax_grid]=visualize_fish_tracking(vidname,trackfile,ss,framenum)
%     [frames,ax_grid]=get_fish_image(fnum,frameind,'vid',vid,'data',data);
end

%% get data stats
f_ind=[frame(seg).ind0(fnum):frame(seg).ind0(fnum+1)-1]; %frame indices for file
f_ind=f_ind(frameind);
fdata=frame(seg).data(f_ind,:);    %frame data for file
ft=frame(seg).t(f_ind,:);  %frame time
eind=[eod(seg).ind0(fnum):eod(seg).ind0(fnum+1)-1];   %eod indices for file
edata=eod(seg).data(eind,:);   %eod data for file
et=eod(seg).t(eind,:);    %eod time
%% get variables
fx=fdata(:,1)*pxl2mm;  fy=fdata(:,2)*pxl2mm;  fa=fdata(:,3);

circle = file(seg).circle*pxl2mm;
rx=circle(3)/2; ry=circle(4)/2;    %ellipse radii
x0=circle(1) + rx; y0=circle(2) + ry; %center point

%get wall egocentric position r,th
phi=atan2((fy-y0),(fx-x0));    %fish azimuth in tank
dphi=[0;diff(phi)];
R1=hypot((fx-x0),(fy-y0)); %distance of fish from center
R2=(rx*ry)./sqrt((ry*cos(phi)).^2 + (rx*sin(phi)).^2); %distance of nearest point from cetner
fr=R2-R1;    %distance of fish to nearest wall
fth=phi-fa;   %egocentric angle of closest wall
fth(fth>pi)=fth(fth>pi)-2*pi;
fth(fth<-pi)=fth(fth<-pi)+2*pi;    

%get wall allocentric position
xw=x0+R2.*cos(phi);  yw=y0+R2.*sin(phi);

% find frames with eod
clear frind Tr idx;
iii=find(et>=ft(1) & et<=ft(end));
clear rind;
rast=cell(1,numel(rastind));
for i=1:numel(rastind)
    rind{i}=find(ismember(eod(seg).raster{rastind(i)}(:,2),eind(iii)) & inrange(eod(seg).raster{rastind(i)}(:,1),[0.001 .02]));
    rast{i}=[rast{i};eod(seg).raster{rastind(i)}(rind{i},1) eod(seg).raster{rastind(i)}(rind{i},2)-eind(iii(1))+1];
end

% rast{1}=sortrows([rast{1};rast{2}],2); %unify rasters
% rast{2}=rast{3};
% rast(3)=[];

e=1; frind=[]; idx=[];
while(e<=numel(iii))
    [m,i]=min(abs(ft-et(iii(e)))); %find frame closest to eod
    frind=[frind i];
    Tr(e,:)=data.EOD.traces(iii(e),:,cvar);
    e=e+1;
end
idx=iii;
edata=edata(idx,:);
et=et(idx);
ex=edata(:,2)*pxl2mm;  ey=edata(:,3)*pxl2mm;  ea=edata(:,4);

% Tr=filtfilt(Bhp,Ahp,Tr')';
% Trx=[1:size(Tr,2)]/ops.samplerate*1e3;
%% get traces
Ix=[30:20*30];
ind=find(ismember(t,et))'; %eod indices in t
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
Tr=amp(Idx);
Tr=filtfilt(Bhp,Ahp,Tr')';
Trx=Ix/ops.samplerate*1e3;
%% get rasters
Sp=nan([size(Tr) numel(rast)]);
for i=1:size(Sp,3)
    for j=1:size(rast{i},1)
        spk=rast{i}(j,:);
        is=round(spk(1)*3e4-30)+[-15:30];
        is(is>size(Sp,2))=[];
        Sp(spk(2),is,i)=Tr(spk(2),is);
    end
end
%% grid 
%allocentric
Rcirc=300;
[X,Y]=meshgrid(ax_grid{1}*pxl2mm,ax_grid{2}*pxl2mm);
mask=(((X-x0).^2/rx^2+(Y-y0).^2/ry^2)<1);
immask=repmat(uint8(((X-x0).^2+(Y-y0).^2)<Rcirc^2)*255,1,1,3);

%% play decoder

F1=figure;
COL=colormap('lines');
COL(4,:)=[236 157 118]/255; %afferent color
COL(5,:)=[64 128 0]/255;%[57 181 74]/255; output color
COL(6,:)=[210 148 37]/255;%COL(7,:);%interneuron
COL(7,:)=[86 124 141]/255;%output color
LFP_t=4.5;
rsymbol=[char(0x2B24) char(0x25AE)];
colormap('parula');
set(F1,'Color',[0 0 0],'Position',[1 41 1920 964]);
A1=axes('XTick',[],'YTick',[]);
set(A1,'Position',[-.1 0.11 0.775 0.815]);
A2=axes('Position',[ 0.52,0.6,0.35,0.3],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0]);
hold on;
% for i=1:numel(rast)
%     dx=0.4/numel(rast);
%     A3(i)=axes('Position',[ 0.584,0.2+dx*(i-1),0.2847,dx],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);    
% end

A3=axes('Position',[ 0.584,0.2,0.2847,0.4],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);
% hold on;
% A3(2)=axes('Position',[ 0.70,0.3029,0.1461,0.2019],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);
% hold on;
A4=axes('Position',[ 0.52,0.2,0.063854166666667,0.4],'XTick',[],'YTick',[],'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1]);
hold on;

N=10; %number of LFP amps to display
H=[];
r=[];
zt=[];

for j=1:numel(frames)

    axes(A1);
    cla(A1);
    hold on;
    imagesc(ax_grid{1}([1 end])*pxl2mm,ax_grid{2}([1 end])*pxl2mm,(frames(j).cdata(:,:,:)));
    imagesc(ax_grid{1}([1 end])*pxl2mm,ax_grid{2}([1 end])*pxl2mm,immask,'AlphaData',~immask(:,:,1))
    view(0,90);    
%     Hwall=plot(xw(j),yw(j),'.');
%     set(Hwall,'MarkerSize',38,'Color',COL(2,:));
    set(gca,'YDir','normal');
    axis('image');
    set(gca,'Xlim',ax_grid{1}([1 end])*pxl2mm,'Ylim',ax_grid{2}([1 end])*pxl2mm);


    %LFP traces
    axes(A2);
    Itr=find(inrange(Trx,[1 LFP_t]));
    if(ismember(j,frind) & ~strcmp(dispmode,'tuning') )
        cla(A2);
        trind=find(frind==j)+[-3:0];
        trind(trind<1)=[];
        H=plot(Trx,Tr(trind,:));
        set(H(1:end-1),'LineWidth',2,'Color',0.1*[1 1 1]);
        set(H(end),'LineWidth',3,'Color',0.75*[1 1 1]); 
        H1=plot(Trx(Itr),Tr(trind(end),Itr));        
        set(H1,'LineWidth',3,'Color',COL(4,:));
        for s=1:size(Sp,3)
            S=plot(Trx,Sp(trind(end),:,s));
            set(S,'LineWidth',3,'Color',COL(4+s,:));
        end
        
        HH=plot([15 16],200*[1 1]);
        set(HH,'LineWidth',2,'Color',[1 1 1]);
        HH=plot([15 15],200+100*[0 1]);
        set(HH,'LineWidth',2,'Color',[1 1 1]);
%         T=text(4.1,-1.1e3,'1ms');
%         set(T,'Color',[1 1 1],'FontSize',20);
%         T=text(3.7,-1e3,'.5mV');
%         set(T,'Color',[1 1 1],'FontSize',20,'Rotation',90);
    end
    set(A2,'Xlim',[min(Trx) max(Trx)],'Ylim',[min(Tr(:)) max(Tr(:))]);
    
    %raster
    axes(A3);
%     cla(A3);            
    for s=1:numel(rast)
%         axes(A3(s));
        if(ismember(j,frind))
%             cla(A3(s));
            ind=find(rast{s}(:,2)<=trind(end));
            TT=text(rast{s}(ind,1)*1e3,rast{s}(ind,2),rsymbol(s),'Color',COL(4+s,:),'FontWeight','bold');
            set(A3,'Xlim',[LFP_t max(Trx)],'Ylim',[0 size(Tr,1)+1]);        
            hold on;
        end
    end

    %lfp image
    axes(A4);
    if(ismember(j,frind))
        cla(A4);
        L=imagesc(Trx(Itr),1:trind(end),Tr(1:trind(end),Itr));
        set(A4,'Xlim',[1 LFP_t],'Ylim',[0 size(Tr,1)+1],'Clim',[-300 600]);  
        colormap('copper');
        
    end
    frames1(j)=getframe(F1);
end
