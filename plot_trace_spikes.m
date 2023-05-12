function [fmed,fpeak,flat,utype] = plot_trace_spikes(ops,segnum,fileind,lfpnum,unitnums,wind,expnum)
%PLOT_TRACE_SPIKES show voltage trace with detected spikes and LFP marked
interactive=str2num(getenv('INTERACTIVE'));

filenum=num2str(ops.seg(segnum).files(fileind));
if(numel(ops.seg(segnum).dates)>1)
    sdate=num2str(ops.seg(segnum).dates(fileind));
else
    sdate=num2str(ops.seg(segnum).dates);
end
sesspath=[ops.datapath,'\',sdate,'\'];
ind=find(ops.outchans);
outchans=zeros(size(ops.outchans));
chans=ops.seg(segnum).LFPgroups{lfpnum};
outchans(ind(chans))=1;
%get voltage trace
[t,amp]=read_bonsai_binary([sesspath,'amp_',filenum],ops.samplerate,ops.chan_num,ops.blocksize,outchans,'adc');
load([sesspath,'data_',filenum],'data');
if(nargin>=6)
    indw=find(inrange(t,wind));
else
    indw=[1:numel(t)];
end
%% trace fig
F=figure;
COL=colormap('lines');
A=axes;
H=plot(t(indw),amp(indw,:));
% set(H,'Color',COL(1,:));
hold on;

% lfp
dt=median(diff(t))*1e3;
trace_B=ceil(ops.traceblnk/dt);
trace_M=ceil(ops.tracesaved/dt);
eind=find(inrange(data.EOD.t+data.FILE.offset,wind));
eodind=round((data.EOD.t(eind)+data.FILE.offset)*ops.samplerate)+trace_B; 
I=eodind'*ones(1,trace_M) + ones(numel(eodind),1)*[0:trace_M-1]; %index matrix
for i=1:numel(chans)
    H=plot(t(I'),data.EOD.traces(eind,:,chans(i))');
    set(H,'Color',[0 0 0]);
end

% spikes in trace
indpre=15;
indpost=45;
if(isfield(data,'RAST')) %ks1 sorting
    rast=data.RAST(inrange(t(data.RAST(:,1)),wind),:);
elseif(isfield(data,'SPIKES'))%my sorting
    rast=[int64((data.SPIKES.raster(:,1)+data.FILE.offset)*ops.samplerate) data.SPIKES.raster(:,2)];
    rast=rast(inrange(t(rast(:,1)),wind),:);
%     rast(:,1)=(rast(:,1)+data.FILE.offset)*ops.samplerate;
end
rast=rast(ismember(rast(:,2),unitnums),:);
k=1;
for i=unique(rast(:,2))'
    i
    ind=rast(find(rast(:,2)==i),1);
    I=double(ind(:))*ones(1,indpre+indpost) + ones(numel(ind),1)*[-indpre:indpost-1]; %index matrix
    S(k)=plot(t(ind),zeros(size(ind)),'^');
    set(S(k),'Color',COL(2+k,:),'MarkerFaceColor',COL(2+k,:)); 
    k=k+1;
end

% save trace
% trace_path='C:\Users\Avner\iCloudDrive\Documents\Mormyrid_Data\Data\traces\'
%trace_path='C:\Users\Ephys\iCloudDrive\Documents\Mormyrid_Data\Data\traces\'
%savefig([trace_path,...
%     num2str(expnum),'_',num2str(segnum),'_',num2str(lfpnum)]);
%% prepare stats fig
F1=figure;
Nu=numel(unitnums);
xcor_ind=[];
for i=1:Nu*(Nu+4)
    A=subplot(Nu,Nu+4,i);    
    P(i,:)=get(A,'Position');    
    if(mod(i,Nu+4)>=5 | mod(i,Nu+4)==0)
        xcor_ind=[xcor_ind,i];
    end
end
p_shape=P(1:(Nu+4):end,:);
pca_xy=[P(2,1) min(P(:,2))];
pca_wh=P(3,[1 2])+P(3,[3 4])-pca_xy;
p_pca=[pca_xy pca_wh];
p_psth=P(4:(Nu+4):end,:);
p_xcor=P(xcor_ind,:);
close(F1);
F2=figure;
%% spikes stats (all file)
indpre=15;
indpost=45;
if(isfield(data,'RAST')) %ks1 sorting
    rast=data.RAST;
elseif(isfield(data,'SPIKES'))%my sorting
    rast=[int64((data.SPIKES.raster(:,1)+data.FILE.offset)*ops.samplerate) data.SPIKES.raster(:,2)];
end
rast=int32(rast(ismember(rast(:,2),unitnums),:));
k=1;
clear R;
cind=[];
for i=unitnums
    i
    ind=rast(find(rast(:,2)==i),1);
    ind=ind(ind>indpre & ind<(size(amp,1)-indpost));
    I=double(ind(:))*ones(1,indpre+indpost) + ones(numel(ind),1)*[-indpre:indpost-1]; %index matrix
    a=[];
    for j=1:size(amp,2)
        b=amp(:,j);
        c=b(I);
        c=c-median(c,2)*ones(1,size(c,2));
        c(:,end)=nan; %create dicontinuity
        a=[a c];
    end
    R{k}=a;
%     R{k}=R{k}-median(R{k},2)*ones(1,size(R{k},2));
    m(k)=quantile(R{k}(:),.01);
    M(k)=quantile(R{k}(:),.99);
    cind=[cind;k*ones(size(R{k},1),1)];
    k=k+1;
end
%% spike stats
for i=1:numel(unitnums)
    ts=[1:size(R{i},2)]/ops.samplerate*1e3;
    ind=randi(size(R{i},1),100,1);
    a_shape(i)=axes('Position',p_shape(i,:));
    H=plot(ts,R{i}(ind,:));
    set(H,'Color',lighter(COL(2+i,:),0.5));
    hold on;
    H=plot(ts,mean(R{i},1));
    set(H,'Color',COL(2+i,:),'LineWidth',5);
    set(gca,'YLim',[min(m) max(M)]*1.5,'Xlim',[ts(1) ts(end)]);
    set(gca,'FontSize',16);
    xlabel('ms');
    ylabel(num2str(unitnums(i)));
%     set(gca,'XTick',[],'YTick',[],'Box','off')
end
%% pca
RR=cell2mat(R');
RR(isnan(RR))=0;
[coefs,score,rsq,explained]=pca(RR);
a_pca=axes('Position',p_pca);
for i=1:numel(R)
    H=plot(score(cind==i,1),score(cind==i,2),'.');
    hold on;
    set(H,'Color',COL(i+2,:));
end
set(gca,'FontSize',16);
xlabel('PCA1');
ylabel('PCA2');
%% EOD triggered psth
eodind=round((data.EOD.t+data.FILE.offset)*ops.samplerate)+trace_B; 
U=unique(rast(:,2));
Nu=numel(U);
Ne=min(5e3,numel(eodind));
erast=cell(Nu,1);
for i=1:Ne
    r=rast(inrange(rast(:,1)-eodind(i),[-ops.rasterpre,ops.rasterpost]*ops.samplerate),:);
    for j=1:Nu        
        ts=(double(r(find(r(:,2)==U(j)),1))-eodind(i))/ops.samplerate;
        erast{j}=[erast{j};(ts) i*ones(size(ts))];
    end
end
for i=1:Nu
    a_psth(i)=axes('Position',p_psth(i,:));
    [a_psth(i),fmed(i),fpeak(i),flat(i)]=rast_psth(erast{i},[-ops.rasterpre,ops.rasterpost],[],[],COL(i+2,:),[1 1 1],6);
    xlabel('T(ms)');
end
%% cov
T=25;
for j=1:Nu
    x{j}=double(rast(find(rast(:,2)==U(j)),1))/ops.samplerate*1e3;    
end

k=1;
for i=1:Nu
    for j=i:numel(U)
        [c,bins]=MyXcor(x{i},x{j},T);
%         subplot(Nu,Nu,k);
        a_xcor(i,j)=axes('Position',p_xcor(k,:));
        H=bar(bins,c);
        if(i==j)
            set(H,'FaceColor',COL(i+2,:),'LineStyle','none');
        else
            set(H,'FaceColor',0.5*[1 1 1],'LineStyle','none');
        end
        
        set(gca,'Xlim',[-T T],'FontSize',16);
        xlabel('T(ms)');
        k=k+1;
    end
    k=k+i;
end
  
%% determine unit type
if(interactive)
    op.WindowStyle='normal';
    for i=1:numel(unitnums)
        prompt{i} = num2str(unitnums(i));
        defaultans{i} = '1';
    end
    dlg_title = 'Set unit type';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans,op);    
    utype=cellfun(@str2num,answer);
else
    utype=nan(1,numel(unitnums));
end
end
