function [ group_struct ] = plot_io_maps_new(eod0,file0,ops,segnum,groupnum,objectseg)
%PLOT_IO_MAPS_NEW plot offset and slope maps due to variables tcols in both lfp
%and units. 
bgcol=[1 1 1];
tcols=[6 1];
T=numel(tcols);
ttitles={'Tail','Rate'};
offset_clim=[-1 1];
offset_c=[0.01 0.99];
slope_clim=[-1 1];
upsamp=1.5;
fsize=20;
[ret,sname]=system('hostname');
sep=filesep;

group_struct=[];
eod=eod0(segnum);
file=file0(segnum);
if(nargin==6)
    if(objectseg==0) %analyze using wall
        file.objects=[];
    else
        file.objects=file0(objectseg).objects;
    end
end
%get cols
lfpcol=find(cellfun(@(x) numel(strfind(x,['lfp',num2str(groupnum)])),eod.fnames));
allunitcols=find(cellfun(@(x) numel(strfind(x,['sc'])),eod.fnames));
indlim_inds=find(cellfun(@(x) numel(strfind(x,['sc',num2str(groupnum)])),eod.fnames(allunitcols)));
unitcols=allunitcols(indlim_inds);
U=numel(unitcols);
ind_lim=ops.seg(segnum).ind_lim(indlim_inds,:);

load([ops.datapath,sep,'fish_images',sep,'straight_image.mat']);
setenv('PXLSIZE','0.423'); %pxl to mm conversion
M=U+1;
N=T+M+2;
figure;
set(gcf,'WindowState','maximized');
[ha, pos] = tight_subplot(M, N, [.0 .0],[.1 .01],[.01 .01]);
ax_ind_remove=[T+2 T+3+[1:U]];
for i=1:numel(ax_ind_remove)
    ha(ax_ind_remove(i)).Visible='off';
end


    axes(ha(1)); %LFP ex
    [Mz,N0]=plot_tuning2d(eod,file,lfpcol,'t_col',tcols(2),'clim',offset_clim,...
        'image',IMGS.cdata,'mfunc','mean','bgcol',bgcol,'ind_lim',[],'upsamp',upsamp);
    colorbar('off');
%     c=quantile(Mz(:),offset_c);
%     set(gca,'Clim',c);
    T1=text(0,0.33,['LFP',eod.fnames{lfpcol}(4:end)],'Units','normalized','Rotation',90,'FontSize',fsize);
    T2=text(0.4,.9,'Offset','Units','normalized','Rotation',0,'FontSize',fsize);
   [group_struct.LFP.Hx,group_struct.LFP.Sx]=entropy(Mz);

    for j=1:T 
        tcol=tcols(j);
        axes(ha(j+1)); %LFP slope
        plot_tuning2d(eod,file,lfpcol,'t_col',tcol,'clim',slope_clim,...
            'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',[],'upsamp',upsamp);
        colorbar('off');
        T2=text(0.4,.9,[ttitles(j)],'Units','normalized','Rotation',0,'FontSize',fsize);

        [ Izt,Izt_xy,Izxy,Izxy_t,izt_xy ] = compute_info(eod,file,lfpcol,'t_col',tcol,'ind_lim',[]); 
        group_struct.LFP.Info(j).Uz_t=Izt; %uncrtainty coef. for variable
        group_struct.LFP.Info(j).Uz_t_xy=Izt_xy; %uncrtainty coef. for variable, given position

    end
    axes(ha(2+T));
    T2=text(0.4,.85,'I/O','Units','normalized','Rotation',0,'FontSize',fsize);

    for i=1:U
        unitcol=unitcols(i);
        i0=i*N;
        axes(ha(i0+1)); %unit ex
        [Mz,N0]=plot_tuning2d(eod,file,unitcol,'clim',offset_clim,...
            'image',IMGS.cdata,'mfunc','mean','t_col',tcols(2),'bgcol',bgcol,'ind_lim',ind_lim(i,:),'upsamp',upsamp);
        colorbar('off');
%         c=quantile(Mz(:),offset_c);
%         set(gca,'Clim',c);
        T1=text(0,0.33,['Unit',eod.fnames{unitcol}(3:end)],'Units','normalized','Rotation',90,'FontSize',fsize);
        [group_struct.Unit(i).Hx,group_struct.Unit(i).Sx]=entropy(Mz);        

        for j=1:T 
            tcol=tcols(j);
            axes(ha(i0+1+j)); %unit slope
            plot_tuning2d(eod,file,unitcol,'t_col',tcol,'clim',slope_clim,...
                'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim(i,:),'upsamp',upsamp);
            [ Izt,Izt_xy,Izxy,Izxy_t,izt_xy ] = compute_info(eod,file,unitcol,'t_col',tcol,'ind_lim',ind_lim(i,:)); 
            group_struct.Unit(i).Info(j).Uz_t=Izt; %uncrtainty coef. for variable
            group_struct.Unit(i).Info(j).Uz_t_xy=Izt_xy; %uncrtainty coef. for variable, given position
            colorbar('off');
        end
        
        axes(ha(i0+2+T)); %I/O
        plot_tuning2d(eod,file,unitcol,'t_col',lfpcol,'clim',slope_clim,...
            'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim(i,:),'upsamp',upsamp);
        colorbar('off');
        [ Izt,Izt_xy,Izxy,Izxy_t,izt_xy ] = compute_info(eod,file,unitcol,'t_col',lfpcol,'ind_lim',ind_lim(i,:)); 
        group_struct.Unit(i).Info(T+1).Uz_t=Izt; %uncrtainty coef. for variable
        group_struct.Unit(i).Info(T+1).Uz_t_xy=Izt_xy; %uncrtainty coef. for variable, given position
        
                
    end    

%rasters
for i=1:U
    unitcol=unitcols(i);
    i0=i*N;
    p=ha(i0+3+T).Position;
    x0=p(1); y0=p(2); w=p(3); h=p(4);
    x0=x0+w/6; h=h*.9; w=w*.8;
    hp(1)=axes('Position',[x0 y0 w h/3]);
    hp(2)=axes('Position',[x0 y0+h/3 w 2*h/3]);
    ha(i0+3+T).Visible='off';

    cmap{1}=repmat(linspace(0.75,0,64)',1,3);
    plot_sorted_raster(eod.raster{indlim_inds(i)},eod.data(:,lfpcol),ops,ind_lim(i,:),1e2,3,.5,cmap);
    close(gcf);
    ftmp=gcf;
    %         set(gca,'YLim',[0 120],'YTick',[],'XTick',[0 20 40]);
    ca=get(gcf,'Children');
    for i=1:numel(ca)
        copyaxes(ca(i),hp(i),1);
    end
    close(ftmp);
end
drawnow;
%LFP traces
chind=ops.seg(segnum).LFPgroups_pca{groupnum};
x=[1:size(eod.avtraces,2)]/ops.samplerate*1e3+ops.traceblnk;
ind=find(x<=6);
axes(ha(T+2));
plot_graded(x(ind),nanmedian(eod.avtraces(:,ind,chind(1),:),4));
ftr=gcf;
% ca=get(ftr,'Children');
copyaxes(gca,ha(T+3),1);
close(ftr);
set(ha(T+3),'XAxisLocation','origin');
P=ha(T+3).Position;
P(1)=P(1)+P(3)/8;
P(3)=0.85*P(3);
set(ha(T+3),'Position',P);

%XCOR
ind=find(ops.seg(segnum).Spkgroups==groupnum);
[H,hc]=plot_all_xcor(file.units.xcor,ind);
% hc=get(h,'Children');
t=1;
for k=1:U
    i0=k*N+T+4; 
    for q=1:U
        if(q<k)
            ha(i0).Visible='off';
        else
            copyaxes(hc(t),ha(i0),1);
            p=ha(i0).Position;
            x0=p(1); y0=p(2); w=p(3); h=p(4);
            x0=x0+w/6; h=h*.9; w=w*.8;
            ha(i0).Position=[x0 y0 w h];
            t=t+1;
        end
        i0=i0+1;            
    end
end
close(H);

axes(ha(T+4));
if(numel(file.objects))
    otype='Object';
else
    otype='Wall';
end
T_header=text(0.1,.65,['Date: ',num2str(ops.seg(1).dates(1)),'\newlineSegment: ',...
    num2str(segnum),'\newlineGroup: ',num2str(groupnum),'\newline',file.title,'\newline',otype],'Units','normalized','Rotation',0,'FontSize',fsize);

% set(gcf,'Position',[360 143 1020 747]);
end

