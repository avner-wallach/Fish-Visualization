function [ ha ] = plot_io_maps(eod,file,lfpcol,unitcols,tcols,object,ops,rastind,ind_lims)
%PLOT_IO_MAPS plot offset and slope maps due to variables tcols in both lfp
%and units. object=0 for wall analysis
if(nargin<7)
    ind_lim=[1 inf];
end
bgcol=[1 1 1];
offset_clim=[-1 1];
offset_c=[0.01 0.99];
slope_clim=[-1 1];
fsize=20;
[ret,sname]=system('hostname');
if(strfind(sname,'Mac'))
    rootdir='/Volumes/aw3057/';
    sep='/';
else
    rootdir='Z:\';
    sep='\';
end
load([rootdir,'mormyrid_data',sep,'fish_images',sep,'straight_image.mat']);
setenv('PXLSIZE','0.423'); %pxl to mm conversion
N=numel(tcols)+3;
M=numel(unitcols)+1;
figure;
[ha, pos] = tight_subplot(M, N, [.0 .0],[.1 .01],[.01 .01]);

if(object==0) %wall analysis
    axes(ha(1)); %LFP ex
    [Mz,N0]=plot_tuning_polar(eod,file,lfpcol,'t_col',tcols(1),'clim',offset_clim,...
        'image',IMGS.cdata,'mfunc','offset','bgcol',bgcol,'ind_lim',ind_lim);
    c=quantile(Mz(:),offset_c);
    set(gca,'Clim',c);
    colorbar('off');
    T1=text(0,0.33,'LFP','Units','normalized','Rotation',90,'FontSize',fsize);
    T2=text(0.4,.9,'Offset','Units','normalized','Rotation',0,'FontSize',fsize);

    for j=1:numel(tcols) 
        tcol=tcols(j);
        axes(ha(j+1)); %LFP slope
        plot_tuning_polar(eod,file,lfpcol,'t_col',tcol,'clim',slope_clim,...
            'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim);
        colorbar('off');
        T2=text(0.4,.9,[num2str(tcol),' Slope'],'Units','normalized','Rotation',0,'FontSize',fsize);
    end
%     if(nargin>6)
%         axes(ha(N));  %rastplot
%         cmap{1}=repmat(linspace(0.5,0,64)',1,3);
%         plot_sorted_raster(eod.raster{rastinf},eod.data(:,lfpcol),ops,ind_lim,1e2,7,.5,cmap);
%     else
        ha(N).Visible='off';
        T2=text(-0.6,.9,'I/O','Units','normalized','Rotation',0,'FontSize',fsize);
%     end
    for i=1:numel(unitcols)
        unitcol=unitcols(i);
        axes(ha(i*N+1)); %unit ex
        [Mz,N0]=plot_tuning_polar(eod,file,unitcol,'clim',offset_clim,...
            'image',IMGS.cdata,'mfunc','mean','bgcol',bgcol,'ind_lim',ind_lim);
        colorbar('off');
        c=quantile(Mz(:),offset_c);
        set(gca,'Clim',c);
        T1=text(0,0.33,'Unit','Units','normalized','Rotation',90,'FontSize',fsize);

        for j=1:numel(tcols) 
            tcol=tcols(j);
            axes(ha(i*N+1+j)); %unit slope
            plot_tuning_polar(eod,file,unitcol,'t_col',tcol,'clim',slope_clim,...
                'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim);
            colorbar('off');
        end

        axes(ha(i*N+N)); %I/O
        plot_tuning_polar(eod,file,unitcol,'t_col',lfpcol,'clim',slope_clim,...
            'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim);
        colorbar('off');

    end
    
    
else %object analysis
    axes(ha(1)); %LFP ex
    [Mz,N0]=plot_tuning2d(eod,file,lfpcol,'t_col',tcols(1),'clim',offset_clim,...
        'image',IMGS.cdata,'mfunc','offset','bgcol',bgcol,'ind_lim',[]);
    colorbar('off');
    c=quantile(Mz(:),offset_c);
    set(gca,'Clim',c);
    T1=text(0,0.33,['LFP',eod.fnames{lfpcol}(4:end)],'Units','normalized','Rotation',90,'FontSize',fsize);
    T2=text(0.4,.9,'Offset','Units','normalized','Rotation',0,'FontSize',fsize);
    

    for j=1:numel(tcols) 
        tcol=tcols(j);
        axes(ha(j+1)); %LFP slope
        plot_tuning2d(eod,file,lfpcol,'t_col',tcol,'clim',slope_clim,...
            'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',[]);
        colorbar('off');
        T2=text(0.4,.9,[num2str(tcol),' Slope'],'Units','normalized','Rotation',0,'FontSize',fsize);
    end
    axes(ha(N-1));
    T2=text(0.4,.9,'I/O','Units','normalized','Rotation',0,'FontSize',fsize);
    ha(N-1).Visible='off';

    axes(ha(N));
    T2=text(0.4,.9,'Raster','Units','normalized','Rotation',0,'FontSize',fsize);
    ha(N).Visible='off';


    for i=1:numel(unitcols)
        unitcol=unitcols(i);
        if(size(ind_lims,1)>1)
            ind_lim=ind_lims(i,:);
        else
            ind_lim=ind_lims;
        end

        axes(ha(i*N+1)); %unit ex
        [Mz,N0]=plot_tuning2d(eod,file,unitcol,'clim',offset_clim,...
            'image',IMGS.cdata,'mfunc','mean','t_col',tcols(1),'bgcol',bgcol,'ind_lim',ind_lim);
        colorbar('off');
        c=quantile(Mz(:),offset_c);
        set(gca,'Clim',c);
        T1=text(0,0.33,['Unit',eod.fnames{unitcol}(3:end)],'Units','normalized','Rotation',90,'FontSize',fsize);

        for j=1:numel(tcols) 
            tcol=tcols(j);
            axes(ha(i*N+1+j)); %unit slope
            plot_tuning2d(eod,file,unitcol,'t_col',tcol,'clim',slope_clim,...
                'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim);
            colorbar('off');
        end
        
        axes(ha(i*N+N-1)); %I/O
        plot_tuning2d(eod,file,unitcol,'t_col',lfpcol,'clim',slope_clim,...
            'image',IMGS.cdata,'mfunc','slope','bgcol',bgcol,'ind_lim',ind_lim);
        colorbar('off');
        
        if(nargin>6)
            axes(ha(N));  %rastplot
            p=ha(i*N+N).Position;
            x0=p(1); y0=p(2); w=p(3); h=p(4);
            x0=x0+w/6; h=h*.9; w=w*.8;
            hp(1)=axes('Position',[x0 y0 w h/3]);
            hp(2)=axes('Position',[x0 y0+h/3 w 2*h/3]);
            ha(i*N+N).Visible='off';

            cmap{1}=repmat(linspace(0.5,0,64)',1,3);
            plot_sorted_raster(eod.raster{rastind(i)},eod.data(:,lfpcol),ops,ind_lim,1e2,7,.5,cmap);
            close(gcf);
            ftmp=gcf;
    %         set(gca,'YLim',[0 120],'YTick',[],'XTick',[0 20 40]);
            ca=get(gcf,'Children');
            for i=1:numel(ca)
                copyaxes(ca(i),hp(i),1);
            end
            close(ftmp);
        end
        
    end
    
end        
set(gcf,'Position',[360 143 1020 747]);
end

