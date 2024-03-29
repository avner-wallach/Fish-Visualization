% plot lfp and units motor modulations 
load('Z:\analysis_data\20190624\output.mat');
load('Z:\mormyrid_data\fish_images\straight_image.mat');
ind_lim=[5e5 inf];
clim=[-.75 .75];
bgcol='w';
mcol2=[74 133 34]/255;
mcol1=[210 149 0]/255;
%% tail effects
t_col=6;
figure, [Mlfp,Nlfp]=plot_tuning_polar(eod(1),file(1),[20],...
    'zfunc',@mean,'clim',clim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1]);
set(gcf,'Color',[1 1 1]);
% colorbar('off');
set(findobj(gca,'Type','Line'),'LineStyle','none')

figure, [Mu1,Nu1]=plot_tuning_polar(eod(1),file(1),[26],...
    'zfunc',@sum,'clim',clim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1]);
set(gcf,'Color',[1 1 1]);
colorbar('off');
set(findobj(gca,'Type','Line'),'LineStyle','none')

figure, [Mu2,Nu2]=plot_tuning_polar(eod(1),file(1),[25],...
    'zfunc',@mean,'clim',clim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1]);
set(gcf,'Color',[1 1 1]);
colorbar('off');
set(findobj(gca,'Type','Line'),'LineStyle','none')

figure;
compare_tunings(Mlfp,Nlfp,Mu1,Nu1,0,mcol1);
compare_tunings(Mlfp,Nlfp,Mu2,Nu2,0,mcol2);
set(gca,'Xlim',[-0.5 0.5],'Ylim',[-0.5 0.5]);
set(gcf,'Position',[616,509,420,420])
%% freq effects
t_col=1;
figure, [Mlfp,Nlfp]=plot_tuning_polar(eod(1),file(1),[20],...
    'zfunc',@mean,'clim',clim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1]);
set(gcf,'Color',[1 1 1]);
colorbar('off');
set(findobj(gca,'Type','Line'),'LineStyle','none')

figure, [Mu1,Nu1]=plot_tuning_polar(eod(1),file(1),[24],...
    'zfunc',@sum,'clim',clim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1]);
set(gcf,'Color',[1 1 1]);
colorbar('off');
set(findobj(gca,'Type','Line'),'LineStyle','none')

figure, [Mu2,Nu2]=plot_tuning_polar(eod(1),file(1),[25],...
    'zfunc',@mean,'clim',clim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1]);
set(gcf,'Color',[1 1 1]);
colorbar('off');
set(findobj(gca,'Type','Line'),'LineStyle','none')

figure;
compare_tunings(Mlfp,Nlfp,Mu1,Nu1,0,mcol1);
compare_tunings(Mlfp,Nlfp,Mu2,Nu2,0,mcol2);
set(gca,'Xlim',[-0.5 0.5],'Ylim',[-.8 0.25]);
set(gcf,'Position',[616,509,420,420])