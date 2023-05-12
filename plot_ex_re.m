function [exp1,exp2,rep1,rep2,exrs,rers,exvr,revr]=plot_ex_re(eod,file,ops,segnum,unitname,lfpname,t_col,IMGS,ind_lim)
% plot ex-afference and re-afference map for LFP and unit
% load('Z:\mormyrid_data\fish_images\straight_image.mat');
exlfpclim=[-1 1];
exuclim=[-.5 .5];
meanuclim=[0 1];
reclim=[-.05 .05];
stdclim=[0 1];
bgcol='w';
% mcol2=[74 133 34]/255;
% mcol1=[210 149 0]/255;
excol=[0.0433    0.7411    0.7394];
recol=[1 0 0];
eod1=eod(segnum);
file1=file(segnum);
seg1=ops.seg(segnum);

unitcol=find(cellfun(@(x) (strcmp(x,unitname)),eod1.fnames));
lfpcol=find(cellfun(@(x) (strcmp(x,lfpname)),eod1.fnames));

if(t_col==6)
    tname='tail';
elseif(t_col==1)
    tname='freq';
else
    tname=eod1.fnames(t_col);
end

if(strfind(file1.title,'brass')) %BRASS OBJECT
    if(t_col==lfpcol)
        F1=figure;
        [Exlfp,Nlfp]=plot_tuning2d(eod1,file1,lfpcol,...
            'clim',exlfpclim,'mfunc','mean','bgcol',bgcol,'image',IMGS.cdata,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A1=findobj(GUI_fig_children,'type','Axes');

        F2=figure;
        [ExU,NU]=plot_tuning2d(eod1,file1,unitcol,...
            'clim',meanuclim,'mfunc','mean','bgcol',bgcol,'image',IMGS.cdata,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A2=findobj(GUI_fig_children,'type','Axes');

        F3=figure;
        [Relfp,Nlfp]=plot_tuning2d(eod1,file1,lfpcol,...
            'clim',stdclim,'mfunc','std','bgcol',bgcol,'image',IMGS.cdata,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A3=findobj(GUI_fig_children,'type','Axes');

        F4=figure;
        [ReU,NU]=plot_tuning2d(eod1,file1,unitcol,...
            'clim',reclim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A4=findobj(GUI_fig_children,'type','Axes');

        F5=figure;
        idx=find(~isnan(Exlfp(:)) & ~isnan(ExU(:)));
        [f,g]=fit(Exlfp(idx),ExU(idx),'poly1');
        exp1=f.p1;
        exp2=f.p2;
        exrs=g.rsquare;
        H=histogram(ReU(idx),[-.5:0.025:.5]);
        set(gca,'FontSize',18)    
    %     xlabel(''); ylabel(''); legend('off');
        rep1=nanstd(ReU(idx));
        rep2=nanmean(ReU(idx));
        rers=nan;
        exvr=nan;
        revr=nan;
        GUI_fig_children=get(gcf,'children');
        A5=findobj(GUI_fig_children,'type','Axes');

    else   

        F1=figure;
        [Exlfp,Nlfp]=plot_tuning2d(eod1,file1,lfpcol,...
            'clim',exlfpclim,'mfunc','offset','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A1=findobj(GUI_fig_children,'type','Axes');

        F2=figure;
        [ExU,NU]=plot_tuning2d(eod1,file1,unitcol,...
            'clim',exuclim,'mfunc','offset','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A2=findobj(GUI_fig_children,'type','Axes');

        F3=figure;
        [Relfp,Nlfp]=plot_tuning2d(eod1,file1,lfpcol,...
            'clim',reclim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A3=findobj(GUI_fig_children,'type','Axes');

        F4=figure;
        [ReU,NU]=plot_tuning2d(eod1,file1,unitcol,...
            'clim',reclim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A4=findobj(GUI_fig_children,'type','Axes');

        F5=figure;
        [exp1,exp2,exrs,exvr]=compare_tunings(Exlfp,Nlfp,ExU,NU,0,excol);
        [rep1,rep2,rers,revr]=compare_tunings(Relfp,Nlfp,ReU,NU,0,recol);
        GUI_fig_children=get(gcf,'children');
        A5=findobj(GUI_fig_children,'type','Axes');

    end
else
    if(t_col==lfpcol)
        F1=figure;
        [Exlfp,Nlfp]=plot_tuning_polar(eod1,file1,lfpcol,...
            'clim',exlfpclim,'mfunc','mean','bgcol',bgcol,'image',IMGS.cdata,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A1=findobj(GUI_fig_children,'type','Axes');

        F2=figure;
        [ExU,NU]=plot_tuning_polar(eod1,file1,unitcol,...
            'clim',meanuclim,'mfunc','mean','bgcol',bgcol,'image',IMGS.cdata,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A2=findobj(GUI_fig_children,'type','Axes');

        F3=figure;
        [Relfp,Nlfp]=plot_tuning_polar(eod1,file1,lfpcol,...
            'clim',stdclim,'mfunc','std','bgcol',bgcol,'image',IMGS.cdata,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A3=findobj(GUI_fig_children,'type','Axes');

        F4=figure;
        [ReU,NU]=plot_tuning_polar(eod1,file1,unitcol,...
            'clim',reclim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A4=findobj(GUI_fig_children,'type','Axes');

        F5=figure;
        idx=find(~isnan(Exlfp(:)) & ~isnan(ExU(:)));
        [f,g]=fit(Exlfp(idx),ExU(idx),'poly1');
        exp1=f.p1;
        exp2=f.p2;
        exrs=g.rsquare;
        H=histogram(ReU(idx),[-.5:0.025:.5]);
        set(gca,'FontSize',18)    
    %     xlabel(''); ylabel(''); legend('off');
        rep1=nanstd(ReU(idx));
        rep2=nanmean(ReU(idx));
        rers=nan;
        exvr=nan;
        revr=nan;
        GUI_fig_children=get(gcf,'children');
        A5=findobj(GUI_fig_children,'type','Axes');

    else   

        F1=figure;
        [Exlfp,Nlfp]=plot_tuning_polar(eod1,file1,lfpcol,...
            'clim',exlfpclim,'mfunc','offset','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A1=findobj(GUI_fig_children,'type','Axes');

        F2=figure;
        [ExU,NU]=plot_tuning_polar(eod1,file1,unitcol,...
            'clim',exuclim,'mfunc','offset','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A2=findobj(GUI_fig_children,'type','Axes');

        F3=figure;
        [Relfp,Nlfp]=plot_tuning_polar(eod1,file1,lfpcol,...
            'clim',reclim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A3=findobj(GUI_fig_children,'type','Axes');

        F4=figure;
        [ReU,NU]=plot_tuning_polar(eod1,file1,unitcol,...
            'clim',reclim,'mfunc','slope','bgcol',bgcol,'image',IMGS.cdata,'t_col',t_col,'ind_lim',ind_lim);
        GUI_fig_children=get(gcf,'children');
        A4=findobj(GUI_fig_children,'type','Axes');

        F5=figure;
        [exp1,exp2,exrs,exvr]=compare_tunings(Exlfp,Nlfp,ExU,NU,0,excol);
        [rep1,rep2,rers,revr]=compare_tunings(Relfp,Nlfp,ReU,NU,0,recol);      
        GUI_fig_children=get(gcf,'children');
        A5=findobj(GUI_fig_children,'type','Axes');

    end
end    
fig=figure;
a1=copyobj(A1,fig);
set(a1,'Position',[0.1300    0.5838    0.2134    0.3412]);
a2=copyobj(A2,fig);
set(a2,'Position',[0.4108    0.5838    0.2134    0.3412]);
a3=copyobj(A3,fig);
set(a3,'Position',[0.1300    0.11      0.2134    0.3412]);
colormap(a3,modified_jet);
a4=copyobj(A4,fig);
set(a4,'Position',[0.4108    0.11      0.2134    0.3412]);
colormap(a4,modified_jet);
a5=copyobj(A5,fig);
set(a5,'Position',[0.6916    0.33      0.2134    0.3412]);
if(t_col~=lfpcol)
    set(a5,'Ylim',[-1 1]);
end
title([tname,';',unitname]);
close(F1,F2,F3,F4,F5)
end