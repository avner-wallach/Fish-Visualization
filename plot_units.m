% plot spike shapes, xcor, raster psth
fsize=12;
imgpath='Z:\mormyrid_data\fish_images\';
load([imgpath,'straight_image.mat']);
str_imgs=IMGS;
%% parameters
expnum=7;
groupnum=1;
daynum=2;
collnum=1;
filenum='30';
figure; COL=colormap('lines'); close(gcf);
%% spike shape, raster, psth, xcor
load(['Z:\mormyrid_data\',experiment(expnum).day(daynum).date,'\data_',filenum]);
unitnum=numel(data.EOD.raster{groupnum});
t=[1:numel(data.SPIKES.trace(1).m)]/samplerate*1e3;
F1=figure;
for n=1:unitnum
    rast_psth(data.EOD.raster{groupnum}{n},[-20e-3 40e-3],[0 1e3],[],COL(n,:));   
    figure(F1);
    my_plotWithConf(t,data.SPIKES.trace(n).m,data.SPIKES.trace(n).s,COL(n,:));
    hold on;
end
ylabel('{\mu}V');
xlabel('ms');

F2=figure;
k=1;
for n=1:unitnum
    for m=1:unitnum
        if(m>=n)
            subplot(unitnum,unitnum,k);
            B=bar(data.SPIKES.xbins,data.SPIKES.xcor(:,n,m));
            if(n==m)
                set(B,'FaceColor',COL(n,:));
            else
                set(B,'FaceColor',0.5*[1 1 1]);
            end
            set(gca,'Xlim',[data.SPIKES.xbins(1) data.SPIKES.xbins(end)]);
            xlabel('lag(ms)');
            ylabel('Xcor');
        end
        k=k+1;
    end
end
%% lfp traces
chans=experiment(expnum).day(daynum).changroups{groupnum};
t=[1:size(eod(collnum).avtraces,2)]/samplerate*1e3+experiment(expnum).day(daynum).trace_b;
for i=1:numel(chans)
    plot_graded(t,eod(collnum).avtraces(:,:,chans(i),str2num(filenum)));
    set(gca,'Xlim',[0 10]);
    ylim=get(gca,'YLim');
    hold on;
    ar=area([0 experiment(expnum).day(daynum).trace_b],[ylim(2) ylim(2)],ylim(1),'FaceColor',.75*[1 1 1]);
    ar.LineStyle='none'    
end
ylabel('{\mu}V');
xlabel('ms');

%% plot tunings
cgname='group2';
for c=8:8%numel(experiment(expnum).collection);
    uind=find(cellfun(@(x) numel(strfind(x,'sc')),peod(c).fnames));
    for u=1:numel(uind)
        unum=uind(u);
%         M=quantile(peod(c).data(:,unum),0.9);
        M=mean(peod(c).data(:,unum))*2;
            switch ffile(c).title
                case 'none' %wall
                    figure;
                    plot_tuning_polar(peod(c),ffile(c),unum,'clim',[0 M],'image',str_imgs.cdata);
                case 'brass' %brass pole
                    figure;
                    plot_tuning2d(peod(c),ffile(c),unum,'clim',[0 M],'image',str_imgs.cdata);
                case 'aqu_brass'
                    figure;
                    plot_tuning2d(peod(c),ffile(c),unum,'clim',[0 M],'headfix','off','xbins',35,'ybins',35);
                    set(gca,'Alim',[-1 1],'YDir','reverse');
                case 'aqu_plastic'
                    figure;
                    plot_tuning2d(peod(c),ffile(c),unum,'clim',[0 M],'headfix','off','xbins',35,'ybins',35);
                    set(gca,'Alim',[-1 1],'YDir','reverse');
                    
                    
            end
    end
    lind=find(cellfun(@(x) numel(strfind(x,cgname)),peod(c).fnames));
    switch ffile(c).title
        case 'none' %wall
            figure;
            plot_tuning_polar(peod(c),ffile(c),lind,'clim',[-1 1],'image',str_imgs.cdata);
        case 'brass' %brass pole
            figure;
            plot_tuning2d(peod(c),ffile(c),lind,'clim',[-1 1],'image',str_imgs.cdata);
        case 'aqu_brass'
            figure;
            plot_tuning2d(peod(c),ffile(c),lind,'clim',[-1 1],'headfix','off','xbins',35,'ybins',35);
            set(gca,'Alim',[-1 1],'YDir','reverse');
        case 'aqu_plastic'
            figure;
            plot_tuning2d(peod(c),ffile(c),lind,'clim',[-1 1],'headfix','off','xbins',35,'ybins',35);
            set(gca,'Alim',[-1 1],'YDir','reverse');
            
    end
    
end
%% analyzed aquacalm data
load('Z:\mormyrid_data\analyzed data\20190131_aquacalm_data.mat');
figure, plot_tuning2d(eod([1 4]),file(1),20,'headfix','off','maxa',0.95,'clim',[-2 1.5],'x_nbins',35,'y_nbins',35);
figure, plot_tuning2d(eod([1 4]),file(1),22,'headfix','off','maxa',0.95,'clim',[0 4],'x_nbins',35,'y_nbins',35);
figure, plot_tuning2d(eod([1 4]),file(1),23,'headfix','off','maxa',0.95,'clim',[0 1.5],'x_nbins',35,'y_nbins',35);

%% analyzed free data
load('Z:\mormyrid_data\20190201\objects_6.mat');
load('Z:\mormyrid_data\fish_images\straight_image.mat');
load('Z:\mormyrid_data\analyzed data\20190131_alldata.mat');
file(3).objects=objects;
figure, plot_tuning_polar(eod([1]),file(1),25,'tank','circ','maxa',0.4,'clim',[-2 1.5],'image',IMGS.cdata);
figure, plot_tuning_polar(eod([1]),file(1),27,'tank','circ','maxa',0.4,'clim',[1 4],'image',IMGS.cdata);
figure, plot_tuning_polar(eod([1]),file(1),28,'tank','circ','maxa',0.4,'clim',[0 0.5],'image',IMGS.cdata);

figure, plot_tuning2d(eod([3]),file(3),22,'maxa',0.2,'clim',[-2 1.5],'image',IMGS.cdata);
figure, plot_tuning2d(eod([3]),file(3),24,'maxa',0.2,'clim',[0 .75],'image',IMGS.cdata);
figure, plot_tuning2d(eod([3]),file(3),25,'maxa',0.2,'clim',[0 3],'image',IMGS.cdata);
figure, plot_tuning2d(eod([3]),file(3),26,'maxa',0.2,'clim',[0 1],'image',IMGS.cdata);

