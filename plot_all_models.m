% pth='/Volumes/sawtell-locker/home/aw3057/'
% pth='z:/';
% load([pth,'mormyrid_data/analyzed data/all_networks.mat']);
% load([pth,'mormyrid_data/fish_images/straight_image.mat']);
str_imgs=IMGS;
c='*+o^';
c='....';
bgcol=1*[1 1 1]; %0=black, 1=white
msize=20;
p=1;
F=figure;
set(gcf,'Color',bgcol);
COL=colormap('lines');
boxcol=lighter(COL(1,:),.5);
COL(1:4,:)=[1;1;1;1]*COL(1,:);
allR=[]; allG=[];
for i=1:numel(wall_nets)
    for j=1:numel(wall_nets(i).net_struct)
        net=wall_nets(i).net_struct(j);
        clear R;
        R=[net.model_loc.rsq net.model_niloc.rsq net.model_lfbloc.rsq net.model_fbloc.rsq net.model_fbniloc.rsq];
        R=[net.model_loc.rsq net.model_niloc.rsq net.model_fbniloc.rsq];
        
%         R=[R R(2)+R(5)-R(1)];
        R=R-(p*R(1));
        allR=[allR;R'];
        allG=[allG;[1:numel(R)]'];
    end    
end
boxplot(allR,allG,'PlotStyle','compact','Colors',ones(5,1)*boxcol,'Symbol','');
A=findobj(gca,'Tag','Median');
delete(A);
A=findobj(gca,'Tag','MedianOuter');
delete(A);
A=findobj(gca,'Tag','Box');
set(A,'LineWidth',10)

hold on;
for i=1:numel(wall_nets)
    for j=1:numel(wall_nets(i).net_struct)
        net=wall_nets(i).net_struct(j);
        clear R;
        R=[net.model_loc.rsq net.model_niloc.rsq net.model_lfbloc.rsq net.model_fbloc.rsq net.model_fbniloc.rsq];
        R=[net.model_loc.rsq net.model_niloc.rsq net.model_fbniloc.rsq];
        
%         R=[R R(2)+R(5)-R(1)];
        R=R-(p*R(1));
        H=plot(R,':');
        hold on;
        set(H,'Color',COL(i,:),'Marker',c(j),'MarkerSize',msize,'LineWidth',2);
    end    
end
set(gca,'Color',bgcol,'XTick',[1:6],'FontSize',18,'XColor',1-bgcol,'YColor',1-bgcol);
set(gca,'XTickLabel',{'no NI','Motor','local FB','global FB','Motor+FB','M+FB'},'XTickLabelRotation',90);


F=figure;
set(gcf,'Color',bgcol);
% COL=colormap('lines');
allR=[]; allG=[];
for i=1:numel(brass_nets)
    for j=1:numel(brass_nets(i).net_struct)
        net=brass_nets(i).net_struct(j);
        clear R;
        R=[net.model_loc.rsq net.model_niloc.rsq net.model_lfbloc.rsq net.model_fbloc.rsq net.model_fbniloc.rsq];
%         R=[R R(2)+R(5)-R(1)];
        R=R-(p*R(1));
        allR=[allR;R'];
        allG=[allG;[1:numel(R)]'];
    end    
end
boxplot(allR,allG,'PlotStyle','compact','Colors',ones(5,1)*boxcol,'Symbol','');
A=findobj(gca,'Tag','Median');
delete(A);
A=findobj(gca,'Tag','MedianOuter');
delete(A);
A=findobj(gca,'Tag','Box');
set(A,'LineWidth',10)

hold on;
for i=1:numel(brass_nets)
    for j=1:numel(brass_nets(i).net_struct)
        net=brass_nets(i).net_struct(j);
        clear R;
        R=[net.model_loc.rsq net.model_niloc.rsq net.model_lfbloc.rsq net.model_fbloc.rsq net.model_fbniloc.rsq];
%         R=[R R(2)+R(5)-R(1)];
        R=R-(p*R(1));
        H=plot(R,':');
        hold on;
        set(H,'Color',COL(i,:),'Marker',c(j),'MarkerSize',msize,'LineWidth',2);
    end    
end
set(gca,'Color',bgcol,'XTick',[1:6],'FontSize',18,'XColor',1-bgcol,'YColor',1-bgcol);
set(gca,'XTickLabel',{'no NI','M','local FB','global FB','MFB','M+FB'},'XTickLabelRotation',90);
%% prsq
k=1
for i=1:numel(wall_nets)
    for j=1:numel(wall_nets(i).net_struct)
        net=wall_nets(i).net_struct(j);
        figure;
        plot_polar_model(net,'model','loc','image',255-str_imgs.cdata,'clim',[-1 1],'r_lim',[0 200]);        
        title(['fish=',num2str(i),'; RF=',num2str(j),'; k=',num2str(k)],'Color',[1 1 1]);
        P(k,:)=net.model_loc_mot.prsq;
        k=k+1;
    end
end
figure;
imagesc(P');

k=1
for i=1:numel(brass_nets)
    for j=1:numel(brass_nets(i).net_struct)
        net=brass_nets(i).net_struct(j);
        figure;
        plot_polar_model(net,'model','loc','image',255-str_imgs.cdata,'r_lim',[0 150]);        
        title(['fish=',num2str(i),'; RF=',num2str(j),'; k=',num2str(k)],'Color',[1 1 1]);
        P(k,:)=net.model_loc_mot.prsq;
        k=k+1;
    end
end
figure;
imagesc(P');