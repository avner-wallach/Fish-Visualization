function [ a  ] = plot_sorted_lfp(lfp,t,inds,avtraces,ops,bgcol,axes_in)
%PLOT_SORTED_LFP plot the effect of variable t on lfp in index subsets
%inds. show scatter, fit and traces
% clim=[-2 2];
clim=[-1 1];
trace_num=7;
markersize=1.667;
nline=1;
wline=2;
% markersize=3;
% nline=.25;
% wline=.5;
fsize=7;
x0=2; y0=-2;
x1=3; y1=-1.5;
p=linspace(0,1,trace_num);
Q=edge2bin(quantile(lfp,linspace(0,1,11)));
N=numel(inds);
if(nargin<7)
    F=figure;
    axes_in=axes;
end
[a1,p1]=tight_subplot(1, 3, [.0 .025],[.0 .00],[.00 .0],[],[],axes_in); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[a11,p11]=tight_subplot(N, 1, [.025 .0],[.0 .00],[.00 .0],[],[],a1(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[a12,p12]=tight_subplot(N, 1, [.025 .0],[.0 .00],[.00 .0],[],[],a1(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
a=[a11;a12;a1(3)];
% COL=colormap('parula'); 
if(bgcol==[0 0 0])
    COL=invert_map(flipud(brewermap(64,'RdBu')));
else
    COL=flipud(brewermap(64,'RdBU'));
end
CN=size(COL,1); 
CC=brewermap(64,'Greys');    CCN=size(CC,1);
if(bgcol==[0 0 0])
    CC=flipud(CC);
end
CC=interp1(1:CCN,CC,linspace(16,64,64));
for i=1:N
    ind=inds{i};
    ind=ind(~isnan(t(ind)) & ~isnan(lfp(ind)));
        
    axes(a(N+i));
    S=scatter(t(ind),lfp(ind),markersize,t(ind),'filled');
    set(S,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%,'MarkerFaceColor',COL(1,:));    
    hold on;    
    f=fit(t(ind),lfp(ind),'poly1');
%     H=plot(f);
%     set(H,'Color',0.5*[1 1 1],'LineWidth',wline);
%     cval=min(max((f(0)-clim(1))/diff(clim)*(CN-1)+1,1),CN);
%     MCOL=interp1(1:CN,COL,cval);
    cval=min(max((f.p1-clim(1))/diff(clim),0),1);
    MCOL=interp1(linspace(0,1,CN),COL,cval);
%     H=plot(0,f(0),'o');
%     set(H,'MarkerFaceColor',MCOL,'MarkerSize',markersize,'MarkerEdgeColor',1-bgcol);
    
    XX=linspace(min(t(ind)),max(t(ind)),1e2);
    p = predint(f,XX,0.95,'functional','on');
    A1=patch([XX fliplr(XX)],[p(:,2)' fliplr(p(:,1)')],MCOL,'LineStyle','none','FaceColor',MCOL,'FaceAlpha',0.5);
    hold on;
    H=plot(XX,f(XX)); set(H,'Color',MCOL,'LineWidth',wline); legend('off');
    H=plot(XX,p(:,1)); set(H,'Color',0.5*[1 1 1],'LineWidth',nline); legend('off');
    H=plot(XX,p(:,2)); set(H,'Color',0.5*[1 1 1],'LineWidth',nline); legend('off');
    ci(1:2,:,i)=confint(f,0.95);
    ci(3,1,i)=f.p1; ci(3,2,i)=f.p2;
    mc(i,:)=MCOL;
    
    legend('off');    
    xlabel('');
    ylabel('');
    if(i==N) %last 
        H=plot(x0+[0 1],y0+[0 0],'k',x0+[0 0],y0+[0 1],'k');
        set(H,'LineWidth',nline);
    end
    set(gca,'Clim',[-3 3],'XAxisLocation','origin','YAxisLocation','origin',...
        'Ylim',[-3 3],'Xlim',[-4 4],'YTick',[[]],'XTick',[],'FontSize',fsize);    
    set(gca,'Color',bgcol,'YColor',1-bgcol,'XColor',1-bgcol,'Colormap',CC);    

    axes(a(i));%=subplot(numel(inds),2,2*i-1);
%     q=quantile(t(ind),p);
    q=linspace(nanmin(t(ind)),nanmax(t(ind)),trace_num);
    y=f(q); %fit value across t range
    y=max(min(y,Q(end)),Q(1));
    for j=1:1%size(avtraces,3)
        traces=interp1(Q,avtraces(:,:,j),y,'linear');
        x=[1:420]/30;       
        plot_graded(x,traces/1e3,CC);
        hold on;
        mtrace=interp1(Q,avtraces(:,:,j),f(0),'linear');        
%         H=plot(x,mtrace/1e3);
%         set(H,'LineWidth',wline,'Color',MCOL);
        set(gca,'Xlim',[x(1) 5],'ylim',[-2 1],'FontSize',16);
        set(gca,'Color','none','YColor','none','XColor','none');
        if(i==N) %last 
            H=plot(x1+[0 1],y1+[0 0],'k',x1+[0 0],y1+[0 1],'k');
            set(H,'LineWidth',nline);
        end
    end                   
end

axes(a(end));
bw=.2;
j=1;
for i=1:N
    H=rectangle('Position',[i-bw/2 ci(1,j,i) bw diff(ci(1:2,j,i))]);
    set(H,'FaceColor',[mc(i,:),0.5],'EdgeColor',0.5*[1 1 1]);
    hold on;
    H=plot(i,ci(3,j,i),'o');
    set(H,'MarkerEdgeColor',1-bgcol,'MarkerFaceColor',mc(i,:),'MarkerSize',5);        
end
set(gca,'Color',bgcol,'YColor',1-bgcol,'XColor',1-bgcol,'Xlim',[.5 i+.5],'Ylim',[-1 1],'XTick',[],'XAxisLocation','origin','Box','off');    


end

