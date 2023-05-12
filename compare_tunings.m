function [p1,p2,rs,vr]=compare_tunings(M1,N1,M2,N2,bgcol,mcol)
C=[1 1 1]-bgcol;
N=max(N1,N2);
ind=find(N>0);
S=scatter(M1(ind),M2(ind),N(ind)/1e2,'o');
set(S,'MarkerEdgeColor',mcol,'MarkerEdgeAlpha',0.5,'MarkerFaceColor',mcol,'MarkerFaceAlpha',0.25);

ind=find(~isnan(M1(:)+M2(:)+N1(:)+N2(:)));
[f,g]=fit(M1(ind(:)),M2(ind(:)),'poly1','Weights',min(N1(ind(:)),N2(ind(:))));
hold on;
H=plot(f);
H.LineWidth=3;
H.Color=mcol;
set(gca,'Color',C,'XColor',1-C,'YColor',1-C,'FontSize',18)
set(gcf,'Color',C);
xlabel(''); ylabel(''); legend('off');
% p = predint(f,M1(ind(:)),0.95,'functional','on');
% plot(M1(ind(:)),p,'m--');
p1=f.p1;
p2=f.p2;
rs=g.rsquare;
vr=var(M1(ind),N1(ind))/var(M2(ind),N2(ind));
end

