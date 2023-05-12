function plot_lfp_preogression(ops,eod)
%PREPCA_INTERFACE go over all segments, check for failed LFP and revise
%   LFPgroups
x=[1:size(eod(1).avtraces,2)]/ops.samplerate*1e3+ops.traceblnk;
Y=[];
D=[];
for i=1:numel(ops.seg)
        Y=[Y;permute(eod(i).avtraces(9,:,:,:),[4 2 3 1])];
        D=[D;size(eod(i).avtraces,4)];
end
DD=[0;cumsum(D)];
for j=1:numel(ops.seg(1).LFPgroups)
    lfgorup=ops.seg(1).LFPgroups{j};
    F=figure;
    for k=1:numel(lfgorup)
        subplot(1,numel(lfgorup),k);
        imagesc(x,1:size(Y,1),Y(:,:,lfgorup(k)));
        hold on;
        plot([x(1) x(end)],DD*[1 1],'k');
        set(gca,'YTick',edge2bin(DD'),'YTickLabel',[1:numel(D)]);
        set(gca,'Clim',[-500 500]);        
    end
    lfpname=num2str(lfgorup);    
end


end
