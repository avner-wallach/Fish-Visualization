function [ lfp_ratio ] = plot_lfp(eod,ops,segnum,lfpnum,fileind)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
trind=ops.seg(segnum).LFPgroups_pca{lfpnum};
x=[1:size(eod(segnum).avtraces,2)]/ops.samplerate*1e3 + ops.traceblnk;
r=[];
for i=1:numel(trind)
    plot_graded(x,eod(segnum).avtraces(:,:,trind(i),fileind));
    [m,mind]=min(eod(segnum).avtraces(:,:,trind(i),fileind)');
    [M,Mind]=max(eod(segnum).avtraces(:,:,trind(i),fileind)');
%     if(mean(Mind)<mean(mind)) %positive first
%         r=[r M./(-m)];
%     else %negative first
    r=[r (M-(abs(m)))./(M+abs(m))];
end
lfp_ratio=mean(r);
end

