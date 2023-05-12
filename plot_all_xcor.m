function [h,a] = plot_all_xcor(xcor,gnum)
%PLOT_ALL_XCOR plot all auto/xcor of spikes
if(~iscell(xcor.val))
    val=xcor.val;
    xcor=rmfield(xcor,'val');
    xcor.val{1}=val;
end
if nargin==1
    gnum=1:numel(xcor.val)
end
k=1;
for i=gnum
    h(k)=figure;
    N=size(xcor.val{i},2);
    ia=1;
    for u1=1:N
        for u2=u1:N
            k=(u1-1)*N+u2;
            subplot(N,N,k);
            B=area(xcor.bins,xcor.val{i}(:,u1,u2));
            B.LineStyle='none';
            a(ia)=gca;
            ia=ia+1;
        end
    end
    k=k+1;
end
end

