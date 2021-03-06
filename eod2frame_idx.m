function [frameidx]=eod2frame_idx(fileidx,eodidx,eod,frame)
    for j=1:numel(fileidx)
        eidx=eodidx{j}+eod.ind0(j)-1; %idx in eod struct
        te=eod.t(eidx);  %times within eod struct
        tf=frame.t(frame.ind0(j):(frame.ind0(j+1)-1)); %times within frame struct
        for k=1:numel(eidx)
            [m,fidx(k)]=min(abs(tf-te(k)));
        end
        frameidx{j}=fidx;
    end
end
