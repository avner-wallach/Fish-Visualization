function [fileidx,frameidx]=get_frame_indices(I,ind0)        
    for j=1:numel(I)
        fidx(j)=find(ind0>I(j),1)-1;
    end
    fileidx=unique(fidx);
    for j=1:numel(fileidx)
        frameidx{j}=I(fidx==fileidx(j))-ind0(fileidx(j))+1;
    end       
                    
end
