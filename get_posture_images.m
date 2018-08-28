function [ imgs ] = get_posture_images(varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% params
params.headfix='on';
params.track='off';
params.bg='remove';
params.imgnum=5;
params.selfile=0;
%% get varins
n=numel(varargin);
baseargs=5;
if(n<baseargs | n>baseargs&(xor(mod(n,2),mod(baseargs,2))))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};  %target variable
ctl_col=varargin{4};    %variables controlled (kept around median)
pctl_ranges=varargin{5};
if(n>baseargs)
    i=baseargs+1;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
%%    
ctl=data_struct.data(:,ctl_col);
ctl_mean=nanmean(ctl,1);
ctl=ctl-ones(size(ctl,1),1)*ctl_mean;
z=data_struct.data(:,z_col);
for i=1:size(pctl_ranges,1)
    ind=1:numel(z);
    if(params.selfile>0)
        ind=ind(find(ind>=data_struct.ind0(params.selfile) & ...
            ind<data_struct.ind0(params.selfile+1)));
    end
    ind=ind(inrange(z(ind),pctl_ranges(i,:)));    
    [vals,idx] = sort(sum(abs(ctl(ind,:)),2),1);
    I=ind(idx(1:params.imgnum)); %image indices
    [fileidx,frameidx]=get_frame_indices(I,data_struct.ind0);    
    if(params.selfile==0)
        [imgs(i,:)] = get_fish_image(fileidx,frameidx,...
            'headfix',params.headfix,'bg',params.bg,'track',params.track);
    else
        [imgs(i,:)] = get_fish_image(fileidx,frameidx,...
            'headfix',params.headfix,'bg',params.bg,'track',params.track,...
            'vid',params.vid,'txt',params.txt,'data',params.data,'num',params.num);
    end
        
end

    function ind=inrange(vec,range)
        [range]=quantile(vec,range);
        ind=find(vec>range(1) & vec<=range(2));
    end

            
            

end