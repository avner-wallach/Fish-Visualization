function [ frames,ax_grid] = get_fish_image(varargin)
%GET_FISH_IMAGE extract image from video file 
%   Detailed explanation goes here
%%env params
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
framescale=str2num(getenv('FRAMESCALE')); %scaling of video file for feature tracking
sesspath=[datapath,'\',sdate,'\'];
%% params
params.headfix='off';
params.track='off';
params.bg='off';
params.objects='off';
params.plotter=[];
params.plotint=100; %in frames
params.cropW=500;
params.cropH=600;
params.vecL=50;
params.bgthreshold=35;
params.bglevel=256;
params.kopen=5;
params.kclose=5;
params.cthreshold=0.8;
params.txt=[];        
params.num=[];

%% vars
x0=0;
y0=0;
a0=0;
F=figure;
COL=colormap('lines');
close(F);
k=1;
%% get varins
n=length(varargin);
switch n
    case 1
        error('too few input arguments!');
    case 2 %get_fish_image(filenum,idx)
        filenums=varargin{1};
        ids=varargin{2};
    otherwise %get_fish_image(filenum,idx,{fname,fvalue})
        filenums=varargin{1};
        ids=varargin{2};
        if(mod(n,2)) %odd number of args
            error('wrong number of input arguments!');
        else
            i=3;
            while(i<n)
                params.(varargin{i})=varargin{i+1};
                i=i+2;
            end
        end
end
for v=1:numel(filenums)
    if(~isfield(params,'vid'))
        % get frame from video
        filenum=num2str(filenums(v));
        vid = VideoReader([sesspath,'video_',filenum,'.avi']);
        % get data file
        datafile=[sesspath,'data_',filenum];
        s=load(datafile);        
        data=s.data;        
        %get tracking file
        if(strcmp(params.track,'features'))
            trackfile=[sesspath,'trackedFeaturesRaw_',filenum];
            [num,txt,raw] = xlsread([trackfile,'.csv']);    
            fnames=setdiff(unique(txt),{'coords','likelihood','scorer','bodyparts','x','y'});
        end        
    else
        vid=params.vid;
        data=params.data;
        if(strcmp(params.track,'features'))
            txt=params.txt;
            fnames=setdiff(unique(txt),{'coords','likelihood','scorer','bodyparts','x','y'});
            num=params.num;
        end
    end
    if(numel(filenums)==1 & ~iscell(ids))
        id=ids;
    else
        id=ids{v};
    end
    
    for j=1:numel(id)
        idx=id(j);
        vid.CurrentTime=(idx-1)/vid.FrameRate;
%         idx=round(data.FRAME.t(id(j))*50);
%         frame = rgb2gray(read(vid, idx));
        frame=rgb2gray(readFrame(vid));
        process_frame;
        frames(k)=getframe;
        k=k+1;
    end
end

    function process_frame       
    %% remove bg if selected
    if(strcmp(params.bg,'remove'))
        mask1=imfuse(data.FILE.BG,frame,'diff');
        mask2=uint8(mask1>params.bgthreshold);
        mask4=(imclose(imopen(mask2,ones(params.kopen)),ones(params.kclose)));
        frame=frame.*mask4+params.bglevel*uint8(~mask4);
    end
    %% if in ego mode or vector tracking mode- get location and heading azimut
    if(strcmp(params.headfix,'on') | ~strcmp(params.track,'off'))
        xind=cellfun(@(x) strcmp(x,'X'),data.FILE.model);
        yind=cellfun(@(x) strcmp(x,'Y'),data.FILE.model);
        aind=cellfun(@(x) strcmp(x,'azim'),data.FILE.model);    
        X=data.FRAME.posture(idx+1,xind);
        Y=data.FRAME.posture(idx+1,yind);
        a=data.FRAME.posture(idx+1,aind);
    end
    %% tracking mode= 'features'
    if(strcmp(params.track,'features'))
        for i=1:numel(fnames)
            ind=find(strcmp(txt(2,:),fnames{i}));    
%             if(~strcmp(txt(1,1),'bodyparts'))
%                 ind=ind+1;
%             end
            if(numel(ind))
                x(i)=num(idx+1,ind(1))*framescale;
                y(i)=num(idx+1,ind(2))*framescale;
                c(i)=num(idx+1,ind(3));        
            end
        end
        x(c<params.cthreshold)=NaN;
        y(c<params.cthreshold)=NaN;
    end
    %%    object mode='on'
    if(strcmp(params.objects,'on'))
        for i=1:numel(data.FILE.objects)
            xobj(i)=data.FILE.objects(i).x;
            yobj(i)=data.FILE.objects(i).y;
        end
    end
    %% ego mode- 'headfix' fish
    if(strcmp(params.headfix,'on'))
        x0=size(frame,2)/2;
        y0=size(frame,1)/2;
        frame=imtranslate(frame,[x0-X,y0-Y],'FillValues',params.bglevel); %center image around fish
        frame=imrotate(frame,rad2deg(a)+90,'bilinear','crop');
        frame=imcrop(frame,[x0-params.cropW/2,y0-params.cropH/3,params.cropW,params.cropH]);        
    %     y0=params.cropH/4;
    %     x0=params.cropW/2;
    %     a0=a+pi/2;
        if(strcmp(params.track,'features'))
            [x,y]=allo2ego(x,y,a,X,Y);
            x=x+params.cropW/2;
            y=-y+params.cropH/3;
        end
        if(strcmp(params.objects,'on'))
            [xobj,yobj]=allo2ego(xobj,yobj,a,X,Y);
%             xobj=xobj+params.cropW/2;
%             yobj=-yobj+params.cropH/3;
        end
        
        ix=([1:params.cropW]-params.cropW/2);
        iy=([1:params.cropH]-2*params.cropH/3);
            
        a=-pi/2;
        X=params.cropW/2;
        Y=params.cropH/3;
    else
            ix=[1:size(frame,2)];
            iy=[1:size(frame,1)];
    end
    %% plotter NOT DONE
    if(numel(params.plotter))
        ind=[max(1,idx-params.plotint):idx];
        t=data.FRAME.t(ind);        
    end
    %% display image ?    
    imagesc(ix,iy,repmat(flipud(frame),[1 1 3]));

    ax_grid={ix,iy};
%     im = imshow(frame);
    hold on;

    if(strcmp(params.track,'vector'))
        U=params.vecL*cos(a);
        V=params.vecL*sin(a);
        H=quiver(X,Y,U,V,'LineWidth',3,'MaxHeadSize',10);    
    elseif(strcmp(params.track,'features'))
        pecind=find(cellfun(@(x) numel(strfind(x,'Pec')),fnames));
        chinind=find(cellfun(@(x) numel(strfind(x,'chin')),fnames));
        for(i=1:numel(x))
            H=plot(x(i),y(i),'.');
            set(H,'Color',COL(1,:),'MarkerSize',20);
            if(ismember(i,pecind))
                set(H,'Color',COL(2,:));
            elseif(ismember(i,chinind))
                set(H,'Color',COL(3,:));
            end
        end
    end
    if(strcmp(params.objects,'on'))
        H=plot(xobj,yobj,'.');
        set(H,'Color',COL(4,:),'MarkerSize',20);        
    end
    hold off;

    set(gca,'YDir','normal');   
    axis('image');
    set(gca,'Xlim',[ix(1) ix(end)],'Ylim',[min(iy) max(iy)]);    
    set(gca,'XTick',[],'YTick',[]);

    end
end

