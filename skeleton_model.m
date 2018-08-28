function [fish,seg,data]=skeleton_model(txt,num)
% get DeepLabCut tracking data and convert to model
% [num,txt,raw] = xlsread('trackedFeaturesRaw.csv');
    fnames=unique(txt);
    for i=1:numel(fnames)
        ind=find(strcmp(txt,fnames{i}))+1;    
        data.(fnames{i}).xy=num(:,(ind(1):ind(2)));
        data.(fnames{i}).c=num(:,(ind(3)));
    end

    %% get midpoint between LED and Trunk1
    data.midpoint.xy=[mean([data.LED.xy(:,1),data.Trunk1.xy(:,1)],2) ...
        mean([data.LED.xy(:,2),data.Trunk1.xy(:,2)],2)];
    data.midpoint.c=min([data.LED.c data.Trunk1.c],[],2);
    
    %% coordinate names
    coornames={'X','Y','Z','azim','b_chin','b_tr2','b_tl1','b_tl2','b_cf',...
        'c_rb','c_rt','c_lb','c_lt'};

    th=0.8;
    %% extract location (I leave Z position for later)
    %X
    fish(:,strcmp(coornames,'X'))=data.LED.xy(:,1);
    fish(data.LED.c<th,strcmp(coornames,'X'))=NaN;

    %Y
    fish(:,strcmp(coornames,'Y'))=data.LED.xy(:,2);
    fish(data.LED.c<th,strcmp(coornames,'X'))=NaN;

    %azim
    mouthv=data.mouth.xy-data.LED.xy;
    azim=atan2(mouthv(:,2),mouthv(:,1));
    fish(:,strcmp(coornames,'azim'))=azim;
    fish(data.LED.c<th | data.mouth.c<th,strcmp(coornames,'azim'))=NaN;
    
    %beta chin
    fish(:,strcmp(coornames,'b_chin'))=get_angle(data.Trunk1,data.mouth,data.chin);

%     %beta trunk1
%     fish(:,strcmp(coornames,'b_tr1'))=get_angle(data.mouth,data.LED,data.Trunk1);

    %beta trunk2
    fish(:,strcmp(coornames,'b_tr2'))=get_angle(data.mouth,data.Trunk1,data.Trunk2);

    %beta tail1
    fish(:,strcmp(coornames,'b_tl1'))=get_angle(data.Trunk1,data.Trunk2,data.Tail1);

    %beta tail2
    fish(:,strcmp(coornames,'b_tl2'))=get_angle(data.Trunk2,data.Tail1,data.Tail2);

    %beta cuadal fork
    fish(:,strcmp(coornames,'b_cf'))=get_angle(data.Tail1,data.Tail2,data.CaudalFork);

    %gamma right pect base
    fish(:,strcmp(coornames,'c_rb'))=get_angle(data.LED,data.midpoint,data.RPecBase);

    %gamma right pect tip
    fish(:,strcmp(coornames,'c_rt'))=get_angle(data.midpoint,data.RPecBase,data.RPecTip);

    %gamma left pect base
    fish(:,strcmp(coornames,'c_lb'))=get_angle(data.LED,data.midpoint,data.LPecBase);

    %gamma left pect tip
    fish(:,strcmp(coornames,'c_lt'))=get_angle(data.midpoint,data.LPecBase,data.LPecTip);
    
    %% extract segment length
    seg(1)=get_segment(data.Trunk1,data.mouth)/2; %LED-mouth
    seg(2)=get_segment(data.mouth,data.chin); %mouth-chin
    seg(3)=get_segment(data.mouth,data.Trunk1)/2; %LED-trunk1
    seg(4)=get_segment(data.Trunk1,data.Trunk2); %trunk1-trunk2
    seg(5)=get_segment(data.Trunk2,data.Tail1); %trunk2-tail1
    seg(6)=get_segment(data.Tail1,data.Tail2); %tail1-tail2
    seg(7)=get_segment(data.Tail2,data.CaudalFork); %tail2-cf
    seg(8)=get_segment(data.LED,data.RPecBase); %LED-rb
    seg(9)=get_segment(data.RPecBase,data.RPecTip); %rpb-rpt
    seg(10)=get_segment(data.LED,data.LPecBase); %LED-lb
    seg(11)=get_segment(data.LPecBase,data.LPecTip); %lpb-lpt

    %% visualize
    F=figure;
    a=axes;
    x0=0;
    y0=0;
    for i=1:0;%size(fish,1)
        x0=0;
        y0=0;
        theta=pi/2;
        x=x0+[0 seg(1)*cos(theta)];
        y=y0+[0 seg(1)*sin(theta)];        
        plot(x,y);
        x0=x(2); y0=y(2);
        theta=theta+fish(i,5);
        hold on;
        x=x0+[0 seg(2)*cos(theta)];
        y=y0+[0 seg(2)*sin(theta)];        
        plot(x,y);
        
        x0=0;
        y0=0;
        theta=-pi/2;
        for j=1:5           
            x=x0+[0 seg(j+2)*cos(theta)];
            y=y0+[0 seg(j+2)*sin(theta)];        
            plot(x,y);
            x0=x(2); y0=y(2);            
            theta=theta+fish(i,j+4);
        end
        
        theta=-pi/2;%+fish(i,5);
        x0=seg(3)/2*cos(theta);
        y0=seg(3)/2*sin(theta);        
%         theta=-pi/2;
        for j=1:2
            theta=theta+fish(i,j+9);
            x=x0+[0 seg(j+7)*cos(theta)];
            y=y0+[0 seg(j+7)*sin(theta)];        
            plot(x,y);
            x0=x(2); y0=y(2);            
        end
        
        theta=-pi/2;%+fish(i,5);
        x0=seg(3)/2*cos(theta);
        y0=seg(3)/2*sin(theta);        
        for j=1:2
            theta=theta+fish(i,j+11);
            x=x0+[0 seg(j+9)*cos(theta)];
            y=y0+[0 seg(j+9)*sin(theta)];        
            plot(x,y);
            x0=x(2); y0=y(2);            
        end
        
        
        hold off;
        set(a,'Xlim',sum(seg(3:7))*[-1 1],'Ylim',[-sum(seg(3:7)) sum(seg(1:2))]);
        drawnow;
        pause(0.1);
            
        
    end
    
    function theta=get_angle(fielda,fieldb,fieldc)
        A=fielda.xy;
        B=fieldb.xy;
        C=fieldc.xy;    
        ca=fielda.c;
        cb=fieldb.c;
        cc=fieldc.c;
        AB=B-A;
        BC=C-B;
%         theta=acos(dot(AB,BC,2)./(norm2(AB).*norm2(BC)));
        theta=atan2(BC(:,2),BC(:,1))-atan2(AB(:,2),AB(:,1));
%         theta=mod(theta+pi,2*pi)-pi;
        theta(ca<th | cb<th | cc<th)=NaN;
    end

    function N=norm2(vec)
        N=hypot(vec(:,1),vec(:,2));
    end
    
    function seg=get_segment(fielda,fieldb)
        A=fielda.xy;
        B=fieldb.xy;
        ca=fielda.c;
        cb=fieldb.c;
        AB=B-A;
        N=norm2(AB);
        N(ca<th | cb<th)=NaN;
        seg=nanmean(N);
    end
        
    
end