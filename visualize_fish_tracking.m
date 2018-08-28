function visualize_fish_tracking(vidname,trackfile,ss,framenum)
    th=0.7;
    sfactor=2;
    [num,txt,raw] = xlsread([trackfile,'.csv']);
    txt=txt(2,2:end);
    [fish,seg,data]=skeleton_model(txt,num);
    fnames1={'chin','mouth','Trunk1','Trunk2','Tail1','Tail2','CaudalFork'};      
    fnames2={'LPecTip','LPecBase','midpoint','RPecBase','RPecTip'};
%     fnames3={'SideView'};
    vid = VideoReader([vidname,'.avi']);
    vidout=VideoWriter('vidout1');
    vidout.FrameRate=25;
    open(vidout);
    
    currentFrame=ss;
    F=figure;
    COL=colormap('lines');
    A1=axes;
    set(A1,'Position',[0.13 0.11 0.5 0.815]);
    A2=axes;
    set(A2,'Position',[0.73 0.11 0.2 0.815]);
    for i=1:framenum
        axes(A1);
        frame = rgb2gray(read(vid, currentFrame));
        im = image(frame, 'CDataMapping', 'scaled');     colormap('gray');
        hold on;
        X=[]; Y=[];        
        for j=1:numel(fnames1)
            xy= data.(fnames1{j}).xy(currentFrame,:);
            c= data.(fnames1{j}).c(currentFrame);
            if(c>th)
              X=[X;xy(1)];
              Y=[Y;xy(2)];
            else
              X=[X;NaN];
              Y=[Y;NaN];
            end
        end
        H=plot(X*sfactor,Y*sfactor);        
        set(H,'Marker','.','Linewidth',2,'MarkerSize',12);    

        X=[]; Y=[];        
        for j=1:numel(fnames2)
            xy= data.(fnames2{j}).xy(currentFrame,:);
            c= data.(fnames2{j}).c(currentFrame);
            if(c>th)
              X=[X;xy(1)];
              Y=[Y;xy(2)];
            else
              X=[X;NaN];
              Y=[Y;NaN];
            end
        end

        H=plot(X*sfactor,Y*sfactor);        
        set(H,'Marker','.','Linewidth',2,'MarkerSize',12);    
        set(A1,'XTick',[],'YTick',[])
        
%         X=[]; Y=[];        
%         for j=1:numel(fnames3)
%             xy= data.(fnames3{j}).xy(currentFrame,:);
%             c= data.(fnames3{j}).c(currentFrame);
%             if(c>th)
%               X=[X;xy(1)];
%               Y=[Y;xy(2)];
%             else
%               X=[X;NaN];
%               Y=[Y;NaN];
%             end
%         end
%         
%         H=plot(X,Y);        
%         set(H,'Marker','.','Linewidth',2,'MarkerSize',20);    
%         set(A1,'XTick',[],'YTick',[])

        hold off;
        axes(A2);
        draw_skeleton();

        drawnow;
        currentFrame=currentFrame+1;
        pause(0.01);
        FF=getframe(F);
        writeVideo(vidout,FF);
    end
    close(vidout);
    
    function draw_skeleton()
        x0=0;
        y0=0;
        theta=pi/2;
        x=x0+[0 seg(1)*cos(theta)];
        y=y0+[0 seg(1)*sin(theta)];        
        H=plot(x,y); set(H,'Color',COL(1,:),'Marker','.','MarkerSize',14,'LineWidth',2);
        x0=x(2); y0=y(2);
        theta=theta+fish(currentFrame,5);
        hold on;
        x=x0+[0 seg(2)*cos(theta)];
        y=y0+[0 seg(2)*sin(theta)];        
        H=plot(x,y); set(H,'Color',COL(3,:),'Marker','.','MarkerSize',14,'LineWidth',2);
        
        x0=0;
        y0=0;
        theta=-pi/2;        
        x=x0+[0 seg(j+2)*cos(theta)];
        y=y0+[0 seg(j+2)*sin(theta)];        
        H=plot(x,y); set(H,'Color',COL(1,:),'Marker','.','MarkerSize',14,'LineWidth',2);
        x0=x(2); y0=y(2);            
        for j=1:5           
            theta=theta+fish(currentFrame,j+4);
            x=x0+[0 seg(j+2)*cos(theta)];
            y=y0+[0 seg(j+2)*sin(theta)];        
            H=plot(x,y); set(H,'Color',COL(1,:),'Marker','.','MarkerSize',14,'LineWidth',2);
            x0=x(2); y0=y(2);                        
        end
        
        theta=-pi/2;
        x0=seg(3)/2*cos(theta);
        y0=seg(3)/2*sin(theta);        
%         theta=-pi/2;
        for j=1:2
            theta=theta+fish(currentFrame,j+9);
            x=x0+[0 seg(j+7)*cos(theta)];
            y=y0+[0 seg(j+7)*sin(theta)];        
            H=plot(x,y); set(H,'Color',COL(2,:),'Marker','.','MarkerSize',14,'LineWidth',2);
            x0=x(2); y0=y(2);            
        end
        
        theta=-pi/2;%+fish(i,5);
        x0=seg(3)/2*cos(theta);
        y0=seg(3)/2*sin(theta);        
        for j=1:2
            theta=theta+fish(currentFrame,j+11);
            x=x0+[0 seg(j+9)*cos(theta)];
            y=y0+[0 seg(j+9)*sin(theta)];        
            H=plot(x,y); set(H,'Color',COL(2,:),'Marker','.','MarkerSize',14,'LineWidth',2);
            x0=x(2); y0=y(2);            
        end
        hold off;
        set(A2,'Xlim',sum(seg(3:7))*[-.55 .55],'Ylim',[-sum(seg(3:7)) sum(seg(1:2))]*1.5,'XDir','reverse');
        set(A2,'XTick',[],'YTick',[])
    end

end

                

