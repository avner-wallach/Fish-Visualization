function [FF,ax_grid]=visualize_fish_tracking(vidname,trackfile,ss,framenum)
    th=0.8;
    sfactor=3;
    [num,txt,raw] = xlsread([trackfile,'.csv']);
%     txt=txt(2,2:end);
    [fish,seg,coornames,data]=get_posture(txt,num);
    fnames1={'chin','mouth','Trunk1','Trunk2','Tail1','Tail2','CaudalFork'};      
    fnames2={'LPecTip','midpoint','RPecTip'};
%     fnames3={'SideView'};
    vid = VideoReader([vidname,'.avi']);
%     vidout=VideoWriter('vidout1');
%     vidout.FrameRate=25;
%     open(vidout);
    
    currentFrame=ss;
    F=figure;
    set(F,'Color',[0 0 0],'Position',[1 41 1920 964]);
    COL=colormap('lines');
    A1=axes;
%     set(A1,'Position',[0.13 0.11 0.5 0.815]);
%     A2=axes;
%     set(A2,'Position',[0.73 0.11 0.2 0.815]);
    for i=1:framenum
        axes(A1);
        frame = rgb2gray(read(vid, currentFrame));
        ix=[1:size(frame,2)];
        iy=[1:size(frame,1)];
        ax_grid={ix,iy};   
        imagesc(ix,iy,repmat((frame),[1 1 3]));
%         im = image(frame, 'CDataMapping', 'scaled');     colormap('gray');
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
        set(H,'Marker','.','Linewidth',2,'MarkerSize',12,'Color',COL(2,:));    
        axis('image');
        set(gca,'Xlim',[ix(1) ix(end)],'Ylim',[min(iy) max(iy)]);    
        set(gca,'XTick',[],'YTick',[]);        
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
%         axes(A2);
%         draw_skeleton();

        drawnow;
        pause(0.01);
        FF(i)=getframe;
        
%         writeVideo(vidout,FF(i));
        currentFrame=currentFrame+1;        
    end
%     close(vidout);
%     save('tracked_frames','FF');
    
    function draw_skeleton()
        x0=0;
        y0=0;
        theta=pi/2;
        x=x0+[0 seg(1)*cos(theta)];
        y=y0+[0 seg(1)*sin(theta)];        
        H=plot(x,y); set(H,'Color',COL(1,:),'Marker','.','MarkerSize',14,'LineWidth',2);
        x0=x(2); y0=y(2);
        theta=theta+fish(currentFrame,4);
        hold on;
        x=x0+[0 seg(2)*cos(theta)];
        y=y0+[0 seg(2)*sin(theta)];        
        H=plot(x,y); set(H,'Color',COL(3,:),'Marker','.','MarkerSize',14,'LineWidth',2);
        
        x0=0;
        y0=0;
        theta=-pi/2;        
        x=x0+[0 seg(3)*cos(theta)];
        y=y0+[0 seg(3)*sin(theta)];        
        H=plot(x,y); set(H,'Color',COL(1,:),'Marker','.','MarkerSize',14,'LineWidth',2);
        x0=x(2); y0=y(2);            
        for j=1:4           
            theta=theta+fish(currentFrame,j+4);
            x=x0+[0 seg(j+3)*cos(theta)];
            y=y0+[0 seg(j+3)*sin(theta)];        
            H=plot(x,y); set(H,'Color',COL(1,:),'Marker','.','MarkerSize',14,'LineWidth',2);
            x0=x(2); y0=y(2);                        
        end
        
%         theta=0;
%         x0=seg(3)/2*cos(theta);
%         y0=seg(3)/2*sin(theta);        
        x0=0;
        y0=0;                
        j=8
%         theta=-pi/2;
%         for j=1:2
            theta=pi+fish(currentFrame,j+1);
            x=x0+[0 seg(j)*cos(theta)];
            y=y0+[0 seg(j)*sin(theta)];        
            H=plot(x,y); set(H,'Color',COL(2,:),'Marker','.','MarkerSize',14,'LineWidth',2);
%             x0=x(2); y0=y(2);            
%         end
        
%         theta=0;%+fish(i,5);
        j=9;
%         x0=seg(3)/2*cos(theta);
%         y0=seg(3)/2*sin(theta);        
%         for j=1:2
            theta=-fish(currentFrame,j+1);
            x=x0+[0 seg(j)*cos(theta)];
            y=y0+[0 seg(j)*sin(theta)];        
            H=plot(x,y); set(H,'Color',COL(2,:),'Marker','.','MarkerSize',14,'LineWidth',2);
%             x0=x(2); y0=y(2);            
%         end
        hold off;
        set(A2,'Xlim',sum(seg(3:7))*[-.55 .55],'Ylim',[-sum(seg(3:7)) sum(seg(1:2))]*1.1,'XDir','reverse');
        set(A2,'XTick',[],'YTick',[])
    end

end

                

