function make_video_from_imgs(imgs,vidname)
    presentpath='Z:\AW\presentations\';    
    rep=1;
%     load([presentpath,filename]);
    
    for i=1:size(imgs,1)
        i
        vidout=VideoWriter([presentpath,vidname,num2str(i)]);
        vidout.FrameRate=2.5;
        open(vidout);
        for k=1:rep
            for j=1:size(imgs,2)            
                writeVideo(vidout,imgs(i,j));
            end
            for j=(size(imgs,2)-1):-1:2
                writeVideo(vidout,imgs(i,j));
            end

        end
        close(vidout);
        F=figure;
        imshowpair(imgs(i,1).cdata,imgs(i,end).cdata);
        
    end

end
   
