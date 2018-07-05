function make_video_from_imgs(imgs)
    presentpath='Z:\AW\presentations\';    
    rep=10;
%     load([presentpath,filename]);
    
    for i=1:size(imgs,1)
        vidout=VideoWriter([presentpath,'vidout',num2str(i)]);
        vidout.FrameRate=1;
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
   
