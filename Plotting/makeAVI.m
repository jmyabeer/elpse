

%  profile = 'MPEG-4'; % not available
  
  vidObj = VideoWriter('ramanBackScatter','Uncompressed AVI');

  vidObj.FrameRate = 6;
  
  open(vidObj);
  
  [~,nframes] = size(Mraman);
  
  for i=1:nframes
      currFrame = Mraman(i);
      writeVideo(vidObj,currFrame);
  end
  
  close(vidObj);