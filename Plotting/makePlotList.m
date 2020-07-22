function success = makePlotList(pltIncl,rayGd)
% MAKEPLOT - to do later
%   

% pltIncl.density = true;
 
 
 if pltIncl.density   % with the straight line
     figure(1);
     clf

     pcolor(rayGd.zSblock,rayGd.rSblock,rayGd.neblock);
     shading interp
     colorbar
     hold on
 
     title('Log10 of electron density (1/cm3)');
     xlabel('z in um');
     ylabel('r in um');          

     %legend('str line');
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",rayGd.time))

     hold off
 end

% pltIncl.temperature = true;
 
 if pltIncl.temperature   % with the straight line
     figure(2)
     clf
     pcolor(rayGd.zSblock,rayGd.rSblock,rayGd.teblock);
     colormap('hot')
     shading interp
     colorbar
     hold on
 
     
     title('electron temperature (eV)');
     xlabel('z in um');
     ylabel('r in um');     
     
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",rayGd.time))
     
     %legend('str line');
     hold off
 end

 if pltIncl.quiverVel
     % Linearly spaced grid
     %
     nskip = 100;     % microns
     rRange = [-rayGd.rSblock(1,end):nskip:rayGd.rSblock(1,end)];
     zRange = [rayGd.zSblock(5,1):nskip:rayGd.zSblock(end,1)];
     % Query points (points to interpolate onto)
     [ZQ,RQ] = meshgrid(zRange,rRange); 
     Pq = [reshape(ZQ,[numel(ZQ),1]) reshape(RQ,[numel(RQ),1])]; 
     [ti,bc] = pointLocation(rayGd.DT,Pq);  
     
     triVals1 = rayGd.valsVz(rayGd.DT(ti,:));
     VqLin1 = dot(bc',triVals1')';

     triVals2 = rayGd.valsVr(rayGd.DT(ti,:));
     VqLin2 = dot(bc',triVals2')';
     
     hold on
     quiver(ZQ,RQ,reshape(VqLin1,size(ZQ)),reshape(VqLin2,size(ZQ)),'k') 
     xt = 600; % microns 
     yt = 1000; 
     text(xt,yt,sprintf("time: %3.2f ns",rayGd.time),'Color','r', ...
          'FontSize', 14)
     hold off
 end

 
 success = 1;
 
end

