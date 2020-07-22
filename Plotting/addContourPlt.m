function success = addContourPlt(rayBundle,rayGd,conType)
% ADDBUNDLETOPLT - to do later
%   

 if ~exist('conType','var')
     conType = 'nc';
 end
     

 switch conType
   case 'nc'
     % add nc
     qcrit = log10(rayBundle.nc);
     contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[qcrit ...
                         qcrit],'k:')     
   case 'nc4'
     % add nc/4 too
     qc4 = log10(rayBundle.nc/4);
     contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[qc4 ...
                         qc4],'k:')
   case 'nc10'
     % add nc/10
     n10 = log10(rayBundle.nc/10);
     contour(rayGd.zSblock,rayGd.rSblock,rayGd.neblock,[n10 ...
                         n10],'k:')     
   otherwise
     success = 0;
 end

 success = 1;
 

end
