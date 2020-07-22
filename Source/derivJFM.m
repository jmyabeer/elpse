 function [dydx,y,pp]=derivJFM(y,x,method)
 % Take derivative of an array
 %   function [dydx,y,pp]=derivJFM(y,x,method)
 %
 % method is 'pp' | 'diff3' | 'diff5' which corresponds to
 % piecewize polynomial (pp), third order or fifth order
 % approximations
     

 if ~exist('method','var')   % Use default method
     method='diff5'
 end

 [nrows,ncols] = size(y);    % y should be a 1-D array
 
 if nrows>1 && ncols>1
   error(' FUNCTION deriv only for vectors at this time')
 elseif nrows>1
   columnvec = true;   % y is a column vector
 else
   columnvec = false;  % y is a row vector
 end

 if ~exist('x','var') || isempty(x); % check to see if x is empty
                                     % or missing
   x=1:length(y);
   if columnvec
     x = x';
   end
 elseif strcmpi(method,'diff5')  % diff5 can only be used on a
                                 % unform grid
   dx=x(2:end)-x(1:end-1);
   if abs(std(dx)/mean(dx))>1e-4 % not evenly spaced
     method='diff3';
     disp('WARNING: deriv diff5 requires evenly spaced x. Using diff3 instead')
   end
 end

 dydx=zeros(nrows,ncols);        % preallocate our derivative vector

 switch lower(method)    % method is 'pp' | 'diff3' | 'diff5' which
                         % corresponds to piecewize polynomial
                         % (pp), third order or fifth order
                         % approximations
   case 'pp'
     pp = interp1(x,y,'pchip','pp');
     dydx(:) = [pp.coefs(:,3); 3*pp.coefs(end,1)+2*pp.coefs(end, ...
                                                       2)+pp.coefs(end,3)];
     y=ppval(pp,x);      % reconstruct y (which is output)
   case 'diff3'
     pp=[];
     if columnvec        % make y and x row vectors
       y=y';
       x=x';
       dydx=dydx';       % dydx is a row vector too
     end
     dx=x(:,2:end)-x(:,1:end-1);  % written for an array but see above
     dy=y(:,2:end)-y(:,1:end-1);
     dydx2=dy./dx;
     dydx(:,2:end-1) = ...
       (dydx2(:,2:end).*dx(:,1:end-1)+dydx2(:,1:end-1).*dx(:,2:end)) ...
       ./(dx(:,1:end-1)+dx(:,2:end));
     dydx(:,1)=2*dydx2(:,1)-dydx(:,2);          % Fix up the endpoints
     dydx(:,end)=dydx2(:,end)-dydx(:,end-1);    % Fis up the endpoints
     if columnvec     % make col vec again for output
       y=y';
       dydx=dydx';
     end

   otherwise
     pp=[];
     if columnvec % make row vec
       y=y';
       x=x';
       dydx=dydx';
     end
     dx=x(2)-x(1);
     dydx(:,3:end-2) = (-y(:,5:end)+8*y(:,4:end-1)-8*y(:,2:end-3)+ ...
                        y(:,1:end-4))/(12*dx);
     dydx(:,2)=(y(:,3)-y(:,1))/(2*dx);
     dydx(:,end-1)=(y(:,end)-y(:,end-2))/(2*dx);
     dydx(:,1)=2*(y(:,2)-y(:,1))/dx - dydx(:,2);
     dydx(:,end)=2*(y(:,end)-y(:,end-1))/dx - dydx(:,end-1);
     if columnvec % make col vec again
       y=y';
       dydx=dydx';
     end
 end

 %Make sure there are no NaNs at beginning/end
 check = isnan(dydx);
 firstGoodValue = find(~check, 1 ,'first');

 if firstGoodValue~=1
    dydx(1:firstGoodValue-1) = dydx(firstGoodValue);
 end
 
 lastGoodValue = find(~check, 1 ,'last');
 
 if lastGoodValue~=length(check)
     dydx(lastGoodValue+1:end) = dydx(lastGoodValue);
 end
 
 end
