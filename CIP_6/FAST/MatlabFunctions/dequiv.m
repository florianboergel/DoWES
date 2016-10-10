function [Req, R] = dequiv(time, data, m)

%
% dequiv.m
%  MATLAB function to determine the damage equivalent loads (Req) from a loads data array (data), 
%    a time vector (time), and Wohler curve exponent for material (m)
%
%    Uses a rainflow counting algorithm published in 'Simple rainflow counting algorithm'
%        by Downing, S.D. and Socie, D.F..  The same algorithm is used in CRUNCH and LIFE2
%
% K. Stol
% 4/22/02: modified from P. Moriarty's MATLAB code
%            Corrected for the case when the first and end points have the same value
%
% for composite material, m ~= 10 to 12
% fiberglass: m = 8
% iron: m=16.8
% steel: m=6.5
% structural steel: m=3

% downloaded by Luis the 25th Jan 2014 
% http://wind.nrel.gov/public/Awright/CART2-data/General%20data/CART_Matlab_Base_Scripts/Utilities/dequiv.m

% Identify peaks/valleys and remove plateaus

peaks = [];
data = [data(end); data; data(1)];    % add end point to start and first point to end to test boundaries
for i = 2:(length(data)-1)   % test interior points
    if (data(i)-data(i-1))*(data(i+1)-data(i)) < 0      % if slope change before and after point
        peaks = [peaks, data(i)];                 %  then point is a peak/valley
    end
end

%' first iteration of rainflow counting' - once through data
% count cycles
start = 1;
i = 1;
k = 1;       % k keeps track of length(v)
R = [];
while (i <= length(peaks))
   %%disp(['start = ',num2str(start)]);
   %disp(['i,k = ',num2str(i), ' ', num2str(k)]);
   if (k <(start+2));
       v(k)=peaks(i);
       %disp(['v = [',num2str(v),']']);
       k = k+ 1;
       i = i+1;
   else   
      v(k)=peaks(i);
      %disp(['v = [',num2str(v),']']);
      Y = abs(v(k-1)-v(k-2));
      X = abs(v(k)-v(k-1));
      %disp(['X = ', num2str(X)]);
      %disp(['Y = ', num2str(Y)]);
      if X < Y
      	k = k + 1;
      	i = i + 1;
      elseif (X == Y)&(k <= (start+2)) %((start==(k-1))|(start==(k-2))))
      	k = k + 1;
      	i = i + 1;
      elseif (X > Y)&(k <= (start+2)) %&((start==(k-1))|(start==(k-2))))
          start = start + 1;
          k = k + 1;
      	i = i + 1;
      elseif (X >= Y)&(k > (start+2)) %&(~((start==(k-1))|(start==(k-2)))))
          %ncycles = ncycles + 1;
          R = [R, Y];
          %disp(['R = [',num2str(R),']']);
          %mcycle(ncycles) = (v(k-1) - v(k-2))/2;
      	%mmean(ncycles) = (v(k-1) + v(k-2))/2;
      	%nfail(ncycles) = a*((momentultimate-mmean(ncycles))/mcycle(ncycles))^b;
      	%lifeused(ncycles) = 1/nfail(ncycles);
          v = [v(1:k-3), v(k)];
          k = k - 2;
      end
   end
end

%'second iteration of rainflow counting' - recycle data
j = 1;
while (j <= start)
   %disp(['j,k = ',num2str(j), ' ', num2str(k)]);
   if (k < (start+2))
      v(k) = v(j);
      k = k + 1;
      j = j + 1;
   else
      v(k) = v(j);
      %disp(['v = [',num2str(v),']']);
      Y = abs(v(k-1)-v(k-2));
      X = abs(v(k)-v(k-1));
      %disp(['X = ', num2str(X)]);
      %disp(['Y = ', num2str(Y)]);
      if X < Y
          k = k + 1;
          j = j + 1;
      elseif X >= Y
          %ncycles = ncycles + 1;
          R = [R, Y];
          %disp(['R = [',num2str(R),']']);
      	%mcycle(ncycles) = (v(k-1)-v(k-2))/2;
      	%mmean(ncycles) = (v(k-1)+v(k-2))/2;
      	%nfail(ncycles)=a*((momentultimate-mmean(ncycles))/mcycle(ncycles))^b;
      	%lifeused(ncycles)=1/nfail(ncycles);
      	v = [v(1:k-3), v(k)];
          k = k - 2;
       end
   end
end  

% count remaining half cycles  - !! can't find a case where this is used!
halfcycles = 0;
if (length(v) > start)
   for i= start:(length(v)-1)
      %ncycles=ncycles + 1;
      halfcycles=halfcycles+1;
      R = [R, abs(v(i+1)-v(i))];
      %mcycle(ncycles) = (v(i+1)-v(i))/2;
      %mmean(ncycles) = (v(i+1)+v(i))/2;
      %if(mcycle(ncycles) ~= 0)
      %   fprintf(1,' value close to zero (?) : mcycle(ncycle):: %e - [%d] \n',mcycle(ncycles),ncycles);
      %end
      %nfail(ncycles)=a*((momentultimate-mmean(ncycles))/mcycle(ncycles))^b;
      %lifeused(ncycles)=0.5/nfail(ncycles);
   end
end

% Calculate damage equivalent load at 1 Hz

Req = (sum(R.^m)/time(end))^(1/m);