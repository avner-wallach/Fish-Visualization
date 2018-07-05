function [x_,y_]=allo2ego(x,y,a,x0,y0);
% convert object position from allocentric to egocentric coordinates
%   x- vector, fish x head position 
%   y- vector, fish y head position
%   a- vector, fish azimouth relative to negative y direction
%   x0- scalar, object x position
%   y0- scalar, object y position

a=a+pi; %rotate angle 180 deg (toward positive y direction)

%translate 
x1=x-x0;
y1=y-y0;

%rotate
x_=(sin(a).*x1-cos(a).*y1);
y_=-(cos(a).*x1+sin(a).*y1);
end
