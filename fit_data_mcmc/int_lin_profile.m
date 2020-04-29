function [integral] = int_lin_profile(x,y)

% Given a series  of points {x} and {y} that define the vertices of a
% piecewise linear profile, this function integrates the profile in
% spherical coordinates.

xk1 = x(1:end-1);
xk2 = x(2:end);
yk1 = y(1:end-1);
yk2 = y(2:end);
slopes = (yk2 - yk1)./(xk2 - xk1);

int_each = yk1.*(xk2.^3 - xk1.^3)/3 + slopes.*(0.25*xk2.^4 - (xk1/3).*xk2.^3 - 0.25*xk1.^4 + (1/3)*xk1.^4);

tot = sum(int_each);

integral = 4*pi*tot;

end

