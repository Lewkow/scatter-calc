function [y] = gauss_sum(x, xx)

a1 = x(1); 
b1 = x(2); 
c1 = x(3); 
a2 = x(4); 
b2 = x(5); 
c2 = x(6); 
a3 = x(7); 
b3 = x(8); 
c3 = x(9); 
C  = x(10); 

y = a1*exp( -((xx-b1)/c1)^2 ) + ...
    a2*exp( -((xx-b2)/c2)^2 ) + ...
    a3*exp( -((xx-b3)/c3)^2 ) + C; 

y = y*180/pi; 

