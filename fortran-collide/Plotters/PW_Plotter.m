
clear all
close all

f = fopen('../AngleFits.dat'); 
d = fscanf(f,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[21,inf]); 

d = d'; 

C = 180/pi; 

pL(1:10) = d(2:11); 
pH(1:10) = d(12:21); 

ff = fopen('../ProbDen.dat'); 
dd = fscanf(ff,'%f %f',[2,inf]); 
dd = dd'; 

Ang = dd(:,1); 
PD  = dd(:,2); 
dx  = Ang(2) - Ang(1); 
dx  = dx/C; 

for i=1:length(Ang)
	tot = 0; 
	for j=1:i
		tot = tot + PD(j)*dx; 
	end
	SumPD(i) = tot; 
end	

r = 0:0.001:1; 

for i=1:length(r)
	yL(i) = pL(1)*exp( -( (r(i)-pL(2))/pL(3))^2) + ...
          pL(4)*exp( -( (r(i)-pL(5))/pL(6))^2) + ...
          pL(7)*exp( -( (r(i)-pL(8))/pL(9))^2) + pL(10);  

	yH(i) = pH(1)*exp( -( (r(i)-pH(2))/pH(3))^2) + ...
          pH(4)*exp( -( (r(i)-pH(5))/pH(6))^2) + ...
          pH(7)*exp( -( (r(i)-pH(8))/pH(9))^2) + pH(10);  
end

figure
plot(r,yL*C,'k',r,yH*C,'b',SumPD,Ang,'r.')
axis([0 1 0 180])


