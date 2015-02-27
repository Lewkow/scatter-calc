
clear all
close all

f1 = fopen('../AngleFits.dat'); 
f2 = fopen('../ProbDen.dat'); 

d1 = fscanf(f1,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[21,1]); 
d2 = fscanf(f2,'%f %f',[2,inf]); 

d1 = d1'; 
d2 = d2'; 

C = pi/180; 

Ang 	 = d2(:,1); 
PD_PD  = d2(:,2); 
dx 		 = (Ang(2)-Ang(1))*C; 

pL(1:10) = d1(2:11);
pH(1:10) = d1(12:21); 

for i=1:length(Ang)
	tot = 0; 
	for j=1:i
		tot = tot + PD_PD(j)*dx; 
	end
	AngProb(i) = tot; 
end		 

FitX = AngProb; 
for i=1:length(FitX)
	if (FitX(i) < 0.8)
		FitY(i) = (pL(1)*exp( - ( (FitX(i) - pL(2))/pL(3))^2) + ...
               pL(4)*exp( - ( (FitX(i) - pL(5))/pL(6))^2) + ...
               pL(7)*exp( - ( (FitX(i) - pL(8))/pL(9))^2) +pL(10))/C;  
	else
		FitY(i) = (pH(1)*exp( - ( (FitX(i) - pH(2))/pH(3))^2) + ...
               pH(4)*exp( - ( (FitX(i) - pH(5))/pH(6))^2) + ...
               pH(7)*exp( - ( (FitX(i) - pH(8))/pH(9))^2) +pH(10))/C;  
	end
	if (Ang(i) == 0) 
		DFF(i) = abs(FitY(i)-Ang(i))*100/FitY(i); 
	else
		DFF(i) = abs(FitY(i)-Ang(i))*100/Ang(i); 
	end
end

T_DFF = sum(DFF)/length(FitX); 
fprintf('Mean Percent Difference %f\n', T_DFF); 


figure
semilogx(Ang,PD_PD,'g.',Ang,PD_PD,'k')

figure
subplot(2,1,1)
semilogy(AngProb,Ang,'k.',FitX,FitY,'r')
subplot(2,1,2)
plot(AngProb,DFF,'g.')

max_r   = 0.5;
max_ang = 0.0;  

while (max_ang < 180.0)
	max_r = max_r + 0.0001; 
	if (max_r < 0.8)
		max_ang = (pL(1)*exp( - ( (max_r - pL(2))/pL(3))^2) + ...
               pL(4)*exp( - ( (max_r - pL(5))/pL(6))^2) + ...
 	             pL(7)*exp( - ( (max_r - pL(8))/pL(9))^2) +pL(10))/C;
	else
		max_ang = (pH(1)*exp( - ( (max_r - pH(2))/pH(3))^2) + ...
               pH(4)*exp( - ( (max_r - pH(5))/pH(6))^2) + ...
 	             pH(7)*exp( - ( (max_r - pH(8))/pH(9))^2) +pH(10))/C;
	end
end

r = rand(100000,1)*max_r; 
for i=1:length(r)
	if (r(i) < 0.8)
		rand_y(i) = (pL(1)*exp( - ( (r(i) - pL(2))/pL(3))^2) + ...
 	            	 pL(4)*exp( - ( (r(i) - pL(5))/pL(6))^2) + ...
 	               pL(7)*exp( - ( (r(i) - pL(8))/pL(9))^2) +pL(10))/C;  
	else
		rand_y(i) = (pH(1)*exp( - ( (r(i) - pH(2))/pH(3))^2) + ...
 	            	 pH(4)*exp( - ( (r(i) - pH(5))/pH(6))^2) + ...
 	               pH(7)*exp( - ( (r(i) - pH(8))/pH(9))^2) +pH(10))/C;  
	end
	if (rand_y(i) > 180.0) 
		rand_y(i) = 180.0; 
	end
end

[histy,histx] = hist(rand_y,1000); 

histy = histy/length(r); 

figure
semilogy(histx,histy,'k.',histx,histy,'r')

fprintf('Mean angle: %f [deg]\n', mean(rand_y)); 



