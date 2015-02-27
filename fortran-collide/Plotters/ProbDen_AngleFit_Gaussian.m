
clear all
close all

f1 = fopen('../AngleFits.dat'); 
f2 = fopen('../ProbDen.dat'); 

d1 = fscanf(f1,'%f %f %f %f',[4,1]); 
d2 = fscanf(f2,'%f %f',[2,inf]); 

d1 = d1'; 
d2 = d2'; 

C = pi/180; 

Ang 	 = d2(:,1); 
PD_PD  = d2(:,2); 
dx 		 = (Ang(2)-Ang(1))*C; 

p(1:3) = d1(2:4);

for i=1:length(Ang)
	tot = 0; 
	for j=1:i
		tot = tot + PD_PD(j)*dx; 
	end
	AngProb(i) = tot; 
end		 

FitX = AngProb; 
for i=1:length(FitX)
	FitY(i) = (p(1)*exp( - ( (FitX(i) - p(2))/p(3))^2))/C;  
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
	max_ang = (p(1)*exp( - ( (max_r - p(2))/p(3))^2))/C;  
end

r = rand(100000,1)*max_r; 
for i=1:length(r)
	rand_y(i) = (p(1)*exp( - ( (r(i) - p(2))/p(3))^2))/C;  
	if (rand_y(i) > 180.0) 
		rand_y(i) = 180.0; 
	end
end

[histy,histx] = hist(rand_y,1000); 

histy = histy/length(r); 

figure
semilogy(histx,histy,'k.',histx,histy,'r')

fprintf('Mean angle: %f [deg]\n', mean(rand_y)); 


