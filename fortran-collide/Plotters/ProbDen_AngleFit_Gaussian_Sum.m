
clear all
close all

f1 = fopen('../AngleFits.dat'); 
f2 = fopen('../ProbDen.dat'); 

d1 = fscanf(f1,'%f %f %f %f %f %f %f %f %f %f %f',[11,1]); 
d2 = fscanf(f2,'%f %f',[2,inf]); 

d1 = d1'; 
d2 = d2'; 

C = pi/180; 

Ang 	 = d2(:,1); 
PD_PD  = d2(:,2); 
dx 		 = (Ang(2)-Ang(1))*C; 

p(1:10) = d1(2:11);
%p(1:9) = d1(2:10);

a1 = p(1); 
b1 = p(2); 
c1 = p(3); 
a2 = p(4); 
b2 = p(5); 
c2 = p(6); 
a3 = p(7); 
b3 = p(8); 
c3 = p(9); 

for i=1:length(Ang)
	tot = 0; 
	for j=1:i
		tot = tot + PD_PD(j)*dx; 
	end
	AngProb(i) = tot; 
end		 

FitX = AngProb; 
for i=1:length(FitX)
	FitY(i) = (a1*exp( - ( (FitX(i) - b1)/c1)^2) + ...
             a2*exp( - ( (FitX(i) - b2)/c2)^2) + ...
             a3*exp( - ( (FitX(i) - b3)/c3)^2) + p(10))/C;  
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
xlabel('Scattering Angle [deg]'); 
ylabel('Probability Density'); 

figure
subplot(2,1,1)
semilogy(AngProb,Ang,'k.',FitX,FitY,'r')
xlabel('Probability [RND]'); 
ylabel('Scattering Angle [deg]'); 
subplot(2,1,2)
plot(AngProb,DFF,'g.')
xlabel('Probability [RND]'); 
ylabel('Percent Error'); 


max_r   = 0.5;
max_ang = 0.0;  

while (max_ang < 180.0)
	max_r = max_r + 0.0001; 
	max_ang = (p(1)*exp( - ( (max_r - p(2))/p(3))^2) + ...
             p(4)*exp( - ( (max_r - p(5))/p(6))^2) + ...
             p(7)*exp( - ( (max_r - p(8))/p(9))^2) + p(10))/C;
end

fprintf('Max r: %f\tMax Ang: %f\n', max_r, max_ang);
 
r = rand(1000,1)*max_r; 
for i=1:length(r)
	rand_y(i) = (p(1)*exp( - ( (r(i) - p(2))/p(3))^2) + ...
             	 p(4)*exp( - ( (r(i) - p(5))/p(6))^2) + ...
               p(7)*exp( - ( (r(i) - p(8))/p(9))^2) + p(10))/C;  
	if (rand_y(i) > 180.0) 
		rand_y(i) = 180.0; 
	end
end

fprintf('Random darts done\n'); 

[histy,histx] = hist(rand_y,1000); 

histy = histy/length(r); 

figure
semilogy(histx,histy,'k.',histx,histy,'r')
xlabel('Scattering Angle [deg]'); 
ylabel('Frequency'); 

fprintf('Mean angle: %f [deg]\n', mean(rand_y)); 


