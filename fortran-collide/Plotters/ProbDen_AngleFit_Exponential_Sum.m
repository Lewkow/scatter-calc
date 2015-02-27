
clear all
close all

f1 = fopen('../AngleFits.dat'); 
f2 = fopen('../ProbDen.dat'); 

d1 = fscanf(f1,'%f %f %f %f %f',[5,1]); 
d2 = fscanf(f2,'%f %f',[2,inf]); 

d1 = d1'; 
d2 = d2'; 

C = pi/180; 

Ang 	 = d2(:,1); 
PD_PD  = d2(:,2); 
dx 		 = (Ang(2)-Ang(1))*C; 

p(1:4) = d1(2:5);

for i=1:length(Ang)
	tot = 0; 
	for j=1:i
		tot = tot + PD_PD(j)*dx; 
	end
	AngProb(i) = tot; 
end		 

FitX = AngProb; 
for i=1:length(FitX)
	FitY(i) = (p(1)*exp(p(2)*FitX(i)) + p(3)*exp(p(4)*FitX(i)))/C; 
	DFF(i)  = abs(FitY(i)-Ang(i)); 
end

T_DFF = sum(DFF)/length(FitX); 
fprintf('Mean Percent Difference %f\n', T_DFF); 


figure
plot(Ang,PD_PD,'k.')

figure
subplot(2,1,1)
plot(AngProb,Ang,'k.',FitX,FitY,'r')
subplot(2,1,2)
plot(AngProb,DFF,'g.')


