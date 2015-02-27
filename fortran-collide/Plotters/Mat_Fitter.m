
clear all
close all
clc

f2 = fopen('../ProbDen.dat'); 

d2 = fscanf(f2,'%f %f',[2,inf]); 

d2 = d2'; 

C = pi/180; 

Y    	 = d2(:,1); 
PD_PD  = d2(:,2); 
dx 		 = (Y(2)-Y(1))*C; 

Lcut = 0; 
Hcut = 0; 

for i=1:length(Y)
	tot = 0; 
	for j=1:i
		tot = tot + PD_PD(j)*dx; 
	end
	X(i,1) = tot; 
	if ( (Lcut == 0) && (X(i,1) > 0.5))
		Lcut = i; 
	end
	if ( (Hcut == 0) && (X(i,1) > 0.8))
		Hcut = i; 
	end
end		 

XL = X(1:Lcut); 
YL = Y(1:Lcut); 
XM = X(Lcut+1:Hcut); 
YM = Y(Lcut+1:Hcut); 
XH = X(Hcut+1:length(X)); 
YH = Y(Hcut+1:length(Y)); 

N = 1e12; 
np = 9; 

Low_Lim(1:np) = -N; 
Upp_Lim(1:np) = N;
SP(1:np)      = 1.0;  

s = fitoptions('Method','NonlinearLeastSquares', ...
							 'Algorithm','Trust-Region', ...
               'Lower',Low_Lim, ...
               'Upper',Upp_Lim, ...
               'Startpoint',SP, ...
               'Display','iter', ...
							 'MaxIter',1000, ...
							 'MaxFunEvals',2000);
f = fittype('a1*exp(-((x-b1)/c1)^n)+a2*exp(-((x-b2)/c2)^n)+a3*exp(-((x-b3)/c3)^n)','problem','n','options',s);

[low_coeff,go_low] = fit(XL,YL,f,'problem',2)

[med_coeff,go_med] = fit(XM,YM,f,'problem',2)

np = 8;
clear SP
clear Low_Lim; 
clear Upp_Lim;  
SP(1:np) = 1.0;
Low_Lim(1:np) = -N; 
Upp_Lim(1:np) = N;
sh = fitoptions('Method','NonlinearLeastSquares', ...
							 'Algorithm','Trust-Region', ...
               'Lower',Low_Lim, ...
               'Upper',Upp_Lim, ...
               'Startpoint',SP, ...
               'Display','iter', ...
							 'MaxIter',1000, ...
							 'MaxFunEvals',5000);
fh = fittype('(a/(b-x))*(c1*exp(d1*x)+c2*exp(d2*x)+c3*exp(d3*x))','options',sh);

[hig_coeff,go_hig] = fit(XH,YH,fh)


figure
h1 = plot(XL,YL,'k.'); 
hold on
h2 = plot(low_coeff,'r'); 
set(h1,'LineWidth',2.5); 
set(h2,'LineWidth',2.5); 
xlabel('Probability')
ylabel('Scattering Angle [deg]')

figure
h6 = plot(XM,YM,'k.'); 
hold on
h7 = plot(med_coeff,'r'); 
set(h6,'LineWidth',2.5); 
set(h7,'LineWidth',2.5); 
xlabel('Probability')
ylabel('Scattering Angle [deg]')

figure
h3 = plot(XH,YH,'k.'); 
hold on
h4 = plot(hig_coeff,'r'); 
set(h3,'LineWidth',2.5); 
set(h4,'LineWidth',2.5); 
xlabel('Probability')
ylabel('Scattering Angle [deg]')



