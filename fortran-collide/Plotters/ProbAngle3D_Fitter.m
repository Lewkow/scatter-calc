
clear all
close all
clc

PLOTS_ON 	   = 1; 
WRITE_COEFFS = 0; 

f = fopen('../AngProb3D.dat'); 
d = fscanf(f,'%f %f %f',[3,inf]); 
d = d'; 
C = pi/180; 

E_A  = d(:,1); 
X_A  = d(:,2); 
Y_A  = d(:,3); 
Enow = E_A(1); 

i = 1; 
while ( E_A(i) == Enow )
	i = i+1; 	
end

i = i - 1;  

N_dp = i; 
N_En = length(X_A)/N_dp; 

fprintf('There are %d data points per energy\n', N_dp); 
fprintf('There are %d total energies\n', N_En); 

Energy = zeros(N_En); 
Z      = zeros(N_En,2,N_dp); 

k = 1; 
for E=1:N_En
	for P=1:N_dp
		Z(E,1,P) = X_A(k); 
		Z(E,2,P) = Y_A(k); 
		k = k+1; 
	end
	Energy(E) = E_A(k-2); 
end

%% Starting energy loop

for E = 1:N_En
	Enow = Energy(E); 
	fprintf('Starting fitting of energy %f eV\n', Enow); 
	
	%% clear data to be reused/resized
	clear X* Y* 
	clear Low_Lim Upp_Lim SP
 
	%% find cut indicies
	Lcut   = 0.5;
	Hcut   = 0.8;
	Lcut_i = 0;
	Hcut_i = 0;

	for i=1:N_dp
		if ( (Lcut_i == 0) && (Z(E,1,i) > Lcut))
			Lcut_i = i; 
		end
		if ( (Hcut_i == 0) && (Z(E,1,i) > Hcut))
			Hcut_i = i; 
		end
	end		 

	XL(1:Lcut_i) 				= Z(E,1,1:Lcut_i); 
	YL(1:Lcut_i) 				= Z(E,2,1:Lcut_i); 
	XM(Lcut_i+1:Hcut_i) = Z(E,1,Lcut_i+1:Hcut_i); 
	YM(Lcut_i+1:Hcut_i) = Z(E,2,Lcut_i+1:Hcut_i); 
	XH(Hcut_i+1:N_dp)   = Z(E,1,Hcut_i+1:N_dp); 
	YH(Hcut_i+1:N_dp)   = Z(E,2,Hcut_i+1:N_dp); 

	XL = XL'; 
	YL = YL'; 
	XM = XM'; 
	YM = YM'; 
	XH = XH'; 
	YH = YH'; 

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
							 'MaxIter',1000, ...
							 'MaxFunEvals',2000);
	f = fittype('a1*exp(-((x-b1)/c1)^n)+a2*exp(-((x-b2)/c2)^n)+a3*exp(-((x-b3)/c3)^n)','problem','n','options',s);

	[low_coeff,go_low] = fit(XL,YL,f,'problem',2)
	[med_coeff,go_med] = fit(XM,YM,f,'problem',2)
	[hig_coeff,go_hig] = fit(XH,YH,f,'problem',2)

	OTHER_FIT = 0; 

	if (OTHER_FIT == 1) 
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
							 'MaxIter',1000, ...
							 'MaxFunEvals',5000);
		fh = fittype('(a/(b-x))*(c1*exp(d1*x)+c2*exp(d2*x)+c3*exp(d3*x))','options',sh);
		[hig_coeff,go_hig] = fit(XH,YH,fh)
	end

	if (WRITE_COEFFS == 1)



	end

	if (PLOTS_ON == 1)

		figure
		subplot(3,1,1)
		h1 = plot(XL,YL,'k.'); 
		hold on
		h2 = plot(low_coeff,'r'); 
		set(h1,'LineWidth',2.5); 
		set(h2,'LineWidth',2.5); 
		xlabel('Probability')
		ylabel('Scattering Angle [deg]')
	
		subplot(3,1,2)
		h6 = plot(XM,YM,'k.'); 
		hold on
		h7 = plot(med_coeff,'r'); 
		set(h6,'LineWidth',2.5); 
		set(h7,'LineWidth',2.5); 
		xlabel('Probability')
		ylabel('Scattering Angle [deg]')

		subplot(3,1,3)
		h3 = plot(XH,YH,'k.'); 
		hold on
		h4 = plot(hig_coeff,'r'); 
		set(h3,'LineWidth',2.5); 
		set(h4,'LineWidth',2.5); 
		xlabel('Probability')
		ylabel('Scattering Angle [deg]')

	end 

end

