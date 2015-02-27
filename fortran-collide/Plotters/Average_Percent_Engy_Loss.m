
clear all
close all

f = fopen('../AvePercentEngyLoss.dat'); 
d = fscanf(f,'%f %f %f %f',[4,inf]); 
d = d'; 

t = d(:,1); 
g = d(:,2); 
p = d(:,3); 
y = d(:,4); 

figure
h = plot(t,y,'k'); 
set(h,'LineWidth',2.5)
xlabel('CM Scattering Angle [deg]')
ylabel('Average Percent Energy Loss')

figure
h1 = plot(t,g,'g'); 
set(h1,'LineWidth',2.5)
xlabel('CM Scattering Angle [deg]')
ylabel('Percent Energy Loss')

figure
h2 = semilogy(t,p,'r'); 
set(h2,'LineWidth',2.5)
xlabel('CM Scattering Angle [deg]')
ylabel('Probability Density')


