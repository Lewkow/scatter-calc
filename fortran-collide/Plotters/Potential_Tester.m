
clear all
close all

fid = fopen('../Pot_Test.dat'); 
dat = fscanf(fid,'%f %f',[2,inf]); 

dat = dat'; 

R   = dat(:,1); 
V   = dat(:,2); 

figure
h = semilogy(R,V,'r.',R,V,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('R [a_{0}]'); 
ylabel('V [eV]'); 

figure
h = plot(R,V,'g.',R,V,'k','LineWidth',2.5); 
set(gca,'FontSize',16); 
xlabel('R [a_{0}]'); 
ylabel('V [eV]'); 



