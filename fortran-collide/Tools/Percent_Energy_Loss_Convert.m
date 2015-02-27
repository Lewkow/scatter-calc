clear all
close all

fprintf('I convert data files in the form (<dE>/E, E) to the form (<dE>, E)\n'); 
file  = input('What file would you like to convert? ', 's'); 
ofile = input('What should I name the output file? ', 's'); 

fid  = fopen(file);
ofid = fopen(ofile,'w');  
dat  = fscanf(fid,'%f %f',[2,inf]); 
dat  = dat'; 

xi = dat(:,1); 
yi = dat(:,2); 

xf = xi; 
for i=1:length(xi)
	yf(i) = yi(i)*xi(i); 
end

for i=1:length(xi)
	fprintf(ofid,'%e %e\n', xf(i), yf(i)); 
end

close all

fprintf('\nAll Done!!\n\n'); 

 
