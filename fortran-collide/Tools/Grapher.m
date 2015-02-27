clear all
close all

file = input('What file would you like to use for plots? ', 's'); 
many = input('How many plots would you like to create? '); 
col = input('How many colums are included in the file? '); 
% col = str2double(co); 

if (col == 1)
    form = '%f'; 
elseif (col == 2)
    form = '%f %f'; 
elseif (col == 3)
    form = '%f %f %f'; 
elseif (col == 4)
    form = '%f %f %f %f'; 
else
    disp('More than 4 columns is not currently supported')
end
    
fid = fopen(file); 
dat = fscanf(fid, form, [col,inf]); 
dat = dat'; 
for i=1:many
    x1 = input('What column would you like to use for the x variable? '); 
    y1 = input('What column would you like to use for the y variable? '); 
    mark = input('What symbol would you like to use for this plot? ','s'); 
    op1 = input('Would you like to state additional options for this plot? Y/N ','s');
    if (op1 == 'Y')
        disp(' '); 
        disp('Axis Options  -> 1');
        disp('Label Options -> 2');
        disp('Both Options  -> 3'); 
        op2 = input('Which options do you choose? '); 
        if (op2 == 1 || op2 == 3)
            disp(' '); 
            disp('Axis Options'); 
            disp('Scale Options  -> 1'); 
            disp('Domain Options -> 2'); 
            disp('Both Options   -> 3'); 
            axop = input('Which options do you choose? '); 
            if (axop == 2 || axop == 3)
                xm = input('xmin = '); 
                xM = input('xmax = '); 
                ym = input('ymin = '); 
                yM = input('ymax = '); 
       
            end
            if (axop == 1 || axop == 3)
                disp(' '); 
                disp('Scale Options'); 
                disp('XLog  -> 1'); 
                disp('YLog  -> 2'); 
                disp('XYLog -> 3'); 
                logop = input('Which options do you choose? '); 
            end
        end
        if (op2 == 2 || op2 == 3)
            disp(' '); 
            disp('Label Options'); 
            tit = input('Plot Title? ','s'); 
            xlab = input('X-axis label? ','s'); 
            ylab = input('Y-axis label? ','s'); 
        end
    else 
        op2 = 0;
		logop = 0; 
		axop = 0; 
    end
    figure
    if (logop == 1)
        semilogx(dat(:,x1),dat(:,y1),mark)
    elseif (logop == 2)
        semilogy(dat(:,x1),dat(:,y1),mark)
    elseif (logop ==3)
        loglog(dat(:,x1),dat(:,y1),mark)
	elseif (op2 == 0) 
        plot(dat(:,x1),dat(:,y1),mark)
    end
    if ((op2 == 1 || op2 == 3) && (axop == 2 || axop == 3))
        axis([xm xM ym yM]); 
    end
    if (op2 == 2 || op2 == 3)
        title(tit)
        xlabel(xlab)
        ylabel(ylab)
    end

    
    
end
    
    
    
