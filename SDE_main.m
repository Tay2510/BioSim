close all
figure

nrows = 2;
ncols = 4;
M = nrows * ncols;


for i = 1:M
    
    
    fprintf('trial %d\r',i);
    
    subplot(nrows,ncols,i);
    SDE
    
    
end