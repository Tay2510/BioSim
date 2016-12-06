close all

nrows = 2;
ncols = 4;
M = nrows * ncols;

SDE_models = {@SDE, @SDE_constant_dW, @SDE_hw_dW, @SDE_test_dW};

for h = 1:length(SDE_models)
    
    figure
    SDE_model = SDE_models(h);
    
    fprintf('model %d of %d\r',h,length(SDE_models));
    
    for i = 1:M

        fprintf('\ttrial %d of %d\r',i,M);

        subplot(nrows,ncols,i);
        SDE

    end
    
    
end