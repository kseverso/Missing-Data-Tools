% Test script of AddMissingness

clear; close all; clc

%create a random dataset
x = rand(100,100);

%create datasets with all five types of missingness
miss_type = 1:5;
miss_type_labels = {'Random','Drop-out','Multi-rate','Censor','Patterned'};

%look at two levels of missingness
miss_level = [0.02, 0.2];

%use spy to see the pattern of the missing data
figure()
for i = 1:length(miss_level)
    for j = 1:length(miss_type)
        
        %note that here, the tolerance is fixed for all missing levels
        [x_m, miss_per] = AddMissingness(x,miss_type(j),miss_level(i),0.005);
        
        subplot(2,5,j + (i-1)*length(miss_type))
        spy(~isnan(x_m))
        fig_title = strcat(miss_type_labels(j),{' '},num2str(miss_per*100),'% missing');
        title(fig_title)
        
    end %end loop through missing types
    
end %end loop through missing levels