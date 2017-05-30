function [ x_miss, miss_per ] = AddMissingness( x, miss_type, level , tol)
%Function to add missingness at a specfied level to the input dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Inputs                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - data set organized such that the rows are samples (or timepoints) and
% the columns are measurements; x may already contain missing elements,
% which should be represented by NaN
% miss_type - an integer specifying the desired type of missingness
%   1 - missing completely at random
%   2 - sensor drop-out: missing correlated in the rows
%   3 - multi-rate: missing due to different sampling rates
%   4 - not missing at random: missing due to particular values
%   5 - patterned: similar to sensor where missing is correlated spatially
%   but less random
% level - percent of missing data, specified by a decimal
% tol - tolerance to specify an acceptble deviation from the desired level
% N.B. there are several factors that impact missingness that are not
% variable in the below implementation. E.g. the maximum number of time
% points for sensor drop-out missingness is 10. Hopefully these factors are
% obvious and you can change them to suit your needs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Outputs                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_miss - the data set with added NaNs to achieve the missing type and
% level
% miss_per - the final result for the level
%
%
% For more information concerning types of missing data, please see
% Rubin, D.B. Inference and missing data. Biometrika 1976, 63, 581-592.
% Little, R.J.A. and Rubin, D.B. Statistical Analysis with Missing Data.
% 2nd ed. John Wiley & Sons, 2014.
%
%
% If you use this software, please cite
% Severson, K.A., Molaro, M.C., and Braatz, R.D. Methods for applying 
% principal component analysis to process datasets wiht missing values.
% Processes 2017, vol, pages.


[n,p] = size(x);
max_iter = 1000; %set maximum number of iterations
iter = 0; %initialize iteration count

%set missingness
if miss_type == 1 %mcar
    x_miss = x;
    base_miss = sum(sum(isnan(x)))/(n*p); %check if the dataset already has missing values
    while (abs(sum(sum(isnan(x_miss)))/(n*p) - level) > tol) && (iter < max_iter)
        x_miss = x;
        wx = random('unif',0,1,size(x)) < level - base_miss;
        x_miss(wx) = NaN;
        iter = iter + 1;
    end
    miss_per = sum(sum(isnan(x_miss)))/(n*p);
end

if miss_type == 2 %sensor drop-out
    x_miss = x;
    c_level = floor(p/4); %number of measurements with missingness
    l_level = 10; %max length of time for missingness
    while (abs(sum(sum(isnan(x_miss)))/(n*p) - level) > tol) && (iter < max_iter)
        x_miss = x;
        cx = randi(p,[c_level,1]); %choose measurements (cols) to have missing measurements
        rx = randi(n,[c_level,1]); %choose time index where missingness starts
        lx = randi(l_level,[c_level,1]); %choose length of time missingness occurs
        for i = 1:c_level
            if rx(i) + lx(i) > n
                x_miss(rx(i):end,cx(i)) = NaN;
            else
                x_miss(rx(i):rx(i) + lx(i),cx(i)) = NaN;
            end
        end;
        if sum(sum(isnan(x_miss)))/(n*p) - level > 0 %too much missing data
            c_level = c_level - 1; %decrease the number of measurements with missingness
        else
            c_level = c_level + 1; %increase the number of measurements with missingness
        end
        iter = iter + 1;
    end
    miss_per = sum(sum(isnan(x_miss)))/(n*p);
end

if miss_type == 3 %multi-rate
    % check if possible
    min_miss = 1/(2*p);
    if min_miss - level > tol
        x_miss = x;
        c_level = 1;
        rx = 2;
        cx = randi(p,[c_level,1]);
        for i = 1:c_level
            x_miss(:,cx(i)) = NaN;
            x_miss(1:rx(i):n,cx(i)) = x(1:rx(i):n,cx(i));
        end
        fprintf('Desired tolerance not possible, missing level set to %0.2f \n',min_miss)
    else
        x_miss = x;
        c_level = floor(p/4); %number of measurements with missingess
        r_level = 5; %max level of subsampling
        while (abs(sum(sum(isnan(x_miss)))/(n*p) - level) > tol) && (iter < max_iter)
            x_miss = x;
            cx = randi(p,[c_level,1]);
            rx = randi(r_level,[c_level,1]);
            for i = 1:c_level
                x_miss(:,cx(i)) = NaN;
                x_miss(1:rx(i):n,cx(i)) = x(1:rx(i):n,cx(i));
            end;
            if sum(sum(isnan(x_miss)))/(n*p) - level > 0 %too much missing data
                c_level = c_level - 1; %decrease number of measurements wiht missingness
            else
                c_level = c_level + 1; %increase number of measurements with missingness
            end
            iter = iter + 1;
        end
    end %end test
    miss_per = sum(sum(isnan(x_miss)))/(n*p);
end

if miss_type == 4 %nmar
    x_miss = x;
    min_c_level = floor(p*level) + 1;
    max_c_level = min(min_c_level + floor(p/4),p);
    c_level = randi([min_c_level,max_c_level],1); %number of measurements with missingness
    meas = randperm(p,c_level);
    std_start = 1.5*nanstd(x(:,meas)); %initialize the threshold as 1.5 standard deviations
    updown = sign(-1 + 2.*rand(1,length(meas)));
    thres = nanmean(x(:,meas)) + std_start.*updown; %randomize up and down thresholds
    while (abs(sum(sum(isnan(x_miss)))/(n*p) - level) > tol) && (iter < max_iter)
        x_miss = x;
        for i = 1:length(meas)
            if updown(i) < 0
                x_miss(x(:,meas(i)) < thres(i),meas(i)) = NaN;
            else
                x_miss(x(:,meas(i)) > thres(i),meas(i)) = NaN;
            end
        end
        if sum(sum(isnan(x_miss)))/(n*p) - level > 0 %too much missing data
            thres = thres + 0.01.*std_start.*updown; 
        else
            thres = thres - 0.05.*std_start.*updown;
        end
        iter = iter + 1;
    end
    miss_per = sum(sum(isnan(x_miss)))/(n*p);
end

if miss_type == 5 %mar
    x_miss = x;
    r_level = 0.25;
    c_level = 0.55;
    while (abs(sum(sum(isnan(x_miss)))/(n*p) - level) > tol) && (iter < max_iter)
        x_miss = x;
        rx = random('unif',0,1,[n,1])<r_level;
        cx = random('unif',0,1,[p,1])<c_level;
        x_miss(rx,cx) = NaN;
        if sum(sum(isnan(x_miss)))/(n*p) - level > 0 %too much missing data
            r_level = r_level/2;
            c_level = c_level/2;
        else
            r_level = r_level*1.1;
            c_level = c_level*1.1;
        end
        iter = iter + 1;
    end
    miss_per = sum(sum(isnan(x_miss)))/(n*p);
end



end

