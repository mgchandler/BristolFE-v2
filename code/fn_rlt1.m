function [output]=fn_rlt1(dat_2_filt)
%RLT the signal
dat_2_filt=dat_2_filt./max(dat_2_filt); % normalises data. ensured the range is from 0 to 1.
use_dat=abs(dat_2_filt);
count_dat=[1:length(use_dat)]';
use_dat=use_dat./max(use_dat);
%first_term=sum(log10(use_dat));
first_term=sum(log(use_dat));
sigma_1=sqrt(cumsum(use_dat.^2).*1./(2.*(count_dat))); % standard deviation
down_sum=flipud(cumsum(flipud(use_dat.^2)));
sigma_2=sqrt(down_sum.*1./(2.*(length(use_dat)-count_dat)));
%output=2.*first_term-2.*count_dat.*log10(sigma_1)-2.*(count_dat(end)-count_dat).*log10(sigma_2)-count_dat(end);
output=2.*first_term - 2.*count_dat.*log(sigma_1) - 2.*(count_dat(end)-count_dat).*log(sigma_2) - count_dat(end);  %% eq3.12
output(end)=output(end-1);
output=output-min(output);
end

% [output]=fn_rlt1(dat_2_filt): This line defines a MATLAB function called fn_rlt1 that takes one input argument dat_2_filt and returns one output output.
% 
% dat_2_filt=dat_2_filt./max(dat_2_filt);: This line normalizes the input data dat_2_filt by dividing it by its maximum value. This step ensures that the data ranges from 0 to 1.
% 
% use_dat=abs(dat_2_filt);: This line takes the absolute value of the normalized data.
% 
% count_dat=[1:length(use_dat)]';: This line creates a vector count_dat containing integers from 1 to the length of use_dat.
% 
% use_dat=use_dat./max(use_dat);: This line normalizes use_dat by dividing it by its maximum value again, ensuring it ranges from 0 to 1.
% 
% first_term=sum(log(use_dat));: This line computes the sum of the natural logarithm of use_dat. It calculates the logarithm of each element in use_dat, sums them up, and stores the result in first_term.
% 
% sigma_1=sqrt(cumsum(use_dat.^2).*1./(2.*(count_dat)));: This line calculates the standard deviation sigma_1 using a formula involving cumulative sum and square of use_dat.
% 
% down_sum=flipud(cumsum(flipud(use_dat.^2)));: This line calculates a reversed cumulative sum of the squares of use_dat.
% 
% sigma_2=sqrt(down_sum.*1./(2.*(length(use_dat)-count_dat)));: This line calculates another standard deviation sigma_2 using down_sum and the length of use_dat.
% 
% output=2.*first_term-2.*count_dat.*log(sigma_1)-2.*(count_dat(end)-count_dat).*log(sigma_2)-count_dat(end);: This line computes the final output using the calculated values. It involves a mathematical formula that combines first_term, count_dat, sigma_1, and sigma_2.
% 
% output(end)=output(end-1);: This line sets the last element of output equal to the second-to-last element, possibly for smoothing purposes.
% 
% output=output-min(output);: This line subtracts the minimum value of output from all elements, effectively shifting the output so that the minimum value becomes zero.
% 
% In summary, the fn_rlt1 function takes an input signal, performs some calculations involving logarithms and standard deviations, and returns a processed output signal.