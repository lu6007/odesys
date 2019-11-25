% function error = error_function(time, y, y_hat)
% Calculate the L2 distance between y and y_hat normalized over time
function error = error_function(time, y, y_hat)
dist2 = (y - y_hat).^2; % both are column vectors
delta_t = diff(time);
num_point = length(delta_t);
% The area of the error function in unit nM. 
% normalize to be independent of time span
error2 = sum(0.5 *delta_t.* (dist2(1:num_point) + dist2(2:num_point+1)));
error = sqrt(error2/time(end)); 
end

