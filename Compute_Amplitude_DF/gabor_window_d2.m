function w=gabor_window_d2( N, L, n )
% w=gabor_window( N, L )
%
% Compute second derivative of Gabor function
%
% 
% INPUT:
% N      : signal length
% L      : window time spread parameter
%
% OUTPUT:
% w    : Gabor function
%
% Author: Q.Legros 
% Date: 

rt=reduced_time_axis(N);
w=  (((rt) ./ (L^2)).^2) .* exp(-((rt).^2 ./ (2*L^2)));
end