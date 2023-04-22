function w=gabor_window_d1( N, L, n )
% w=gabor_window( N, L )
%
% Compute first derivative of Gabor function
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
w=  -((rt) ./ (L^2)) .* exp(-((rt).^2 ./ (2*L^2)));
end