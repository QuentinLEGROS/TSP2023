function [F_mat]=comp_Fc(M,L)
%
% Compute convolutional Gabor kernel. The Gabor
% kernel is approximated using a Gaussian distribution.
%
% INPUT:
% M             : Number of frequential bin
% L             : analysis window size (in bin)
% 
% 
% OUTPUT:
% F_mat         : Postulated obervation model

% Compute the sliding tfr window
val = transpose(Fh(-(M/4)+1:((3*M)/4), M, L ));

% Generate a 2D array for later use. Accelerate further computations
F_mat = zeros(M,M/2);
F_mat(:,1) = val;
for i = 1:(M/2)-1
    F_mat(:,i+1) = [F_mat(end,i);F_mat(1:end-1,i)];
end
% F_mat = F_mat((M/4)+1:((3*M)/4),:); % Truncation to the same lenght than the data
F_mat = F_mat((M/4):((3*M)/4)-1,:); % Truncation to the same lenght than the data
F_mat = F_mat./sum(F_mat);% Normalization


