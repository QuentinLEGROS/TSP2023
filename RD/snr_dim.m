function Y = snr_dim(s,M,dim)
%SNR_DIM compute the SNR along the dim-th dimention
%   [Y] = SNR_DIM(s, M)
%   [Y] = SNR_DIM(s, M, dim)
%
% INPUTS:
%   s        : signal
%   M        : matrix with variants of s
%   dim      : dimonsion in M containing variants of s
%              default is 1.
%
% OUTPUTS:
%   Y        : matrix of output SNRs between s and (s - M)

if nargin == 2
    dim = 1;
end

if dim == 2
    s = s(:).';
else
    s = s(:);
end

Y_sz = size(M);
Y_sz(dim) = 1;
Y = squeeze(zeros(Y_sz));
for idx = 1:numel(Y)
    I_Y = cell(1, ndims(Y));
    [I_Y{:}] = ind2sub(size(Y),idx);
    I_s = [I_Y(1:dim-1), {':'}, I_Y(dim:end)];

    Y(I_Y{:}) = snr(s, s - squeeze(M(I_s{:})));
end

end

