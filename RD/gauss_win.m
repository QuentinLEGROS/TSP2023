function [g, Lh] = gauss_win(L, sigma_w, varargin)
%GAUSS_WIN Gaussian window
%   [g, Lh] = GAUSS_WIN(L, sigma_w)
%
% INPUTS:
%   L        : samples per second
%   sigma_w  : time spreading parameter
%
%   --  List of possible name-value pair argument
%   'prec'   : precision at which the window shoud be truncated
%              default is 1E-3.
%
% OUTPUTS:
%   g        : the gaussian window
%   norm2h   : the L2 norm of g

defaultPrec = 10^(-3);
p = inputParser;
addRequired(p,'L');
addRequired(p,'sigma_w');
addOptional(p,'prec',defaultPrec);
parse(p,L,sigma_w,varargin{:});
prec = p.Results.prec;

Lh = floor(L*sigma_w*sqrt(-log(prec)/pi))+1;
Lg = 2*Lh + 1;

t=(1:Lg)'/L - (Lh + 1)/L;
g = exp(-(t/sigma_w).^2 * pi);

end