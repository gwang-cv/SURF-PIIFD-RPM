function [E, G] = RPM_cost(param, X, Y, K, U, lambda, sigma2, is_grad,r)
% ========================================================================
% Robust Point Matching For Multimodal Retinal Image Registration
% Copyright(c) 2015 Gang Wang
% All Rights Reserved.
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereQ
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
% For more information, please refer to
% Gang Wang et al., "Robust Point Matching Method for Multimodal Retinal
% Image Registration", Biomedical Signal Processing and Control, 2015, 
% Vol. 19, pp. 68-76.
%----------------------------------------------------------------------
%Author: Gang Wang
%Contact: gwang.cv@gmail.com
%Version: 2014-11-30 
% ========================================================================
[N, D] = size(X);
M = size(U, 2);
C = reshape(param, [M D]);
E0 = 1/( 2^D * pi^(D/2) * (sigma2*(((r+1)/2)^2))^(D/4) );
E1 = E0 + lambda/2 * trace(C'*K*C);
a = -2 / N / (2*pi*sigma2)^(D/2) / ((r+1)/2)^D;
V = U*C - Y;
Va=sum(V,2);
F0 = exp(-sum(V.^2, 2) / (2*sigma2));
F1 = exp(-sum(V.^2, 2) / (2*sigma2*r*r));
F=((Va<=0).*F0+(Va>0).*F1);
E = E1 + a * sum(F);
G = [];
if is_grad
    Vs=repmat(sum(V,2),[1 D]);
    FF=repmat(F, [1 D]);
    Vb= V .* FF;
    Vc=(Vs<=0).*Vb/sigma2+(Vs>0).*Vb/(sigma2*r*r);
    G = - a * U' * ( Vc ) + lambda * K * C;
    G = G(:);
end