function  [Y_reg,idt]=RPM(X,Y)
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
%
%Input : (1) X: featuer points extracted from the reference image
%        (2) Y: featuer points extracted from the second image
%
%Output: the index of the correct matches
%-----------------------------------------------------------------------
%Author: Gang Wang 
%Contact: gwang.cv@gmail.com
%Version: 2014-11-30
% ========================================================================

% initialize
normalize = 1;
n_iter=55;
k=1;
s=1;
Xk=X;
Yk=Y;
[N1, D] = size(X);
[N2, D] = size(Y);
idt_mat={};
v_mat=[];
min_value_mat=[];
sigma2=(N1*trace(X'*X)+N2*trace(Y'*Y)-2*sum(X)*sum(Y)')/(N1*N2*D);
r=1;

while s
    X2 = Xk;
    Y2=Yk;
    normal.xm=0; normal.ym=0;
    normal.xscale=1; normal.yscale=1;
    if normalize, [nX, nY, normal]=norm2s(X2,Y2);
        if k<15
            [idt, V, param,min_value] = L2E_AGM_Aff(nX, nY-nX, 0.93, 1.0,sigma2);
        else
            [idt, V, param,min_value] = L2E_AGM_Nonrigid(nX, nY-nX, 0.93, r,sigma2);
        end
        min_value_mat=[min_value_mat,min_value];
        idt_mat{k}=idt;
        sigma2= sigma2*0.75;
        V=(V+nX)*normal.yscale+repmat(normal.ym,size(Y2,1),1);
        v_mat=[v_mat,V];
    end
    Xk = V;
    if k==n_iter
        s=0;
    else
        k=k+1;
    end
end


[cc,inx]=min(min_value_mat);
V=v_mat(:,2*inx-1:2*inx);
idt=idt_mat{inx};
Y_reg=V;

%% functions
    function [idt,Vs,param,min_value]= L2E_AGM_Aff(X,Y,thresh,r,sigma2)
        
        [N,D]=size(X);
        lambda = 0.0005;
        % Optimal
        param = [1 0 0; 1 0 0; 0 0 1] ;
        param(9) =1;
        options = optimset(  'display','off', 'MaxIter', 1000,'Algorithm','interior-point');
        options = optimset(options, 'GradObj', 'off');
        X_=X;
        Y_=Y;
        X_(:,3)=ones(N,1);
        Y_(:,3)=ones(N,1);
        % constraint minimizing
        Lb = [-Inf; -Inf; -Inf; -Inf; -Inf; -Inf];
        Ub = [Inf; Inf; Inf; Inf; Inf; Inf];
        Lb(7:9,:) = [0;0;1];
        Ub(7:9,:) = [0;0;1];
        [param, min_value] = fmincon(@RPM_cost_Affine, param, [],[], [],[], Lb, Ub, [], options, X_, Y_, sigma2, r,lambda);
        theta = reshape(param, [3 3]);
        V0=X_*theta;
        Vs=V0(:,1:2);
        Pb =  exp(-sum((Y-Vs).^2, 2) / (2*sigma2)) ;
        idt = find(Pb > thresh);
    end

    function [E, G] = RPM_cost_Affine(param, X, Y, sigma2,r,lambda)
        
        [M, dd] = size(X);
        D=dd-1;
        theta = reshape(param, [dd dd]);
        V0=X*theta;
        V=Y-V0;
        E0 = 1/( 2^D * pi^(D/2) * (sigma2*(((r+1)/2)^2))^(D/4) );
        E1 = E0 + 0.5*lambda*trace(V0'*V0);
        a = -2 / M / (2*pi*sigma2)^(D/2) / ((r+1)/2)^D;
        Va=sum(V(:,1:2),2); %n-by-1
        F0 = exp(-diag(V*V') / (2*sigma2));
        F1 = exp(-diag(V*V') / (2*sigma2*r*r));
        F=((Va>=0).*F0+(Va<0).*F1);
        E = E1 + a * sum(F);
    end

    function [idt,V,param,min_value]= L2E_AGM_Nonrigid(X,Y,thresh,r,sigma2)
        
        [N1, D] = size(X);
        [N2, D] = size(Y);
        beta =10;
        lambda = 0.1;
        is_grad = 1;
        n_ker = 15;
        [Q,K]=RPM_Lowrank(X, beta, n_ker, 1);
        U=Q*K;
        x0 = zeros(n_ker*D, 1);
        options = optimset( 'display','off', 'MaxIter', 1000,'Algorithm','interior-point');%,
        if is_grad
            options = optimset(options, 'GradObj', 'on');
        end
        [param,min_value] = fminunc(@(x)RPM_cost(x, X, Y, K, U, lambda, sigma2, is_grad,r), x0, options);
        C = reshape(param, [n_ker D]);
        V=U*C;
        
        if N1~=N2
            V0 = Y- [V;zeros(abs(N2-N1),2)] ;
        else
            V0 = Y- V;
        end
        Pb =  exp(-sum((Y-V).^2, 2) / (2*sigma2)) ;
        idt = find(Pb > thresh);
    end
end