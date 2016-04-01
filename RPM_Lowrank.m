function [Q,S]=RPM_Lowrank(Y, beta, numeig, eigfgt)
[M,D]=size(Y);
hsigma=sqrt(2)*beta;
OPTS.issym=1;
OPTS.isreal=1;
OPTS.disp=0;
if ~eigfgt
    G=con_K(Y,Y,beta);
    [Q,S]=eigs(G,numeig,'lm',OPTS);
    return;
end
e          = 10;
K          = round(min([sqrt(M) numeig]));
p          = 8;
[Q,S]=eigs(@grbf,M,numeig,'lm',OPTS);
    function y=grbf(x,beta)
        [xc , A_k] = fgt_model(Y' , x', hsigma, e,K,p);
        y = fgt_predict(Y' , xc , A_k , hsigma,e);
    end

end