function [pred_w c1 err logprob]=run_glLogisticR(xx,ll,groupind,group0,rho)
%help
opts=[];
% Starting point

opts.init=2;        % starting from a zero point
opts.tFlag=5;
opts.maxIter=1000;   % maximum number of iterations
opts.nFlag=0;       % without normalization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
opts.G=group0;

%opts.gWeight=(groupind(2,:)-groupind(1,:)+1).^2;

%w=[ind(1:end-1)+1; ind(2:end) ; ones(1,length(ind)-1) ]; %??
opts.ind=groupind;       % set the group indice
opts.q=2;           % set the value for q
opts.sWeight=[1,1]; % set the weight for positive and negative samples

tic;
[pred_w, c1, funVal1, ValueL1]= overlapping_LogisticR(xx, ll, rho, opts);
%[x, c, funVal, ValueL]= overlapping_LogisticR(A, y, z, opts);
etime=toc;
fprintf(1,'%.3f ',etime);
logprob=1;
if (nargout>=3)
    prob=exp(xx*pred_w+c1)./(exp(xx*pred_w+c1)+exp(-xx*pred_w-c1));
    llpred=(prob>0.5)*2-1;
    err=sum(llpred~=ll)/length(ll);
    
    logprob=sum(log(abs((ll+1)/2-prob)));
end;

