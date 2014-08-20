function [err logprob]=run_glLogisticR_pred(xx,ll,pred_w,c1)
    prob=exp(xx*pred_w+c1)./(exp(xx*pred_w+c1)+exp(-xx*pred_w-c1));
    llpred=(prob>0.5)*2-1;
    err=sum(llpred~=ll)/length(ll);
    logprob=sum(log(abs((1-ll)/2-prob)));
    

