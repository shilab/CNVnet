function [mus errors logprobs]=cross_validation(sim_feature,sim_label,groupind,groupmap,murange)
%take input of original sim_feature and sim_labels
%take untouched label={0,1,2}
foldidx=randperm(length(sim_label));
oldmu=0;
nsample=length(sim_label);
%tuning on single population labeling 
for targetlabel=0:2
  mus=[];
  errors=[];
  logprobs=[];
  
  
  for mu=murange
    sim_err=0;
    avelogprob=0;
    for fold=1:3
      trainidx=ceil(fold/3*nsample) : floor( (fold+2)/3*nsample );
      testidx=ceil( (fold+2)/3*nsample) : floor( (fold+3)/3*nsample );
      %      keyboard;
      trainidx=mod(trainidx,length(sim_label))+1;
      testidx=mod(testidx,length(sim_label))+1;
      trainidx=foldidx(trainidx);
      testidx=foldidx(testidx);
      %   tic;
      %      keyboard;
      [sim_w sim_c]=run_glLogisticR(sim_feature(trainidx,:),(sim_label(trainidx)==targetlabel)*2-1, ...
					groupind',groupmap,[mu  mu]);
      %toc
      [err logprob]=run_glLogisticR_pred(sim_feature(testidx,:),(sim_label(testidx)==targetlabel)*2-1, ...
      			    sim_w,sim_c);

      sim_err=sim_err+err;
      avelogprob=avelogprob+logprob;
    end

  errors=[errors sim_err/3];
  mus=[mus mu];
  logprobs=[logprobs avelogprob];
 
  end
  %  save(sprintf('cv_res_%d',targetlabel));
%  if(sim_err>0.05*9)
%    fprintf(1,'mu=%f\n',mu);
%    break;
%  end;
end
mu=oldmu;
