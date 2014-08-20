function [EmpiricalDistr zscore pvalue nzgroup weight]=run_gllassopvalue(xx,ll,groupind, ...
						  group0,rho,selectedgroup,estype)
zscore=0;
pvalue=0;
%!!! SET NREPEAT VALUE FOR YOUR APPLICATION. 1000 IS RECOMMENDED. !!!
nrepeat=10;
%nrepeat=1;
%groupind indicates all begining and ending point in extended features
%matrix.
EmpiricalDistr=zeros(length(groupind), nrepeat);
%permutate labels and do the group lasso to compute empirical distribution
for i=1:nrepeat 
  mp_ll_rand=ll(randperm(length(ll)));
  %compute the weights with right labels.
  %compute each group total weight.
  if(i==1)mp_ll_rand=ll;end;
  [x1 c1]=run_glLogisticR(xx,mp_ll_rand,groupind,group0,rho);
  if(i==1)
    [x1 c1 err]=run_glLogisticR(xx,mp_ll_rand,groupind,group0,rho);
    weight=x1;
    nzgroup=findnzgroup(x1,groupind',group0);
    err
  end;
%return ;
  absx1=abs(x1);
  [aa,ii]=sort(absx1,'descend');
  %tie breaking.
  zz_idx=find(aa==0);
  zzperm=ii(zz_idx);
  zzperm=zzperm(randperm(length(zz_idx)));
  ii(zz_idx)=zzperm;
  invii=ii;
  for h=1:length(ii);
    invii(ii(h))=h;
  end
  nn=length(x1);
  %compute every enrichemnt score 
  for j=1:nn;
    x1corr(j)= corr2(mp_ll_rand, xx(:,j));
  end
  for j=selectedgroup
    %ks=EnrichScore(x1corr,aa,invii,groupind,group0,j,0);
    %ks=groupnorm(x1,aa,invii,groupind',group0,j); %nzgroup have small
                                                  %p-value, so are many
                                                  %other groups
 %nzgroup have small
    switch estype
      case{''}
      ks=groupnorm(x1,aa,invii,groupind',group0,j);%by default.
     case {'groupnorm'}
      ks=groupnorm(x1,aa,invii,groupind',group0,j);
     case{'groupnormRank'}
      ks=groupnormRank(x1,aa,invii,groupind',group0,j);
     case{'esp0'}
      ks=EnrichScore(x1,aa,invii,groupind,group0,j,0);
     case{'corresp0'}
      ks=EnrichScoreCorr(x1,x1corr,aa,invii,groupind,group0,j,0);
     case{'esp1'}
      ks=EnrichScore(x1,aa,invii,groupind,group0,j,1);
     case{'corresp1'}
      ks=EnrichScoreCorr(x1,x1corr,aa,invii,groupind,group0,j,1);
     case{'esp2'}
      ks=EnrichScore(x1,aa,invii,groupind,group0,j,2);
    end
    EmpiricalDistr(j,i)=ks;
  end    
end
[null ksrank]=sort(EmpiricalDistr(:,i),'descend');
invKSrank= zeros(length(ksrank), 1);
for j = 1:length(ksrank)
  invKSrank(ksrank(j)) = j;
end
%EmpiricalDistr(:,i)=invKSrank;
%showx=[selectedgroup' EmpiricalDistr(selectedgroup,1:i)]

zscore=zeros(length(groupind), 1);
pvalue=-ones(length(groupind), 1);
for j=1:length(groupind)
  if( ~isnan(EmpiricalDistr(j,1)) )
    %zscore(j)=EmpiricalDistr(j,1);
    zindex=find(~isnan(EmpiricalDistr(j,1:nrepeat)));
    if(length(zindex)<=1) 
      continue; 
    end;       
    %EmpiricalDistr(j,zindex)=0;
    %zindex=2:nrepeat;
    %[hh hhcc]=hist(EmpiricalDistr(j,2:end),100);
    pvalue(j)=length(find(EmpiricalDistr(j,zindex)<= EmpiricalDistr(j,1)))/length(zindex);
    %zmean=mean(EmpiricalDistr(j,zindex));
    %zstd=std(EmpiricalDistr(j,zindex));
    %zsore(j)=(zscore(j)-zmean)/zstd;
    %pvalue(j)=cdf('norm',zscore(j),0,1);
  end;
end;

%----------------------function EnrichScore--------
function [ks]=EnrichScore(x1,aa,invii,groupind,group0,j,p)
%input the weight of each feature, aa:sorted weight, invii:rank
%information, j which f:ature to be computed.

  %sum=norm(x1(groupind(j):groupind(j+1)-1))^2;
  %compute K-S sumedit
  absx1=abs(x1);
  nn=length(absx1);
  ks=-1e20;
  kssum=0;
  nh=groupind(2,j)-groupind(1,j)+1;
  if(p~=0)
    nr=norm(x1( group0(groupind(1,j):groupind(2,j)) ) ,p) ; 
  else
   nr=groupind(2,j)-groupind(1,j)+1; 
  end
  %p=1, as in paper of kai wang, mingyao li, maja bucan
  geneset=group0(groupind(1,j):groupind(2,j));
  
  %	[ggss ggii]=sort(absx1(geneset),'descend');
  %  %find all the rank of geneset members in all x1
  gsrank=invii(geneset);
  %rank of each element in geneset
  %the global rank of the ith gene in geneset
  [tmp ggii]=sort(gsrank);
  
  for(iks=1:length(geneset))
    if(p~=0)
      partnr=norm(x1( geneset(ggii([1:iks])) ),p )  ;
    else
      partnr=iks;
    end
%    kssum=iks/nr-(gsrank(ggii(iks))-iks)/(nn-nh);
    kssum=partnr/nr-((ggii(iks))-iks)/(nn-nh);
    %	  if(ismember(ii(iks),geneset))
    %	    kssum=kssum+1/nr;
    %	  else
    %	    kssum=kssum-1/(nn-nh);
    %	  end;
    if(ks ==-1e-20)
      ks=ksum;
    end;
    if(ks<kssum)
      ks=kssum;
    end;
  end
  %       EmpiricalDistr(j,i)=ks;% use the max ksum
  %compute the rank of ES score.

%----------------
function [ks]=EnrichScoreCorr(x1,x1corr,aa,invii,groupind,group0,j,p)
%input the weight of each feature, aa:sorted weight, invii:rank
%information, j which f:ature to beecomputed.

  %sum=norm(x1(groupind(j):groupind(j+1)-1))^2;
  %compute K-S sumedit
  absx1=abs(x1);
  nn=length(absx1);
  ks=-1e20;
  kssum=0;
  nh=groupind(2,j)-groupind(1,j)+1;
  %        nr=sum( absx1( group0(groupind(1,j):groupind(2,j)) ) ); 
%  nr=groupind(2,j)-groupind(1,j)+1; 
  if(p~=0)
    nr=norm(x1( group0(groupind(1,j):groupind(2,j)) ) ,p) ; 
  else
   nr=groupind(2,j)-groupind(1,j)+1; 
  end
%  nr=norm(x1( group0(groupind(1,j):groupind(2,j)) ) ,p) ; 
  %p=1, as in paper of kai wang, mingyao li, maja bucan
  geneset=group0(groupind(1,j):groupind(2,j));
  
  %	[ggss ggii]=sort(absx1(geneset),'descend');
  ggss=aa(invii(geneset));
  ggii=invii(geneset);
  %  %find all the rank of geneset members in all x1
  gsrank=invii(geneset);
  %rank of each element in geneset
  %the global rank of the ith gene in geneset
  [tmp ggii]=sort(gsrank);
  
  for(iks=1:length(geneset))
    if(p~=0)
      partnr=norm(x1( geneset(ggii([1:iks])) ),p )  ;
    else
      partnr=iks;
    end
%    partnr=norm(x1( geneset(ggii([1:iks])) ),p )  ;
    kssum=partnr/nr-((ggii(iks))-iks)/(nn-nh);
%    partnr=sum(absx1( gsrank(ggii([1:iks])) ));
%    kssum=iks/nr-(gsrank(ggii(iks))-iks)/(nn-nh);
    %	  if(ismember(ii(iks),geneset))
    %	    kssum=kssum+1/nr;
    %	  else
    %	    kssum=kssum-1/(nn-nh);
    %	  end;
    if(ks ==-1e-20)
      ks=ksum;
    end;
    if(ks<kssum)
      ks=kssum;
    end;
  end
  %       EmpiricalDistr(j,i)=ks;% use the max ksum
  %compute the rank of ES score.


%---------only norm---------------
function [ks]=groupnorm(x1,aa,invii,groupind,group0,j)
ks=norm(x1(group0(groupind(j,1):groupind(j,2))));

function [ks]=groupnormRank(x1,aa,invii,groupind,group0,j)
[nzgroupid nzgroupnorm allgnorm]=findnzgroup(x1,groupind,group0);
ks=length(find(allgnorm>allgnorm(j)));

function [ks]=groupcorr(x1,aa,invii,groupind,group0,j)
ks=norm(x1(group0(groupind(j,1):groupind(j,2))));
