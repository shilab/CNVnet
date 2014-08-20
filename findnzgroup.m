function [res gnorm allgnorm]=findnzgroup(xx,groupind,groupmap)
res=[];
gnorm=[];
allgnorm=[];
for ii=1:length(groupind);
%    ss=length( find((xx(groupmap([groupind(ii,1):groupind(ii,2)])))~=0) );
    ss=norm(xx(groupmap([groupind(ii,1):groupind(ii,2)])));
    allgnorm=[allgnorm ss];
    if(ss~=0)
        res=[res ii];
        gnorm=[gnorm ss];
    end
end

