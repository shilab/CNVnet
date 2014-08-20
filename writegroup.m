function writegroup(groupidx,groupind3,groupmap,fname,inwebIdx);
% writegroup(groupidx,groupind3,groupmap,fname);
fid = fopen(fname, 'w');
%load('/home-nfs/xing/genonetwork/svn/data/156Data_inweb.mat','inwebIdx');

for i=1:length(groupidx)
    pos=groupidx(i);
    %norm(x(groupind(i)+1:groupind(i+1)));
%    if(y>0) 
        s=sprintf('%d',pos); % group index of this group,st0
        for(j=groupind3(pos,1):groupind3(pos,2));
            %            s=sprintf('%s,%d',s,groupmap(j) );
            s=sprintf('%s,%d',s,inwebIdx(groupmap(j)) );
        end;
        fprintf(fid, '%s\n', s);
%    end
end
fclose(fid);
