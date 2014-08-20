function run_cv_pvalue_group_iter_maxoccu(datafile,groupIdxFile,groupMapFile)
%[mus1 errors1 logprobs1]=cross_validation(matrix_pilot_1_feature, ...
                                           %				  matrix_pilot_1_label, groupIndFull,groupMapFull,[1e-3:1e-3:1e-2]);%get a group partition
                                           %for gsize=3:2:15		
                                           %    for g_rep=0:2	
                                           %gsize=3;%getenv('GSIZE');
                                           %g_rep=2;%getenv('G_REP')%;
addpath(genpath('../bin/'));
addpath(genpath('./'));
outputPrefix='groupTest';
%load ../data/156DataNew;
%use feature_normalized, label, maxoccu, 
%load('/home-nfs/xing/genonetwork/svn/data/156Data_inweb.mat','inwebIdx');

%keyboard;
%read feature_normalized
%read label
allsample=dlmread(datafile);
feature_normalized=allsample(:,2:end);
label=allsample(:,1);
matrix_pilot_1_feature=feature_normalized;
matrix_pilot_1_label=label;
%keyboard

%!!! SET THE MU VALUES FOR YOUR APPLICATION. !!!
muRange=[1e-3 1e-2 0.1];

if(~exist('maxoccu'))
maxoccu='-1';
end

if(exist('inwebIdx.txt','file'))
    inwebIdx=dlmread('inwebIdx.txt');
else
    inwebIdx=1:size(feature_normalized,2);
end


ReInwebIdx=zeros(1,length(matrix_pilot_1_feature));
for i=1:length(inwebIdx)
    ReInwebIdx(inwebIdx(i))=i;
end;
%keyboard;
matrix_pilot_1_feature=matrix_pilot_1_feature(:,inwebIdx);

%        gidx=dlmread(sprintf('../NodeEdgeMapPilot1-101510/group_iter_lowredun.%s.idx.txt',maxoccu));
        gidx=dlmread(groupIdxFile);
%        gmap=dlmread(sprintf('../NodeEdgeMapPilot1-101510/group_iter_lowredun.%s.map.txt',maxoccu));
        gmap=dlmread(groupMapFile);
        gmap=gmap+1; % !!! caution about extract all group CNV infor....
        gmap=ReInwebIdx(gmap);
        gmap=reshape(gmap,1,length(gmap));
        gidx(:,3)=gidx(:,2)-gidx(:,1)+1;
        if(gidx(end,3)<=0)
            gidx=gidx(1:end-1,:);
        end
%        GeneCNVidx=unique(sort(gmap));
        featlen=length(matrix_pilot_1_feature);
%        NotInwebCNVidx=setdiff([1:featlen],GeneCNVidx);
%        PsdoGroupMap=NotInwebCNVidx;
 %        PsdoGroupMap=PsdoGroupMap(randperm(length(PsdoGroupMap)));
        
        lenRG=ceil(mean(gidx(:,3))*3); %size of random groups
%        PsdoGroupInd1=[1:lenRG:(length(NotInwebCNVidx)-lenRG)];
%        PsdoGroupInd2=PsdoGroupInd1+lenRG-1;
%        PsdoGroupInd2(end)=length(NotInwebCNVidx);
%        PsdoGroupInd3=ones(1,length(PsdoGroupInd1))*lenRG;
%        PsdoGroupInd=[PsdoGroupInd1' PsdoGroupInd2' PsdoGroupInd3'];
%        groupMapFull=[gmap PsdoGroupMap];
        groupMapFull=gmap ;


%        PsdoGroupInd1=PsdoGroupInd1+length(gmap);
%        PsdoGroupInd2=PsdoGroupInd2+length(gmap);
%        PsdoGroupInd=[PsdoGroupInd1' PsdoGroupInd2' PsdoGroupInd3'];
%        groupIndFull=[gidx; PsdoGroupInd];
        groupIndFull=gidx;

        %[mus errors logprobs]=cross_validation(matrix_pilot_1_feature,matrix_pilot_1_label,gidx,gmap,[1e-4 5e-4 1e-3:1e-3:1e-2]);
        %add random groups
        %                return;
        %        keyboard;
        %        keyboard;
        %        [mus errors
        %        logprobs]=cross_validation(matrix_pilot_1_feature,matrix_pilot_1_label,groupIndFull,groupMapFull,[1e-4
        %        5e-4 1e-3:1e-3:1e-2]);
        %if(strcmp(maxoccu,'nomax'))
        %[mus errors logprobs]=cross_validation(matrix_pilot_1_feature, ...
        %                                       matrix_pilot_1_label, ...
        %                                       groupIndFull,groupMapFull,[1e-5 ...
        %                    :1e-5:1e-4]);
        %else
        [mus errors logprobs]=cross_validation(matrix_pilot_1_feature, ...
                                               matrix_pilot_1_label, ...
                                               groupIndFull,groupMapFull,muRange);
        %end
        % ...
        %                    5e-4 1e-3:3e-3:1e-2 0.03:0.03:0.1]);
        
        %        errors
        %mus
        %save(sprintf('cv_group_iter_%s_rep',maxoccu),'mus','errors','logprobs','groupIndFull','groupMapFull');
        %        return;
        for  tlabel=0:2
            errorsindx=find(errors<=0.02); %error cut
            e=errorsindx(end);
            mu=mus(e);
            %fprintf(1,'Trgt lble %d mu %f', tlabel,mu);
            [em zs pv nzgroup x1]=run_gllassopvalue(matrix_pilot_1_feature, ...
                                                    (matrix_pilot_1_label==tlabel)*2-1, ...
                                                    groupIndFull',groupMapFull,[mu mu],[1:length(groupIndFull)],'groupnorm');
            %   save(sprintf('pvalue_group_iter_tlabel_%d_%s',tlabel,maxoccu),'em','zs','pv','nzgroup','x1');
            writegroup(nzgroup,groupIndFull,groupMapFull,sprintf('%s.%d.nzgroup.csv',outputPrefix,tlabel),inwebIdx);
            %            dlmwrite(sprintf('%s.%d.nzgroup.csv',outputPrefix,tlabel), ...
            %                     nzgroup,' ');
            dlmwrite(sprintf('%s.%d.pvalue.csv',outputPrefix,tlabel), ...
                     pv,' ');
            dlmwrite(sprintf('%s.%d.weight.csv',outputPrefix,tlabel), ...
                     x1,' ');
        end
