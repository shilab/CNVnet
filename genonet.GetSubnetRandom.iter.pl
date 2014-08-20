#!/usr/bin/perl -w
#This program samples out subnetwork of 5 edges.
die "\nUsage: GetSubnetRandom.pl [Nodes.list] [Edge list with weight] [maps from seqs to genes(maybe more than nodes.list)] [Genome Seqs list] [length of path] [outputfile prefix]\n" if (@ARGV<5);

use List::Util qw(first max maxstr min minstr reduce shuffle sum );
use sort 'stable'; 

# for i in `seq 3 2 15` ; do for j in 0 1 2 ; do ./genonet.GetSubnetRandom.iter.pl Nodes.txt Edges.txt Map Seqs $i group$i.$j ; done; done
$DEBUG=0;
$VVV=0;
$pathlen=$ARGV[4];
$prefix="";
if(defined $ARGV[5])
{
    $prefix=$ARGV[5];
}
my $max_occu=$ARGV[6]; # the max number of occurrence of each node
#if(!defined $max_occu){
    $max_occu=60;
#}
#read nodes
open FN,"<$ARGV[0]";
my @allnodes=<FN>;
close FN;
my %nodeIdx;
for(my $i=0;$i<scalar(@allnodes);$i++)
{
    chomp($allnodes[$i]);
    $nodeIdx{$allnodes[$i]}=$i;
}

#Reado edges

open FE,"<$ARGV[1]";
my @linklist;
while(<FE>)
{
    chomp;
    my @p=split/\s+/;
    my $n1=$nodeIdx{$p[0]};
    my $n2=$nodeIdx{$p[1]};
    my $w=$p[2];
    #print "$n1,$n2,$w\n";next;
    my @list;
    if(defined($linklist[$n1])){
	@list=@{$linklist[$n1]};
    }
    $list[$n2]=$w;
    $linklist[$n1]=\@list;

    my @list1;
    if(defined($linklist[$n2])){
	@list1=@{$linklist[$n2]};
    }
    $list1[$n1]=$w;
    $linklist[$n2]=\@list1;
    
}
close FE;

#normalize linklist
for(my $i=0;$i<scalar(@allnodes);$i++)
{
    my  $sum=0;
    for(my $j=0;$j<scalar(@allnodes);$j++)
    {
	if(defined $linklist[$i][$j]){
	    $sum=$sum+$linklist[$i][$j];
	}
	else
	{
	    $linklist[$i][$j]=0;
	}
    }
}

my %SeqIdx;
open Fseq,"<$ARGV[3]";
my @allseqs=<Fseq>;
for(my $i=0;$i<scalar(@allseqs);$i++){
    chomp($allseqs[$i]);
    $SeqIdx{$allseqs[$i]}=$i;
}
close Fseq;

my @GeneMapSeq;
open FM,"<$ARGV[2]";
while(<FM>)
{
    chomp;
    my @p=split/\s+/;
    my @seqlist;
    if(defined $nodeIdx{$p[1]} && defined $GeneMapSeq[$nodeIdx{$p[1]}])
    {#the $p[1] gene is in the nodes list AND we already have infor about the genome sequence of this gene.
	@seqlist=@{$GeneMapSeq[$nodeIdx{$p[1]}]};
    }
    push @seqlist,$SeqIdx{$p[0]};#add the sequence idx to the genes
    if(defined $nodeIdx{$p[1]}){
	$GeneMapSeq[$nodeIdx{$p[1]}]=\@seqlist;#store the seq ids into list.
    };
}
close FM;



#begin iteration

my $klchange=10000;
my @allgroups;
@oldAllGroups=qw();
@allGroups=qw();
$maxrep=30;
$nrep=0;
my $eps=0.01;
our @oldDistCoverage=qw();
my @KLchange=qw(1 1 1);
#while($klchange>$eps  && $nrep<$maxrep){
while(max(@KLchange) > $eps  && $nrep<$maxrep){
    $nrep++;
    print  "Iteration $nrep : KL divergence changed $KLchange[2] from last iteration\n";
#keep an copy of current idx map
    push @oldAllGroups,@allGroups;
    print  "Number of groups in last iteration: ",scalar(@oldAllGroups),"\n";
    PrintSetState(@oldAllGroups);
    @allGroups=qw();
    #filtering out all groups with  max occurance
    my $n1=scalar(@oldAllGroups);
    if($max_occu>=0){
	@oldAllGroups=FilterbyOccu(\@oldAllGroups, $max_occu);
#	my $n2=scalar(@oldAllGroups);
#	print STDERR $n1-$n2," filtered ", $n2," remained\n";
#	$n1=scalar(@oldAllGroups);
#	@oldAllGroups=FilterbyOccu(\@oldAllGroups, $max_occu);
#	my $n2=scalar(@oldAllGroups);
#	print STDERR $n1-$n2," filtered ", $n2," remained\n";
    }
    my $n2=scalar(@oldAllGroups);
    print  $n1-$n2," filtered ", $n2," remained\n" if($VVV);
    #end of filtering
    my $begin=1;
    for(my $st=0;$st<scalar(@allnodes);$st++)
    {
	if($DEBUG==1){print $st,," ;";}
	my $nextnode=$st;
	my @color;#self avoidning, but maybe overlap with other paths.
	my @path;
	$color[$nextnode]=1;
	push @path,$st;
	my $nonext;
	for(my $j=0;$j<$pathlen;$j++)
	{
	    my $sumprob=0;
	    my $nextprob=rand(); #toss a dye;
	    my $totalprob=0;
	    $nonext=1;
	    for(my $nn=0;$nn<scalar(@allnodes);$nn++)
	    {
		if(!defined $color[$nn]){ # not visited
		    $totalprob=$totalprob+$linklist[$nextnode][$nn];
		};
	    }
	    $nextprob=$totalprob*$nextprob;
	    for(my $nn=0;$nn<scalar(@allnodes);$nn++)
	    {
		if(defined $color[$nn]){ # visited
		    next;
		}
		if($sumprob<$nextprob && $nextprob<=$sumprob+$linklist[$nextnode][$nn]){
		    if($linklist[$nextnode][$nn]==0){print STDERR "0 weight edges\n";}
		    push @path,$nn; # add the node to the path
		    $color[$nn]=1;
		    $nextnode=$nn;
		    $nonext=0;
		    last;
		}
		$sumprob=$sumprob+$linklist[$nextnode][$nn];
	    }
	    if($nonext==1)
	    {
		last;
	    }
	}
	if(@path<2){
	    #$st=$st-1;
	    next;} # repeat if find a path <= 2;

	if($nonext==0 || 1 )    { # allow all paths with length < maxlength= nHops
	    my @seqgroup;
	    for(my $t=0;$t<scalar(@path);$t++)	{
		my @seqs=@{$GeneMapSeq[$path[$t]]};
		push @seqgroup, @seqs;
	    }
	    for(my $t=0;$t<scalar(@path);$t++)	{
	    }
	    for(my $si=0;$si<@seqgroup;$si++){
	    }
	    $strseqgroup=join(" ",@seqgroup);
	    if(!defined $checkrepeat{$strseqgroup}){
		$checkrepeat{$strseqgroup}=1;
		push @allGroups,\@seqgroup;
	    }
	}
    }

    #print FOUT $begin;
    print  "Number of candidate groups in this iteration ",scalar(@allGroups),"\n";
    $klchange=CheckConverge(\@oldAllGroups,\@allGroups);
    push @KLchange,$klchange;
    shift @KLchange;
    print  "KL divergence changes $klchange in this iteration\n";
}
SaveSets("$prefix.idx.txt","$prefix.map.txt",\@oldAllGroups);


sub CheckConverge{
    my @oldAllGroups=@{$_[0]};
    my @allgroups=@{$_[1]};
    #@oldDistCoverage;
    return 1 if(@oldAllGroups+0==0); #there is no any groups sampled.
    if(@oldDistCoverage+0==0){
	@oldDistCoverage=GroupCoverage(@_);
	return 1;
    }
    @dist=GroupCoverage(@_);
    print   "length of distribution ",@oldDistCoverage+0," ",@dist+0,"\n" if $VVV;

    #compute distirbution of ExpectedCOverage; 
    print  join(" ",@dist),"\n" if $VVV;
#    print  join(" ",@oldDistCoverage),"\n";
    my     $kl=KLdiv(\@oldDistCoverage,\@dist);
    @oldDistCoverage=@dist;
    return $kl;

}

sub KLdiv{
    my @x=@{$_[0]};
    my @y=@{$_[1]};
    my $res=0;
#    print STDERR join(" ",@x),"\n";
#    print STDERR join(" ",@y),"\n";

    for(my $i=0;$i<@x;$i++){
	if($x[$i]>0){
	    $y[$i]=0.001 if($y[$i]==0);
	    $res=$res+$x[$i]*log($x[$i])-$x[$i]*log($y[$i]);
	    print  $x[$i]," ",$y[$i],"\n" if $VVV;
	}
    }
    return $res;
}

sub GroupSimi{
    #
    my @g1=@{$_[0]};
    my @g2=@{$_[1]};
    my @g12=@g1;
    push @g12,@g2; # union of g1 and g2

    my %mg1=map{$_=>1}@g1;
    my @inter=grep{$mg1{$_}}@g2;
    my %inter=map{$_=>1}@inter;
    @inter=keys %inter;
#    return scalar(@inter)/max(@g1+0,@g2+0);
    
    my $res= scalar(@inter)/max(@g1+0,@g2+0);
    die join(" ",@inter)."\n".join(" ",@g1) if($res>1);
    return $res;
}

sub GroupCoverage{
my @grps=@{$_[0]};
my @testgrps=@{$_[1]};
print scalar(@grps)," ",scalar(@testgrps);

my @ecvg; #all group similarity
#return qw() if(@testgrps==0)
for(my $i=0;$i<@testgrps;$i++){
    my @r;
    $max_sim=0;
    for(my $j=0;$j<@grps;$j++){
	my $gsimi=GroupSimi($testgrps[$i],$grps[$j]);
	if($max_sim<$gsimi){
	 $max_sim=$gsimi;      
	}
	last if($max_sim>0.999);
    }
    push @ecvg, $max_sim;
}
#print STDERR join(" ",@ecvg),"\n";
return Histogram(@ecvg);
}

sub Histogram{
    my @bin;
    for( my $i=0;$i<20;$i++){
	push @bin, $i/20;
    }
    #input a vector
    my @x=@_;
    my $i;
    my @dist;
    for( my $i=0;$i<=20;$i++){
	$dist[$i]=0;
    }
#    die "input size 0" if (@x+0==0);
    my $sum=0;
    for $i(@x){
	$dist[int($i*20)]++;
	$sum=$sum+1;
    }

    for( my $i=0;$i<=20;$i++){
	print  "$i " if(!defined $dist[$i]);
	$dist[$i]=$dist[$i]/$sum;
    }
    #output a distribution on a given binning 
    return @dist;
}

sub SamplingSets
{
}

sub SaveSets
{
    open FOUT,">$_[0]";
    open FMAP,">$_[1]";
    my @allgroups=@{$_[2]};
    my $begin=1;
    for(my $i=0;$i<@allgroups;$i++){
	my @seqgroup=@{$allgroups[$i]};
	print FMAP join(" ",@seqgroup)," ";
	print FOUT $begin," ",$begin+scalar(@seqgroup)-1,"\n";
	$begin=$begin+scalar(@seqgroup);
    }
    close FOUT;
    close FMAP;
    
}

sub LoadSets
{
}

#compute the distribution of coverage.

sub FilterbyOccu
{
   
    my @allgroups=@{$_[0]};
    my $occu_cutoff=$_[1];
    return if(@allgroups<10);
    #1. first sort all groups by increasing representativeness
    #2. remove in this order one by one 

    #1. finding all element with occu > max_occu
    #2. remove the group with largest occur first.
    my @groupOccu;
    my %featOccu;

    #build @featOccu
    for(my $i=0;$i<@allgroups;$i++)    {
	my @onegroup=@{$allgroups[$i]};
	for(my $j=0;$j<@onegroup;$j++)	{
	    $featOccu{$onegroup[$j]}++ if(defined $featOccu{$onegroup[$j]});
	    $featOccu{$onegroup[$j]}=1 if(!defined $featOccu{$onegroup[$j]});
	}
    }
    for(my $i=0;$i<@allgroups;$i++)    {
	$groupOccu[$i]=0;
	my @onegroup=@{$allgroups[$i]};
	for(my $j=0;$j<@onegroup;$j++)	{
	    $groupOccu[$i]=$groupOccu[$i]+$featOccu{$onegroup[$j]};
	}
    }
    #end of build @featOccu, @groupOccu
    @groupOccu=sort{$a <=> $b}(@groupOccu);
    #debug
    @featOccu=sort{$a <=> $b }(values %featOccu);
    print  "Occurrence ",join(" ",@featOccu[@featOccu-10..@featOccu-1]),"\n" if(scalar(@featOccu)>=10);
#remove the group with large groupOccu repeatly , until no featOccu > occu_cutoff
    #sort once, remove many times 
    #remove
    my $i;
    my $max_featoccu=0;
    my @selected=qw();
    for($i=@allgroups-1;$i>=0;$i--)    {
	my @onegroup=@{$allgroups[$i]};
	my $flagone=0; #the group has at least one node with occu = 1 . We should keep this group.
	for(my $j=0;$j<@onegroup;$j++)	{
	    $flagone=1 if ($featOccu{$onegroup[$j]}==1);
	}
	if($flagone==0){
	}
	else{
	    push @selected,$i; # keep group $i
	    next;
	}

	for(my $j=0;$j<@onegroup;$j++)	{
	    $featOccu{$onegroup[$j]}-- if($featOccu{$onegroup[$j]}>0);
	    #rebuild the heap of feat
	    $max_featoccu=max(values %featOccu);
	}
	if($max_featoccu < $occu_cutoff){
	    last;
	}
    }

    my @selgroups=@allgroups[0..$i];
    for(my $j=0;$j<@selected;$j++){
	push @selgroups,$allgroups[$selected[$j]];
    }
	
    @allgroups=@selgroups;
    return @allgroups;
}

sub PrintSetState{
    #representative
    #redundance
    #occurrance
    
}
