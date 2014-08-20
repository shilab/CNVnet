#!/usr/bin/perl -w

# checkgroup.pl 


#read nodes


open FN,"<Nodes.txt";
my @allnodes=<FN>;
close FN;
my %nodeIdx;
for(my $i=0;$i<scalar(@allnodes);$i++)
{
    chomp($allnodes[$i]);
    $nodeIdx{$allnodes[$i]}=$i;
}


#read edges
open FE,"<Edges.txt";
my %edgeWeight=qw();
my @linklist;
while(<FE>)
{
    chomp;
    my @p=split/\s+/;
    my $n1=$nodeIdx{$p[0]};
    my $n2=$nodeIdx{$p[1]};
    my $w=$p[2];
    #print "$n1,$n2,$w\n";next;
    $edgeWeight{$p[0]." ".$p[1]}=$w;
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
#print join(" ",keys %edgeWeight),"\n";die;
#normalize linklist
for(my $i=0;$i<scalar(@allnodes);$i++)
{
    my    $sum=0;
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
open Fseq,"<Seqs";
my @allseqs=<Fseq>;
for(my $i=0;$i<scalar(@allseqs);$i++){
    chomp($allseqs[$i]);
    $SeqIdx{$allseqs[$i]}=$i;
}
close Fseq;

my %CNVGeneId;
if(0 && -e "cnv2geneid.txt"){
open FHgene,"<cnv2geneid.txt" or die $!;
my @CNVid;
my @allgeneid=<FHgene>;
for(my $i=0;$i<scalar(@allgeneid);$i++)
{
    chomp($allgeneid[$i]);
    my @p=split /\s+/, $allgeneid[$i];
    $CNVGeneId{$p[0]}=$CNVGeneId{$p[0]}." ".$p[1] if(defined $CNVGeneId{$p[0]});
    $CNVGeneId{$p[0]}=$p[1] if(!defined $CNVGeneId{$p[0]});
    
    push @CNVid, $p[0];
}
close FHgene;
}
my %CNVGeneName;
if(0 && -e "cnv2genename.txt"){
open FHgenename,"<cnv2genename.txt" or die $!; 
# from ens gene name to ref gene name
my @allgenename=<FHgenename>;
for(my $i=0;$i<scalar(@allgenename);$i++)
{
    chomp($allgenename[$i]);
    my @p=split /\s+/, $allgenename[$i];
    $CNVGeneName{$p[0]}=$p[1];
}
close FHgenename;
}


my @GeneMapSeq;
my @SeqMapGene;
open FM,"<Map";
while(<FM>)
{
    chomp;
    my @p=split/\s+/;
    my @seqlist;
    if(defined $nodeIdx{$p[1]} && defined $GeneMapSeq[$nodeIdx{$p[1]}])
    {
	@seqlist=@{$GeneMapSeq[$nodeIdx{$p[1]}]};
    }
    push @seqlist,$SeqIdx{$p[0]};
    if(defined $nodeIdx{$p[1]}){
	$GeneMapSeq[$nodeIdx{$p[1]}]=\@seqlist;
    };
    my @genelist;
    if(defined $SeqIdx{$p[0]} && defined $SeqMapGene[$SeqIdx{$p[0]}])
    {
	@genelist=@{$SeqMapGene[$SeqIdx{$p[0]}]};
    }
#    push @genelist,$CNVGeneName{$p[1]};
    push @genelist,$p[1];
    if(defined $SeqIdx{$p[0]}){
	$SeqMapGene[$SeqIdx{$p[0]}]=\@genelist;
    };
}
close FM;

if(0){
open FFST,"<../data/GenotypeGenicCNVE.Pilot1.Fst.onePopVsOthers.dat";
@allfst=qw();
while(<FFST>)
{
    chomp;
    my @p=split/\s+/;
    push @allfst,\@p;
}
#print STDOUT "",scalar(@allfst);
close FFST;

open FHPHA,"<../data/all.3260.cnvs.phastcons44way.degrees.txt";
@allphast=qw();
#$null=<FHPHA>;
while(<FHPHA>)
{
    chomp;
    my @p=split/\s+/;
    @p=@p[@p-2..@p-1];
    push @allphast,\@p;
}
print STDOUT "phastcon ",scalar(@allfst)," read\n";
close FHPHA;
}

open FR,"<$ARGV[0]" or die;
my @allgroups=<FR>;
for(my $k=0;$k<scalar(@allgroups);$k++){
    chomp($allgroups[$k]);
    my @resultseqs=split/,/,$allgroups[$k];
    @resultseqs=@resultseqs[1..@resultseqs-1];#first is the group number
    print "Group $k: \n";
    my $groupfst=0;
    my $grouppha=0;
    my $groupdeg=0;
    my @genes_g=qw();
    for $ii(@resultseqs)
    {
	chomp($ii);
	$ii=$ii-2;
	my @t=qw(); # all genes 
	@t=@{$SeqMapGene[$ii+1]} if(defined $SeqMapGene[$ii+1]);
	my $flagIsNode=0;
	for($i=0;$i<scalar(@t);$i++)
	{
	    if(defined $nodeIdx{$t[$i]}){$flagIsNode=1;}
#	    $t[$i]=$CNVGeeName{$t[$i]} if (defined $CNVGeneName{$t[$i]});
	}
	if($flagIsNode==0)
	{
	    die "$ii+1,$allseqs[$ii+1], Get some gene not in graph!";
	}
	if(scalar(@t)==0)
	{
	    print "NA\n";
	}
	else{
	    push @genes_g, @t;
	    print $ii+1,",".$allseqs[$ii+1]." ,";
#	    $groupfst=$groupfst+$allfst[$ii+1][0];
#	    $grouppha=$grouppha+$allphast[$ii+1][0];
#	    $groupdeg=$groupdeg+$allphast[$ii+1][1];
	    print join(" ",@t),",";
	    print $CNVGeneId{$allseqs[$ii+1]} if(defined $CNVGeneId{$allseqs[$ii+1]});
	    
	    print "\n";
	}
	#output a separate file for the subgraph
    }
#    print "groupfst ".$groupfst."\n";
#    print "groupPhastCon ".$grouppha."\n";
#    print "groupDegree ".$groupdeg."\n";
    open FSG,">$ARGV[0].group.$k.edges.csv" or die;
    for(my $j=0;$j<@genes_g;$j++)
    {
	for(my $j2=$j+1;$j2<@genes_g;$j2++)
	{
	    my $g_key=$genes_g[$j]." ".$genes_g[$j2];
	    print FSG $g_key." ".$edgeWeight{$g_key}."\n" if(defined $edgeWeight{$g_key});
	    $g_key=$genes_g[$j2]." ".$genes_g[$j];
	    print FSG $g_key." ".$edgeWeight{$g_key}."\n" if(defined $edgeWeight{$g_key});
	}	
    }
    close FSG;
#    @res=ExtractSubnet(\@genes_g,\@linklist);
#    print join("\n",@res),"\n";
#    exit;
#    %hashmap=map {$_ => 1} @genes_g;
#    print join(" ",keys %hashmap );

#    print "\n";
}
close FR;


sub ExtractSubnet
{
    my @res;
    my @subNodes=@{$_[0]};
#    my @Edges=@{$_[1]};
    print(join(" ",@subNodes));
#    print(join(" ",$linklist[0][0][1][1]));
    for(my $i=0; $i<@subNodes;$i++)
    {
	my $aa=$nodeIdx{$subNodes[$i]};
#	next if(!defined $aa);
	for(my $j=$i+1; $j<@subNodes;$j++)
	{
	    my $bb=$nodeIdx{$subNodes[$j]};
	#    next if(!defined $bb);
	    if(defined $linklist[$aa][$bb])
	    {
		push @res,"$aa $bb ";#.$Edges[$aa][$bb];
	    }
	}
    }
    return @res;
}
