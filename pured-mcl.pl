#!/usr/bin/perl -w

$| = 1;

use Getopt::Long;
use LWP::Simple;
#use XML::Validate;
use XML::LibXML::Reader;

$version = localtime((stat($0))[9]);

$usage = "# PuReD-MCL - clustering of related PubMed documents
# Theodosiou T. and Darzentas N.
# $version
# http://bat.infspire.org | bat\@infspire.org

welcome

before you run PuReD-MCL, be sure to:
- have Stijn van Dongen's MCL installed and in your path
- have installed the following Perl modules:
\tGetopt::Long
\tLWP::Simple
\tXML::LibXML::Reader

pured-mcl.pl

    either
--query     query you want to perform in PubMed, enclose it in double quotes
    or
--xml       results from PubMed in XML format

    other (optional) arguments
--min       define the minimum inflation for MCL (smaller the value, larger/fewer the clusters) [default=2.0]
--max       define the maximum inflation for MCL (larger the value, smaller/more the clusters) [default=2.0]
--step      define the increment step to reach the maximun inflation from the minimum [default=1.0]
--max_doc   define the maximum number of documents you want to handle [default=10,000]

";

my $inflation     = 2.0;
my $max_inflation = 2.0;
my $step          = 1.0;
my $max_documents = 10000;
my $pwd           = '.';
my $list          = 0;

GetOptions(
    "query:s"   =>\$query,
    "xml:s"     =>\$xml_filename,
    "min:f"     =>\$inflation,
    "max:f"     =>\$max_inflation,
    "step:f"    =>\$step,
    "max_doc:f" =>\$max_documents,
    "pwd:s"     =>\$pwd
    );

if (!defined($query) && !defined($xml_filename) || (defined($query) && defined($xml_filename))) {
    die $usage;
}

$randomer = int(rand(1001)) . $$;
$randomer = sprintf "%010s", $randomer;

$step = '1.0' if $step == 0;

open STDERR, '>/dev/null';

if (defined($query) =~ /\w/) {
    # Define library for the 'get' function used in the next section.
    # $utils contains route for the utilities.
    # $db, $query, and $report may be supplied by the user when prompted;
    # if not answered, default values, will be assigned as shown below.

    print "(?) querying PubMed with: $query\n";

    my $utils   = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
    my $db      = "Pubmed";
    my $report  = "xml";
    my $esearch = "$utils/esearch.fcgi?" .
                  "db=$db&usehistory=y&term=";
    my $retstart;
    my $retmax=1000;
    my $results = $pwd."/pured-mcl.query$randomer.xml";
    my $esearch_result = get($esearch . $query);

    $esearch_result =~
      m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

    my $count    = $1;
    my $queryKey = $2;
    my $webEnv   = $3;

    if ($count > $max_documents){ $count = $max_documents; print "(!) we found $count documents but only the first $max_documents will be downloaded\n";
    } elsif ($count > 1) {                                 print "(=) we found $count documents - downloading\n";
    } else {                                               print "(=) we found <2 documents, maybe resubmit in case PubMed did not respond - exiting\n"; exit; }

    # this area defines a loop which will display $retmax citation results from
    # Efetch each time the the Enter Key is pressed, after a prompt.
    for($retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $left = $count - $retstart;
        print " $left\n" unless $retstart == 0;
        my $efetch = "$utils/efetch.fcgi?" .
        "retmode=$report&retstart=$retstart&retmax=$retmax&" .
        "db=$db&query_key=$queryKey&WebEnv=$webEnv";

        my $efetch_result = get($efetch);
        $efetch_results .= $efetch_result;
        sleep(2);
    }
    print " ok\n" if $retstart > $retmax;
    @split_efetch_results = (split /\n/, $efetch_results);
    foreach $line (@split_efetch_results) {
        unless ($line =~ /(^\<\?xml version|^\<\!DOCTYPE|^\<PubmedArticleSet|^\<\/PubmedArticleSet)/) {
            $cleaned_efetch_results .= "$line\n";
        }
    }
    $xml_results = <<XMLCODE;
<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2017//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_170101.dtd">
<PubmedArticleSet>
$cleaned_efetch_results
</PubmedArticleSet>
XMLCODE
    open(RES,">$results") or die "Cannot create file query_res.xml: $!\n";
    binmode(RES, ":utf8");
    print RES $xml_results;
    close RES;
    $xml_filename = $results;
    print "(?) download complete and XML saved\n";
} elsif (defined($xml_filename)) {
    print "(?) loading PubMed XML\n";
}

my $mcli_filename = $xml_filename . ".mcli";

if ($inflation == $max_inflation) { print "(?) will run PuReD-MCL with inflation $inflation\n";
} else {                            print "(?) will run PuReD-MCL with inflation from $inflation to $max_inflation in $step steps\n"; }

$original_inflation = $inflation;
 main_maketab($xml_filename,$list);
 main_mcl_evaluation($mcli_filename,$xml_filename,$inflation,$step,$max_inflation);
$inflation = $original_inflation;
 main_select_mesh($inflation,$step,$max_inflation);

print "\n$html\n" if $pwd ne '.';
print "(=) done\n";

exit(0);


## main sub

sub main_maketab {
    my $pmids;
    my $pairs;

    my $min_related;
    my $max_related;

    my $xml_filename = $_[0];
    my $list = $_[1];

    $pmids = read_pmids($xml_filename,$list);
    ($pairs,$min_related,$max_related) = create_pairs($pmids,$xml_filename);
    create_mcl_graph($pairs,$min_related,$max_related,$xml_filename);
}


## secondary subs

sub main_select_mesh {
    my ($inflation,$step,$max_inflation) = @_;
    while($inflation <= $max_inflation) {
        select_mesh("$xml_filename.mcl$inflation.clustering.infl_cluster_pmid_mesh");
        $inflation = sprintf "%.1f", $inflation + $step;
    }
}


sub read_pmids {
    my $file_name = $_[0];
    my $list = $_[1];
    my $seen = 0;
    my %pmids;
    $count = 0;

    goto LIST if $list == 1;
    my $reader = new XML::LibXML::Reader(location => $xml_filename) or
        die "Cannot open file $xml_filename for reading\n";

    while ($reader->read) {
        if($reader->name eq "PMID" and $reader->depth == 3 and $reader->nodeType == 1){
            $reader->read;
            ++$count;
            if ($count == $max_documents+1){
                print "(!) we have reached the maximum number of documents ($max_documents), so we will stop reading any more\n";
                goto ENOUGH;
            } else {
                $pmids{$reader->value} = 1;
                $pmid = $reader->value;
            }
        }
        elsif ($reader->name eq "ArticleTitle") {
            $reader->read;
            if($reader->hasValue){
                $title = $reader->value;
                $title =~ s/"//g;
                $titles{$pmid} = "$title";
            } else {
                $titles{$pmid} = "no title";
            }
            $reader->skipSiblings;
        }
        elsif($reader->name eq "DateRevised"){
            $reader->read;
            $reader->read;
        	if($reader->name eq "Year"){
	            $reader->read;
	            $pub_years{$pmid}= $reader->value;
	        }
        }
        elsif ($reader->name eq "PublicationType" and $reader->nodeType == 1) {
            $reader->read;
            unless ($reader->value =~ /Corrected|Duplicate|English Abstract|Journal Article|Research Support/i) {
                $pubtypes{$pmid} .= $reader->value . ";";
            }
        }
    }
    ENOUGH:
    $reader->finish;
    goto TELOS;
    LIST:
    open (XML, "<$file_name") || die "\nProblem reading XML or PMID list file: $!\n\n";
    while(<XML>){
        chomp;
        $pmids{$_} = 1;
    }
    close XML;
    TELOS:

    print "(?) we have " . scalar(keys %pmids) . " documents\n";
    return (\%pmids);
}


sub create_pairs {

    my %pmids = %{ $_[0] };
    my $xml_filename = $_[1];
    my $count = scalar(keys %pmids);

    print "(?) getting related documents and their scores for $count documents in chunks of 500 - please wait\n";

    my $max_related = 0;
    my $min_related = 0;
    my %pairs;

    my @tmp_pmid;
    my $pmid;

    my %chunks_500 = %{ divide_500(\%pmids) };

    foreach my $chunk (sort { $a <=> $b } keys %chunks_500){
        #print " $chunk\n" if keys(%chunks_500) > 1;
        my $file_name = "$xml_filename.chunk$chunk.related";
        $pmid = join("&id=",@ { $chunks_500{$chunk} });
        if(-e $file_name) {
            #print "(?) skipping downloading chunk number: $chunk\n";
            goto PAIRS;
        }
        my $query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=$pmid&cmd=neighbor_score";
        my $related_results = getstore($query,$file_name);
        #print "$related_results\n";
        PAIRS:
        open(RELATED,"<$file_name") or die "Cannot open file $file_name for reading: $!\n";
        my $main = 0;
        my $related_id = 0;
        my $first_time = 0;

        while (<RELATED>){
            chomp;
            if(/<LinkSet>/){
                $main = 1;
            }
            elsif (/<Id>(\d+)<\/Id>/ and $main == 1){
                if(!defined($pairs{$pmid}) and $first_time == 1){
                  $pairs{$pmid}{$pmid} = 0;
                }
                $pmid = $1;
                $main = 0;
                $first_time = 1;
                next;
            }
            elsif(/<Id>(\d+)<\/Id>/ and $main == 0){
                $related_id = $1;
            }
            elsif(/<Score>(\d+)<\/Score>/){
                if(defined($pmids{$related_id}) && $pmid != $related_id){
                    $pairs{$pmid}{$related_id} = $1;
                    if ($1 < $min_related or $min_related == 0) { $min_related = $1; }
                    if ($1 > $max_related or $max_related == 0) { $max_related = $1; }
                }
            }
        }
        if(!defined($pairs{$pmid}) and $first_time == 1){
            $pairs{$pmid}{$pmid} = 0;
        }
    }
    if ($min_related == 0 && $max_related == 0) { print "(!) there was an error getting related documents and their scores - exiting\n"; exit;
    } else {
        #print " ok\n" if keys(%chunks_500) > 1;;
        return (\%pairs,$min_related,$max_related);
    }
}


sub create_mcl_graph{
    my %pairs = %{ $_[0] };
    my $min_related = $_[1];
    my $max_related = $_[2];
    my $file_name = $_[3];

    print "(?) the minimum related article score is $min_related and the maximum is $max_related\n";
    print "(?) we have " . scalar(keys %pairs) . " unique MeSH terms\n";

    open (BLO, ">$file_name.mcli") || die "\nCannot create $file_name.mcli\n\n"; #the file describing the graph in MCL format(space instead of tab)

    foreach my $pmid (sort {$a<=>$b} keys %pairs) {
        foreach my $linkpmid (sort {$a<=>$b} keys %{$pairs{$pmid}}) {
            if ($pmid ne $linkpmid) {
                my $related_score = (($pairs{$pmid}{$linkpmid}-$min_related)/($max_related-$min_related))*100;
                print BLO "\"$pmid\" \"$linkpmid\" $related_score\n";
            } else {
                print BLO "\"$pmid\" \"$linkpmid\" 0\n";
            }
        }
    }
    close BLO;
}


sub divide_500 {
    my %pairs = %{ shift() };
    my @tmp_pmids = keys %pairs;
    my $count = 0;
    my %chunks_500;
    my $i = 1;

    while (@tmp_pmids) {
        $count = scalar(@tmp_pmids);
        if ($count <= 500) {
            $chunks_500{$i} = [ splice(@tmp_pmids,0,scalar(@tmp_pmids)) ];
            return (\%chunks_500);
        }
        $chunks_500{$i} = [ splice(@tmp_pmids,0,500) ];
        $i++;
    }
}


## helper subs

sub main_mcl_evaluation {
    my $xml_filename = $_[1];
    my $mcli_filename = $_[0];
    my $inflation = $_[2];
    my $step = $_[3];
    my $max_inflation = $_[4];

    mcl($mcli_filename,$xml_filename,$inflation,$step,$max_inflation);
    #evaluation($mcli_filename,$inflation,$step,$max_inflation);
}


sub filename {
    my $filename = shift;
    @_ =  split(/\./,$filename);
    $filename = $_[0];
    return $filename;
}


sub mcl {
    my $ori_mcli_filename = shift;
    my $mcli_filename = $ori_mcli_filename;
    my $xml_filename = shift;
    my $inflation = shift;
    my $step = shift;
    my $max_inflation = shift;

    sleep 1;

    #The file must be in space seperated format in the form of NODE_1 NODE_2 Weight

    print "(?) running MCL with inflation $inflation\n";
    `mcl $ori_mcli_filename --abc -scheme 6 -I $inflation -o $xml_filename.mcl$inflation.clustering -write-graph $mcli_filename.mclgraph 2> $xml_filename.mcl$inflation.clustering.log`;
    #  `clmformat -dump label_cluster$inflation -icl $mcli_filename.cluster$inflation -tab $mcli_filename.map`;

    if (-e "$xml_filename.mcl$inflation.clustering.log") {
        mesh_chemical($inflation,"$xml_filename.mcl$inflation.clustering",$xml_filename);
    } else { print "(!) running MCL failed, is mcl installed? - exiting\n"; exit; }

    $inflation = sprintf "%.1f", $inflation + $step;
    while ($inflation <= $max_inflation) {
        print "(?) running MCL with inflation $inflation\n";
        `mcl $ori_mcli_filename --abc -scheme 6 -I $inflation -o $xml_filename.mcl$inflation.clustering 2> $xml_filename.mcl$inflation.clustering.log`;
        #    `clmformat -dump label_cluster$inflation -icl $mcli_filename.cluster$inflation -tab $mcli_filename.map`;
        mesh_chemical($inflation,"$xml_filename.mcl$inflation.clustering",$xml_filename);
        #    $inflation =sprintf "%.1f", $inflation + $step;
        $inflation = sprintf "%.1f", $inflation + $step;
    }
    return;
}


#sub evaluation {
#    my $ori_mcli_filename = shift;
#    my $mcli_filename = $ori_mcli_filename;
#    my $inflation = shift;
#    my $step = shift;
#    my $max_inflation = shift;
#
#    print "(?) performing evaluation tests\n";
#    open(EVALRES,">$xml_filename.mcl$inflation.clustering.evaluation") or die "Cannot write the evaluation results: $!\n";
#    sleep 1;
#
#    my $i = 0;
#    my $maxclusters = 0;
#    my $maxinflation = 0;
#    my $maxefficiency = 0;
#
#    while($inflation <= $max_inflation){
#        my $clminfo = `clminfo $mcli_filename.mclgraph $xml_filename.mcl$inflation.clustering > $xml_filename.mcl$inflation.evaluation.log`;
#        $clminfo =~ /clusters=(\d+)/;
#        my $current_cluster = $1;
#        $clminfo =~ /efficiency=(\d+\.\d+)/;
#        my $current_efficiency = $1;
#        print EVALRES "Inflation_$inflation\tNo_Clusters: $current_cluster\tEfficiency: $current_efficiency\n";
#        if($current_efficiency > $maxefficiency){
#            $maxefficiency = $current_efficiency;
#            $maxinflation = $inflation;
#            $maxclusters = $current_cluster;
#        }
#        $inflation = sprintf "%.1f", $inflation + $step;
#    }
#
#    print EVALRES "\n(=) maximum efficiency is $maxefficiency at inflation $maxinflation, which corresponds to $maxclusters number of clusters\n";
#    close EVALRES;
#
#    return;
#}


sub select_mesh {
    $| = 1;

    my %clust_pmid_mesh = (); #holds the info from the input file (Inflation\tCluster_id\tPMID\tMeSH)
    my %clust_mesh_pmid = (); # data structure to hold the pmids foreach MeSH term in a cluster
    my %cluster_size = (); # the pmids for each cluster
    my %mesh_in_cluster = (); # the number of documents each mesh term belongs inside a cluster
    my %mesh_across_clusters = (); # the number of documents containing each MeSH term across all clusters
    my %mesh_scores = (); # the scores for each MeSH term

    my $file_name = shift;
    $file_name =~ /\.mcl(\d+\.*\d*)\./;
    $inflation = $1;
    my $mesh_selected_file = "$xml_filename.mcl$inflation.out";
    my $total_pmids = 0; # total number of articles
    my $no_clusters = 0; # total number of clusters

    open (FILE,"<$file_name") or die "Cannot open file $file_name for reading: $!\n";
    while (<FILE>){
        chomp;
        @_ = split(/\t/,$_);
        $clust_pmid_mesh{$_[1]}{$_[2]}{$_[3]} = 1;
    }
    close FILE;

    print "(=) at inflation $inflation MCL returned " . scalar(keys %clust_pmid_mesh) . " clusters\n";

    $singletons = 0;
    foreach my $cluster (sort keys %clust_pmid_mesh){
        $cluster_size{$cluster} = scalar(keys %{ $clust_pmid_mesh{$cluster} });
        $total_pmids += $cluster_size{$cluster};
        $no_clusters++;
        foreach my $pmid (keys %{ $clust_pmid_mesh{$cluster} }){
            foreach my $mesh (keys %{ $clust_pmid_mesh{$cluster}{$pmid} }) {
                $mesh_in_cluster{$cluster}{$mesh}++;
                $mesh_across_clusters{$mesh}++;
                $clust_mesh_pmid{$cluster}{$mesh}{$pmid}++;
            }
        }
        if ($cluster_size{$cluster} == 1) {
            ++$singletons;
        } else {
            print "\t(=) cluster $cluster has $cluster_size{$cluster} documents\n";
        }
    }
    print "\t(=) we also found $singletons clusters with one document each\n" if $singletons > 0;
    print "\t(=) in total, we have $total_pmids documents with " . scalar(keys %mesh_across_clusters) . " MeSH terms\n";

    foreach my $cluster (keys %cluster_size) {
        foreach my $mesh (keys %{ $mesh_in_cluster{$cluster} }) {
            if($mesh =~ /no MeSH terms/){
                $mesh_scores{$cluster}{0}{$mesh} = 1;
                next;
            }
        # if the pmid does not belong to a class. It is the only article in a cluster.
            my $prob_in_cluster = $mesh_in_cluster{$cluster}{$mesh} / $cluster_size{$cluster};  # frequency of MeSH term inside  cluster
            my $prob_across_clusters = -log2($mesh_across_clusters{$mesh} / $total_pmids);      # frequency of MeSH term across clusters
            my $score = $prob_in_cluster * $prob_across_clusters;
            $mesh_scores{$cluster}{$score}{$mesh} = 1;
        }
    }
    my $prev_mesh = undef;
    my $pmids_seen = 0;

    @split_xml_filename = (split /\//, $xml_filename);
    $xml_filename4html = pop(@split_xml_filename);

    open (SELECT, ">$mesh_selected_file.html") or die "Cannot open file for writing selected MeSH terms: $!\n";
    print SELECT <<HTMLCODE;
<HTML>
<HEAD>
<TITLE>PuReD-MCL @ the BAT cave | results</TITLE>
</HEAD>
<BODY>
<pre>

# PuReD-MCL - clustering of related PubMed documents
# Theodosiou T. and Darzentas N.
# $version
# http://bat.infspire.org | bat\@infspire.org
#
# query | $query
# inflation | $inflation
#
# MeSH terms ordered by their score (freq in cluster * freq across clusters),
# followed by any PMIDs they add to the cluster

HTMLCODE
    foreach my $cluster (sort keys %mesh_scores){
        my %mem = ();
        my %tmpmem = ();
        my $cur_pmids = 0;
        my $previous_pmids = 0;
        my $diff = 0;

        $mem4print = "$cluster\n";
        foreach my $score (sort {$b<=>$a} keys %{ $mesh_scores{$cluster} }){
            $score2print = sprintf "%.2f", $score;
            foreach my $mesh (sort keys %{ $mesh_scores{$cluster}{$score} }){
                unless ($mesh eq 'no MeSH terms') {
                    $mem4print .= qq(\t<a href=http://www.ncbi.nlm.nih.gov/mesh?term=$mesh>$mesh</a>  {$score2print}\n);
                } else {
                    $mem4print .= qq(\t$mesh\n);                        
                }
                foreach my $pmid (keys %{ $clust_mesh_pmid{$cluster}{$mesh} }){
                    $tmpmem{$pmid} = 1;
                }
                $cur_pmids = scalar(keys %tmpmem);
                    $diff = $cur_pmids - $previous_pmids;
                if ($diff >= 1) {
                    foreach my $pmid (sort keys %{ $clust_mesh_pmid{$cluster}{$mesh} }){
                        $pubtypes{$pmid} =~ s/;$//;
                        $mem4print .= qq(\t\t<a href="http://www.ncbi.nlm.nih.gov/pubmed?term=$pmid">$pmid</a>\t$pub_years{$pmid}\t$titles{$pmid}\t$pubtypes{$pmid}\n);
                        $mem{$pmid} = 1;
                    }
                    $previous_pmids = $cur_pmids;
                } else {
                    %tmpmem = %mem;
                }
            }
        }
        print SELECT "$mem4print\n<hr>\n";
    }
    print SELECT <<HTMLCODE;
</pre>
</BODY>
</HTML>
HTMLCODE
    close SELECT;

    $html .= "\t<a href=\"\/pured-mcl\_results\/$xml_filename4html.mcl$inflation.out.html\">[HTML] clusters, their MeSH terms, and their PubMed IDs</a>\n";
}


sub log2 {
    my $n = shift;
    return log($n)/log(2);
}


sub mesh_chemical {
    my ($inflation_no,$clust_filename,$mesh_filename) = @_;

    open (CLUSTERS, "$clust_filename") or die "Cannot open file $clust_filename for reading. $!\n";

    my %clusters;
    my %mesh_terms;
    my %tfdf;
    my $count = 1;
    my $inflation = "mcl$inflation_no";

    while (<CLUSTERS>){
        chomp;
        $count2print = sprintf "%06d", $count;
        my $cur_cluster = "PuReD-MCL-$inflation-$count2print";
        my $members = $_;
        $members =~ s/"//g;
        my @tmp = split /\t/,$members;
        foreach my $label (@tmp) {
            $clusters{$cur_cluster}{$label} = 1;
        }
        $count++;
    }

    close CLUSTERS;

    my $pmid = "";

    my $reader = new XML::LibXML::Reader(location => $mesh_filename) or
    die "Cannot open file $xml_filename for reading\n";

    while($reader->read){
        if($reader->name eq "PMID" and $reader->depth == 3 and $reader->nodeType == 1) {
            $reader->read;
            if(!exists $mesh_terms{$pmid}){
                $mesh_terms{$pmid}{"no MeSH terms"} = 1;
            }
            $pmid = $reader->value;
        }
        elsif ($reader->name eq "DescriptorName" or $reader->name eq "QualifierName" or $reader->name eq "NameOfSubstance") {
            $reader->read;
            if($reader->hasValue){
                if ($reader->value =~ /\w/) {
                    $mesh_terms{$pmid}{$reader->value}++;
                    if ($mesh_terms{$pmid}{$reader->value} < 2) {
                        $tfdf{$reader->value}++;
                    }
                }
            }
        }
    }

    $reader->finish;

    open(PM,">$clust_filename.infl_cluster_pmid_mesh") or die "Cannot create file $clust_filename.infl_cluster_pmid_mesh. $!\n";
    foreach my $cluster (sort keys %clusters){
        foreach $pmid (sort keys %{ $clusters{$cluster} }) {
            foreach my $mesh_term (sort keys %{ $mesh_terms{$pmid} }) {
                print PM "$inflation\t$cluster\t$pmid\t$mesh_term\n";
            }
        }
    }
    close PM;
}
