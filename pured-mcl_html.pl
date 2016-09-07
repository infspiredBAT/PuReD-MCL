#!/usr/bin/perl -w
$| = 1;

use CGI;
use LWP::Simple;

$randomer = int(rand(1001)) . $$;
$randomer = sprintf "%010s", $randomer;

$names{official} = 'PuReD-MCL';
$names{unix} = 'pured-mcl';
$names{description} = 'clustering of related PubMed documents';
#$safe_filename_characters = "a-zA-Z0-9_.-";  
$upload_dir = "/labs/bat/www/tools/results/$names{unix}";
$results_filename = "$upload_dir/out$randomer";

undef($/); # read whole file into a variable    

$max_doc = 2000;
 
$cgi = new CGI;
$cgi->cgi_error and error( $cgi, "error transferring file: " . $cgi->cgi_error );
$query              = $cgi->param( "query" ); # || error( $cgi, "no query FASTA received" );
unless ($query) {
    error( $cgi, "you have provided no input" );
}
if ($query =~ /href ?=/) {
    error( $cgi, "if you are not a spammer, remove href tag and have another go" );    
}

if ($cgi->param( "count" )) {
    my $utils   = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
    my $db      = "Pubmed";
    my $report  = "xml";
    my $esearch = "$utils/esearch.fcgi?" .
                  "db=$db&usehistory=y&term=";
    my $esearch_result = get($esearch . $query);
    $esearch_result =~ 
      m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;    
    my $count    = $1;
    my $queryKey = $2;
    my $webEnv   = $3;
print "Content-type: text/html\n\n"; #needed for CGI scripts
$index_html = do {local (@ARGV,$/) = "/labs/bat/www/tools/www/pured-mcl/index.html"; <>}; # file slurp!
# before
# <input type="text" size="80" name="query" /> <input type="submit" name="count" value="see count before submit" />
# after
# <input type="text" size="80" name="query" value="$query" /> <input type="submit" name="count" value="see count before submit" /> = $count
$index_html =~ s/<input type="text" size="80" name="query" \/> <input type="submit" name="count" value="see count before submit" \/>/<input type="text" size="80" name="query" value="$query" \/> <input type="submit" name="count" value="see count before submit" \/> = $count/;
print <<HTML;
$index_html
HTML
    exit;
}

print "Content-type: text/html; charset=iso-8859-1\n\n";
print <<HTML;
<html>
<head>

<pre>
<link rel="stylesheet" href="/css/blueprint/screen.css" type="text/css" media="screen, projection, print"/>
<link rel="stylesheet" href="/css/css4tools.css" type="text/css" />
</pre>

<title>$names{official} @ the BAT cave | results</title>
</head>
<body>

<div class="container" style="opacity:0.92;filter:alpha(opacity=92)">

<div id="header" class="span-24 last" style="background-color:white">
<h1 style="text-align:center">$names{official}</h1>
<h2 style="text-align:center">$names{description}</h2>
</div>

<div id="content" class="span-24 last">
<fieldset>
<p style="text-align:left">

HTML

$inflations         = $cgi->param( "inflations" );
@split_inflations   = (split /-/, $inflations);
$min_inflation      = $split_inflations[0]; 
$max_inflation      = $split_inflations[1]; 
$step               = $split_inflations[2]; 

print "<label>running $names{official} - please wait...</label></br>";

$log = "/labs/bat/www/tools/logs/$names{unix}.log";
$max_time = 600;
$max_mem = 500000000;

`perl /labs/bat/www/tools/cgi-bin/limitresources.pl --command "perl pured-mcl.pl --query '$query' --min $min_inflation --max $max_inflation --step $step --max_doc $max_doc --pwd $upload_dir" --log $log --seconds $max_time --memory $max_mem > $results_filename &`;

print "\n<pre>\n";
$progress = '';
$done = 0;
$naptime = 0.1;
while ($progress !~ /done/){
    select(undef, undef, undef, $naptime);
    ++$time;
    open (OUT, "<$results_filename") or die "cannot read output: $!\n";
    $progress = <OUT>;
    @split_progress = (split /\n/, $progress);
    foreach $line (@split_progress) {
        unless (defined($printed{$line})) {
            if ($line =~ /done/) {
#                @split_line = (split / /, $line);
#                $gff = $split_line[7];
#                $input = (split /\./, $gff)[0];
                $done = 1;
            } elsif ($line =~ /ERROR/ || $line =~ /exiting/) {
                print "$line</br>";
                goto MOVEON;
            } elsif ($line =~ /^html/) {
                $line =~ s/^html //g;
                $files .= "$line</br>";
            } elsif ($line !~ /[A-Za-z]/) {
                print "$line";
            } else {
                print "$line</br>";
            }
            $printed{$line} = 1;
        }
    }
    close OUT;
    if ($time*$naptime > $max_time*3) {
        print "(!) this run has timed out</br>";
        goto MOVEON;
    }
}
MOVEON:
print "\n</pre>\n";
if ($done) {
    print
    "
    $files
    ";
}

print <<HTML;
</p>

</fieldset>

<p class="small quiet terms" style="text-align:right">
hosted at the <a href="http://bat.infspire.org/" target="_blank"> Bioinformatics Analysis Team / BAT</a>
</p>

</div>
</div>

</body>
</html>

HTML

exit;

sub error {
    my( $cgi, $reason ) = @_;
    
    print $cgi->header( "text/html" ),
          $cgi->start_html( "error @ the BAT cave" ),
          $cgi->h1( "error @ the BAT cave" ),
          $cgi->p( "your upload was not processed because the following error occured: " ),
          $cgi->p( $cgi->i( $reason ) ),
          $cgi->end_html;
    exit;
}
