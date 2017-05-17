# PuReD-MCL
clustering of related PubMed documents

</br>
server: http://bat.infspire.org/tools/pured-mcl/</br>
citation: http://www.ncbi.nlm.nih.gov/pubmed/18593717</br>
code: we provide the main Perl script, and the server's CGI script and HTML page - you'll need MCL and few Perl modules</br>
usage:
<pre>
before you run PuReD-MCL, be sure to:
- have Stijn van Dongen's MCL installed and in your path
- have installed the following Perl modules:
	Getopt::Long
	LWP::Simple
	XML::LibXML::Reader

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
</pre>
