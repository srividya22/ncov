#!/usr/bin/env perl -w
use strict;

#initialling variables:
my $fh = undef;
my $seqs_per_file = 10000;
my $file_number = 1;
my $curseq = undef;
my $filename = join('',$ARGV[1],"/summary");
my $name=undef;
#open files:
open (my $summary, ">", $filename) || die "Could not open file '$filename: $!\n";

#Fasta file:
open(FILE1, $ARGV[0]) || die "Could not open file '$ARGV[0]: $!\n";
my $count = 0;

while( <FILE1> ){
    chomp $_;
       if($_=~ /^>/g ) 
       {
             $count++; 
    }
}

print STDERR "Total number of sequences : $count\n";
seek FILE1, 0, 0;
my $seq_ctr = 0;
## If the number of sequneces is less than 100 then split to individual contigs
if($count < $seqs_per_file ) {
  while (<FILE1>) #reads in line after line
  {
  	if (/>(\S+)/) #If fasta header:
        {
        $curseq = $1;
        print STDERR "\tSplitting $curseq\n";
        $name=join('',$ARGV[1],"/$curseq.fa");
        open $fh, ">",$name or die "Cant open $curseq.fa ($!)\n";
        print $fh ">$curseq\n"; #print heaer in new file
        print $summary "$name\n"; #adds the name of the file to the summary file
        }  
        else{
           print $fh $_;	
        }
  } 
 }
else{ #If the number of sequneces is greater than 100 then split to chunks of 100 sequences
   while (<FILE1>) #reads in line after line
   {
        if (/>(\S+)/) {  #If fasta header:
 
        if ( $seq_ctr++ % $seqs_per_file == 0) {
        my @arr=split /(\d+)/, $1;
        my $seq_base=$arr[0];
        my $end_contig = $seq_ctr + 99;
        if ( $end_contig  > $count ){ $end_contig = $count; }  
        $curseq = $seq_base .$seq_ctr. "-" . $end_contig;
        $name=join('',$ARGV[1],"/$curseq.fa");
        #close($fh);
    	open $fh, ">",$name or die "Cant open $curseq.fa ($!)\n";

    	#print $fh ">$curseq\n"; #print heaer in new file
	
    	print $summary "$name\n"; #adds the name of the file to the summary file
        }
        }
        print $fh $_;
    }	
} 
