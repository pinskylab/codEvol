#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

dummyGT ANGD.vcf

Description:

Adds dummy GTs

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;

my $file = shift;
die $usage unless $file;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

while (<$IN>) {
    chomp;
    if($_ =~ /^\#/){
	print "$_\n";
    }
    else{
	my @line  = split /\t/, $_;
    
	$line[8] = "GT:DP:AD:GP:GL";

	print join "\t", @line[0..8];

	foreach my $g (@line[9..scalar(@line)-1]){
	    print"\t1|1:$g"
	}
	print "\n";
    }
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

