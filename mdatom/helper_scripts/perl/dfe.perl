#!/usr/bin/perl

#set up some timing stuff
#require 'sys/#syscall.ph';

# open in- and output file
open(INPUT, "< $ARGV[0]") || die("can't open ",$ARGV[0]);

# set commands for execution
$MDATOM="/root/mdatom/C/greg/promdac.exe";
$MDORIG="/root/mdatom/C/orig/promdac.exe";


# set standard output file, overridden when series is made
$label="output";
$series=0;
$stopwatch=0.0;

# make a list for sorting the parameters
$rank{'Title'}=0;
$rank{'NAT'}=1;
$rank{'NTXI'}=2;
$rank{'NTXO'}=3;
$rank{'BOX[1]'}=4;
$rank{'BOX[2]'}=5;
$rank{'BOX[3]'}=6;
$rank{'NBOX[1]'}=7;
$rank{'NBOX[2]'}=8;
$rank{'NBOX[3]'}=9;
$rank{'DX'}=10;
$rank{'IG'}=11;
$rank{'TEMPI'}=12;
$rank{'NTT'}=13;
$rank{'TEMP0'}=14;
$rank{'TAUT'}=15;
$rank{'BOLTZ'}=16;
$rank{'NSTLIM'}=17;
$rank{'T'}=18;
$rank{'DT'}=19;
$rank{'AMAS'}=20;
$rank{'EPSLJ'}=21;
$rank{'SIGLJ'}=22;
$rank{'RCUTF'}=23;
$rank{'NTPR'}=24;
$rank{'NTWX'}=25;
$rank{'NTWXM'}=26;
$rank{'NTPW'}=27;
$rank{'NGR'}=28;
$rank{'RCUTG'}=29;        
$rank{'UF'}=30;
# go through every line of the input file
while(<INPUT>){

# match anything beginning with a "=" and ending with a ";" 
        /(.*)=(.*);/;
# $<digits> corresponds to what matched the () in regexp
	$param=$1;
	$value=$2;
# check if a series shall be calculated: match xxx:xxx:xxx
	if($value =~ /(.*):(.*):(.*)/ ) {
	  $varparam=$param;
	  $from=$1;
	  $step=$2;
	  $endval=$3;
	  $label="var";
	  $series=1;
	}
	$file_list[$rank{$param}]=$value;
}

#make prefix if given as argument
if ($#ARGV==1){
  $label=$ARGV[1]."/".$label;
}

#save old stdout
open(OLDOUT,">&STDOUT");

#open timing file
open(TIMING,">$label.timing") || die("can't open timing file $label.timing");
open(OTIMING,">$label.otiming") || die("can't open otiming file $label.timing");

#do the iteration if series else just execute once
if ($series){
  for($i=$from;$i<=$endval;$i+=$step){
    open(STDOUT,">&OLDOUT");
    printf(STDOUT "Running with %s=%f\n",$varparam,$i);
    $file_list[$rank{$varparam}]=$i;
    open(DATAFILE, "> mdtemp") || die("can't open temporary file mdtemp");
    printf(DATAFILE "%s\n%d %d %d\n%f %f %f\n%d %d %d %f %d %f\n%d %f %f %f\n%d %f %f\n%f %f %f %f\n%d %d %d %d\n%d %f\n%d", @file_list);
    $outputfilename=$label."-".$varparam."=".$i;    
    close(DATAFILE);

    open(STDOUT, "> $outputfilename") || die("can't redirect stdout");

#do the modified version
    @command=($MDATOM,"mdtemp");
 
    $TIMEVAL_T = "LL";
 
    $done = $start = pack($TIMEVAL_T, ());


    ##syscall( &SYS_gettimeofday, $start, 0) != -1
    #  or die "gettimeofday: $!";

    system(@command);

    #syscall( &SYS_gettimeofday, $done, 0) != -1
      #or die "gettimeofday: $!";

    @start = unpack($TIMEVAL_T, $start);
    @done  = unpack($TIMEVAL_T, $done);

    # fix microseconds
    for ($done[1], $start[1]) { $_ /= 1_000_000 }
    
    $delta_time = sprintf "%.4f", ($done[0]  + $done[1]  )
      -
	($start[0] + $start[1] );
    printf(TIMING "%f\t%f\n", $i,$delta_time);

#do the original version

    open(STDOUT, "> $outputfilename.orig") || die("can't redirect stdout");
    @command=($MDORIG,"mdtemp");
 
    $TIMEVAL_T = "LL";
 
    $done = $start = pack($TIMEVAL_T, ());

    #syscall( &SYS_gettimeofday, $start, 0) != -1
     # or die "gettimeofday: $!";

    system(@command);

    #syscall( &SYS_gettimeofday, $done, 0) != -1
    #  or die "gettimeofday: $!";

    @start = unpack($TIMEVAL_T, $start);
    @done  = unpack($TIMEVAL_T, $done);

    # fix microseconds
    for ($done[1], $start[1]) { $_ /= 1_000_000 }
    
    $delta_time = sprintf "%.4f", ($done[0]  + $done[1]  )
      -
	($start[0] + $start[1] );

    printf(OTIMING "%f\t%f\n", $i,$delta_time);

  }
}

else{
    open(STDOUT,">&OLDOUT");
    open(DATAFILE, "> mdtemp") || die("can't open temporary file mdtemp");
    printf(DATAFILE "%s\n%d %d %d\n%f %f %f\n%d %d %d %f %d %f\n%d %f %f %f\n%d %f %f\n%f %f %f %f\n%d %d %d %d\n%d %f\n%d", @file_list);
    close(DATAFILE);
    close(STDOUT);
    open(STDOUT, "> $label") || die("can't redirect stdout");
    @command=($MDATOM,"mdtemp");
    print @command;
    system(@command);
}
close(DATAFILE);
close(INPUT);

