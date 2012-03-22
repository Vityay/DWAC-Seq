#!/usr/bin/perl -w

use strict;
use IO::Handle;
use Getopt::Long;
use List::Util qw/ min max sum /;
use Math::Complex;

# Last mod 20110409

# Parameters section
my $nOfReadsInWindow = 250; # number of reads in a window
my $bootstraps       = 1000; # the number of bootstraps (the more the better)
my $window           = 100; #size of the window (of windows) for segmentation
my $step             = 9; #step of slicing window (the less the better)
my $conf_level       = 0.95; #confidence level in change points finding
my $cnv_lower_ratio  = 0.8; #Do not fine-tune if ratio is between _lower_ and _higher_ 
my $cnv_higher_ratio = 1.2;
my $read_length      = 48; # Define read length of NGS read

my $usage = "options:--test test_bam
	--ref ref_bam
	--chr chromosome_name
	--amb ambiguity (optional, default 10)
	--scale number_of_reads_in_window (optional 1000 default)
	--no_fine_tuning (the results are not fine tuned if present)";

my $nOfParameters = @ARGV;
die($usage) unless $nOfParameters >= 3;
STDOUT->autoflush(1);

my $no_fine_tuning = 0; # by default do fine tuning
my $ambiguity      = 1; #default ambiguity

my ($test_bam, $ref_bam, $target_chr);

GetOptions ( 'no_fine_tuning+' => \$no_fine_tuning,
             'test=s'          => \$test_bam,
             'ref=s'           => \$ref_bam,
             'chr=s'           => \$target_chr,
             'amb:i'           => \$ambiguity,
             'scale:i'         => \$nOfReadsInWindow
           );

die("no chromosome specified") unless $target_chr;

# extract chromosome $target_chr from the BAM files and create output file


my $test = $test_bam;
$test =~ s/\.bam$//;
my $ref = $ref_bam;
$ref =~ s/\.bam$//;

my $output_name = join( '_', 'tst', $test,
                             'ref', $ref,
                             'amb', $ambiguity,
                             'win', $nOfReadsInWindow,
                             'chr', $target_chr,
#                             'TEST'
                      );


my @ref  = ();
warn "Reading REF\n";
#open (F, "samtools view $ref_bam $target_chr:1-10000000 |") or die "No test BAM given";
open (F, "samtools view $ref_bam $target_chr |") or die "No test BAM given";
my $i = 0;
while (<F>){
	next unless m/X0\:i\:(\d+)/;
	next unless $1 <= $ambiguity;
	my @values = split /\t/;
	push @ref, $values[3];
	print "\r",$i/1_000_000,"M reads " if ++$i % 1_000_000 == 0;
}
print "\rTotal: $i reads\n";
close F;

my @ref_copy = @ref;

die("No reference reads on this chromosome or not sorted BAM file") unless $#ref > 0;
my @test = ();
warn "Reading TEST\n";
$i = 0;
#open (F, "samtools view $test_bam $target_chr:1-10000000 |") or die "No reference BAM given";
open (F, "samtools view $test_bam $target_chr |") or die "No reference BAM given";
while (<F>){
	next unless m/X0\:i\:(\d+)/;
	next unless $1 <= $ambiguity;
	my @values = split /\t/;
	push @test, $values[3];
	print "\r",$i/1_000_000,"M reads " if ++$i % 1_000_000 == 0;
}
print "\rTotal: $i reads\n";
close F;

my @test_copy = @test;

die("No test reads on this chromosome or not sorted BAM file") unless $#test > 0;

# step 1 - counting number of reads per window

my $window_start = $ref[0];
my $window_end   = -1;
my $read;
my $test_reads = 0;
my $ref_reads  = 0;


# Setting windows' boundaries
warn "Setting boundaries of dynamic windows\n";
my %windows = ();
foreach my $ref_pos ( @ref ) {
	die "Reference BAM is not sorted" if ($ref_pos + 1) < $window_start; 
	$ref_reads++;
	if ( $ref_reads >= $nOfReadsInWindow ) {
	    $window_end = $ref_pos;
	    $windows{$window_start} = $window_end;
	    $window_start = $ref_pos + 1;
	    $ref_reads = 0;
	}
}


warn "Counting reads in windows for test set\n";
my %test_win  = ();
my @starts = sort {$a<=>$b} keys %windows;
foreach my $test_pos ( @test ) {
    last unless @starts;
    if ( $test_pos < $starts[0] ) {
        ;
    } 
    elsif ( $test_pos <= $windows{$starts[0]} ) {
        $test_win{$starts[0]}++;
    }
    else {
        shift @starts;
        redo;
    }
}

# Counting number of reads in a reference window
warn "Counting reads in windows for reference set\n";
my %ref_win = ();
my @full = ();
@starts = sort {$a<=>$b} keys %windows;
foreach my $ref_pos ( @ref ) {
    next unless @starts;
    if ( $ref_pos < $starts[0] ) {
        ;
    } 
    elsif ( $ref_pos <= $windows{$starts[0]} ) {
        $ref_win{$starts[0]}++;
    }
    else {
        $test_win{$starts[0]} = 0 unless exists($test_win{$starts[0]}); 
        push @full, $test_win{$starts[0]} / $ref_win{$starts[0]};
        shift @starts;
        redo;
    }
}

# Step 2 - segmentation
warn "Segmenting coverage profile\n";
open F, '>', 'breakpoints.'.$output_name;
my $data_size = scalar keys %windows;
@starts = sort {$a<=>$b} keys %windows;
my %change_points = ();
my $progress_prev = 0;

for ( my $counter = 0; $counter < $data_size - 1; $counter += $step ) {
    print "\rwindow ",$counter+1,' / ', $data_size, ' ';
    my $s = $counter;
    my $e = $counter + $window - 1;
    next if $e > $#starts;
    my ( $point, $conf ) = find_cp( $s, $e );
    if ( $conf >= $conf_level and ( !exists($change_points{ $point }) or $change_points{ $point } < $conf ) ) { # change point is located after window with number $point
        $change_points{ $point } = $conf;
        my $start = $starts[$point];
        print "Found breakpoint at $windows{$start} confidence $conf\n";
        print F $windows{$start}, "\t", $conf, "\n";
    }
}
close F;

## add last point
$change_points{ $data_size-1} = 1;

# Averaging of copy-number for each segment
print "Averaging ratio between breakpoints\n";
my @readsPerWAveraged = ();
my $nOfWindows = scalar keys %windows;
my $previousChangePoint = -1; # before window 0
foreach  my $currentCP ( sort { $a <=> $b } keys %change_points ) {
    # average from previous change point
    my @values = ();
    foreach my $k ( ($previousChangePoint + 1) .. $currentCP ) {
        die "Attempt to get over end of chromosome" if $k > $#starts;
        $test_win{$starts[$k]} = 0 unless exists($test_win{$starts[$k]});
        push @values, $test_win{$starts[$k]} / $ref_win{$starts[$k]};
    }
    my $median = median(\@values);
    # save averaged data
    for my $k ( ( $previousChangePoint + 1) .. $currentCP ) {
        $readsPerWAveraged[$k] = $median;
    }
    $previousChangePoint = $currentCP;
}

my $med = median(\@readsPerWAveraged);
print "Correcting ratios using median ratio ($med)\n";
my @readsNormalized = ();
foreach ( @readsPerWAveraged ) { 
    my $rounded =  sprintf ( "%.3f", $_/$med );
    push ( @readsNormalized, $rounded );
}

print "Merging adjacent regions that have similar copy number\n";
my %data_out = ();
my $r = $readsNormalized[0];
my $be = 0;
my $e;
foreach  my $counter ( 1 .. $nOfWindows - 1 ) {
    next unless abs( $readsNormalized[$counter] - $r ) > 0.1;
    $e = $counter;
    $data_out{ $be."\t".($e - 1) } = $r > 0 ? $r : 0.00001; #otherwise logn will fail
    $r = $readsNormalized[$counter];
    $be = $counter;
}
# ... and the last segment
$data_out{ $be."\t". $#starts} = $r > 0 ? $r : 0.00001;

print "Saving results\n";

my %data_out_sorted;

foreach (keys %data_out) {
    my ( $s, $e ) = split /\t/;
    $data_out_sorted{$s} = join("\t", $data_out{ $s."\t".$e}, $starts[$s], $windows{$starts[$e]});
}

open F,  '>', 'segmented.'.$output_name;
foreach (sort{$a<=>$b} keys %data_out_sorted){
    print F "$data_out_sorted{$_}\n";
}
close F;

open F, '>', 'windows.'.$output_name;
open F1, '>', 'avg4hilbert.'.$output_name;
foreach my $i( 0 .. $nOfWindows-1 ) {
     print F  join( "\t", $starts[$i], $windows{$starts[$i]}, $test_win{$starts[$i]}, $ref_win{$starts[$i]} ),"\n";
     print F1 join( "\t", $starts[$i], $windows{$starts[$i]}, $readsPerWAveraged[$i]),"\n";
}
close F;
close F1;

# Do fine-tuning when requested
exit if $no_fine_tuning;
warn "Fine-tuning CNV boundaries\n";

my %seen = ();

## process windows sizes 
open F, '>', 'tuned.'.$output_name;
my $nOfSegments = scalar keys %data_out;
foreach my $range ( sort { abs(logn($data_out{$b},2)) <=> abs(logn($data_out{$a},2)) } keys %data_out ) {
    my ( $swin, $ewin ) = split /\t/, $range;
    my $r = $data_out{$range}; # Starting ratio 

#    my ( $r_best, $s_best, $e_best, $refbest, $test_best ) = ( 1, $starts[$swin], $windows{$starts[$ewin]}, 0, 0 );

    my $s1 = $swin > 0 ? $starts[ $swin - 1 ] : $starts[0]; 
    my $e1 = $windows{ $starts[ $swin ] }; 
    my $s2 = $starts[ $ewin ]; 
    my $e2 = $ewin == $#starts ? $windows{ $starts[ $ewin ] } : $windows{ $starts[ $ewin + 1 ] };

#    warn "Range $s1 - $e1 with $s2 - $e2, ratio $r\n";
    next unless $r <= $cnv_lower_ratio or $r >= $cnv_higher_ratio; 
    warn "Fine-tuning $s1 - $e1 with $s2 - $e2 (was $r)\n";
    my $mid_point = int( $s1/2 +  $e2/2 );

    print "\rPreparing finetuning for test set, left side, pass 1   ";
    my %left_test = ();
    my $cnt = 0;
    foreach my $pos ( sort {$b<=>$a} ( grep { $s1 <= $_ and $_ <= $mid_point } @test ) ) {
        $left_test{$pos} = ++$cnt;
    } 

    print "\rPreparing finetuning for test set, left side, pass 2   ";
    my $curr = 0;
    foreach my $pos ( sort {$b<=>$a} ( grep { $s1 <= $_ and $_ <= $mid_point } ( @test, @ref ) ) ) {
        if ( exists($left_test{$pos} ) ) {
            $curr = $left_test{$pos}
        }
        else {
            $left_test{$pos} = $curr;
        }
    } 

    print "\rPreparing finetuning for test set, right side, pass 1   ";
    my %right_test = ();
    $cnt = 0;
    foreach my $pos ( sort {$a<=>$b} ( grep { $mid_point < $_ and $_ <= $e2 } @test ) ) {
        $right_test{$pos} = ++$cnt;
    } 

    print "\rPreparing finetuning for test set, right side, pass 2   ";
    $curr = 0;
    foreach my $pos ( sort {$a<=>$b} ( grep { $mid_point < $_ and $_ <= $e2 }  ( @test, @ref ) ) ) {
        if ( exists($right_test{$pos} ) ) {
            $curr = $right_test{$pos}
        }
        else {
            $right_test{$pos} = $curr;
        }
    } 

    print "\rPreparing finetuning for ref set, left side, pass 1   ";
    my %left_ref = ();
    $cnt = 0;
    foreach my $pos ( sort {$b<=>$a} ( grep { $s1 <= $_ and $_ <= $mid_point } @ref ) ) {
        $left_ref{$pos} = ++$cnt;
    } 

    print "\rPreparing finetuning for ref set, left side, pass 2   ";
    $curr = 0;
    foreach my $pos ( sort {$b<=>$a} ( grep { $s1 <= $_ and $_ <= $mid_point } ( @test, @ref ) ) ) {
        if ( exists($left_ref{$pos} ) ) {
            $curr = $left_ref{$pos};
        }
        else {
            $left_ref{$pos} = $curr;
        }
    } 

    print "\rPreparing finetuning for ref set, right side, pass 1   ";
    my %right_ref = ();
    $cnt = 0;
    foreach my $pos ( sort {$a<=>$b} ( grep { $mid_point < $_ and $_ <= $e2 } @ref ) ) {
        $right_ref{$pos} = ++$cnt;
    } 

    print "\rPreparing finetuning for ref set, right side, pass 2   ";
    $curr = 0;
    foreach my $pos ( sort {$a<=>$b} ( grep { $mid_point < $_ and $_ <= $e2 }  ( @test, @ref ) ) ) {
        if ( exists($right_ref{$pos} ) ) {
            $curr = $right_ref{$pos};
        }
        else {
            $right_ref{$pos} = $curr;
        }
    } 

    print "\rFinding best range                                     ";
    my $s3 = $windows{$starts[$swin]};
    my $e3 = $windows{$starts[$ewin]};
    my ( $best_score, $best_s, $best_e, $best_size, $best_r ) = ( 0, $s3, $e3, $e3-$s3+1, $r );

    my $test_count = $s3  > $mid_point ? ( $right_test{ $e3 } - $right_test{ $s3 } ) : 
                     $e3 <= $mid_point ? ( $left_test{  $s3 } - $left_test{  $e3 } ) :
                     ( $left_test{ $s3 } + $right_test{ $e3 } );
    my $ref_count  = $s3  > $mid_point ? ( $right_ref{ $e3 } - $right_ref{ $s3 } ) : 
                     $e3 <= $mid_point ? ( $left_ref{  $s3 } - $left_ref{  $e3 } ) :
                     ( $left_ref{ $s3 } + $right_ref{ $e3 } );
    $test_count = 1 unless $test_count;
    $ref_count  = 1 unless $ref_count;
    my $init_score = $r > 1 ? ( $test_count/10000 + $test_count / $ref_count ) : ( $ref_count/10000 + $ref_count / $test_count );

    my @s1s = grep { $s1 <= $_ and $_ <= $e1 } ( @test, @ref );
    my @s2s = grep { $s2 <= $_ and $_ <= $e2 } ( @test, @ref );
    foreach my $s ( @s1s ) {
        foreach my $e ( @s2s ) {
            next if $s >= $e;
            $test_count = $s  > $mid_point ? ( $right_test{ $e } - $right_test{ $s } ) : 
                          $e <= $mid_point ? ( $left_test{  $s } - $left_test{  $e } ) :
                          ( $left_test{ $s } + $right_test{ $e } );
            $ref_count  = $s  > $mid_point ? ( $right_ref{ $e } - $right_ref{ $s } ) : 
                          $e <= $mid_point ? ( $left_ref{  $s } - $left_ref{  $e } ) :
                          ( $left_ref{ $s } + $right_ref{ $e } );
            $test_count = 1 unless $test_count;
            $ref_count  = 1 unless $ref_count;
            # Avoid sub-window
            next if $r > 1 and $test_count < $nOfReadsInWindow;
            next if $r < 1 and $ref_count < $nOfReadsInWindow;
            my $score = $r > 1 ? ( $test_count/10000 + $test_count / $ref_count ) : ( $ref_count/10000 + $ref_count / $test_count );
            if ( !$best_score or $best_score < $score or ( $best_score == $score and $best_size < ( $e - $s +1 ) ) ) {
                ( $best_score, $best_s, $best_e, $best_size, $best_r ) = ( $score, $s, $e, $e - $s + 1, $test_count / $ref_count ) ;
            }
        } 
    }
    $best_r = sprintf ( "%.3f", $best_r/$med);
    $test_count = $best_s  > $mid_point ? ( $right_test{ $best_e } - $right_test{ $best_s } ) : 
                  $best_e <= $mid_point ? ( $left_test{  $best_s } - $left_test{  $best_e } ) :
                  ( $left_test{ $best_s } + $right_test{ $best_e } );
    $ref_count  = $best_s  > $mid_point ? ( $right_ref{ $best_e } - $right_ref{ $best_s } ) : 
                  $best_e <= $mid_point ? ( $left_ref{  $best_s } - $left_ref{  $best_e } ) :
                  ( $left_ref{ $best_s } + $right_ref{ $best_e } );
    my $test_count2 = scalar grep { $best_s <= $_ and $_ <= $best_e } @test;
    my $ref_count2  = scalar grep { $best_s <= $_ and $_ <= $best_e } @ref;
    print join( "\t", "\nTuned: ", $best_s, $best_e+$read_length-1, $best_r, $test_count, $ref_count ), "\n";
    print F join( "\t", $best_s, $best_e+$read_length-1, $best_r, $test_count, $ref_count ), "\n";
    
}
close F;


### SUBROUTINES
###


# Change Point candidate
sub find_cp {
    my $s1 = shift;
    my $e1 = shift; 
    my $size = $e1 - $s1 + 1; 
    return if $size < 3;
    my $m = 0; 
    map { $m += $_ } @full[$s1..$e1];
    $m /= $size; # m = average ratio

    # compute cumulative sums
    my @S = ( $full[$s1] - $m ) ;
    foreach my $counter ( $s1+1 .. $e1 ) {
        push @S, $S[-1] + $full[$counter] - $m;
    }

    map {$_ = abs($_)} @S;
    my $maximum = max(@S);
    my $Sd = $maximum - min(@S);

    my $point = 0; # max point
    while ( $S[$point] < $maximum ) { $point++ } 
    die "Did not reach max" if $S[$point] != $maximum;

    # Bootstrap
    my @B = ();
    foreach  my $strap ( 1 .. $bootstraps ) {
        my @arr = @full[$s1..$e1];
        foreach my $i (reverse ( 1 .. $#arr ) ) {
            my $j = int( rand( $i + 1 ) );
            @arr[$i,$j] = @arr[$j,$i];
        }
        # compute cumulative sums
        @S = ( $arr[0] - $m );
        for my $counter ( 1 .. $#arr ) {
            push @S, $S[-1] + $arr[$counter] - $m;
        }
        map {$_ = abs($_)} @S;
        push @B, (max(@S) - min(@S));
    }

    my $x = 0;
    foreach (@B) { if ($_ < $Sd) {$x++} }

    my $conf = $x/$bootstraps;
    return ($s1 + $point, $conf);
}
# end of Change Point 

sub median { 
    @_ == 1 or die ('Sub usage: $median = median(\@array);'); 
    my ($array_ref) = @_; 
    my $count = scalar @$array_ref; 
    # Sort a COPY of the array, leaving the original untouched 
    my @array = sort { $a <=> $b } @$array_ref; 
    if ( $count % 2 ) { # odd 
        return $array[int($count/2)]; 
    }
    else { # even = avg of 2 middle elements
        return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
    } 
} 

