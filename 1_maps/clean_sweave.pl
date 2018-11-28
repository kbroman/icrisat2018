#!/usr/bin/perl
######################################################################
# clean_sweave.pl
#
# Clean up the results produced by Sweave, so that the R input is
# exactly as it appears in the .Rnw file.
#
# (Sweave reformats the R input and sometimes the result runs into
# the margins; here we read the .Rnw file and find all of the code
# chunks, then read the .tex file and replace the input chunks with
# the code as it appeared in the .Rnw file.)
######################################################################

(@ARGV > 0) or die("Give file stem.");
$rnwfile = "$ARGV[0]\.Rnw";
$texfile = "$ARGV[0]\.tex";
$bakfile = "$ARGV[0]\.tex.bak";

$flag = 0;
if(-e $rnwfile and -e $texfile) {
    open(IN, $rnwfile) or die("Cannot read from $rnwfile.");
    while($line = <IN>) {
    chomp($line);
    if($line =~ /begin\{Sinput\}/) {
        while($line = <IN>) {
        if($line =~ /end\{Sinput\}/) {
            last;
        }
        if(substr($line, 0, 1) eq " ") {
            $thecode[@thecode - 1] .= "+ " . $line;
        }
        else {
            push(@thecode, "> " . $line);
        }
        }
    }
    if($line =~ /^<<(.*)>>/) {
        @v = split(/,/, $1);

        if($1 =~ /echo=FALSE/) {
        $dontprint = 1;
        }
        else {
        $dontprint = 0;
        }

        if($v[0] =~ /=/) {
        $dontsave = 1;
        }
        else {
        $dontsave = 0;
        $codename = $v[0];
        }

        while($line = <IN>) {
        if($line =~ /^@/) {
            last;
        }
        unless($dontsave) {
            # <<something>> to be replaced with previous code
            if($line =~ /^<<(.+)>>/) {
            unless($dontprint) {
                push(@thecode, @{$code{$1}});
            }
            }
            else {
            $temp = substr($line, 0, 1);
            if($temp eq " " or $temp eq "}") {
                $line = "+ " . $line;
                $n = @{$code{$codename}};
                $code{$codename}[$n-1] .= $line;
                unless($dontprint) {
                $flag = 1;
                $thecode[@thecode-1] = $code{$codename}[$n-1];
                }
            }
            else {
                push(@{$code{$codename}}, "> " . $line);
                unless($dontprint) {
                $flag = 1;
                push(@thecode, "> " . $line);
                }
            }
            }
        }
        }
    }
    }
    close(IN);
}
else {
    print(" --No files to clean.\n");
}


if($flag) {

    if(-e $bakfile) {
    system("\\rm $bakfile");
    }
    system("\\mv $texfile $bakfile");

    open(IN, $bakfile) or die("Cannot read from $bakfile");
    open(OUT, ">$texfile") or die("Cannot write to $texfile");
    $n = 0;
    while($line = <IN>) {
    if($line =~ /begin\{Sinput\}/) {
        $line = "\\vspace*{-6pt}\n" . $line;
    }
    if($line =~ /begin\{Soutput\}/) {
        $line = "\\vspace*{-6pt}\n" . $line;
    }
    if($line =~ /end\{Soutput\}/) {
        $line = $line . "\\vspace*{-6pt}\n";
    }
    print OUT ("$line");
    if($line =~ /begin\{Sinput\}/) {
        while($line = <IN>) {
        if($line =~ /end\{Sinput\}/) {
            $line = $line . "\\vspace*{-6pt}\n";
            print OUT ("$line");
            last;
        }
        if($line =~ /^> /) {

            $thecode[$n] =~ s/^> > /> /;
            print OUT $thecode[$n];
            $n++;
        }
        }
    }
    }
    close(IN);
}
else {
    print(" --No need for cleaning.\n");
}
