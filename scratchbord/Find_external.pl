#! /usr/bin/perl -w

use warnings;
use strict;
use File::Spec::Functions 'catfile';

#setMethod('as_cellexalvrR',

my $path = $ARGV[0];
$path = './' unless ( -d $path );

my $translate = $ARGV[1];
my $replaceWith;

my $trivial = { map { $_ => 1 } 'function', 'c', 'nrow', 'ncol', 'DataFrame' ,'IN', 'in',
        'any', 'apply', 'attributes', 'cbind', 'rbind', 'vector' ,'ceiling',
        'class', 'basename', 'col' , 'colFun', 'colMeans', 'dim', 'dirname',
        'eval', 'exists', 'base', 'factor', 'cellexalvrr' ,'unique', 'exp', 'for', 'fun',
        'unlist', 'get','getwd', 'grep', 'gsub' , 'hc', 'hclustfun'};

unless ( -f $translate ) {
    warn
"You could give me a second file with the definitions of the external function - I would replace them for you.\n"
      . "on each line: funtion<TAB>package\n";
}
else {
    open( IN, "<$translate" ) or die $!;
    my @line;
    while (<IN>) {
        chomp();
        @line = split( /\s/, $_ );
        next unless ( $line[1] =~ m/\w/ );
	if ( $line[1] eq "base" ){
		#die "Added $line[0] to tzrivial\n";
		$trivial->{$line[0]} = 1;
	}else{
	        $replaceWith->{ $line[0] } = "$line[1]::$line[0]";
	}
    }
    close(IN);
}

opendir( DIR, $path );
my ( $in, $data, @tmp );

foreach my $file ( readdir(DIR) ) {
    next unless ( $file =~ m/.R$/ );
    open( $in, "<" . catfile( $path, $file ) )
      or die "I could not open the file $file\n$!\n";
    while (<$in>) {

        #setMethod('as_cellexalvrR',
        next if ( $_ =~ m/^\s*#/ );
        chomp();
        if ( $_ =~ m/setMethod\s*\('([\w_\.]+)',/ ) {
            $data->{$1} = 1;
	    #warn "found internal function $1\n";
        }elsif ( $_ =~ m/^\s*([\w_\.]*)\s*[<=]-?\s*function\s*\(/ ){
	    $data->{$1} = 1;
	    #warn "found internal function $1\n";
    	}

    }
    close($in);
}

closedir(DIR);

# Now I have all internal functions.

opendir( DIR, $path );
my ( @mod, @missing, $modified, $noStrings, $lineChanged );

foreach my $file ( readdir(DIR) ) {
    next unless ( $file =~ m/.R$/ );
    open( $in, "<" . catfile( $path, $file ) )
      or die "I could not open the file $file\n$!\n";
    $modified = 0;
    @mod = undef;
    while (<$in>) {
        my $new = $_;
        #setMethod('as_cellexalvrR',
        ## get all function calls not preficed with package::
	$lineChanged =0;
        unless ( $_ =~ m/^\s*#/ ) {
            
	    $noStrings = $_;
	    $noStrings =~ s/".*?"/---/g;
	    $noStrings =~ s/'.*?'/---/g;
            foreach my $function ( $noStrings =~ m/[,\s]([\w_\.]+)\s*\(/g ) {
                next if ( $data->{$function} );
                next if ( $trivial->{$function} );

                #		 warn "found external function $function\n";
                if ( defined $replaceWith->{$function} ) {
		    $lineChanged = 1;
                    $new =~ s/$function*\(/$replaceWith->{$function}(/;
                    $modified = 1;
                }
                else {
                    push( @missing, $function );
                }
            }
	    #print "line changed from/to\n\t$_\t$new\n" if( $lineChanged ) ;
        }
        push( @mod, $new );
	#print "I opushed this line: $new";
    }
    close($in);
    if ( $modified == 1 ) {
	    print "I have a modified file $file\n";
	open ( OUT, ">". catfile( $path, $file) ) or die "I could not open the file $file\n$!\n";
	print OUT join("", @mod);
	close ( OUT );
    }
}

closedir(DIR);

if ( @missing > 0 ) {
    my $functs;
    map { $functs->{$_} = 1 } @missing;
    warn "I lack knowledge about some external functions ("
      . scalar( keys %$functs )
      . " - please fill the file 'function_definitions.tab'\n";

    open( OUT, ">function_definitions.tab" ) or die $!;
    print OUT map { "$_\t\n" } sort keys %{$functs};
    warn "rdocumentation R function " . join( "\nrdocumentation R function ", sort keys %{$functs} );
    close(OUT);
}

print "Done\n"
