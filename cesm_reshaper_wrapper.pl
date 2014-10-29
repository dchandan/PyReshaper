#!/bin/env perl

#require XML::Lite;
require XML::LibXML;

# KMP: The following to be replaced when Alice's XML parser is finished...
#
#require ConfigCase;
#
#my $model = $ARGV[0];
#
#my $time_series_env = ConfigCase->new("config_definition.xml","config_tseries.xml");
#
#Set the variables
#my $pio_type = $time_series_env->getresolved('FILE_TYPE');
#my $i_directory = $time_series_env->getresolved('TSERIES_INPUT_DIRECTORY');
#my $o_directory = $time_series_env->getresolved('TSERIES_OUTPUT_DIRECTORY');
#my $in_format = $time_series_env->getresolved('TSERIES_IN_TIME_FORMAT');
#my $out_format = $time_series_env->getresolved('TSERIES_OUT_TIME_FORMAT');

my $pio_type = "netcdf4";
my $i_directory = "/foo/bar/hist";
my $o_directory = "/foo/bar/tseries";
my $in_format = "YYYY-MM-DD";
my $out_format = "YYYYMMDD";
my $start_time = "0001-01-01";
my $end_time = "0010-12-01";
my $root_name = "my.run";
my $i_time_period = "mon";
my $o_time_option = "nmonths"; #nsteps,nstep,nseconds,nsecond,nminutes,nminute,nhours,nhour,ndays,nday,nmonths,nmonth,nyears,nyear
my $o_time_n = "2";
my $compset = "B1850C5";
my $grid = "ne30_g16";
my $tseries_exclusions = "foo";
my @meta_vars = ("time", "lat", "lon", "time_bounds");


# Create a new xml document
my $doc = XML::LibXML::Document->new( "1.0", "UTF-8" );
my $root = $doc->createElement( "pyreshaper" );
$doc->setDocumentElement( $root );

my $timestamp = time;
my $config = $doc->createElement( "config" );
$root->addChild( $config );
$config->addChild( $doc->createAttribute( id => $timestamp ) );

my $node = $doc->createElement( "file_type" );
$node->addChild( $doc->createTextNode( $pio_type ) );
$config->addChild( $node );

my $node = $doc->createElement( "input_type" );
$node->addChild( $doc->createTextNode( "timeslice" ) );
$config->addChild( $node );

my $node = $doc->createElement( "output_type" );
$node->addChild( $doc->createTextNode( "timeseries" ) );
$config->addChild( $node );

my $node = $doc->createElement( "input_directory" );
$node->addChild( $doc->createTextNode( $i_directory ) );
$config->addChild( $node );

my $node = $doc->createElement( "output_directory" );
$node->addChild( $doc->createTextNode( $o_directory ) );
$config->addChild( $node );

my $node = $doc->createElement( "input_time_format" );
$node->addChild( $doc->createTextNode( $in_format ) );
$config->addChild( $node );

my $node = $doc->createElement( "output_time_format" );
$node->addChild( $doc->createTextNode( $out_format ) );
$config->addChild( $node );

my $node = $doc->createElement( "start_time" );
$node->addChild( $doc->createTextNode( $start_time ) );
$config->addChild( $node );

my $node = $doc->createElement( "end_time" );
$node->addChild( $doc->createTextNode( $end_time ) );
$config->addChild( $node );

my $node = $doc->createElement( "root_name" );
$node->addChild( $doc->createTextNode( $root_name ) );
$config->addChild( $node );

my $node = $doc->createElement( "input_time_period" );
$node->addChild( $doc->createTextNode( $i_time_period ) );
$config->addChild( $node );

my $node = $doc->createElement( "output_time_option" );
$node->addChild( $doc->createTextNode( $o_time_option ) );
$config->addChild( $node );

my $node = $doc->createElement( "output_time_n" );
$node->addChild( $doc->createTextNode( $o_time_n ) );
$config->addChild( $node );

my $node = $doc->createElement( "compset" );
$node->addChild( $doc->createTextNode( $compset ) );
$config->addChild( $node );

my $node = $doc->createElement( "grid" );
$node->addChild( $doc->createTextNode( $grid ) );
$config->addChild( $node );

my $node = $doc->createElement( "tseries_exclusions" );
$node->addChild( $doc->createTextNode( $tseries_exclusions ) );
$config->addChild( $node );

for (my $i = 0; $i < @meta_vars; $i++) {
  my $node = $doc->createElement( "meta_vars" );
  $node->addChild( $doc->createTextNode( $meta_vars[$i] ) );
  $config->addChild( $node ); 
}

$doc->toFile("pyreshaper_config2.xml",1);

#end

