## wrangle.R - generates a summary table that ranks dark kinases on a priority scale
##
## by Artem Sokolov, et al.

library( dplyr )
library( synapseClient )

synapseLogin()

## Resolves a filename by downloading the file if it's a synapse ID
## Returns a filename that can be directly used for loading by, e.g., read.delim
resolve.filename <- function( fn )
{
    if( substr( fn, 0, 3 ) == "syn" )
    {
        s <- synGet( fn, downloadLocation = "." )
        return( s@filePath )
    }
    return( fn )
}

## Computes the overall score
totalScore <- function( v )
{
    val3 <- as.numeric( v["SmMolLig_Spec_Norm"] )

    s1 <- mean(v[c("RecombK","Struct3D","Antibody","ThermoFisher")] == "Y")
    s2 <- ((v["Disease"] == "Y") + as.numeric(v["AveGenomAberr_Norm"])) / 2
    s3 <- ((v["SmMolLig"] == "Y") + ifelse( is.na(val3), 0, val3 )) / 2
    s4 <- ifelse( is.na( v["BiolActivity_Norm"] ), 0, as.numeric( v["BiolActivity_Norm"] ) )

    s1 + 2*s2 + s3 + s4 + 5*as.integer( v["inLSPdatasets"] == "Y" )
}

## Miscellaneous functions
binarize <- function(v) { ifelse( v == "", "", "Y" ) }
wipe.NA <- function(v) { v[is.na(v)] <- ""; v }
valnorm <- function( v )
{ (v - min( v, na.rm=TRUE )) / diff( range( v, na.rm=TRUE ) ) }

## ################## ##
## Real-valued fields ##
## ################## ##

## Specificity of Small Molecule Ligand
X1 <- read.csv( resolve.filename("syn8377981"), check.names=FALSE, as.is=TRUE ) %>%
    select( Kinase = symbol, SmMolLig_Spec = 2 )

## Genomic Aberration Score
X2 <- read.csv( resolve.filename("syn8377980"), check.names=FALSE, as.is=TRUE ) %>%
    select( Kinase, AveGenomAberr = 3 )

## Biological Activity
X3 <- read.csv( resolve.filename("syn8377984"), check.names=FALSE, as.is=TRUE, header=FALSE ) %>%
    rename( Kinase = V1, BiolActivity = V2 )

## ############# ##
## Binary Fields ##
## ############# ##

## Expression in HMS Celllines
X4 <- scan( resolve.filename("syn8377975"), what=character() ) %>%
    data.frame( Kinase = ., inLSPdatasets = "Y", stringsAsFactors = FALSE )

## 3D Structure, Antibody, Recombinant Kinase, Known Link to Disease
X5 <- read.csv( resolve.filename("syn8377978"), check.names=FALSE, as.is=TRUE ) %>%
    select( Kinase = 2, RecombK = 4, Disease, Struct3D = 18, Antibody = 22 ) %>%
    mutate_each( funs( binarize ), -Kinase ) %>%
    group_by( Kinase ) %>%
    summarize_each( funs( unique ) ) %>%
    as.data.frame()

## ThermoFisher knock-out cell lines
X6 <- read.csv( resolve.filename("syn8377983"), check.names=FALSE, as.is=TRUE ) %>%
    select( Kinase = 3 ) %>%
    mutate( Kinase = gsub( ".$", "", Kinase ), ThermoFisher = "Y" )

## Presence of Small Molecular Ligand
X7 <- read.csv( resolve.filename("syn8377981"), check.names=FALSE, as.is=TRUE ) %>%
    select( Kinase = symbol ) %>% mutate( SmMolLig = "Y" )

## ############ ##
## Master table ##
## ############ ##

XX <- full_join( X1, X2 ) %>% full_join( X3 ) %>% full_join( X4 ) %>%
    full_join( X5 ) %>% full_join( X6 ) %>% full_join( X7 ) %>%
    mutate_each( funs( wipe.NA ), -(1:4) ) %>%
    mutate_each( funs( Norm = valnorm ), 2:4 ) %>%
    mutate( OverallScore = apply( ., 1, totalScore ) ) %>%
    slice( order( OverallScore, decreasing=TRUE ) )

write.table( XX, file="summary.tsv", quote=FALSE, sep="\t", row.names=FALSE )
