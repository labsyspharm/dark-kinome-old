## plot.R - plots the results in the summary.tsv table
##
## by Artem Sokolov, et al.

library( dplyr )
library( reshape2 )
library( ggplot2 )

## Converts character { "", "Y" } to numeric {NA, 1}
yesno2int <- function(v) { ifelse( v == "Y", 1, NA ) }

## Maps Source to Panel
source2panel <- function(v)
{
    c( inLSPdatasets = "Expr", RecombK = "Info", Disease = "Disease", Struct3D = "Info",
      Antibody = "Info", ThermoFisher = "Info", SmMolLig = "Info", AveGenomAberr = "Disease",
      BiolActivity = "Expr", OverallScore = "Total" )[v]
}

## Scales the values to be between 0 and 1
valnorm <- function( v )
{ (v - min( v, na.rm=TRUE )) / diff( range( v, na.rm=TRUE ) ) }

## Load the summary table
X <- read.delim( "summary.tsv", as.is=TRUE ) %>%
    select( -(2:4), -12 ) %>%
    mutate_each( funs( yesno2int ), 2:8 ) %>%
    mutate( OverallScore = valnorm(OverallScore),
           Kinase = factor(Kinase, levels=Kinase) ) %>%
    rename( AveGenomAberr = AveGenomAberr_Norm,
           BiolActivity = BiolActivity_Norm ) %>%
    melt( "Kinase", variable.name = "Source", value.name = "Value" ) %>%
    na.omit() %>%
    mutate( Panel = source2panel(Source),
           Type = ifelse( Source %in% c("AveGenomAberr","BiolActivity","OverallScore"), "RV", "Binary" ) )

    

## Make the summary plot
g <- ggplot( X, aes( x = Kinase, y = Source ) ) +
    geom_point( data = subset(X, Type=="Binary"), color = "tomato" ) +
    geom_tile( data = subset(X, Type=="RV"), aes( fill = Value ) ) +
    facet_grid( Panel ~ ., scales = "free", space = "free" ) +
    scale_x_discrete( limits = levels(X$Kinase) ) +
    scale_fill_gradient2( high = "tomato", low = "steelblue", midpoint = 0.5, mid = "white" ) +
    theme( axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(angle=45, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank() )

ggsave( "summary.pdf", g, width=20, height=3 )

