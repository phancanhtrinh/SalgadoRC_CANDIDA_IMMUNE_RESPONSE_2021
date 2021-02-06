
#Install GOplot
install.packages("GOplot")
install.packages("RColorBrewer")
install.packages("viridis")
install.packages("colorspace")
install.packages("scales")
install.packages('unikn')
library('unikn') 

#Access the library
library(GOplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(colorspace)


################################################################
#GoChord
#choose one of the three colours below:

#n. of colours (length) should be the same as GOs.

#1 BREWER PAL:  
# BROWN BLUE https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/
col = brewer.pal (7, 'BrBG')
col = usecol(pal = pal_unikn, n = 25)

# BLUES (https://rdrr.io/cran/viridis/f/vignettes/intro-to-viridis.Rmd)
col = gradient_n_pal(brewer_pal(type="seq")(7))(seq(0, 1, length=7)) 

# YELLOW GREEN BLUE (https://rdrr.io/cran/viridis/f/vignettes/intro-to-viridis.Rmd)
col = gradient_n_pal(brewer_pal(type="seq", palette = 'YlGnBu')(7))(seq(0, 1, length=7)) 

#2 VIRIDIS:   (https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
col=rev(viridis(17)) 
#3 COLORSPACE (https://www.groundai.com/project/colorspace-a-toolbox-for-manipulating-and-assessing-colors-and-palettes/1)
col = sequential_hcl(17, palette = "Terrain 2")
col2 = sequential_hcl(17, palette = "Purple-Blue")


library(circlize)
install.packages("circlize")
# Create the GOChord 

GOChord(test.goplot, 
        space = 0.02, 
        gene.order = 'logFC', 
        gene.space = 0.25, 
        gene.size = 20,
        ribbon.col = col,
        lfc.col = c('black', 'black', 'black'),
        process.label = 20)


if ('logFC' > 2) { 
        lfc.col = c("coral2")
} else {lfc.col = c("red") }
############################################################
#GOheat
#to show the genes in the heatmap
if ('logFC' > 2) { 
        lfc.col = c("coral2")
} else {lfc.col = c("red") }

GOHeat_fix <- function (data, nlfc, fill.col) 
{
        x <- y <- z <- NULL
        if (missing(nlfc)) 
                nlfc <- 0
        else nlfc <- nlfc
        if (missing(fill.col)) 
                fill.col <- c("firebrick", "white", "dodgerblue")
        else fill.col <- fill.col
        distance <- dist(data)
        cluster <- hclust(distance)
        M <- dim(data)[2]
        nterm <- M - nlfc
        if (nlfc == 0) {
                s <- rowSums(data[, 1:nterm])
                tmp <- NULL
                for (r in 1:nrow(data)) {
                        tmp <- c(tmp, as.numeric(gsub(1, s[r], data[r, 1:nterm])))
                }
        }
        else {
                tmp <- NULL
                for (r in 1:nrow(data)) {
                        tmp <- c(tmp, as.numeric(gsub(1, data[r, (nterm + 
                                                                          1)], data[r, 1:nterm])))
                }
        }
        df <- data.frame(x = factor(rep(cluster$order, each = nterm)), y = rep(colnames(data[, 
                                                                                             1:nterm]), length(rownames(data))), z = tmp, lab = rep(rownames(data), 
                                                                                                                                                    each = nterm))
        df_o <- df[order(df$x), ]
        g <- ggplot() +
                geom_tile(data = df_o, aes(x = x, y = y, fill = z)) +
                scale_x_discrete(breaks = 1:length(unique(df_o$x)), labels = unique(df_o$lab)) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(), 
                      axis.text.y = element_text(size = 14),
                      panel.background = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                coord_fixed()
        if (nlfc == 0) {
                g + scale_fill_gradient2("Count", space = "Lab", low = fill.col[2], 
                                         mid = fill.col[3], high = fill.col[1])
        }
        else {
                g + scale_fill_gradient2("logFC", space = "Lab", low = fill.col[3], 
                                         mid = fill.col[2], high = fill.col[1])
        }
}

#plot the GOHeat
GOHeat_fix(chord, nlfc = 1, fill.col = c('blue', 'black', 'yellow'))

