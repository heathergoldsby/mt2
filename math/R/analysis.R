library(stringr)
library(ggplot2)
library(plyr)

find.files <- function(file_pattern=".*", dir_pattern=".*", path=".") {
	grep(dir_pattern, 
		list.files(recursive=TRUE, pattern=file_pattern, path=path, full.names=TRUE),
		value=TRUE)
}

filename.tags <- function(filename) {
	m = str_match(filename,"([-[:alnum:]]+)/([[:alnum:]]+)_([[:digit:]]+)/[[:alnum:]\\._]+$")
	return(list(filename=m[[1]],
		expr=as.factor(m[[2]]),
		treatment=as.factor(m[[3]]),
		trial = as.factor(m[[4]])))
}

load.files <- function(file_pattern, dir_pattern=".*", path="./", tag=TRUE, ...) {
	files = find.files(file_pattern, dir_pattern, path)
	data.frame(do.call(rbind,lapply(files, 
		function(x) {
			if(length(grep("\\.gz$",x))) {
				y=read.table(gzfile(x), header=TRUE, ...)
			} else {
				y=read.table(file(x), header=TRUE, ...)
			}
			if(tag) {
				t = filename.tags(x)
				y$filename = x
				y$expr = t$expr
				y$treatment = t$treatment
				y$trial = t$trial
			} else {
				return(y)
			}
			return(y)
		}
	)))
}

dominant <- function(D) {
	D[which.max(D$max_fitness),]
}

STYLE="draft" # or "final"
WIDTH=6
HEIGHT=3.75

quick_theme <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(colour = "black",size=0.75), legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(1,0), text=element_text(size=14))


# Show figure x.
#
showfig <- function(x) {
	quartz(width=WIDTH, height=HEIGHT)
	print(x)
}

# Save figure x in in PDF format.
#
savefig <- function(x, name, width=WIDTH, height=HEIGHT) {
	f = paste(figpath,name,".pdf",sep="")
	pdf(file=f, width=width, height=height, family="Helvetica")
	print(x)
	dev.off()
}
