source("~/.Rprofile")

install.packages(c("bibtex",
	"cluster",
	"data.table",
	"DBI",
	"devtools",
	"dplyr",
	"gdata",
	"ggplot2",
	"gplots",
	"igraph",
	"knitcitations",
	"knitr",
	"knitrBootstrap",
	"lattice",
	"locfit",
	"pamr",
	"pander",
	"plyr",
	"RColorBrewer",
	"Rcpp",
	"RcppArmadillo",
	"RCurl",
	"RefManageR",
	"RJSONIO",
	"RSQLite",
	"stringr",
	"survival",
	"XML",
	"xtable",
	"yaml"
	)
)

require(devtools)
install_github("rCharts", "ramnathv")