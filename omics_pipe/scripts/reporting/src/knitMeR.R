### created by Tobias Meissner
###
### get arguments from the command line
### read in configs from the omics_pipe parameters file
### kick of the report using knitr bootstrap
library(knitr)
library(knitrBootstrap)
library(yaml)
library(reshape)

inFile <- commandArgs(trailingOnly=TRUE)[1]
inConfig <- commandArgs(trailingOnly=TRUE)[2]
patientID <- commandArgs(trailingOnly=TRUE)[3]
rootDir <- commandArgs(trailingOnly=TRUE)[4]

config <- yaml.load_file(inConfig)
outFile <- paste(config$REPORT_RESULTS, '/', patientID, '/', patientID, '.html', sep='')

if (!file.exists(rootDir)){
  dir.create(rootDir)
} 

opts_knit$set(root.dir = rootDir) 

# modify knitr bootstrap to write custom css to .html header...
my.knit_bootstrap_md <- function(input, output = NULL, boot_style=NULL, code_style=NULL, chooser=NULL,
                                 text = NULL, thumbsize=3, show_code=FALSE,
                                 show_output=TRUE, show_figure=TRUE,
                                 markdown_options=c('mathjax', 'base64_images', 'use_xhtml'),
                                 graphics = getOption("menu.graphics"), ...) {
  
  header = create_header(boot_style=boot_style, code_style=code_style,
                         chooser=chooser,
                         thumbsize=thumbsize,
                         show_code=show_code, show_output=show_output,
                         show_figure=show_figure, graphics=graphics)
  
  if(is.null(output))
    output <- sub_ext(input, 'html')
  
  if (is.null(text)) {
    markdown::markdownToHTML(input, header=header, stylesheet='custom.css',
                             options=markdown_options, output = output, ...)
  }
  else {
    markdown::markdownToHTML(text = input, header=header, stylesheet='custom.css',
                             options=markdown_options, output = output, ...)
  }
  invisible(output)
}

my.knit_boostrap <- function(input, output = NULL, boot_style=NULL, code_style=NULL, chooser=NULL,
                             thumbsize=3, show_code=FALSE, show_output=TRUE, show_figure=TRUE,
                             markdown_options=c('mathjax', 'base64_images', 'use_xhtml'),
                             ..., envir = parent.frame(), text = NULL,
                             quiet = FALSE, encoding = getOption('encoding'),
                             graphics = getOption("menu.graphics")) {
  
  knitr::render_html()
  opts_chunk$set(tidy=FALSE, highlight=FALSE)
  md_file =
    knit(input, NULL, text = text, envir = envir,
         encoding = encoding, quiet = quiet)
  
  my.knit_bootstrap_md(md_file, output, boot_style=boot_style,
                    code_style=code_style, chooser=chooser,
                    markdown_options = markdown_options,
                    thumbsize=thumbsize,
                    show_code=show_code, show_output=show_output,
                    show_figure=show_figure, ..., graphics=graphics)
  invisible(output)
}

# create the report...
my.knit_boostrap(inFile, output=outFile, chooser=c('boot','code'))
