library(extrafont)
#Only once!
# font_import()
loadfonts()

###  2 fonctions
figPath <- ""
savePdf <- function(name, width = 8, height = 8){
  pdf(paste0(figPath, name, "_tmp.pdf"), width =  width,
      height = height, family = "CM Roman")
}
embedPdf <- function(name){
  embed_fonts(paste0(figPath, name, "_tmp.pdf"),
              outfile = paste0(figPath, name, ".pdf"))
  file.remove(paste0(figPath, name, "_tmp.pdf"))
}


# Procedure



