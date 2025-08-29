

box::use(
  ./mods/gtf[...],
  ./mods/sims[...]
)


GR <- read_gtf("gencode.v44.annotation.gtf") 
GR2 <- tidyr::unnest(GR[,'trans'], trans)
GR3 <- as_genomicRange_(GR2)
