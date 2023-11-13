

box::use(
  ./gtf
)

expected_gencode_file <- "gencode.v44.annotation.gtf"
gencode_file <- box::file(expected_gencode_file)
if (!file.exists(gencode_file))
  stop(paste0('file: `',expected_gencode_file,'` was not found in path `', box::file(), '`'), call. = F)

if (!file.exists(rds_file <- box::file('genomicRanges.rds'))) {
  GR <- gtf$read_gtf2GR(gencode_file)
  saveRDS(GR, file = rds_file)
} else {
  GR <- readRDS(rds_file)
}

box::export(GR)


