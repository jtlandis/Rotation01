

box::use(
  ./gtf,
  rlang[abort]
)

expected_gencode_file <- "gencode.v44.annotation.gtf"
gencode_file <- box::file(expected_gencode_file)
if (!file.exists(gencode_file))
  abort(
    c(sprintf("file: `%s` was not found in path `%s`", expected_gencode_file, box::file()),
      "i" = "module `gr` requires a gtf file reference which likely hasnt been downloaded or is not named correctly.",
      "i" = "download from: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz'"))


if (!file.exists(rds_file <- box::file('genomicRanges.rds'))) {
  GR <- gtf$read_gtf2GR(gencode_file)
  saveRDS(GR, file = rds_file)
} else {
  GR <- readRDS(rds_file)
}

box::export(GR)


