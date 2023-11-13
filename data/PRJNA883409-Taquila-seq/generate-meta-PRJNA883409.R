#! /usr/bin/Rscript

box::use(readr[read_csv,write_csv], dplyr[...])

data <- readr::read_csv("data/PRJNA883409-Taquila-seq/meta_PRJA.txt") |>
  filter(nanopore_library_type=="TEQUILA-seq") |>
  head(n=18)

data |>
  mutate(
    sample = c(
      rep_len("brain", 3),
      rep_len("neuroblastoma", 9),
      rep_len("breast-cancer-basal-A", 6)
    ),
    cell_line = Cell_Line,
    reps = c(
      rep_len(1:3, 12), rep_len(1:2, 6)
    ), 
    time_points = c(
      rep_len(NA_character_, 3),
      rep(c("4Hr", "8Hr", "48Hr"), each = 3),
      rep_len(NA_character_, 6)
    ),
    .before = Run
  ) |>
  select(-c(`Assay Type`, starts_with("DATASTORE"),Cell_Line:version)) |>
  write_csv(file = "data/PRJNA883409-Taquila-seq/meta_PRJA-filtered.csv")

# data <- read.csv(commandArgs(trailingOnly = T)[1])
# data <- read.csv("data/PRJNA883409-Taquila-seq/meta-data.csv")
# run_num_acc <- 678811L:678921L
# 
# template_url <- "https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR21%i&display=data-access"
# 
# urls <- sprintf(template_url, run_num_acc)
# 
# 
# rD <- RSelenium::rsDriver(port = netstat::free_port(), browser = "firefox", chromever = NULL)
# rD$server$log()
# clnt <- rD[["client"]]
# clnt$navigate(urls[2])
