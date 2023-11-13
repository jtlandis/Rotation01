
if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}

tryCatch(invisible(find.package("box")),
         error = function(cnd) install.packages("box"),
         finally = {function(){
           p <- Sys.getenv("R_BOX_PATH") |>
             strsplit(split = .Platform$path.sep) |>
             _[[1L]]
           #set box mods in path
           path <- paste(c(box::file(), p), collapse = .Platform$path.sep)
           Sys.setenv(R_BOX_PATH = path)
           }}())

