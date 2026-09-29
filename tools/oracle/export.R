# Export RData arrays without a decimal roundtrip. Source files remain untouched.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], mustWork = TRUE)
output <- args[2]
files <- list.files(root, pattern = "\\.RData$", recursive = TRUE, full.names = TRUE)
for (file in files) {
  environment <- new.env(parent = emptyenv())
  objects <- load(file, envir = environment)
  relative <- substring(file, nchar(root) + 2)
  for (name in objects) {
    value <- get(name, envir = environment)
    if (!is.atomic(value) || is.factor(value)) stop(paste("Unsupported object", file, name))
    folder <- file.path(output, relative, name)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    metadata <- list(type = typeof(value), length = length(value), dim = dim(value),
                     names = names(value), dimnames = dimnames(value), order = "F",
                     endian = "little")
    if (typeof(value) == "double") {
      writeBin(as.vector(value), file.path(folder, "data.bin"), size = 8, endian = "little")
    } else if (typeof(value) %in% c("integer", "logical")) {
      writeBin(as.integer(value), file.path(folder, "data.bin"), size = 4, endian = "little")
    } else if (typeof(value) == "character") {
      jsonlite::write_json(as.vector(value), file.path(folder, "values.json"), na = "null")
    } else stop(paste("Unsupported type", typeof(value)))
    # NA and NaN are distinct in R and need explicit masks beside binary values.
    writeBin(as.integer(is.na(value)), file.path(folder, "missing.bin"), size = 4, endian = "little")
    writeBin(as.integer(is.nan(value)), file.path(folder, "nan.bin"), size = 4, endian = "little")
    jsonlite::write_json(metadata, file.path(folder, "metadata.json"), auto_unbox = TRUE,
                         null = "null", pretty = TRUE)
  }
}
