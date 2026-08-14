# Integration checks for the Shiny app's non-interactive helpers and exports.
# Run from the repository root:
#   Rscript R/shiny/tests/test_dada2_parameter_selection_app.R

app_dir <- normalizePath(file.path("R", "shiny"), mustWork = TRUE)
previous_directory <- getwd()
on.exit(setwd(previous_directory), add = TRUE)
setwd(app_dir)

app_environment <- new.env(parent = globalenv())
sys.source("dada2_parameter_selection_app.R", envir = app_environment)

tests_run <- 0L
tests_failed <- 0L
check <- function(condition, message) {
  tests_run <<- tests_run + 1L
  if (isTRUE(condition)) {
    cat(sprintf("  [ok]   %s\n", message))
  } else {
    tests_failed <<- tests_failed + 1L
    cat(sprintf("  [FAIL] %s\n", message))
  }
}

cat("== App integration ==\n")
check(!exists("DADA2_DEFAULT_TRUNCQ", envir = app_environment, inherits = FALSE) &&
        identical(app_environment$RETENTION_DEFAULT_TRUNCQ, 2),
      "the engine owns the single truncQ default")

selection <- app_environment$select_validation_samples(
  c(missing = NA_real_, low = 20, high = 80), 2L
)
check(setequal(selection$name, c("low", "high")),
      "missing retention does not displace finite validation candidates")

workbook_path <- tempfile("app-report-", fileext = ".xlsx")
table_data <- data.frame(Parameter = c("count", "label"),
                         Value = c("12", "text"), stringsAsFactors = FALSE)
app_environment$add_sheet_to_excel(workbook_path, "Info", table_data, overwrite = TRUE)
app_environment$fix_excel_numeric_typed_cells(
  workbook_path, "Info", table_data, c(count = 12)
)
round_trip_numeric <- openxlsx::read.xlsx(workbook_path, rows = 2, cols = 2,
                                          colNames = FALSE)
check(is.numeric(round_trip_numeric[[1]]) && identical(round_trip_numeric[[1]][1], 12),
      "numeric report cells survive the atomic workbook update")
unlink(workbook_path)

app_source <- paste(readLines("dada2_parameter_selection_app.R", warn = FALSE), collapse = "\n")
check(!grepl("msg.text;", app_source, fixed = TRUE) &&
        grepl("createTextNode(String(msg.text))", app_source, fixed = TRUE),
      "console messages are inserted as text rather than executable HTML")

server_started <- tryCatch({
  shiny::testServer(app_environment$server, {
    session$setInputs(main_nav = "visualizer")
  })
  TRUE
}, error = function(e) {
  message("Server startup error: ", conditionMessage(e))
  FALSE
})
check(server_started, "the Shiny server initializes with its background tasks")

cat(sprintf("\n%d checks run, %d failed.\n", tests_run, tests_failed))
if (tests_failed > 0L) quit(status = 1L)
cat("ALL APP INTEGRATION TESTS PASSED\n")
