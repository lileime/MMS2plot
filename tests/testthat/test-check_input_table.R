context("check input table")

test_that("check input table and path", {
    setwd("../../")
    mqpar_filepath = "inst/extdata/mqpar_batch.txt"
    par_xml_path = "inst/extdata/modifications.xml"
    id_table_path = "inst/extdata/silac/msms_SILAC.txt"
    mqpar_files<-data.table::fread(mqpar_filepath, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)
    mqpar_ppm <- data.table::rbindlist(apply(mqpar_files, 1, readMQPar_ppm))
    input_table <- data.table::fread(id_table_path, na.strings = "NA",
        sep = "\t", fill = TRUE, header = TRUE)
    # test for check_input_table
    input_table_check <- check_input_table(input_table, id_table_path, mqpar_ppm)
    input_table$`Modified sequence` <-
        gsub("_","", input_table$`Modified sequence`)
    expect_equal(input_table, input_table_check)
    # test for add_mod_aa

    MS2FileName<-unique(input_table$`Raw file`)[1]
    list_aaMwModTable_ppm_ama<-add_mod_aa(par_xml_path, basename(MS2FileName),
        mms2plot::aa_mw_table, mqpar_ppm)
    list_aaMwModTable_ppm=readRDS(file = "tests/testthat/list_aaMwModTable_ppm.rds")
    expect_equal(list_aaMwModTable_ppm, list_aaMwModTable_ppm_ama)
})
