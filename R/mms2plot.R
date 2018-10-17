#' @author Lei Li
#' @title mms2plot
#' @description Visualization of multiple MS/MSs for (un)modified peptides
#' @param id_table_path File path name of a table that contains MS2 information
#'        of identified (un)modified peptides plus a group labelling. The format
#'        of MS2 information is referred to as the output file ms2.txt from
#'        the Maxquant search software.
#' @param par_xml_path Xml file path of parameters for modifications and
#'        labelling. The file format is referred to as modifications.xml in
#'        Maxquant.
#' @param mqpar_filepath File path name that includes a list of parameter files
#'        for search engines. The parameter file format is referred to as
#'        mqpar.xml in Maxquant.
#' @param min_intensity_ratio minimal percentage threshold of MS2 intensity,
#'        compared with the highest intensity. (default=0.01).
#' @param pdf_width The width of a single PSM figure area in inches. The area
#'        includes both plot region and outer margin area. The default is 3.35
#'        for a single column. The width is 7 for double-column.
#' @param pdf_height The height of a single PSM figure area in inches. The
#'        default is pdf_width/2.4.
#' @param xmai margin of the figure in number of inches for x axis.
#'        (default=pdf_width*0.15/3.35).
#' @param ymai margin of the figure in number of inches for y axis.
#'        (default=pdf_width*0.15/3.35).
#' @param ppm The threshold of mass error in parts per million(ppm):
#'        (exactMass-accurateMass)/exactMass*1E6. (default=20).
#' @param y_ion_col y ion color. The default is red.
#' @param b_ion_col b ion color. The default is blue.
#' @param peaks_col color of all M/Z peaks. The default is grey.
#' @param ymax The height of the plot region relative to the default.
#'        (default=1.6).
#' @param info_height The height of MS2 annotation. (default=1.5)
#' @param peptide_height The height of peptide sequence annotation in the plot
#'        region, relative to the default. (default=1.3)
#' @param mod_height The height of modification annotation relative to the
#'        location where peptide sequence is annotation. (default=0.07).
#' @param len_annoSpace The length of b/y ion annotation segments. (default=0.1).
#' @param lwd line width relative to the default. (default=pdf_width/3.35).
#' @param cex A numerical value giving the amount by which plotting text and
#'        symbols is magnified relative to the default. (default=pdf_width/3.35).
#' @param show_letterBY Logical: should "b"/"y" characters are shown on the peak
#'        annotation? The default is FALSE.
#'
#' @return NULL
#' @import xml2
#' @importFrom MSnbase readMSData
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines par plot segments text
#' @importFrom data.table fread
#'
#rm(list=ls())
#.libPaths( c( .libPaths(), "D:/Rpackages_tmp") )
#library(xml2)
#library(MSnbase)
#library(data.table)
#library(DescTools)  # MixColor

# load aa_mw and atom_mw files
#aa_mw_path = "data/AA_MW.txt"
#atom_mw_path = "data/atom_MW.txt"
#aa_mw_table <-   data.table::fread(aa_mw_path,   sep = "\t", check.names = FALSE, fill = TRUE, header = TRUE)
#atom_mw_table <- data.table::fread(atom_mw_path, sep = "\t", check.names = FALSE, fill = TRUE, header = TRUE)

#id_table_path = "D:/R_packages/MMS2plot/ext/msms_SILAC.txt"
#id_table_path = "D:/R_packages/MMS2plot/ext/msms_TMT.txt"
#id_table_path = "D:/R_packages/MMS2plot/ext/msms_dim.txt"
#par_xml_path = "D:/R_packages/MMS2plot/ext/modifications.xml"
#parfile = "ext/mqpar.txt"

#source("R/plot_mirror_or_group.R")
#source("R/add_mod_aa.R")
#source("R/psm_calculation.R")
#source("R/plot_components.R")
#' @export
mms2plot <-function(id_table_path, #="ext/msms_test.txt",
                    par_xml_path, #"ext/modifications.xml",
                    mqpar_filepath,  #
                    min_intensity_ratio=0.01, # mininum peak intensity percentage
                    pdf_width=3.35, # one column  7  # two column
                    pdf_height=pdf_width/2.4,
                    xmai = 0.15*pdf_width/3.35, #margin for x axis
                    ymai = 0.3*pdf_width/3.35, #margin for y axis
                    ppm=20,
                    y_ion_col="red",
                    b_ion_col="blue",
                    peaks_col = "grey", #color of all MS2 peaks
                    ymax = 1.6,
                    info_height = 1.5,   # height of MS2 information (e.g. Title, Gene name, Charge.)
                    peptide_height = 1.3, # y value for labelling PSMs(e.g. peptides)
                    mod_height = 0.07, # y value for labelling PSMs(e.g. peptides)
                    len_annoSpace = 0.1,   # the length of space for annotation (i.e. empty space between peaks and peakLabel or ion length for PSManno )
                    lwd=1*pdf_width/3.35,
                    cex=1*pdf_width/3.35,
                    show_letterBY=FALSE){
  srt = 0
  # read a batch of mqpar.xml files and extract modifications and label information stored in mqpar
  mqpar_files=readLines(mqpar_filepath, warn=FALSE)
  mqpar = data.table::rbindlist(lapply(mqpar_files, readMQPar))

  input_table <- data.table::fread(id_table_path, na.strings = "NA", sep = "\t",
                                   check.names = FALSE, fill = TRUE, header = TRUE, stringsAsFactors = FALSE)

  input_table$base_rawFile = basename(input_table$`Raw file`)
  #browser()
  input_table = check_input_table(input_table, id_table_path, mqpar)

  lapply(unique(input_table$`Raw file`), drawms2plot_samerawfile, input_table, par_xml_path, mqpar, min_intensity_ratio, pdf_width, pdf_height,
          xmai, ymai, ppm, y_ion_col, b_ion_col, peaks_col, ymax, peptide_height, info_height, mod_height, len_annoSpace, lwd, cex, show_letterBY, srt) # call for individual raw_files
  invisible(gc())
}

#mqpar_filepath = "mqpar_batch.txt"
#mm2plot(id_table_path=id_table_path, mqpar_filepath=mqpar_filepath)

