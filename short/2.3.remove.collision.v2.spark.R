library(ISAnalytics)
library(dplyr)
library(reshape2)

ProjectID <- "SPARK"
proj_folder <- "/Users/cipriani.carlo/Dropbox (HSR Global)/carlo/AAV-Short/SPARK/FINAL"
report_folder <- fs::path(proj_folder, "reports")
if (!fs::dir_exists(report_folder)) {
  fs::dir_create(report_folder)
}
af_path <- "/Users/cipriani.carlo/Dropbox (HSR Global)/carlo/AAV-Short/SPARK/AF.allpools.SPARK.modifAF1.tsv"
## AF IMPORT
af <- import_association_file(path = af_path, root = NULL, separator = )
af$ReplicateNumber<-1
samples_to_filter <- c()
new_mand_vars <- tibble::tribble(
  ~names, ~types, ~transform, ~flag, ~tag,
  "chr", "char", ~stringr::str_replace_all(.x, "chr", ""), "required",
  "chromosome",
  "integration_locus", "int", NULL, "required", "locus",
  "strand", "char", NULL, "required", "is_strand",
  "gap", "int", NULL, "required", NA_character_,
  "junction", "int", NULL, "required", NA_character_,
  "aav_alignments_start_end", "int", NULL, "required", NA_character_
)
set_mandatory_IS_vars(new_mand_vars)
#isa_dataframe <- import_single_Vispa2Matrix(path_to_file = "/home/andrea/ISA/AAV-Taylor/ShsCount_Taylor_gap_junction_new_groupd_diff.no0.annotated.tsv")
isa_dataframe <- read.csv(file = "/Users/cipriani.carlo/Dropbox (HSR Global)/carlo/AAV-Short/SPARK/FINAL/SeqCount_from_ShsCount_SPARK.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1.tsv", header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c('', "NA"))
isa_dataframe_no1 <- cbind(isa_dataframe, 
                           data.frame("SUM" = apply(isa_dataframe[setdiff(colnames(isa_dataframe), c(mandatory_IS_vars(), annotation_IS_vars()))], 1, function(x) {sum(x, na.rm = T)})))
isa_dataframe_no1 <- isa_dataframe_no1[which(isa_dataframe_no1$SUM > 1), setdiff(colnames(isa_dataframe_no1), c("SUM", samples_to_filter))]
isa_dataframe_no1_molten <- melt(data = isa_dataframe_no1, id.vars = c(mandatory_IS_vars(), annotation_IS_vars()), variable.name = "CompleteAmplificationID", na.rm = T, value.name = "seqCount")


# isa_dataframe_molten <- melt(data = isa_dataframe, id.vars = c(mandatory_IS_vars(), annotation_IS_vars(), "gap", "junction"), variable.name = "CompleteAmplificationID", na.rm = T, value.name = "SSC")
# isa_dataframe_molten_no1 <- isa_dataframe_molten %>%
#   # distinct(SubjectID, CellMarker, Tissue, TimePoint, PCRMethod) %>%
#   group_by(across(c(mandatory_IS_vars(), annotation_IS_vars(), "gap", "junction"))) %>%
#   summarise(sum = sum(SSC), count = n(), .groups = "drop") %>%
#   filter(sum > 1)

## COLLISIONS
coll <- remove_collisions(tibble::as_tibble(isa_dataframe_no1_molten), af, 
                          quant_cols = c(seqCount="seqCount"), report_path = report_folder)
# coll_sparse <- dcast(data = coll, formula = )
coll_sparse <- coll %>% tidyr::pivot_wider(id_cols = c(mandatory_IS_vars(), annotation_IS_vars()), names_from = "CompleteAmplificationID", values_from = "seqCount")
# coll_sparse <- as.data.frame(as_sparse_matrix(x = coll, seqCount = "SSC"))

coll_sparse_no1 <- cbind(coll_sparse, 
                         data.frame("SUM" = apply(coll_sparse[setdiff(colnames(coll_sparse), c(mandatory_IS_vars(), annotation_IS_vars()))], 1, function(x) {sum(x, na.rm = T)})))
# coll_sparse_no1_test <- coll_sparse_no1[which(coll_sparse_no1$SUM > 1), ]
coll_sparse_no1 <- coll_sparse_no1[which(coll_sparse_no1$SUM > 1), setdiff(colnames(coll_sparse_no1), c("SUM", samples_to_filter))]


write.table(x = coll_sparse_no1, file = "/Users/cipriani.carlo/Dropbox (HSR Global)/carlo/AAV-Short/SPARK/FINAL/SeqCount_from_ShsCount_SPARK.minfvalue30.gapabs120.aggregAF1.window8_20.cleaned1.no0.annotated.removed1.nocoll.tsv", sep = "\t", col.names = T, row.names = F, quote = F, na = '')
