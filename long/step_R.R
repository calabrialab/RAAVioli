packages <- c("optparse", "tools","sqldf")

install.packages(setdiff(packages, rownames(installed.packages())))
if("GenomicAlignments" %in% rownames(installed.packages()) == FALSE)
{
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("GenomicAlignments",update=FALSE, ask = FALSE)
}

library(optparse)
library(tools)
option_list = list(
    make_option(c("-o", "--outdir"), type="character", default=NULL, 
                help="output directory. Results will be in the subdir results/", metavar="character"),
    make_option(c("-c", "--pipedir"), type="character", default=NULL, 
                help="output directory. Results will be in the subdir results/", metavar="character"),
    make_option(c("-i", "--input_labels"), type="character", default=NULL, 
                help="The file .tsv with the paths to fastq files as last column", metavar="character"),
    make_option(c("-s", "--summary"), type="character", default=NULL, 
                help="The file _Summary.*.sorted.annotated.strandness.tsv", metavar="character"),
    make_option(c("-m", "--min_cigar_alm_width"), type="integer", default=NULL, 
                help="The file _Summary.*.sorted.annotated.strandness.tsv", metavar="integer"),
    make_option(c("-n", "--min_alm_size"), type="integer", default=NULL, 
                help="The file _Summary.*.sorted.annotated.strandness.tsv", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd(opt$outdir)
source(paste0(opt$pipedir,"/isa_utils_functions.R"))

#library(openxlsx)
sample_labels <- opt$input_labels
sample_df <- read.csv(file = sample_labels, header=TRUE, fill=T, sep='\t', 
                      check.names = FALSE)

# maybe read this from file?
sample_label_to_esclude <- c("Test", "Negative")
summary_tsv <- opt$summary
basename_summary <- file_path_sans_ext(basename(summary_tsv))

mixed_aav_hg19 <- read.csv(file = summary_tsv,
                           header=TRUE, fill=T, check.names = FALSE, sep = '\t')

cols_to_include <- c("start", "end", "gene_chr", "gene_annotationsource", 
                     "gene_elementtype", "gene_start", "gene_end", "gene_strand",
                     "gene_details", "distance_to_gene", "sample", "sourcefile", 
                     "score")


reads_all_alm_by_cigar_chimera <- parseCigarForChimera(df_alm = mixed_aav_hg19,
                                                       quiet = T, cols_to_include = cols_to_include, 
                                                       min_cigar_alm_width = opt$min_cigar_alm_width)


reads_all_alm_by_cigar_chimera_ext <- merge(x = reads_all_alm_by_cigar_chimera,
                                            y = sample_df, by = c("sample"), 
                                            all.x = T)
dest_dir <-"final_results/"
cigar_parsed_out <- paste0(dest_dir,basename_summary,"CIGAR_parsed.tsv")
write.table(x = reads_all_alm_by_cigar_chimera_ext, file = cigar_parsed_out, sep = "\t", col.names = T, row.names = F, na = '')

reads_all_alm_by_cigar_chimera_details <- sqldf("select name, cigar, chr, start,
                                        end, strand, score, integration_locus,
                                        vector_junction_locus, sample, 
                                        sourcefile, sample_label, sample_name,
                                        read_len, min(read_alm_start) as read_alm_start,
                                        max(read_alm_end) as read_alm_end,
                                        min(query_start) as query_start,
                                        max(query_end) as query_end, 
                                        gene_chr, gene_start, gene_end, 
                                        gene_strand, gene_details, 
                                        distance_to_gene, chimera, junction,
                                        target_genome_position_wrt_junction, 
                                        read_alm_type, target_genome_multimapping
                                        from reads_all_alm_by_cigar_chimera_ext 
                                        where 1 
                                        group by name, cigar")


cigar_parsed_grouped_out <- paste0(dest_dir,basename_summary, ".CIGAR_parsed.grouped.tsv")
write.table(x = reads_all_alm_by_cigar_chimera_details, file = cigar_parsed_grouped_out, 
            sep = "\t", col.names = T, row.names = F, na = '')
print(paste(cigar_parsed_grouped_out, "done."))

reads_all_alm_by_cigar_chimera_details_geneid <- as.data.frame(
    t(as.data.frame(lapply(
        strsplit(as.character(reads_all_alm_by_cigar_chimera_details$gene_details), 
                 ';', fixed = T), function(x) {c(x)}))) )

reads_all_alm_by_cigar_chimera_details_geneid_2 <- as.data.frame(t(as.data.frame(
    lapply(
        strsplit(as.character(reads_all_alm_by_cigar_chimera_details_geneid$V1),' ', fixed = T),
        function(x) {c(x)}))))

# ??? why geneIDcode is gene_id for all
names(reads_all_alm_by_cigar_chimera_details_geneid_2) <- c("GeneIDcode", "GeneName")

reads_all_alm_by_cigar_chimera_details <- cbind(reads_all_alm_by_cigar_chimera_details, reads_all_alm_by_cigar_chimera_details_geneid_2)
reads_all_alm_by_cigar_chimera_details$alm_size <- abs(reads_all_alm_by_cigar_chimera_details$end - reads_all_alm_by_cigar_chimera_details$start)

cigar_parsed_grouped_ext <- paste0(dest_dir,basename_summary, ".CIGAR_parsed.grouped.ext.tsv")

write.table(x = reads_all_alm_by_cigar_chimera_details, file = cigar_parsed_grouped_ext , sep = "\t", col.names = T, row.names = F, na = '')
print(paste(cigar_parsed_grouped_ext, "done."))





# 1) filter data**********************************************************************************************************************
# remove too short hg19 alignements because other potential genomes or situations that I cannot manage so far
# A: too short -> mm9
# B: multi chr mapping
# C: un-manageble situation LRLR*
# D: contraddictions: chimera but then map only in human (not chrV)
reads_all_alm_by_cigar_chimera_details_toremove <- unique(
    as.character(reads_all_alm_by_cigar_chimera_details[which(
        ( reads_all_alm_by_cigar_chimera_details$alm_size < opt$alm_size &
              !(reads_all_alm_by_cigar_chimera_details$chr %in% c("chrV")) ) |
            ( reads_all_alm_by_cigar_chimera_details$target_genome_multimapping == TRUE ) |
            ( nchar(reads_all_alm_by_cigar_chimera_details$target_genome_position_wrt_junction) > 3 ) |
            ( reads_all_alm_by_cigar_chimera_details$chimera == FALSE & 
                  !(reads_all_alm_by_cigar_chimera_details$chr %in% c("chrV"))  )
    ), c("name")]))

reads_all_alm_by_cigar_chimera_details_rel <- reads_all_alm_by_cigar_chimera_details[which( 
    !(reads_all_alm_by_cigar_chimera_details$name %in% reads_all_alm_by_cigar_chimera_details_toremove) ), ]

#write.csv(reads_all_alm_by_cigar_chimera_details_rel,"joins.csv",row.names = FALSE)
cigar_parsed_grouped_ext_cleaned<- paste0(dest_dir,basename_summary, ".CIGAR_parsed.grouped.ext.CLEANED.tsv")

write.table(x = reads_all_alm_by_cigar_chimera_details_rel, file = cigar_parsed_grouped_ext_cleaned , sep = "\t", col.names = T, row.names = F, na = '')
print(paste(cigar_parsed_grouped_ext_cleaned, "done."))

q0_reads_with_aav <- sqldf("select sample_label, count(distinct name) as nReads_withAAV
                            from reads_all_alm_by_cigar_chimera_details_rel 
                            where 1
                            --integration_locus > 0 and chimera like TRUE 
                            group by sample_label 
                           ")

q0_reads_with_aav_len <- sqldf("select sample_label, name, read_len
                            from reads_all_alm_by_cigar_chimera_details_rel 
                            where sample_label not in ('Test')
                            --integration_locus > 0 and chimera like TRUE 
                            group by name
                           ")

q0_reads_with_aav_chimera <- sqldf("select sample_label, count(distinct name) as nReads_Chimeric, count(distinct (chr||':'||integration_locus) ) as nIS
                            from reads_all_alm_by_cigar_chimera_details_rel 
                            where chimera == 1 
                              and integration_locus > 0
                            -- chr in ('chrV')
                            -- integration_locus > 0 and chimera like TRUE 
                            group by sample_label 
                           ")

q0_reads_stats <- merge(x = q0_reads_with_aav, 
                        y = q0_reads_with_aav_chimera,
                        by = c("sample_label"), all.x = T)


write.table(x = q0_reads_stats, file = paste0(dest_dir, "q0_reads_stats.tsv"), sep = "\t", col.names = T, row.names = F, na = '', quote = F)
#write.xlsx(x = q0_reads_stats, file = paste0(dest_dir, "q0_reads_stats.xlsx"))
print(paste0(dest_dir, "q0_reads_stats.tsv","done."))

