###############################################################
#' @title Parse read CIGAR and get info of chimera
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-10-05
#'
#' @rdname parseCigarForChimera
#' @docType methods
#' @aliases parseCigarForChimera
#'
#' @param df_alm input df with alignment results
#' @param min_cigar_alm_width minimum alignment subread widt for each CIGAR string. default = 5. This value will be applied to ALL tags (indels included)
#'
#' @return a df of sureads results 
#' @usage TODO
#' @description parse CIGAR string and get data for AAV studies (and general chimera). Definition of chimera (as subread): the ONLY portion of the read (identified by CIGAR) that flanks target genome alignment.
#' @note : 
#'
###############################################################
parseCigarForChimera <- function(df_alm, 
                                 col_read_name = "name",
                                 cols_to_include = c("start", "end", "gene_chr", "gene_annotationsource", "gene_elementtype", "gene_start", "gene_end", "gene_strand", "gene_details", "distance_to_gene", "sample", "sourcefile", "score"),
                                 key_cols = c("name", "cigar"),
                                 cigar_opts_to_keep = c('M', 'D', 'I'),
                                 vector_chr_string = "chrV",
                                 min_cigar_alm_width = 5, 
                                 quiet = T
                                 ){
  require(sqldf)
  require(GenomicAlignments)
  
  message(paste("[AP]\tParse CIGAR string to find query alignments (re-arrangements) and Chimera"))
  reads_id <- levels(factor(df_alm[, col_read_name]))
  df_query_summary <- NULL
  cigar_opts_to_keep <- paste0("('", paste(cigar_opts_to_keep, collapse = "', '"), "')")
  for (read_id in reads_id) {
    if (!(quiet)) {
      message(paste0("[AP]\t\t", read_id))
    }
    df_ir <- NULL
    slice_t <- df_alm[which(df_alm$name == read_id),]
    for (cigar_string in as.character(slice_t[,c("cigar")])) {
      ir_t <- cigarRangesAlongQuerySpace(cigar_string, 
                                         flag=NULL, before.hard.clipping=T, after.soft.clipping=F, 
                                         ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE, with.ops=T)
      ir_t_len <- as.numeric(width(cigarRangesAlongQuerySpace(cigar_string, flag=NULL, before.hard.clipping=T, after.soft.clipping=F, ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=T, with.ops=T)))
      ir_t_chr <- slice_t[which(slice_t$cigar == cigar_string), "chr"]
      ir_t_strand <- slice_t[which(slice_t$cigar == cigar_string), "strand"]
      cigar_df <- as.data.frame(ir_t[[1]])
      names(cigar_df) <- c("read_alm_start", "read_alm_end", "read_alm_width", "cigar_flag")
      read_info <- data.frame( "name" = rep(read_id, nrow(cigar_df)),
                               "cigar" = rep(cigar_string, nrow(cigar_df)),
                               "chr" = ir_t_chr,
                               "strand" = ir_t_strand,
                               "read_len" = ir_t_len
      )
      bind_df <- cbind(read_info, cigar_df)
      # if read is reverse, invert iranges 
      if (levels(factor(bind_df$strand)) == '-') {
        bind_df$query_start <- apply( bind_df[c("read_len", "read_alm_start", "read_alm_end", "read_alm_width")], 1, function(x) {ifelse( x[4] > 0, 
                                                                                                                                          ir_t_len - x[3] + 1, 
                                                                                                                                          ir_t_len - x[2] + 1
        )
        }
        )
        bind_df$query_end <- ir_t_len - bind_df$read_alm_start + 1
      } else {
        bind_df$query_start <- bind_df$read_alm_start
        bind_df$query_end <- bind_df$read_alm_end
      }
      bind_df <- merge(x = bind_df, y = slice_t[which(slice_t$cigar == cigar_string), c(key_cols, cols_to_include)], by = key_cols, all.x = T)
      
      # merge df with all other cigar strings
      if (length(df_ir) == 0) {
        df_ir <- bind_df
      } else {
        df_ir <- rbind(df_ir, bind_df)
      } # if (length(df_ir) == 0)
    } # for (cigar_string in as.character(slice_t[,c("cigar")])) 
    
    # summary 
    read_query_summary <- sqldf( paste0("select * from df_ir where cigar_flag in ", cigar_opts_to_keep, " and read_alm_width >= ", min_cigar_alm_width, " order by query_start") )
    
    # create a string of read composition
    read_query_string_alm <- sqldf( paste0("select cigar, chr from df_ir where cigar_flag in ", cigar_opts_to_keep, " and read_alm_width >= ", min_cigar_alm_width, " group by cigar order by query_start") )
    read_query_string_alm_type <- gsub("chr", "", paste0(as.character(read_query_string_alm[,"chr"]), collapse = '-'))
    # do we have multi-mapping on target genome?
    read_query_string_alm_multimapping <- ifelse( (length(levels(factor( (read_query_string_alm[which( !(read_query_string_alm$chr %in% c(vector_chr_string)) ), "chr"]) ))) > 1), TRUE, FALSE)
    
    # ---- find chimera(s) ----
    # init chimera
    read_query_summary$chimera <- FALSE
    read_query_summary$junction <- FALSE
    read_query_summary$target_genome_position_wrt_junction <- NA
    read_query_summary$read_alm_type <- read_query_string_alm_type
    read_query_summary$target_genome_multimapping <- read_query_string_alm_multimapping
    
    vector_subread_chimera_cigar <- NULL
    alignment_on_target_genome <- F
    target_genome_position_wrt_junction <- NA
    target_genome_position_wrt_junction_found <- F
    if (nrow(read_query_summary) > 1) {
      for (i in seq(1, (nrow(read_query_summary)-1))) {
        # case A: V+G+
        # case B: G+V+
        # case C: G+V+G+ -> this case is included in A and B
        
        # case A: V+G+
        if (read_query_summary[i,"chr"] == vector_chr_string & read_query_summary[(i+1),"chr"] != vector_chr_string) {
          vector_subread_chimera_cigar <- c(vector_subread_chimera_cigar, as.character(read_query_summary[i,"cigar"]))
          if ( !target_genome_position_wrt_junction_found ) {
            target_genome_position_wrt_junction <- "R"
            target_genome_position_wrt_junction_found <- T
          } else {
            target_genome_position_wrt_junction <- paste0(target_genome_position_wrt_junction, "R")
          }
        } 
        # case B: G+V+
        if (read_query_summary[i,"chr"] != vector_chr_string & read_query_summary[(i+1),"chr"] == vector_chr_string) {
          vector_subread_chimera_cigar <- c(vector_subread_chimera_cigar, as.character(read_query_summary[(i+1),"cigar"]))
          if ( !target_genome_position_wrt_junction_found ) {
            target_genome_position_wrt_junction <- "L"
            target_genome_position_wrt_junction_found <- T
          } else {
            target_genome_position_wrt_junction <- paste0(target_genome_position_wrt_junction, "L")
          }
        } 
        # end if
      } # for (i in seq(1, (nrow(read_query_summary)-1)))   
    } else {
      vector_subread_chimera_cigar <- c()
    } # if (nrow(read_query_summary) > 1)
    read_query_summary$target_genome_position_wrt_junction <- target_genome_position_wrt_junction
    read_query_summary$junction <- ifelse(read_query_summary$cigar %in% vector_subread_chimera_cigar, TRUE, FALSE)
    # read_query_summary$chimera <- ifelse("TRUE" %in% as.character(read_query_summary$junction), TRUE, FALSE)
    read_query_summary$chimera <- ifelse(length(vector_subread_chimera_cigar) > 0, TRUE, FALSE)
    # now find the aav junction side, to best characterize the aav juction - integration. complex cases are here NOT managed (TODO)
    read_query_summary$vector_junction_locus <- ifelse(read_query_summary$junction == T, 
                                                       ifelse(read_query_summary$strand =='+',
                                                              ifelse(read_query_summary$target_genome_position_wrt_junction %in% c('L', 'LR'), 
                                                                     read_query_summary$start,
                                                                     read_query_summary$end), 
                                                              ifelse(read_query_summary$target_genome_position_wrt_junction %in% c('L', 'LR'), 
                                                                     read_query_summary$end,
                                                                     read_query_summary$start)),
                                                       NA)
    read_query_summary$query_order <- seq(1, nrow(read_query_summary), 1)
    # find locus
    locus_found_and_reported <- F
    indels <- NULL
    
    # read_query_summary$integration_locus <- apply(read_query_summary[c("chr", "strand", "target_genome_position_wrt_junction", "start", "end", "chimera")], 1, function(x) {
    #         .locus <- NA    
    #         if (x[1] != vector_chr_string & x[6] == "TRUE" & !is.na(x[3])) { 
    #           # simple cases L, R with - +
    #           if (x[2] == "+" & x[3] == "L") { .locus <- x[5]}
    #           if (x[2] == "+" & x[3] == "R") { .locus <- x[4]}
    #           if (x[2] == "-" & x[3] == "L") { .locus <- x[4]}
    #           if (x[2] == "-" & x[3] == "R") { .locus <- x[5]}
    #           # complex cases LR with +-
    #           if (!(locus_found_and_reported) & x[3] == "LR" & x[2] == "+") { 
    #             .locus <- x[5]
    #             locus_found_and_reported <- T
    #           } # LR+
    #           if (!(locus_found_and_reported) & x[3] == "LR" & x[2] == "-") { 
    #             .locus <- x[4]
    #             locus_found_and_reported <- T
    #           } # LR-
    #           
    #           # # just look at indels at AAV site with complex cases
    #           # if (locus_found_and_reported & x[3] == "LR" & x[2] == "+") { 
    #           #   .locus <- x[5]
    #           #   locus_found_and_reported <- T
    #           # } # LR+
    #           
    #         } # if (x[1] != vector_chr_string) 
    #         .locus
    #       } # function
    #     ) # apply
    
    integration_locus_tobind <- NA
    # slice_data <- sqldf(paste0('select distinct chr, strand, target_genome_position_wrt_junction, start, end, chimera from read_query_summary where chr not like "', vector_chr_string, '"'))
    slice_data <- sqldf(paste0('select distinct chr, strand, target_genome_position_wrt_junction, start, end, chimera from read_query_summary where 1'))
    locus_found_and_reported <- FALSE
    for (elem in seq(1, nrow(slice_data)) ) {
      .locus <- integration_locus_tobind
      # message(paste(slice_data[elem, "name"], slice_data[elem, "start"], slice_data[elem, "end"], slice_data[elem, "target_genome_position_wrt_junction"], "\n"))
      # message(paste("--> e", elem, as.character(locus_found_and_reported), "\n"))
      if ( slice_data[elem, "chr"] != vector_chr_string & slice_data[elem, "chimera"] == "TRUE" & !is.na(slice_data[elem, "target_genome_position_wrt_junction"]) ) { 
        # simple cases L, R with - +
        if (slice_data[elem, "strand"] == "+" & slice_data[elem, "target_genome_position_wrt_junction"] == "L") { .locus <- slice_data[elem, "end"]}
        if (slice_data[elem, "strand"] == "+" & slice_data[elem, "target_genome_position_wrt_junction"] == "R") { .locus <- slice_data[elem, "start"]}
        if (slice_data[elem, "strand"] == "-" & slice_data[elem, "target_genome_position_wrt_junction"] == "L") { .locus <- slice_data[elem, "start"] }
        if (slice_data[elem, "strand"] == "-" & slice_data[elem, "target_genome_position_wrt_junction"] == "R") { .locus <- slice_data[elem, "end"] }
        # complex cases LR with +-
        if (!(locus_found_and_reported) & slice_data[elem, "target_genome_position_wrt_junction"] == "LR" & slice_data[elem, "strand"] == "+") { 
          .locus <- slice_data[elem, "end"]
          locus_found_and_reported <- T
          # message(paste("\tFOUND:", as.numeric(.locus), as.character(locus_found_and_reported), "\n"))
        } # LR+
        if (!(locus_found_and_reported) & slice_data[elem, "target_genome_position_wrt_junction"] == "LR" & slice_data[elem, "strand"] == "-") { 
          .locus <- slice_data[elem, "start"]
          locus_found_and_reported <- T
          # message(paste("\tFOUND:", as.numeric(.locus), as.character(locus_found_and_reported), "\n"))
        } # LR-
        
        # # just look at indels at AAV site with complex cases
        # if (locus_found_and_reported & x[3] == "LR" & x[2] == "+") { 
        #   .locus <- x[5]
        #   locus_found_and_reported <- T
        # } # LR+
        integration_locus_tobind <- as.numeric(.locus)    
      } # if (x[1] != vector_chr_string) 
    } # for (elem in nrow(slice_data))
    
    read_query_summary$integration_locus <- ifelse(!(read_query_summary$chr %in% vector_chr_string), integration_locus_tobind, NA)
    
    # merge df with all other cigar strings
    if (length(df_query_summary) == 0) {
      df_query_summary <- read_query_summary
    } else {
      df_query_summary <- rbind(df_query_summary, read_query_summary)
    } # if (length(df_ir) == 0)
  } # for (header in reads_id) 
  
  return (df_query_summary)
}


plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, ...) {
  height <- 1
  # if (is(xlim, "IntegerRanges"))
  xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}

extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}



###############################################################
#' @title Plot chimera in CIRCOS
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-11-19
#'
#' @rdname plot_circos_rearrangement_imp
#' @docType methods
#' @aliases plot_circos_rearrangement_imp
#'
#' @param df_alm input df with alignment results
#' @param min_cigar_alm_width minimum alignment subread widt for each CIGAR string. default = 5. This value will be applied to ALL tags (indels included)
#'
#' @return plot
#' @description parse CIGAR string and get data for AAV studies (and general chimera). Definition of chimera (as subread): the ONLY portion of the read (identified by CIGAR) that flanks target genome alignment.
#' @note : 
#' @usage 
#' plot_circos_rearrangement_imp(df = alm_reads_for_rearrangements_withstatsbyread, sample_read = "m64047_200620_005306/19269857/ccs", 
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.circos.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), vector_chr = "chrV", zoomed_chr_index = 23, vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband")
###############################################################
plot_circos_rearrangement_imp <- function(df, sample_read, outfile_png, outfile_pdf,
                                          cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), 
                                          vector_chr = "chrV", 
                                          zoomed_chr_index = 23, 
                                          vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
                                          species = "hg19",
                                          bp_res= 300,
                                          color_vector_fwd="orange",
                                          color_vector_rev="green",
                                          color_target_rev="violet",
                                          color_target_fwd="blue") {
  if (FALSE %in% (cols_to_search %in% colnames(df)) ) {
    message(paste0("[AP]\tERROR: Columns in df are different from expected: ", paste(cols_to_search, collapse = ',')))
    return (NULL)
  } else {
    # slice data
    t <- sample_read
    slice_t <- df[which(df$name == t),]
    slice_t <- slice_t[order(slice_t$query_start),]
    bed_1 <- slice_t[(1:(nrow(slice_t)-1)), cols_to_search]
    names(bed_1) <- c("chr", "start", "end", "value1")
    bed_2 <- slice_t[(2:(nrow(slice_t))), cols_to_search]
    names(bed_2) <- c("chr", "start", "end", "value1")
    
    # orientations
    bed_all_stranded <- slice_t[c("chr", "start", "end", "strand")]
    bed_all_stranded$strand <- ifelse(bed_all_stranded$strand == '+', 1, 0)
    names(bed_all_stranded) <- c("chr", "start", "end", "value1")
    
    # get start-end of each read
    bed_all_long <- slice_t[c("chr", "start", "end", "strand", "query_start", "query_end")]
    outbed_line_sorted_start <- NULL
    if (nrow(bed_all_long) > 1) {
      for (i in seq(1, (nrow(bed_all_long)-1))) {
        if (bed_all_long[i,"strand"] == '+') {
          from_chr <- bed_all_long[i,c("chr")]
          from_pos <- bed_all_long[i,c("end")]
          if (bed_all_long[i+1,"strand"] == '+') {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("start")]
          } else {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("end")]
          } # if (bed_all_long[i+1,"strand"] == '+')
        } else {
          from_chr <- bed_all_long[i,c("chr")]
          from_pos <- bed_all_long[i,c("start")]
          if (bed_all_long[i+1,"strand"] == '+') {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("start")]
          } else {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("end")]
          } # if (bed_all_long[i+1,"strand"] == '+')
        } # if (bed_all_long[i,"strand"] == '+')
        if (length(outbed_line_sorted_start) > 0) {
          outbed_line_sorted_start <- rbind(outbed_line_sorted_start, data.frame("from_chr" = from_chr,
                                                                                 "from_pos" = from_pos,
                                                                                 "to_chr" = to_chr,
                                                                                 "to_pos" = to_pos))  
        } else {
          outbed_line_sorted_start <- data.frame("from_chr" = from_chr,
                                                 "from_pos" = from_pos,
                                                 "to_chr" = to_chr,
                                                 "to_pos" = to_pos)
        } # if (length(outbed_line_sorted_start) > 0)
        
      } # for (i in seq(1, nrow(bed_all_long)))
    } # if (nrow(bed_all_long) > 1)
    outbed_line_sorted_start$value1 <- 1
    bed_1_line <- outbed_line_sorted_start[grep("from", colnames(outbed_line_sorted_start))]
    names(bed_1_line) <- c("chr", "start")
    bed_1_line$end <- bed_1_line$start
    bed_1_line$value1 <- 1
    bed_2_line <- outbed_line_sorted_start[grep("to", colnames(outbed_line_sorted_start))]
    names(bed_2_line) <- c("chr", "start")
    bed_2_line$end <- bed_2_line$start
    bed_2_line$value1 <- 1
    
    target_chr <- gsub("chr", "", setdiff( levels(factor(slice_t$chr)), vector_chr))
    # if 0 -> only AAV (chrV), else target too
    if (length(target_chr) > 0) {
      # go into genomics
      human_cytoband <- read.cytoband(species = species)$df
      vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
      cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
      cytoband <- read.cytoband(cytoband_rbind)
      
      cytoband_df = cytoband$df
      chromosome = cytoband$chromosome
      # acquire target chr
      if (target_chr %in% c("X", "Y")) {
        if (target_chr == "X") {target_chr <- grep("X", chromosome)}
        if (target_chr == "Y") {target_chr <- grep("Y", chromosome)}
      } else {
        target_chr <- as.numeric(target_chr)
      }
      # determine range of zoom
      xrange = c(cytoband$chr.len, cytoband$chr.len[vector_chr])
      normal_chr_index = target_chr
      zoomed_chr_index = zoomed_chr_index
      sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                       xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index])) 
      # PNG
      png(file = outfile_png, height=5, width=5, units = "in", res = 600)
      extended <-extend_chromosomes(cytoband_df, vector_chr)
      chr_target_ff<- bed_all_stranded[bed_all_stranded$chr!=vector_chr,]$chr
      s <- bed_all_stranded[bed_all_stranded$chr==chr_target_ff,]$start
      e <- bed_all_stranded[bed_all_stranded$chr==chr_target_ff,]$end
      m <- (s+e)/2
      r <- extended[(extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m) ,]
      increase<-slice_t$read_len[1]+1000
      r$V2 <- m-(increase/2)
      r$V3 <- m+(increase/2)
      extended[(extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m) ,] <- r
      itr5 <- extended[extended$V4=="ITR.5p" & extended$V1==vector_chr,]$V2
      itr3 <- extended[extended$V4=="ITR.3p" & extended$V1==vector_chr,]$V3
      extended2 <- extended[((extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m)) | (extended$V1!=chr_target_ff & extended$V1!=vector_chr)
                            | ((extended$V1==vector_chr & (extended$V2>=itr5 & extended$V3<=itr3))),]
      extended2$V6 <-''
      extended2[extended2$V1==vector_chr,]$V6 <- extended2[extended2$V1==vector_chr,]$V4
        
      circos.initializeWithIdeogram(extended2, 
                                    chromosome.index = c(chromosome[target_chr], vector_chr),
                                    sector.width = sector.width, tickLabelsStartFromZero = FALSE, major.by = bp_res)
    
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      l <- seq(m-5000,m+5000,600)
      for (subread in seq(1, nrow(bed_all_stranded))) {
          circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                              track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                              panel.fun = function(region, value, ...) {
                                                  x1=region$start
                                                  y1=region$end
                                                  if(value==1){
                                                    v<-value
                                                    if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                      v<-v+2
                                                      circos.arrow(x1,y1, col = f(v), 
                                                                   border = 1, arrow.head.length = cm_x(0.2) ) 
                                                  }
                                                  else{
                                                    v<-value
                                                    if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                      v<-v+2
                                                      circos.arrow(x1, y1, col = f(v), 
                                                                   border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                                  }
                                                 
                                                  
                                                  i = getI(...)
                                                  cell.xlim = get.cell.meta.data("cell.xlim")
                                                  # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                              })
      

      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
      # PDF
      pdf(file = outfile_pdf, height=5, width=5)
      circos.initializeWithIdeogram(extended2, 
                                    chromosome.index = c(chromosome[target_chr], vector_chr),
                                    sector.width = sector.width, tickLabelsStartFromZero = FALSE, major.by = bp_res)
    
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                              track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                              panel.fun = function(region, value, ...) {
                                                x1=region$start
                                                y1=region$end
                                                if(value==1){
                                                  v<-value
                                                  if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                    v<-v+2
                                                  circos.arrow(x1,y1, col = f(v), 
                                                               border = 1, arrow.head.length = cm_x(0.2) ) 
                                                }
                                                else{
                                                  v<-value
                                                  if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                    v<-v+2
                                                  circos.arrow(x1, y1, col = f(v), 
                                                               border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                                }
                                                
                                                
                                                i = getI(...)
                                                cell.xlim = get.cell.meta.data("cell.xlim")
                                                # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                              })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
      
    } else {
      vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
      # cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
      cytoband <- read.cytoband(vector_cytoband)
      
      cytoband_df = cytoband$df
      chromosome = cytoband$chromosome
      # PNG
      png(file = outfile_png, height=5, width=5, units = "in", res = 600)
      circos.par("start.degree" = 90)
      
      itr5 <- cytoband_df[cytoband_df$V4=="ITR.5p" & cytoband_df$V1==vector_chr,]$V2
      
      itr3 <- cytoband_df[cytoband_df$V4=="ITR.3p" & cytoband_df$V1==vector_chr,]$V3 
      
      cytoband_df <- cytoband_df[(cytoband_df$V2>=itr5 & cytoband_df$V3<=itr3),]
      
      circos.initializeWithIdeogram(cytoband = cytoband_df,
                                    chromosome.index = c(vector_chr), 
                                    tickLabelsStartFromZero = FALSE, major.by = bp_res)
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
      # PDF
      pdf(file = outfile_pdf, height=5, width=5)
      circos.par("start.degree" = 90)
      circos.initializeWithIdeogram(cytoband = cytoband_df, chromosome.index = c(vector_chr),
                                    tickLabelsStartFromZero = FALSE, major.by = bp_res)
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
    }  # if (length(target_chr) > 0) 
  } # if (FALSE %in% (cols_to_search %in% colnames(df)) )
  
} # function



###############################################################
#' @title Plot rearrangements linearized
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_rearrangement
#' @docType methods
#' @aliases plot_linear_rearrangement
#'
#' @param df_alm input df with alignment results
#' @param reads_to_plot vector with names of the reads to plot
#' @param cytoband_file path to cytoband file NOTE: the file must have the following
#'                      header: chr	start	end	name	gieStain 
#' @param vector_length vector genome length
#' @return 
#' @description Plot the vector portions of the read
#' @note : 
#' @usage 
#' plot_linear_rearrangement(df = alm_reads_for_rearrangements_withstatsbyread, reads_to_plot = reads_to_plot<-c("m64047_200324_183710/83362987/ccs" ,"m64047_200324_183710/102237007/ccs"),
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' vector_chr = "chrV", cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
#' annotation_file="source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.bed")
plot_linear_rearrangement <- function(df_alm, reads_to_plot, outfile_png, outfile_pdf,
                                      cytoband_file, annotation_file, vector_chr = "chrV",
                                      vector_length=6000,
                                      height=8, width=8, res = 300,
                                      offset_arrow=110
                                      ){
  library(karyoploteR)
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 15
  plot.params$data1height <- 10
  plot.params$data2height<-500
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  increase1<-0.2/length(reads_to_plot)
  increase2<-0.22/length(reads_to_plot)
  increase3<-0.3/length(reads_to_plot)
  
  i_r0<-0.1
  i_r1<-0.15
  y0<-0.1
  y1<-y0+increase1
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                     'target_genome_position_wrt_junction')]
    
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    for (i in 1:nrow(x))
    {
      if (x[i,]$strand=='+')
      {       
        i_col<-"orange"
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-"green"
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.05,lwd=15,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.05,lwd=13,lend='butt',ljoin='mitre', data.panel = 2)
      print(i_r0)
      print(i_r1)
      
      i_r0<-i_r0+increase1
      i_r1<-i_r1+increase1
    }
    
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    if (!is.na(x[1,]$vector_junction_locus))
    {
      j<-x[1,]$vector_junction_locus
     # if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
      #  j<-j+offset_arrow
      #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
        #j<-j-offset_arrow
      kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
    }
    if (nrow(x) > 1)
    {
      for (i in 2:nrow(x))
      {
        
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y0+0.0025*length(reads_to_plot) ,y1=y1,length=0.1, lty=2, data.panel = 2)
        y0<-y0+increase2
        y0_line<-y0_line+increase2
        y1<-y1+increase2
        if (!is.na(x[i,]$vector_junction_locus))
        {
          j<-x[i,]$vector_junction_locus
          #if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
           # j<-j+offset_arrow
          #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
           # j<-j-offset_arrow
          print(x[i,]$vector_junction_locus)
          kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
        }
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end-(offset_arrow/2) else x[i,]$start+(offset_arrow/2)
        
      }
    }
    y0<-y0+increase1
    y1<-y1+increase1
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(y0,y0),lty=2,data.panel = 2,col='grey',lwd=2)
    i_r0<-i_r0+increase3
    i_r1<-i_r1+increase3
    y0<-y0+increase3
    y1<-y1+increase3
    y_text<-y0-increase1
    plotted_segs<-plotted_segs+1
    y0_line<-y0_line+increase1+increase3+increase1/((30/length(reads_to_plot))*plotted_segs)
  }
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 15
  plot.params$data1height <- 10
  plot.params$data2height<-500
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  increase1<-0.2/length(reads_to_plot)
  increase2<-0.22/length(reads_to_plot)
  increase3<-0.3/length(reads_to_plot)
  
  i_r0<-0.1
  i_r1<-0.15
  y0<-0.1
  y1<-y0+increase1
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                         'target_genome_position_wrt_junction')]
    
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    for (i in 1:nrow(x))
    {
      if (x[i,]$strand=='+')
      {       
        i_col<-"orange"
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-"green"
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.05,lwd=15,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.05,lwd=13,lend='butt',ljoin='mitre', data.panel = 2)
      print(i_r0)
      print(i_r1)
      
      i_r0<-i_r0+increase1
      i_r1<-i_r1+increase1
    }
    
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    if (!is.na(x[1,]$vector_junction_locus))
    {
      j<-x[1,]$vector_junction_locus
      # if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
      #  j<-j+offset_arrow
      #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
      #j<-j-offset_arrow
      kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
    }
    if (nrow(x) > 1)
    {
      for (i in 2:nrow(x))
      {
        
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y0+0.0025*length(reads_to_plot) ,y1=y1,length=0.1, lty=2, data.panel = 2)
        y0<-y0+increase2
        y0_line<-y0_line+increase2
        y1<-y1+increase2
        if (!is.na(x[i,]$vector_junction_locus))
        {
          j<-x[i,]$vector_junction_locus
          #if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
          # j<-j+offset_arrow
          #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
          # j<-j-offset_arrow
          print(x[i,]$vector_junction_locus)
          kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
        }
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end-(offset_arrow/2) else x[i,]$start+(offset_arrow/2)
        
      }
    }
    y0<-y0+increase1
    y1<-y1+increase1
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(y0,y0),lty=2,data.panel = 2,col='grey',lwd=2)
    i_r0<-i_r0+increase3
    i_r1<-i_r1+increase3
    y0<-y0+increase3
    y1<-y1+increase3
    y_text<-y0-increase1
    plotted_segs<-plotted_segs+1
    y0_line<-y0_line+increase1+increase3+increase1/((30/length(reads_to_plot))*plotted_segs)
  }
  dev.off()
}




###############################################################
#' @title Plot Read linearized
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_rearrangement
#' @docType methods
#' @aliases plot_linear_read
#'
#' @param df_alm input df with alignment results
#' @param read_to_plot name of the read to plot
#' 
#' @return 
#' @description Plot the different alignments against the read
#' @note : 
#' @usage 
#' plot_linear_read(df,read_to_plot = "m64047_200324_183710/102237007/ccs",
#'        outfile_png = "read.png",outfile_pdf = "read.pdf")
plot_linear_read <- function(df_alm, read_to_plot, outfile_png, outfile_pdf,
                                      vector_chr = "chrV",
                                      height=8, width=8, res = 300,
                                      offset_arrow=110,
                                      color_vector_fwd="orange",
                                      color_vector_rev="green", 
                                      color_target_fwd="blue",
                                      color_target_rev="violet"
){
  library(karyoploteR)
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction')
  i_r0<-1
  i_r1<-0.95
  for (i in 1:nrow(x))
  {
    if (x[i,]$strand=='+')
    {       
      ar_start<-x[i,]$start
      ar_end<-x[i,]$end-offset_arrow
      if (ar_start>=ar_end)
        ar_start<-ar_end-1
      if (x[i,]$genome==vector_chr)
        i_col<-color_vector_fwd
      else
        i_col<-color_target_fwd
    }
    else
    {        
      ar_start<-x[i,]$end
      ar_end<-x[i,]$start+offset_arrow
      if (ar_start<=ar_end)
        ar_start<-ar_end+1
      if (x[i,]$genome=='chrV')
        i_col<-color_vector_rev
      else
        i_col<-color_target_rev
    }
    kpArrows(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.1,lwd=15,lend='butt',ljoin='mitre')
    kpArrows(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.1,lwd=13,lend='butt',ljoin='mitre')
    if (x[i,]$junction==TRUE)
    {
        if(x[i,]$target_genome_position_wrt_junction=='R')
          j<-x[i,]$end
        else
          j<-x[i,]$start
     
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0,i_r0-2*(i_r0-i_r1)),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-0.1
    i_r1<-i_r1-0.1
  }
  legend(x = "bottomleft", fill = c(color_vector_fwd, color_vector_rev,color_target_fwd,color_target_rev), title=c('Strand and Genome'),legend = legend(x = "bottomleft", fill = c("orange", "green","blue","violet"), title=c('Strand and Genome'),legend = c("Vector (+)", "Vector (-)", "Target (+)","Target (-)")))
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction')
  i_r0<-1
  i_r1<-0.95
  for (i in 1:nrow(x))
  {
    if (x[i,]$strand=='+')
    {       
      ar_start<-x[i,]$start
      ar_end<-x[i,]$end-offset_arrow
      if (ar_start>=ar_end)
        ar_start<-ar_end-1
      if (x[i,]$genome==vector_chr)
        i_col<-"orange"
      else
        i_col<-"blue"
    }
    else
    {        
      ar_start<-x[i,]$end
      ar_end<-x[i,]$start+offset_arrow
      if (ar_start<=ar_end)
        ar_start<-ar_end+1
      if (x[i,]$genome=='chrV')
        i_col<-"green"
      else
        i_col<-"violet"
    }
    kpArrows(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.1,lwd=15,lend='butt',ljoin='mitre')
    kpArrows(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.1,lwd=13,lend='butt',ljoin='mitre')
    if (x[i,]$junction==TRUE)
    {
      if(x[i,]$target_genome_position_wrt_junction=='R')
        j<-x[i,]$end
      else
        j<-x[i,]$start
      
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0,i_r0-2*(i_r0-i_r1)),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-0.1
    i_r1<-i_r1-0.1
  }
  legend(x = "bottomleft", fill = c(color_vector_fwd, color_vector_rev,color_target_fwd,color_target_rev), title=c('Strand and Genome'),legend = legend(x = "bottomleft", fill = c("orange", "green","blue","violet"), title=c('Strand and Genome'),legend = c("Vector (+)", "Vector (-)", "Target (+)","Target (-)")))
  dev.off()
}




###############################################################
#' @title Plot rearrangements linearized
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_rearrangement
#' @docType methods
#' @aliases plot_linear_rearrangement
#'
#' @param df_alm input df with alignment results
#' @param reads_to_plot vector with names of the reads to plot
#' @param cytoband_file path to cytoband file NOTE: the file must have the following
#'                      header: chr	start	end	name	gieStain 
#' @param vector_length vector genome length
#' @param annotation_file bed file with annotation
#' @param color_bed_file bed file with rgb colors, include all intervals, see example one
#' @param lwd_black width of the external arrow (this is used to make the arrow border, should be always >= lwd_color)
#' @param lwd_color width of the arrow
#' @param offset_arrow once you have set the two lwd pars you can use this parameter to be sure the arrow end in the right position
#' to do that change this parameter according to red lines which indicates the junction locus. You probably need to change this parameters for different
#' set of reads to plot
#' @return 
#' @description Plot the vector portions of the read
#' @note : 
#' @usage 
#' plot_linear_rearrangement(df = alm_reads_for_rearrangements_withstatsbyread, reads_to_plot = reads_to_plot<-c("m64047_200324_183710/83362987/ccs" ,"m64047_200324_183710/102237007/ccs"),
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' outfile_pdf = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.pdf", sep = ""), 
#' vector_chr = "chrV", offset_arrow=110, lwd_black = 160, lwd_color=148,
#' cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
#' annotation_file="source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.bed",
#' )
plot_linear_rearrangement_impr <- function(df_alm, reads_to_plot, outfile_png, outfile_pdf,
                                      cytoband_file, annotation_file,color_bed_file, vector_chr = "chrV",
                                      vector_length=6000,
                                      height=8, width=8, res = 300,
                                      offset_arrow=110, lwd_black = 160, lwd_color=148,
                                      color_vector_fwd="orange", color_vector_rev="green",
                                      ideogram_height = 5, aav_colored_thickness = 0.3
){
  library(karyoploteR)
  library(bezier)
  library(ggplot2)
  all_slice<-df_alm[is.element(df_alm$name, reads_to_plot),]
  all_slice<-all_slice[all_slice$chr==vector_chr,]
  lwd_black<-lwd_black/nrow(all_slice)
  lwd_color<-lwd_color/nrow(all_slice)
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  
  
 
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- ideogram_height
  plot.params$data1height <- 15
  plot.params$data2height<-200
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  bed_color <- read.table(file = color_bed_file, sep = '\t', header = FALSE)
  bed_color['rgb']<-gsub(",","",bed_color$V9)
  for(i in seq(1:(nrow(bed_color)-1)))
  {
    rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
    r<-rgb_string[1]
    g<-rgb_string[2]
    b<-rgb_string[3]
    kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i+1,]$V2,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  }
  i<-i+1
  rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
  r<-rgb_string[1]
  g<-rgb_string[2]
  b<-rgb_string[3]
  kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i,]$V3,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  
  increase<-0.65/nrow(all_slice)
  
  i_r0<-0.1
  i_r1<-i_r0+increase
  y0<-0.1
  y1<-y0+increase
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  y_r0_ar<-i_r0
  y_r1_ar<-i_r1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                         'target_genome_position_wrt_junction')]
    
    #y_r0_ar<-y_r0_ar+(0.12/length(reads_to_plot))*(plotted_segs-1)+0.01
    #y_r1_ar<-y_r1_ar+0.007*(plotted_segs-1)+0.01
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    for (i in 1:nrow(x))
    {
      
      if (x[i,]$strand=='+')
      {       
        i_col<-color_vector_fwd
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-color_vector_rev
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      #draw segments
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', length=0.05,lwd=lwd_black,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col,length=0.05,lwd=lwd_color,lend='butt',ljoin='mitre', data.panel = 2)
      if (!is.na(x[i,]$vector_junction_locus))
      {
        j<-x[i,]$vector_junction_locus
        kpLines(kp,chr = vector_chr,x=c(j,j), y=c(i_r0-increase/3,i_r0+increase/3),col='red', lwd=5, data.panel = 2)
      }
      #draw connecting curves
      if (i>1){
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        t<-seq(0,1,length=100)
        bz<-bezier(t,p=list(c(arrow_start,arrow_start,arrow_start,arrow_end,arrow_end,arrow_end),seq(i_r0-increase,i_r0, length=6)))
        kpLines(kp,chr="chrV",x=bz[,1],y=bz[,2],data.panel = 2, lty=1)
        #kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y_r0_ar ,y1=y_r0_ar+increase1-0.01,length=0.1, lty=2, data.panel = 2)
        
        
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end else x[i,]$start
        
        
      }
      i_r0<-i_r0+increase
      i_r1<-i_r1+increase
      
    }
    
    
    
    
    
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(i_r0,i_r0),lty=2,data.panel = 2,col='grey',lwd=2)
    y_text<-i_r0+increase/2
    i_r0<-i_r0+increase*1.5
    i_r1<-i_r1+increase*1.5
    
    
    
  }
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- ideogram_height
  plot.params$data1height <- 15
  plot.params$data2height<-200
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  bed_color <- read.table(file = color_bed_file, sep = '\t', header = FALSE)
  bed_color['rgb']<-gsub(",","",bed_color$V9)
  for(i in seq(1:(nrow(bed_color)-1)))
  {
    rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
    r<-rgb_string[1]
    g<-rgb_string[2]
    b<-rgb_string[3]
    kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i+1,]$V2,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  }
  i<-i+1
  rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
  r<-rgb_string[1]
  g<-rgb_string[2]
  b<-rgb_string[3]
  kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i,]$V3,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  
  increase<-0.65/nrow(all_slice)
  
  i_r0<-0.1
  i_r1<-i_r0+increase
  y0<-0.1
  y1<-y0+increase
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  y_r0_ar<-i_r0
  y_r1_ar<-i_r1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                         'target_genome_position_wrt_junction')]
    
    #y_r0_ar<-y_r0_ar+(0.12/length(reads_to_plot))*(plotted_segs-1)+0.01
    #y_r1_ar<-y_r1_ar+0.007*(plotted_segs-1)+0.01
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    for (i in 1:nrow(x))
    {
      
      if (x[i,]$strand=='+')
      {       
        i_col<-color_vector_fwd
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-color_vector_rev
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      #draw segments
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', length=0.05,lwd=lwd_black,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col,length=0.05,lwd=lwd_color,lend='butt',ljoin='mitre', data.panel = 2)
      if (!is.na(x[i,]$vector_junction_locus))
      {
        j<-x[i,]$vector_junction_locus
        kpLines(kp,chr = vector_chr,x=c(j,j), y=c(i_r0-increase/3,i_r0+increase/3),col='red', lwd=5, data.panel = 2)
      }
      #draw connecting curves
      if (i>1){
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        t<-seq(0,1,length=100)
        bz<-bezier(t,p=list(c(arrow_start,arrow_start,arrow_start,arrow_end,arrow_end,arrow_end),seq(i_r0-increase,i_r0, length=6)))
        kpLines(kp,chr="chrV",x=bz[,1],y=bz[,2],data.panel = 2, lty=1)
        #kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y_r0_ar ,y1=y_r0_ar+increase1-0.01,length=0.1, lty=2, data.panel = 2)
        
        
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end else x[i,]$start
        
        
      }
      i_r0<-i_r0+increase
      i_r1<-i_r1+increase
      
    }
    
    
    
    
    
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(i_r0,i_r0),lty=2,data.panel = 2,col='grey',lwd=2)
    y_text<-i_r0+increase/2
    i_r0<-i_r0+increase*1.5
    i_r1<-i_r1+increase*1.5
    
    
    
  }
dev.off()
}

###############################################################
#' @title Plot Read linearized with bed colors
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_read_colorized
#' @docType methods
#' @aliases plot_linear_read
#'
#' @param df_alm input df with alignment results
#' @param read_to_plot name of the read to plot
#' @param space_btw_rects space between two rects, default = 0.05
#' @param rect_thickness thickness of a rect, default = 0.05
#' @param cex_legend scaling factor of legend, default = 0.5 
#' @return 
#' @description Plot the different alignments against the read using the value in the bed as colors
#' @note : 
#' @usage 
#' plot_linear_read(df,read_to_plot = "m64047_200324_183710/102237007/ccs",
#'        outfile_png = "read.png",outfile_pdf = "read.pdf")
plot_linear_read_colorized <- function(df_alm, read_to_plot, outfile_png, outfile_pdf,
                             vector_chr = "chrV",
                             height=8, width=8, res = 300,
                             offset_arrow=110,
                             color_bed_file,
                             color_target_fwd="blue",
                             color_target_rev="violet",
                             space_btw_rects=0.05,
                             rect_thickness = 0.05,
                             cex_legend = 0.5
){
  library(karyoploteR)
  bed_color <- read.table(file = color_bed_file, sep = '\t', header = FALSE)
  sr<-bed_color[1,]
  sr$V1<-"start"
  sr$V2<-0
  sr$V3<-0
  se<-bed_color[nrow(bed_color),]
  se$V1<-"end"
  se$V2<-se$V3+1
  se$V3<-se$V3+1
  bed_color<-rbind(sr,bed_color,se)
  bed_color['rgb']<-gsub(",","",bed_color$V9)
  
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction','start','end')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction','orig_start','orig_end')
  space <- space_btw_rects
  thickness <- rect_thickness
  i_r0<-1
  i_r1<-i_r0 - thickness
  for (i in 1:nrow(x))
  {
    ar_start<-x[i,]$start
    ar_end<-x[i,]$end
    ltype = 1
    if(x[i,]$genome==vector_chr)
    {
      if(x[i,]$strand=='-')
      {
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start
      }
      
      range_start<-x[i,]$orig_start
      range_new_end<-x[i,]$orig_start
      range_end<-x[i,]$orig_end
      for (jj in seq(1:(nrow(bed_color)-1))){
        if  (range_start>=bed_color[jj,]$V2 & range_start<=bed_color[jj+1,]$V2)
        {
          rgb_string<-strsplit(bed_color[jj,]$rgb," ",fixed=TRUE)[[1]]
          r<-rgb_string[1]
          g<-rgb_string[2]
          b<-rgb_string[3]
          i_col<-rgb(r,g,b,maxColorValue = 255)
          range_new_end<-bed_color[jj+1,]$V2-1
          span<-range_new_end-range_start
          if (range_new_end>=range_end)
          {
            
            range_new_end<-range_end
            span<-range_new_end-range_start
            kpRect(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
            ar_start<-ar_start+span
            break
          }
          else
          {
            
            if(x[i,]$strand=='+')
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start+span, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
              ar_start<-ar_start+span
            }
            
            else
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start-span, y0=i_r0 ,y1=i_r0-thickness, col=i_col, lty=ltype)
              ar_start<-ar_start-span
            }
            
          }
          
          range_start<-bed_color[jj+1,]$V2
          
        }
      }
    }
    else
    {
      if (x[i,]$strand=='+')
      {
        i_col = color_target_fwd
        ltype=1
      }
      else
      {
        i_col = color_target_rev
        ltype=2
      }
      kpRect(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
    }
    if (x[i,]$junction==TRUE)
    {
      if(x[i,]$target_genome_position_wrt_junction=='R')
        j<-x[i,]$end
      else
        j<-x[i,]$start
      
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0, i_r0-thickness),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-thickness
    i_r1<-i_r1-thickness
    i_r0<-i_r0-space
    i_r1<-i_r1-space
  }
  legend(x = "bottomleft", lty = c(1,2), legend=c("+", "-"), title=c('Target Strand'), cex=cex_legend)
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction','start','end')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction','orig_start','orig_end')
  space <- space_btw_rects
  thickness <- rect_thickness
  i_r0<-1
  i_r1<-i_r0 - thickness
  for (i in 1:nrow(x))
  {
    ar_start<-x[i,]$start
    ar_end<-x[i,]$end
    ltype = 1
    if(x[i,]$genome==vector_chr)
    {
      if(x[i,]$strand=='-')
      {
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start
      }
      
      range_start<-x[i,]$orig_start
      range_new_end<-x[i,]$orig_start
      range_end<-x[i,]$orig_end
      for (jj in seq(1:(nrow(bed_color)-1))){
        if  (range_start>=bed_color[jj,]$V2 & range_start<=bed_color[jj+1,]$V2)
        {
          rgb_string<-strsplit(bed_color[jj,]$rgb," ",fixed=TRUE)[[1]]
          r<-rgb_string[1]
          g<-rgb_string[2]
          b<-rgb_string[3]
          i_col<-rgb(r,g,b,maxColorValue = 255)
          range_new_end<-bed_color[jj+1,]$V2-1
          span<-range_new_end-range_start
          if (range_new_end>=range_end)
          {
            
            range_new_end<-range_end
            span<-range_new_end-range_start
            kpRect(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
            ar_start<-ar_start+span
            break
          }
          else
          {
            
            if(x[i,]$strand=='+')
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start+span, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
              ar_start<-ar_start+span
            }
            
            else
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start-span, y0=i_r0 ,y1=i_r0-thickness, col=i_col, lty=ltype)
              ar_start<-ar_start-span
            }
            
          }
          
          range_start<-bed_color[jj+1,]$V2
          
        }
      }
    }
    else
    {
      if (x[i,]$strand=='+')
      {
        i_col = color_target_fwd
        ltype=1
      }
      else
      {
        i_col = color_target_rev
        ltype=2
      }
      kpRect(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
    }
    if (x[i,]$junction==TRUE)
    {
      if(x[i,]$target_genome_position_wrt_junction=='R')
        j<-x[i,]$end
      else
        j<-x[i,]$start
      
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0, i_r0-thickness),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-thickness
    i_r1<-i_r1-thickness
    i_r0<-i_r0-space
    i_r1<-i_r1-space
  }
  legend(x = "bottomleft", lty = c(1,2), legend=c("+", "-"), title=c('Target Strand'), cex=cex_legend)
  dev.off()
  
}




