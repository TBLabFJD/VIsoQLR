library(shiny)
# library(ggplot2)
library(vctrs)
library(plotly)
library(DT)
library(htmlwidgets)
# library(shinyFiles)
library(reticulate)
library(RColorBrewer)

# install.packages('reticulate')
# reticulate::install_miniconda()
# reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
# reticulate::use_miniconda('r-reticulate')

rm(list=ls()) 

options(shiny.maxRequestSize=30000*1024^2)
reticulate::py_run_string("import sys")


#############################################################################################
###   FUNCTIONS   ###########################################################################
#############################################################################################


##########################
# Data loading functions #
##########################

gff3_data_loading <- function(input_path){
  gff3 = read.delim(input_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  gff3 = gff3[gff3$V3 == "exon",c("V1", "V4","V5","V9", "V6")]
  colnames(gff3) = c("gene", "start", "end", "id", "score")
  gff3$score = as.numeric(gff3$score)
  gff3$id=gsub("^.*Name=", "", gff3$id, perl = TRUE)
  gff3$id=gsub(";.*$", "", gff3$id, perl = TRUE)
  return(gff3)
}


bed_data_loading <- function(input_path){
  bed6 = read.delim(input_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  colnames(bed6) = c("gene", "start", "end", "id", "score", "orientation")
  bed6$score = as.numeric(bed6$score)
  bed6$start = bed6$start + 1
  # bed6 = bed6[bed6$score == 60,]
  return(bed6)
}


transcript_data_loading <- function(input_path,  inputformat = "GTF"){
  if (inputformat == "GTF") {separation = " "}
  if (inputformat == "GFF3") {separation = "="}
  
  gtf = read.delim(input_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  
  transcripts = gtf[gtf$V3 == "transcript","V9"]
  if (length(transcripts) == 0){transcripts = gtf[,"V9"]}
  transcriptinfo = gsub("; ", ";", transcripts)
  transcriptinfo = do.call("rbind", strsplit(x = transcriptinfo, split = ";"))
  cnames = do.call("rbind", strsplit(x = transcriptinfo[1,], split = separation))[,1]
  transcriptinfo = gsub("^.* ", "", transcriptinfo, perl = T)
  colnames(transcriptinfo) = cnames
  transcriptinfo = transcriptinfo[!duplicated(transcriptinfo[,"transcript_id"]),]
    

  gtf = gtf[gtf$V3 == "exon", c("V1", "V4","V5","V9", "V6")]
  colnames(gtf) = c("gene", "start", "end", "id", "score")
  gtf$score = as.numeric(gtf$score)
  gtf$id=gsub("^.*transcript_id ", "", gtf$id, perl = TRUE)
  gtf$id=gsub(";.*$", "", gtf$id, perl = TRUE)
  
  gtf = merge(gtf, transcriptinfo, by.x = "id", by.y = "transcript_id")
  
  # Size annotation
  gtf_split = split(gtf, gtf$id)
  gtf_lengths = unlist(lapply(gtf_split, function(x) sum(x$end - x$start + 1)))
  gtf$size = paste0(gtf_lengths[gtf$id], "bp")
  
  gtf = type.convert(gtf, as.is = T)
  return(gtf)
}


#########################
# Break point functions #
#########################
break_point_calculation  <- function(raw_exons, num_reads_post_trimming, breakpoint_freq_threshold=3, very_close_bp=3, breakpoint_padding=5){
  
  # Break points frequency
  start_count = table(raw_exons$start)
  end_count = table(raw_exons$end)
  
  
  # Most frequent breakpoints
  # num_reads_post_trimming = length(unique(raw_exons$id))
  
  start_count = start_count[start_count>(num_reads_post_trimming*breakpoint_freq_threshold/100)]
  end_count = end_count[end_count>(num_reads_post_trimming*breakpoint_freq_threshold/100)]
  
  start_position = data.frame(breakp = as.numeric(names(start_count)), freq = as.numeric(start_count), stringsAsFactors = F)
  end_position = data.frame(breakp = as.numeric(names(end_count)), freq = as.numeric(end_count), stringsAsFactors = F)
  
  
  # very_close_bp
  # If two breakpoints from the same type (start/end) are closer than 3bp (default) and one of them is more than twice as frequent than the other, the less frequent is removed
  if (very_close_bp > 0){
    if (breakpoint_padding >= very_close_bp){
      
      start_pos_to_remove=c()
      for (i in 1:(nrow(start_position))){
        
        if (is.element(TRUE, start_position[i, "breakp"] + very_close_bp >= start_position[-i, "breakp"] & 
                       start_position[i, "breakp"] - very_close_bp <= start_position[-i, "breakp"] &
                       start_position[i, "freq"] < start_position[-i, "freq"])) {start_pos_to_remove = c(start_pos_to_remove, i)}
                       # start_position[i, "freq"] * 2 <= start_position[, "freq"])) {start_pos_to_remove = c(start_pos_to_remove, i)}
      }
      if (!is.null(start_pos_to_remove)) {start_position = start_position[-start_pos_to_remove,]}
      
      end_pos_to_remove=c()
      for (i in 1:(nrow(end_position))){
        if (is.element(TRUE, end_position[i, "breakp"] + very_close_bp >= end_position[-i, "breakp"] & 
                       end_position[i, "breakp"] - very_close_bp <= end_position[-i, "breakp"] &
                       end_position[i, "freq"] < end_position[-i, "freq"])) {end_pos_to_remove = c(end_pos_to_remove, i)}
                       # end_position[i, "freq"] * 2 <= end_position[, "freq"])) {end_pos_to_remove = c(end_pos_to_remove, i)}
      }
      if (!is.null(end_pos_to_remove)) {end_position = end_position[-end_pos_to_remove,]}
      
    } else{
      stop("-v, --very_close_bp cannot be bigger than -p, --padding")
    }
  }
  
  
  # Break point boundaries
  start_position$left = start_position$breakp - breakpoint_padding
  start_position$right = start_position$breakp + breakpoint_padding
  end_position$left = end_position$breakp - breakpoint_padding
  end_position$right = end_position$breakp + breakpoint_padding
  
  
  
  for (i in 1:(nrow(start_position)-1)){
    if (start_position[i, "right"] > start_position[i+1, "left"]){
      media_pos = (start_position[i, "right"] + start_position[i+1, "left"])/2
      if (media_pos%%1==0) { # Check if the number is decimal or interger
        start_position[i, "right"] = media_pos - 1
        start_position[i+1, "left"] = media_pos + 1
      }else{
        start_position[i, "right"] = floor(media_pos)
        start_position[i+1, "left"] = ceiling(media_pos)
      }
    }
  }
  
  
  
  for (i in 1:(nrow(end_position)-1)){
    if (end_position[i, "right"] > end_position[i+1, "left"]){
      media_pos = (end_position[i, "right"] + end_position[i+1, "left"])/2
      if (media_pos%%1==0) { # Check if the number is decimal or interger
        end_position[i, "right"] = media_pos - 1
        end_position[i+1, "left"] = media_pos + 1
      }else{
        end_position[i, "right"] = floor(media_pos)
        end_position[i+1, "left"] = ceiling(media_pos)
      }
    }
  }
  
  start_position = start_position[,c("breakp", "left", "right")]
  end_position = end_position[,c("breakp", "left", "right")]
  
  return(list(start_position, end_position))
}




####################################
# Known break points replace/merge #
####################################
known_sites_replace <- function(known_sites_file, gene){

  known_sites = read.delim(known_sites_file, header = F, stringsAsFactors = F, comment.char = "#")[,1:5]
  colnames(known_sites) = c("tipo", "gene", "breakp", "left", "right")
  known_sites = known_sites[known_sites$gene == gene,]
  known_sites$tipo = tolower(known_sites$tipo)
  
  start_position = known_sites[known_sites$tipo %in% c("start"), c("breakp", "left", "right")]
  end_position = known_sites[known_sites$tipo %in% c("end", "stop"), c("breakp", "left", "right")]
  

  return(list(start_position, end_position))
}



known_sites_merge <- function(known_sites_file, break_points_list, gene){
  start_position = break_points_list[[1]]
  end_position = break_points_list[[2]]
  
  known_sites = read.delim(known_sites_file, header = F, stringsAsFactors = F, comment.char = "#")[,1:5]
  colnames(known_sites) = c("tipo", "gene", "site", "lower", "upper")
  known_sites = known_sites[known_sites$gene == gene,]
  known_sites$tipo = tolower(known_sites$tipo)
  
  for (i in 1:nrow(known_sites)){
    if (known_sites[i, "tipo"] == "start"){
      start_position = start_position[!(start_position$breakp >= known_sites[i, "lower"] & start_position$breakp <= known_sites[i, "upper"]), ]
      start_position[start_position$breakp < known_sites[i, "lower"] & start_position$right >= known_sites[i, "lower"], "right"] = known_sites[i, "lower"] - 1
      start_position[start_position$breakp > known_sites[i, "lower"] & start_position$left <= known_sites[i, "upper"], "left"] = known_sites[i, "upper"] + 1
      start_position = rbind(start_position, c(known_sites[i, "site"],
                                               known_sites[i, "lower"],
                                               known_sites[i, "upper"]))
    } else if (known_sites[i, "tipo"] %in% c("end", "stop")) {
      end_position = end_position[!(end_position$breakp >= known_sites[i, "lower"] & end_position$breakp <= known_sites[i, "upper"]), ]
      end_position[end_position$breakp < known_sites[i, "lower"] & end_position$right >= known_sites[i, "lower"], "right"] = known_sites[i, "lower"] - 1
      end_position[end_position$breakp > known_sites[i, "lower"] & end_position$left <= known_sites[i, "upper"], "left"] = known_sites[i, "upper"] + 1
      end_position = rbind(end_position, c(known_sites[i, "site"],
                                           known_sites[i, "lower"],
                                           known_sites[i, "upper"]))
    }
  }
  
  return(list(start_position, end_position))
}




#########################
# Break point asignment #
#########################
break_point_asignment <- function(raw_exons, break_points_list){
  start_position = break_points_list[[1]]
  end_position = break_points_list[[2]]
  
  bp_asigned_exons <- raw_exons
  
  bp_asigned_exons$start_tag = NA
  for (i in 1:nrow(start_position)){
    bp_asigned_exons[bp_asigned_exons$start >= start_position[i, "left"]  & bp_asigned_exons$start <= start_position[i, "right"], "start_tag"] = start_position[i, "breakp"]
  }
  
  bp_asigned_exons$end_tag = NA
  for (i in 1:nrow(end_position)){
    bp_asigned_exons[bp_asigned_exons$end >= end_position[i, "left"]  & bp_asigned_exons$end <= end_position[i, "right"], "end_tag"] = end_position[i, "breakp"]
  }
 
  return(bp_asigned_exons)
}



###########################
# Break point information #
###########################
# break_points_list = rv$break_points; num_reads_post_trimming = rv$num_reads_post_trimming; bp_asigned_exons = rv$bp_asigned_exons


break_point_info_calculator <- function(break_points_list, num_reads_post_trimming, bp_asigned_exons){
  start_position_in = break_points_list[[1]]
  end_position_in = break_points_list[[2]]

  start_position_in$type = "start"
  end_position_in$type = "end"

  # print("1")
  start_position = merge(start_position_in, data.frame(table(bp_asigned_exons$start)), by.x = "breakp", by.y = "Var1", all.x = T, sort = F)
  end_position = merge(end_position_in, data.frame(table(bp_asigned_exons$end)), by.x = "breakp", by.y = "Var1", all.x = T, sort = F)

  # print("2")
  start_position = merge(start_position, data.frame(table(bp_asigned_exons$start_tag)), by.x = "breakp", by.y = "Var1", all.x = T, sort = F)
  end_position = merge(end_position, data.frame(table(bp_asigned_exons$end_tag)), by.x = "breakp", by.y = "Var1", all.x = T, sort = F)
  
  # print("3")
  start_end_position = rbind(start_position, end_position)
  start_end_position$relative_freq = round(start_end_position$Freq.x/num_reads_post_trimming*100, 2)
  start_end_position$relative_freq_padding = round(start_end_position$Freq.y/num_reads_post_trimming*100, 2)

  # print("4")
  start_end_position$increase_percentage = round((start_end_position$Freq.y - start_end_position$Freq.x)/start_end_position$Freq.x*100, 2)
  start_end_position$gene = bp_asigned_exons$gene[1]
  start_end_position = start_end_position[,c("type", "gene", "breakp", "left", "right",  "Freq.x", "relative_freq", "Freq.y", "relative_freq_padding", "increase_percentage")]
  colnames(start_end_position) = c("#type", "gene", "breakpoint", "lower_limit", "upper_limit", "num_of_reads_breakpoint", "percent_of_reads_breakpoint", "number_of_reads_interval", "percent_of_reads_interval", "num_of_reads_increase_percentage_interval")
  
  # print("5")
  start_end_position2transform = start_end_position
  start_end_position2transform$`num_of_reads (%)` = paste0(start_end_position$num_of_reads_breakpoint, " (", start_end_position$percent_of_reads_breakpoint,"%)")
  start_end_position2transform$`total_num_of_reads (%)` = paste0(start_end_position$number_of_reads_interval, " (", start_end_position$percent_of_reads_interval,"%)")

  start_output = start_end_position2transform[start_end_position2transform$`#type` == "start", c("breakpoint", "lower_limit", "upper_limit", "num_of_reads (%)", "total_num_of_reads (%)")]
  end_output = start_end_position2transform[start_end_position2transform$`#type` == "end", c("breakpoint", "lower_limit", "upper_limit","num_of_reads (%)",  "total_num_of_reads (%)")]
  # print(start_output)
  # print(start_position)
  return(list(start_output, end_output, start_end_position))
}




###################
# Exon difinition #
###################
exon_definition <- function(bp_asigned_exons){
  # Read filtering if start/end do not match with with the estimated break points
  exon_id_filter1 = bp_asigned_exons[rowSums(is.na(bp_asigned_exons)) > 0, "id"]
  exon_id_filter2 = bp_asigned_exons[bp_asigned_exons$start_tag > bp_asigned_exons$end_tag, "id"]
  dup1 = duplicated(paste(bp_asigned_exons$id, "start", bp_asigned_exons$start_tag, sep = "_"))
  dup2 = duplicated(paste(bp_asigned_exons$id, "end", bp_asigned_exons$end_tag, sep = "_"))
  exon_id_filter3 = bp_asigned_exons[dup1 | dup2, "id"]
  
  exon_id_filter = unique(c(exon_id_filter1, exon_id_filter2, exon_id_filter3))
  defined_exons = bp_asigned_exons[!bp_asigned_exons$id %in% exon_id_filter,]
  
  # max_with = max(floor(log10(c(defined_exons$start_tag, defined_exons$end_tag))) + 1)
  # start_tag_formated = formatC(defined_exons$start_tag, width = max_with, flag = "0", digits = max_with)
  # end_tag_formated = formatC(defined_exons$end_tag, width = max_with, flag = "0", digits = max_with)
  # defined_exons$coordinates = paste(start_tag_formated, end_tag_formated, sep = "-")
  defined_exons$coordinates = paste(defined_exons$start_tag, defined_exons$end_tag, sep = "-")
  
  return(list(defined_exons,exon_id_filter))
}




# Exon information
exon_info_fun <- function(defined_exons){
  exon_info_df = as.data.frame(table(defined_exons$coordinates))
  colnames(exon_info_df) = c("exon", "number_of_reads")
  exon_info_df$percent_of_reads = round(exon_info_df$number_of_reads / length(unique(defined_exons$id)) * 100, 1)
  exon_info_df = exon_info_df[order(exon_info_df$exon),]
  
  exon_info_df$Size = unlist(lapply(gsub("_", "+", exon_info_df[,"exon"]), function(x) -eval(parse(text = x))+1))
  exon_info_df = exon_info_df[,c("exon", "Size", "number_of_reads", "percent_of_reads")]
  
  return(exon_info_df)
}





######################
# Isoform definition #
######################
isoform_definition <- function(defined_exons){
  # Isoform definition
  read_list = split(defined_exons, defined_exons$id) ### *** Tarda un poco (no es inmediato)
  read_isoform = unlist(lapply(read_list, function(x) paste(sort(x$coordinates), collapse = "_")))
  
  # frequency calculation
  if(length(unique(read_isoform)) == 1){
    isoform_frequencies = data.frame(read_isoform = unique(unique(read_isoform)), Freq = length(read_isoform))
    
  }else{
    isoform_frequencies = data.frame(sort(table(read_isoform), decreasing = TRUE))
  }
  isoform_frequencies$read_isoform = as.character(isoform_frequencies$read_isoform)
  isoform_frequencies$perc = round(isoform_frequencies$Freq * 100 / sum(isoform_frequencies$Freq), 1)

  return(list(isoform_frequencies, read_isoform))
}


# Isoform full-length filter
isoform_full_length_filter = function(isoforms, break_points){
  isoform_info = isoforms[[1]]
  isoform_read = isoforms[[2]]
  start = min(break_points[[1]][,1])
  end = max(break_points[[2]][,1])
  
  positions_split = strsplit(isoform_info$read_isoform, split = "-|_", perl = T)
  keep = c()
  discard = c()
  
  for (i in 1:nrow(isoform_info)){
    if (start %in% positions_split[[i]] & end %in% positions_split[[i]]) {
      keep = c(keep, i)
    } else { 
      discard = c(discard, i)
    }
  }
  
  keep_isoform_info = isoform_info[keep,]
  keep_isoform_info$perc = round(keep_isoform_info$Freq/sum(keep_isoform_info$Freq)*100, digits = 2)
  keep_isoform_read = isoform_read[isoform_read %in% keep_isoform_info$read_isoform]
  partial_reads = isoform_read[isoform_read %in% isoform_info[discard,"read_isoform"]]
  return(list(keep_isoform_info, keep_isoform_read, partial_reads))
}



# Isoform information
isoform_info_fun <- function(isoform_frequencies, n_inicial, n_post_trim, n_final, n_not_full_lenght = NULL){
  
  if (is.null(n_not_full_lenght)){
    n_vector_reads = n_inicial - n_post_trim
    n_no_consensous_breakpoint_reads = n_post_trim - n_final
    
    other_groups_df = data.frame(isoform_id = c("Only vector reads", "No consensous breakpoint reads"),
                                 read_isoform = c("-", "-"), 
                                 Freq = c(n_vector_reads, n_no_consensous_breakpoint_reads),
                                 perc = c("-", "-"), 
                                 Size = c(0, 0),
                                 stringsAsFactors = FALSE)
  } else {
    n_vector_reads = n_inicial - n_post_trim
    n_no_consensous_breakpoint_reads = n_post_trim - n_final - n_not_full_lenght
    
    other_groups_df = data.frame(isoform_id = c("Only vector reads", "No consensous breakpoint reads", "Partial length reads"),
                                 read_isoform = c("-", "-", "-"), 
                                 Freq = c(n_vector_reads, n_no_consensous_breakpoint_reads, n_not_full_lenght),
                                 perc = c("-", "-", "-"), 
                                 Size = c(0, 0, 0),
                                 stringsAsFactors = FALSE)
  } 
  
  isoform_id = paste("Iso", 1:nrow(isoform_frequencies), sep = "")
  isoform_frequencies = cbind(isoform_id, isoform_frequencies)
  isoform_frequencies$Size = unlist(lapply(gsub("_", "-1+", isoform_frequencies[,2]), function(x) -eval(parse(text = x))+1))
  
  all_groups = rbind(other_groups_df, isoform_frequencies)
  all_groups$perc_total = round(all_groups$Freq*100/sum(all_groups$Freq), 2)
  colnames(all_groups) =  c("Isoform_id", "Isoform", "N_of_reads", "Prerc_partial", "Size", "Perc_total")
  
  all_groups = all_groups[,c("Isoform_id", "Isoform", "Size", "N_of_reads", "Prerc_partial", "Perc_total")]
  return(all_groups)
}






# Read classification
read_clasification_fun <- function(all_ids, no_consensus_ids, classified_ids, not_full_length = NA, iso_dict){
  classified_df <- data.frame(read_id = names(classified_ids), type = as.character(classified_ids))
  iso_dict = iso_dict[grep("^Iso[0-9]+$",iso_dict$Isoform_id), c("Isoform_id", "Isoform")]
  rownames(iso_dict) = iso_dict$Isoform
  classified_df$type = iso_dict[classified_df$type,"Isoform_id"]
  
  only_vector_ids = setdiff(unique(all_ids), c(unique(no_consensus_ids), names(classified_ids), names(not_full_length)))
  if (length(only_vector_ids) != 0){
    trimmed_df <- data.frame(read_id = only_vector_ids, type = "Only_vector_reads")
    classified_df = rbind(classified_df, trimmed_df)
  }
  
  if (length(unique(no_consensus_ids)) != 0){
    no_consensus_df <- data.frame(read_id = unique(no_consensus_ids), type = "No_consensous_breakpoint_reads")
    classified_df = rbind(classified_df, no_consensus_df)
  }
  
  if (!is.na(not_full_length[1])){
    not_full_length_df <- data.frame(read_id = names(not_full_length), type = "Partial_length_reads")
    classified_df = rbind(classified_df, not_full_length_df)
  }
  
  return(classified_df)
}






################
# Data to plot #
################
isoform_df_plot_fun <- function(isoform_frequencies){
  isoforms = as.character(isoform_frequencies$read_isoform)
  isoform_num = 1
  x_pos_start = c()
  x_pos_end = c()
  tipo = c()
  exon = c()
  y_pos = c()
  isoform_id_all = c()
  perc = c()
  star_stop  =c()
  
  for (i in 1:nrow(isoform_frequencies)) {
    
    #isoform_id = paste("Iso", isoform_num, ": ", isoform_frequencies[i,2] ," reads (", isoform_frequencies[i,3], "%)", sep = "")
    # isoform_id = paste("Iso", isoform_num, " (", isoform_frequencies[i,3], "%)", sep = "")
    isosize = paste0(-eval(parse(text = gsub("_", "-1+", isoform_frequencies[i,1])))+1, "bp")
    isoform_id = paste("Iso", isoform_num, " | ", isosize, " | ", isoform_frequencies[i,3], "%", sep = "")
    isoform_id_all = c(isoform_id_all, isoform_id)
    isoform_num = isoform_num + 1
    
    positions = strsplit(isoform_frequencies[i,1], split = "_")[[1]]

    for (j in positions){
      
      start = as.numeric(gsub("\\*", "", strsplit(j, split = "-")[[1]][1]))
      end = as.numeric(gsub("\\*", "", strsplit(j, split = "-")[[1]][2]))
      
      x_pos_start = c(x_pos_start, start)
      x_pos_end = c(x_pos_end, end)
      y_pos = c(y_pos, isoform_id)
      perc = c(perc, isoform_frequencies[i,3])
      star_stop = c(star_stop, j)
      
    }
  }
  
  x_pos_start = as.integer(x_pos_start)
  x_pos_end = as.integer(x_pos_end)
  df_pos = data.frame(x_pos_start, x_pos_end, y_pos, perc, star_stop, stringsAsFactors = F)
  
  df_pos$y_pos = ordered(df_pos$y_pos, rev(isoform_id_all))
  df_pos$grupo = "VisoQLR detected isoforms"
  
  return(df_pos)
}



# Break point superplot
break_point_freq_fun <- function(raw_exons){
  
  num_reads_initial = length(unique(raw_exons$id))
  
  start_df = data.frame(table(raw_exons[, "start"]))
  start_df$tipo = "start"
  end_df = data.frame(table(raw_exons[, "end"]))
  end_df$tipo = "end"
  
  df_superplot = rbind(start_df, end_df)
  
  df_superplot$Var1 = as.numeric(as.character(df_superplot$Var1))
  df_superplot$Freq = as.numeric(as.character(df_superplot$Freq))
  df_superplot$perc = df_superplot$Freq/num_reads_initial*100

  return(df_superplot)
}








##################
# Plot functions #
##################

# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }


isoform_plot_fun <- function(df_pos, ylabel, pal = NULL, start_end_levels, showlegend){
  df_pos$star_stop = factor(df_pos$star_stop, levels = start_end_levels)
  
  df_pos %>%
    plot_ly(colors = pal) %>%
    add_lines(
      x = ~x_pos_start, y = ~y_pos,
      showlegend = FALSE, split = ~y_pos, color = I("blue"), size = I(0.8)
    )%>%
    add_segments(
      x = ~x_pos_start, y = ~y_pos,
      xend = ~x_pos_end, yend = ~y_pos, 
      split = ~star_stop, size = I(10), color = ~star_stop, legendgroup = ~star_stop, showlegend = showlegend
    )%>% layout(xaxis = list(title = "Sequence position", showgrid = FALSE),
                yaxis = list(title = ylabel, fixedrange=TRUE, showgrid = FALSE),
                legend = list(title = list(text="Exon coordinates")))
}












break_point_freq_plot_fun <- function(df_superplot, break_points, start_dt, end_dt){
  
  # Assign colors
  df_superplot$current_color[df_superplot$tipo == "start"] = "#2166AC"
  df_superplot$current_color[df_superplot$tipo == "end"] = "#B2182B"
  
  # Text box
  start_bp = break_points[[1]]
  end_bp = break_points[[2]]
  start_bp$tipo = "start"
  end_bp$tipo = "end"
  # start_bp$current_color = "#2166AC"
  # end_bp$current_color = "#B2182B"
  bp_join = rbind(start_bp, end_bp)
  
  # # Change color of the selected rows
  # if (length(start_dt)) {
  #   selected_starts = break_points[[1]][start_dt,]
  #   df_superplot$current_color[df_superplot$tipo == "start" & df_superplot$Var1 %in% selected_starts$breakp] = "green"
  #   # bp_join$current_color[bp_join$tipo == "start" & bp_join$breakp %in% selected_starts$breakp] = "green"
  # }
  # 
  # if (length(end_dt)) {
  #   selected_ends = break_points[[2]][end_dt,]
  #   df_superplot$current_color[df_superplot$tipo == "end" & df_superplot$Var1 %in% selected_ends$breakp] = "orange"
  #   # bp_join$current_color[bp_join$tipo == "end" & bp_join$breakp %in% selected_starts$breakp] = "orange"
  # }
  
  colnames(bp_join)[1] = "Var1"
  bp_join = merge(bp_join, df_superplot, by = c("Var1", "tipo"))
  
  plot_ly(
    df_superplot,
    x = ~Var1,
    y = ~perc,
    marker = list(color = ~current_color),
    # color = ~tipo,
    # colors = c("#B2182B", "#2166AC"), #gg_color_hue(2),
    type = "bar",
    showlegend = FALSE
  ) %>% layout(xaxis = list(title = "Sequence position"),
               yaxis  = list(title = "Percentage of reads (%)", fixedrange=TRUE),
               legend = list(title = list(text="Exon coordinates", showlegend = FALSE)),
               barmode = "stack"
  )%>% add_markers(text = ~paste(tipo, "<br />", Var1, "(", round(perc,1), "%)"), hoverinfo = "text", 
                   marker = list(color = ~current_color), data = bp_join)

}

# 
# plot_ly(
#   df_superplot,
#   x = ~Var1,
#   y = ~perc,
#   marker = list(color = ~current_color),
#   # color = ~current_color,
#   type = "bar"
# ) %>% layout(xaxis = list(title = "Sequence position"),
#              yaxis  = list(title = "Percentage of reads (%)", fixedrange=TRUE),
#              legend = list(title = list(text="Exon coordinates")),
#              barmode = "stack") %>%
#   add_polygons(x = c(600,600,1200, 1200),
#                y = c(0, 100, 100, 0), showlegend = F,
#                color = I("#7F98C6"), opacity = 0.5)



multiplot_fun1 <- function(isoform_plot, break_point_freq_plot, num_iso, iso_pixels, histo_pixels){
  subplot(isoform_plot, break_point_freq_plot, 
          nrows = 2, margin = 0.03, 
          shareX = TRUE, titleX = T, titleY = T,
          heights = c(num_iso*iso_pixels/(num_iso*iso_pixels+histo_pixels), histo_pixels/(num_iso*iso_pixels+histo_pixels)))
}

multiplot_fun2 <- function(isoform_plot1, isoform_plot2, break_point_freq_plot, num_iso1, num_iso2, iso_pixels, histo_pixels){
  subplot(isoform_plot1, isoform_plot2, break_point_freq_plot, 
          nrows = 3, margin = 0.03, 
          shareX = TRUE, titleX = T, titleY = T,
          heights = c(num_iso1*iso_pixels/(num_iso1*iso_pixels+num_iso2*iso_pixels+histo_pixels), num_iso2*iso_pixels/(num_iso1*iso_pixels+num_iso2*iso_pixels+histo_pixels),histo_pixels/(num_iso1*iso_pixels+num_iso2*iso_pixels+histo_pixels)))
}


# 
# library(seqinr)
# 
# seq_path="/home/gonzalo/tblab/home/gonzalo/pax6/references/pSPL3_PAX6_Ex5-7_v2.fa"
# seq = read.fasta(seq_path)
# 
# aa = data.frame(y="Sequence", x=1:length(as.character(seq[["Secuencia"]])), seq=as.character(seq[["Secuencia"]]))
# p55  = ggplot(aa,aes(x=x,y=y)) +
#   geom_text(aes(label = seq)) +
#   theme_minimal() +
#   theme(panel.grid = element_blank()) +
#   theme(axis.title.y=element_blank(),
#         axis.title.x=element_blank(),
#         # axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) +
#   theme(plot.background = element_rect(fill = "white", colour = "white"))
# 
# 
# 
# p55 <- plot_ly(aa, x = ~x, y = ~y, text = ~seq)
# p55 <- p55 %>% add_text()
# p55 <- p55 %>% layout(margin = list(b = 0,  t = 0,  pad = 0))
# 
# # p44 <- plot_ly(
# #   df_superplot,
# #   x = ~Var1,
# #   y = ~perc,
# #   color = ~tipo,
# #   colors = gg_color_hue(2),
# #   type = "bar"
# # ) %>% layout(xaxis = list(title = "Break point position"),
# #              yaxis  = list(title = "Percentage of reads", fixedrange=TRUE),
# #              legend = list(title = list(text="Exon coordinates")),
# #              barmode = "stack")%>% add_annotations(text = aa$seq,
# #                                                    x = aa$x,
# #                                                    y = -1,
# #                                                    showarrow = FALSE)
# 



shinyInput <- function(FUN, n, id, ses, ...) {
  as.character(FUN(paste0(id, n), ...))
}

getRemoveButton <- function(n, idS = "", lab = "Pit") {
  # if (n == 0) { return() }
  if (stringr::str_length(idS) > 0) idS <- paste0(idS, "-")
  ret <- shinyInput(actionButton, n,
                    paste0(runif(1),'button_'), label = "Remove", style = "color: red;background-color: white",
                    onclick = sprintf('Shiny.onInputChange(\"%sremove_button_%s\",  this.id)' ,idS, lab))
  return (ret)
}


# getSelectButton <- function(n, idS = "", lab = "Pit") {
#   # if (n == 0) { return() }
#   if (stringr::str_length(idS) > 0) idS <- paste0(idS, "-")
#   ret <- shinyInput(actionButton, n,
#                     paste0(runif(1),'button_'), label = "Remove", style = "color: red;background-color: white",
#                     onclick = sprintf('Shiny.onInputChange(\"%sremove_button_%s\",  this.id)' ,idS, lab))
#   return (ret)
# }



isoinfo2gtf <- function(isoforms_information, exon_information, selected_gene){
  isoforms_information = isoforms_information[grepl("^Iso[0-9]+$", isoforms_information$Isoform_id, perl = T),]
  rownames(exon_information) = exon_information$exon
  
  isolist = list()
  listpos = 1
  for (i in 1:nrow(isoforms_information)){
    iso_row = isoforms_information[i,]
    exoncoord = do.call("rbind", strsplit(strsplit(iso_row$Isoform, split = "_")[[1]], split = "-"))
    colnames(exoncoord) = c("start", "end")
    exoncoord = data.frame(type.convert(exoncoord, as.is = T))
    exoncoord = exoncoord[order(exoncoord$start),]
    trascript_row = c(selected_gene, "VIsoQLR", "transcript", min(exoncoord$start), max(exoncoord$end), ".", "+", ".", 
                      paste0("gene_id ", selected_gene, "; ",
                             "transcript_id ", iso_row$Isoform_id, "; ",
                             "size ", iso_row$Size,"; ",
                             "n_of_reads ", iso_row$N_of_reads,"; ",
                             "prerc_partial ", iso_row$Prerc_partial,"%; ",
                             "perc_total ", iso_row$Perc_total, "%"
                      ))
    
    isolist[[listpos]] <- trascript_row
    listpos = listpos + 1
    for (j in 1:nrow(exoncoord)){
      exon_id = paste(exoncoord[j,"start"], exoncoord[j,"end"], sep = "-")
      exon_row = c(selected_gene, "VIsoQLR", "exon", exoncoord[j,"start"], exoncoord[j,"end"], ".", "+", ".", 
                   paste0("gene_id ", selected_gene, "; ",
                          "transcript_id ", iso_row$Isoform_id, "; ",
                          "exon_number ", j, "; ",
                          "size ", exon_information[exon_id, "Size"], "; ",
                          "number_of_reads ", exon_information[exon_id, "number_of_reads"], "; ",
                          "percent_of_reads ", exon_information[exon_id, "percent_of_reads"], "%"
                   ))
      isolist[[listpos]] <- exon_row
      listpos = listpos + 1
    }
  }
  
  isotable = do.call("rbind", isolist)
  return(isotable)
}

























#############################################################################################
###   USER INTERFACE   ######################################################################
#############################################################################################

ui <- fluidPage(
  navbarPage(
    title = "VIsoQLR",
    selected = "Isoform analysis",
    
    tabPanel("Mapping",
      
       mainPanel(
       width = 3, class = "well",
       h3("GMAP reference index building"),
       fileInput("fasta2index", "Reference sequence(s) in FASTA format"),
       uiOutput("wait_gmap_index_download")
      ),
      
      mainPanel(
        width = 3, class = "well",
        h3("GMAP mapping"),
        radioButtons("output_format_gmap", "Output format", c("GFF3", "BED6", "BAM"), selected="GFF3"),
        fileInput("reference_gmap", "Reference index in ZIP format (file generated in the first column)"),
        fileInput("rawreads_gmap", "Raw reads in FASTQ format"),
        numericInput("gmap_threads", "Number of threads", min = 1, value = 4),
        actionButton("run_gmap_alignment", "Run mapping"),
        h4("Download aligned reads"),
        uiOutput("wait_gmap_mapping_download"),
        h4("Download bam index (.bai)"),
        uiOutput("wait_gmap_mapping_download_bai")
      ),
      
      mainPanel(
        width = 3, class = "well",
        h3("Minimap2 mapping"),
        radioButtons("output_format_minimap", "Output format", c("BED6", "BAM"), selected="BED6"),
        fileInput("reference_minimap", "Reference sequence(s) in FASTA format"),
        fileInput("rawreads_minimap", "Raw reads in FASTQ format"),
        numericInput("minimap_threads", "Number of threads", min = 1, value = 4),
        actionButton("run_minimap_alignment", "Run mapping"),
        h4("Download aligned reads"),
        uiOutput("wait_minimap_mapping_download"),
        h4("Download bam index (.bai)"),
        uiOutput("wait_minimap_mapping_download_bai")
      ) 
    ),
    
    
    tabPanel("Isoform analysis",
      sidebarLayout(
        sidebarPanel(
          width=2,
          
          wellPanel(
            h3("Input"),
            radioButtons("inputtype", "Input format", c("GFF3", "BED6", "BAM"), selected="GFF3"),
            fileInput("inputfile", "Input file")
            # fileInput("reference", "Reference sequence")
          ),
          uiOutput("input_wait_sidebarpanel"),
        ),
        
        mainPanel(
          width = 10,
          uiOutput("input_wait_mainpanel"),
          ))
      )
))









#############################################################################################
###   SERVER   ##############################################################################
#############################################################################################

server <- function(input, output) {
  
  rv <- reactiveValues(data = NULL, orig=NULL)
  
  #=========#
  # Mapping #
  #=========#


  # GMAP index building
  observeEvent(ignoreInit=T, c( 
    input$fasta2index
  ),{ 
    
    rv$index2download <- NULL
    
    file <- input$fasta2index
    ext <- tools::file_ext(file$datapath)
    
    print("GMAP index building (prerequisite)")
    req(file)
    
    withProgress(message = 'GMAP index building', value = 1/2, {

      validate(need(ext %in% c("fasta", "fa"), "Please upload a FASTA file"))
      
      prefix=gsub("(.fasta|.fa)$", "", file$name, perl = T)
      system(paste0("gmap_build -D ./ -d ", prefix, " " , file$datapath))
      system(paste0("zip ", prefix, ".zip ./", prefix, "/*"))
      rv$index2download = paste0(prefix, ".zip")

    })
    print("GMAP index building")
  })
  
  
  output$wait_gmap_index_download<-renderUI({
    req(rv$index2download)
    tagList(
      downloadButton("download_gmap_index", "Download GMAP index (.zip)")
    )
  })
  
  
  output$download_gmap_index <- downloadHandler(
    filename = function(){rv$index2download},
    content = function(file){file.copy(paste0("./", rv$index2download), file)}
  )
  
  
  
  
  # GMAP mapping
  
  observeEvent(ignoreInit=T, c( 
    input$run_gmap_alignment
  ),{ 
    
    if (is.null(input$gmap_threads)) gmap_threads = 4 else gmap_threads = input$gmap_threads
    
    rv$gmap2download <- NULL
    rv$gmap2download_bai <- NULL
    
    index_file <- input$reference_gmap
    index_ext <- tools::file_ext(index_file$datapath)
    
    reads_file <- input$rawreads_gmap
    reads_ext <- tools::file_ext(reads_file$datapath)
    
    print("GMAP mapping (prerequisite)")
    req(index_file)
    req(reads_file)
    req(input$output_format_gmap)
    
    withProgress(message = 'GMAP mapping', value = 1/2, {

      validate(need(index_ext %in% c("zip"), "Please upload a ZIP file"))
      validate(need(reads_ext %in% c("fastq", "fq", "gz", "fastq.gz", "fq.gz"), "Please upload a FASTQ file"))
      
      system(paste0("unzip -o " , index_file$datapath))
      index_prefix = system("ls -rt | tail -n 1", intern = T)

      if (grepl(".gz$", reads_file$name, perl = T)) { cat_type = "zcat " } else { cat_type = "cat " }
      sample_prefix = gsub("(.fastq|.fq|.fastq.gz|.fq.gz)$", "", reads_file$name, perl = T)
      
      if (input$output_format_gmap == "GFF3"){
        rv$gmap2download = paste0(sample_prefix, ".gff3")
        print(paste0(cat_type, reads_file$datapath, " | gmap -n1 -t ", gmap_threads, " --cross-species --gff3-add-separators=0 -f 2 -z auto -D ./ -d ", index_prefix, " > ", rv$gmap2download, " 2> log.err"))
        system(paste0(cat_type, reads_file$datapath, " | gmap -n1 -t ", gmap_threads, " --cross-species --gff3-add-separators=0 -f 2 -z auto -D ./ -d ", index_prefix, " > ", rv$gmap2download, " 2> log.err"))
      
      } else if (input$output_format_gmap == "BAM"){
        rv$gmap2download = paste0(sample_prefix, ".bam")
        print(paste0(cat_type, reads_file$datapath, " | gmap -n1 -t ", gmap_threads, " --cross-species --gff3-add-separators=0 -f samse -z auto -D ./ -d ", index_prefix, "  2> log.err | samtools view -Su | samtools sort > ", rv$gmap2download))
        system(paste0(cat_type, reads_file$datapath, " | gmap -n1 -t ", gmap_threads, " --cross-species --gff3-add-separators=0 -f samse -z auto -D ./ -d ", index_prefix, "  2> log.err | samtools view -Su | samtools sort > ", rv$gmap2download))

        rv$gmap2download_bai = paste0(sample_prefix, ".bam.bai")
        print(paste0("samtools index ", rv$gmap2download))
        system(paste0("samtools index ", rv$gmap2download))
        
      } else if (input$output_format_gmap == "BED6"){
        rv$gmap2download = paste0(sample_prefix, ".bed6")
        print(paste0(cat_type, reads_file$datapath, " | gmap -n1 -t ", gmap_threads, " --cross-species --gff3-add-separators=0 -f samse -z auto -D ./ -d ", index_prefix, "  2> log.err | samtools view -Su | samtools sort | bamToBed -split -i stdin > ", rv$gmap2download))
        system(paste0(cat_type, reads_file$datapath, " | gmap -n1 -t ", gmap_threads, " --cross-species --gff3-add-separators=0 -f samse -z auto -D ./ -d ", index_prefix, "  2> log.err | samtools view -Su | samtools sort | bamToBed -split -i stdin > ", rv$gmap2download))
        
      }
        
    })
    print("GMAP mapping")
  })
  
 
  output$wait_gmap_mapping_download<-renderUI({
    req(rv$gmap2download)
    tagList(
      downloadButton("download_gmap_reads", label = "Download aligned reads")
    )
  })
 
  
  output$download_gmap_reads <- downloadHandler(
    filename = function(){rv$gmap2download},
    content = function(file){file.copy(paste0("./", rv$gmap2download), file)}
  )
  
  
  
  output$wait_gmap_mapping_download_bai<-renderUI({
    req(rv$gmap2download_bai)
    tagList(
      downloadButton("downloadgmap_reads_bai", label = "Download bam index (.bai)")
    )
    
    output$downloadgmap_reads_bai <- downloadHandler(
      filename = function(){rv$gmap2download_bai},
      content = function(file){file.copy(paste0("./", rv$gmap2download_bai), file)}
    )
  })
  
  
  
  
  
  
  
  # Minimap2 mapping
  
  observeEvent(ignoreInit=T, c( 
    input$run_minimap_alignment
  ),{ 
    
    if (is.null(input$minimap_threads)) minimap_threads = 4 else minimap_threads = input$minimap_threads
    
    rv$minimap2download <- NULL
    rv$minimap2download_bai <- NULL
    
    index_file <- input$reference_minimap
    index_ext <- tools::file_ext(index_file$datapath)
    
    reads_file <- input$rawreads_minimap
    reads_ext <- tools::file_ext(reads_file$datapath)
    
    print("minimap mapping (prerequisite)")
    req(index_file)
    req(reads_file)
    req(input$output_format_minimap)
    
    withProgress(message = 'Minimap2 mapping', value = 1/2, {

      validate(need(index_ext %in% c("fa", "fasta"), "Please upload a ZIP file"))
      validate(need(reads_ext %in% c("fastq", "fq", "gz", "fastq.gz", "fq.gz"), "Please upload a FASTQ file"))
      
      
      if (grepl(".gz$", reads_file$name, perl = T)) { cat_type = "zcat " } else { cat_type = "cat " }
      sample_prefix = gsub("(.fastq|.fq|.fastq.gz|.fq.gz)$", "", reads_file$name, perl = T)
      
      if (input$output_format_minimap == "BAM"){
        rv$minimap2download = paste0(sample_prefix, ".bam")
        print(paste0("minimap2 -ax splice --secondary=no -t ", minimap_threads, " ", index_file$datapath, " ", reads_file$datapath, " > ", rv$minimap2download ))
        system(paste0("minimap2 -ax splice --secondary=no -t ", minimap_threads, " ", index_file$datapath, " ", reads_file$datapath, " > ", rv$minimap2download ))
        
        rv$minimap2download_bai = paste0(sample_prefix, ".bam.bai")
        print(paste0("samtools index ", rv$minimap2download))
        system(paste0("samtools index ", rv$minimap2download))
        
      } else if (input$output_format_minimap == "BED6"){
        rv$minimap2download = paste0(sample_prefix, ".bed6")
        print(paste0("minimap2 -ax splice --secondary=no -t ", minimap_threads, " ", index_file$datapath, " ", reads_file$datapath, " | bamToBed -split -i stdin > ", rv$minimap2download ))
        system(paste0("minimap2 -ax splice --secondary=no -t ", minimap_threads, " ", index_file$datapath, " ", reads_file$datapath, " | bamToBed -split -i stdin > ", rv$minimap2download ))
      }
      
    })
    print("Minimap2 mapping")
  })
  
  output$wait_minimap_mapping_download<-renderUI({
    req(rv$minimap2download)
    tagList(
      downloadButton("downloadminimap_reads", label = "Download aligned reads")
    )
    
    output$downloadminimap_reads <- downloadHandler(
      filename = function(){rv$minimap2download},
      content = function(file){file.copy(paste0("./", rv$minimap2download), file)}
    )
  })
  
  output$wait_minimap_mapping_download_bai<-renderUI({
    req(rv$minimap2download_bai)
    tagList(
      downloadButton("downloadminimap_reads_bai", label = "Download bam index (.bai)")
    )
    
    output$downloadminimap_reads_bai <- downloadHandler(
      filename = function(){rv$minimap2download_bai},
      content = function(file){file.copy(paste0("./", rv$minimap2download_bai), file)}
    )
  })
  
  
  
  
  
  
  
  #===================#
  # Isoform detection #  
  #===================#

    output$input_wait_sidebarpanel<-renderUI({
    req(input$inputfile)
    tagList(
      
      wellPanel(
        h3("Analysis bounding"),
        selectInput(inputId='selectgene', label='Select gene', choices = unique(rv$orig$gene)),
        uiOutput("studied_range"),
        checkboxInput("full_reads", "Select only complete PCR sequences")
      ),
      
      wellPanel(
        h3("Exon coordinates"),
        wellPanel(
          h4("Automatic detection"),
          numericInput("peak_threshold", "Read threshold (%)", min = 0, max = 100, value = 2),
          numericInput("padding", "Padding (# of bases)", min = 0, step = 1, value = 5),
          numericInput("very_close_bp", "Merge close splice sites (# of bases)", min = 0, step = 1, value = 3),
          actionButton("peak_detection_apply", "Apply")
        ),
        wellPanel(
          h4("Custom coordinates"),
          radioButtons("known_sites_merge", "Replace or merge with the existing coordinates", c("Replace", "Merge"), selected="Replace"),
          fileInput("known_sites", "File with exon coordinates")
        )
      ),
      
      wellPanel(
        h3("Display options"),
        numericInput("abundance", "Isoform abundance threshold to display (%)", min = 0, max = 100, value = 1),
        sliderInput("iso_pixels", "Isoform separation (px)", min = 10, max = 40, value = 20, step = 1, round = TRUE),
        sliderInput("histo_pixels", "Barplot height (px)", min = 100, max = 500, value = 250, step = 10, round = TRUE),
        actionButton("abundance_apply", "Apply")
      ),
      
      wellPanel(
        h3("Defined isoforms for comparison"),
        radioButtons("transcriptinputtype", "Input format", c("GTF", "GFF3"), selected="GTF"),
        fileInput("transcriptinputfile", "Input file"),
        uiOutput("transcripts_data_wait")
      ),
      
      
      wellPanel(
        # h3("Download"), 
        # shinyDirButton('output', 'Download directory', 'Please select a folder', FALSE),
        # textOutput("output_dir"),
        textInput(inputId='prefix_output', label='Output prefix', 
                  value = gsub("(.gff3|.gff|.bed|.bed6|.bed12)$", "", input$inputfile$name, perl = T)),
        # actionButton("save_apply", "Save"),
        
      ),
      
    )
  })
 #   download_plot_html, download_plot_pdf, download_bp_info, download_exon_info, download_iso_information_tsv, download_iso_information_gtf, download_read_clasification

  
  output$input_wait_mainpanel<-renderUI({
    req(input$inputfile)
    tagList(
      fluidRow(
        column(12, class = "well", 
               fluidRow(uiOutput('ui_plot')),
               
               fluidRow(
                 column(width = 2, downloadButton("download_plot_html", "Download dynamic plot (.html)")),
                 
                 column(width = 2, downloadButton("download_plot_pdf", "Download static plot")),
                 
                 tags$head(
                   tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: top; }
                          #inline .form-group { display: table-row; }
                          #inline .selectize-control.single div.item {padding-right: 15px; }
                          #inline label.control-label {padding-right: 10px; }")
                 ),
                 column(width = 2, tags$div(id = "inline", selectInput(inputId='plot_format', label='Format', choices = c("png", "jpg", "jpeg", "webp", "svg", "pdf"), selected = "pdf"))),
                 column(width = 2, tags$div(id = "inline", numericInput(inputId="plot_width", label="width (px)", min = 100, value = 1000))),
                 column(width = 2, uiOutput("download_plot_wait"))
               )
               
        )
      ),
      
      
      fluidRow(class = "well",
               column(6, class = "well", h3("Exonic starting points"), DTOutput('start_dto'), actionButton("start_bp_add", "Add row")),
               column(6, class = "well", h3("Exonic ending points"), DTOutput('end_dto'), actionButton("end_bp_add", "Add row")),
               fluidRow(column(12, align = "center",actionButton("break_point_apply", "Apply changes"))),
               fluidRow(column(12, align = "left", downloadButton("download_bp_info", "Download exon coordinates (.tsv)")))
      ),
      
      fluidRow(class = "well",
               column(7, class = "well", h3("Isoform information"), 
                      DTOutput('isoform_information'),
                      downloadButton("download_iso_information_tsv", "Download isoform information (.tsv)"),
                      downloadButton("download_iso_information_gtf", "Download isoform information (.gtf)"),
                      downloadButton("download_read_clasification", "Download read-isoform clasification (.tsv)")), #isoform_df_plot
               
               column(5, class = "well", h3("Exon information"), 
                      DTOutput('exon_information'),
                      downloadButton("download_exon_info", "Download exon information (.tsv)")),
               
               
      ),
    )
  })
  
  
  
  
  
  
  
  
  output$download_plot_wait<-renderUI({
    req(rv$multiplot$height)
    tagList(
          # tags$head(
          #   tags$style(type="text/css", "#inline label{ display: table-cell; text-align: center; vertical-align: top; } 
          #            #inline .form-group { display: table-row; }
          #            #inline .selectize-control.single div.item {padding-right: 15px; }
          #            #inline label.control-label {padding-right: 10px; }")
          # ),
          tags$div(id = "inline", numericInput(inputId="plot_height", label="height (px)", min = 100, value = rv$multiplot$height))
        )

  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ######################
  # Input data loading #
  ######################
  observeEvent(input$inputfile, {
    file <- input$inputfile
    ext <- tools::file_ext(file$datapath)
    
    print("Loading data (prerequisite)")
    req(file)
    
    withProgress(message = 'Loading data', value = 1, {
      if (input$inputtype == "GFF3"){
        validate(need(ext %in% c("gff3", "gff"), "Please upload a GFF3 file"))
        rv$orig <- gff3_data_loading(file$datapath)
      } else if(input$inputtype == "BED6"){
        validate(need(ext %in% c("bed", "bed6"), "Please upload a 6 column BED file"))
        rv$orig <- bed_data_loading(file$datapath)
      } else if(input$inputtype == "BAM"){
        validate(need(ext %in% c("bam"), "Please upload a BAM file"))
        system(paste0("bamToBed -split -i ", file$datapath, " > ", gsub(".bam$", "", file$name, perl = T), ".bed"))
        rv$orig <- bed_data_loading(paste0(gsub(".bam$", "", file$name, perl = T), ".bed"))
      }
    })
    print("Loading data")
  })
  
  
  
  
  #############################
  # Gene selection to analyce #
  #############################

  
  observeEvent(ignoreInit=T, c( 
    input$selectgene
    ),{ 
      print("Gene selection (prerequisite rv$orig)")
      req(rv$orig)
      print("Gene selection (prerequisite input$selectgene)")
      print(input$selectgene)
      req(input$selectgene)
      
      withProgress(message = 'Selecting gene', value = 1, {
        rv$data <- rv$orig[rv$orig$gene == input$selectgene,]
        rv$num_reads_initial = length(unique(rv$data$id))
      })
      print("Gene selection")
  })
  
  
  
  ############
  # Trimming #
  ############
  output$studied_range<-renderUI({
    req(rv$data)
    tagList(
      numericInput("numericInput_beginning", "Start", min = 0, max = (max(rv$data$end) + 1), value = 0),
      numericInput("numericInput_final", "End", min = 0, max = (max(rv$data$end) + 1), value = (max(rv$data$end) + 1)),
      actionButton("analysis_range_apply", "Apply")
    )
  })
  
  observeEvent(ignoreInit=T, c(
    rv$data,
    input$analysis_range_apply
    ),{ 
      print("Trimming (prerequisites rv$data)")
      req(rv$data)

      
      if (is.null(input$numericInput_beginning)) numericInput_beginning = 0 else numericInput_beginning = input$numericInput_beginning
      if (is.null(input$numericInput_final)) numericInput_final = max(rv$data$end) + 1 else numericInput_final = input$numericInput_final

      withProgress(message = 'Trimming', value = 1, {
        rv$trimmed_data = rv$data[rv$data$end >= numericInput_beginning & rv$data$start <= numericInput_final, ]
        rv$num_reads_post_trimming = length(unique(rv$trimmed_data$id))
      })
      print("Trimming")
    })
  
  
  
  ##########################
  # Break point definition #
  ##########################
  
  # Automatic break points

  observeEvent(ignoreInit=T, c(
    rv$trimmed_data,
    input$peak_detection_apply
    ), {
      print("Automatic break point calculation (prerequisites)")
      req(rv$trimmed_data)
      req(input$peak_threshold)
      req(input$very_close_bp)
      req(input$padding)
      validate(need(input$very_close_bp < input$padding, "Very close break points should be smaller than the padding"))
      print("Automatic break point calculation (prerequisites)")
      
      withProgress(message = 'Automatic break point calculation', value = 1, {
        rv$break_points = break_point_calculation(raw_exons=rv$trimmed_data,
                                                  num_reads_post_trimming=rv$num_reads_post_trimming,
                                                  breakpoint_freq_threshold=input$peak_threshold,
                                                  very_close_bp=input$very_close_bp,
                                                  breakpoint_padding=input$padding)
        rv$break_points_original = rv$break_points
      })
      print("Automatic break point calculation")

    })


  
  
  #######################
  # Break point edition #
  #######################
  
  # Cut break points
  proxy_start <- dataTableProxy('start_dto')
  output$start_dto <-renderDT({
    if (is.null(rv$bp_to_show[[1]])) return(NULL)
    cbind(merge(rv$break_points[[1]], rv$bp_to_show[[1]], by.x = c("breakp", "left", "right"), by.y = c("breakpoint", "lower_limit", "upper_limit"), all.x = T, sort = F),
          Remove = unlist(lapply(1:nrow(rv$break_points[[1]]), function(x) getRemoveButton(x, idS = "", lab = "Tabstart"))))},
    editable = list(target = "cell", disable = list(columns = c(3,4,5))),
    colnames = c("Breakpoint", "Lower limit", "Upper limit", "# reads at breakpoint (%)", "# of reads in the interval  (%)", ""),
    rownames = FALSE,
    escape   = FALSE
    )
  
  # output$start_dto <-renderDT({rv$break_points[[1]]}, editable = TRUE, rownames = FALSE)
  # Edit row
  observeEvent(input$start_dto_cell_edit, {
    info = input$start_dto_cell_edit
    i = info$row
    j = info$col + 1
    v = info$value
    rv$break_points[[1]][i, j] <<- DT::coerceValue(v, rv$break_points[[1]][i, j])
    replaceData(proxy_start, rv$break_points[[1]], resetPaging = FALSE, rownames = FALSE)
    print("start modification")
  })
  
  # Remove row
  observeEvent(input$remove_button_Tabstart, {
    myTable <- rv$break_points[[1]]
    s <- as.numeric(strsplit(input$remove_button_Tabstart, "_")[[1]][2])
    myTable <- myTable[-s,]
    replaceData(proxy_start, myTable, resetPaging = FALSE, rownames = FALSE)
    rv$break_points[[1]] <- myTable
  })

  # Add row
  observeEvent(input$start_bp_add, {
    myTable <- rv$break_points[[1]]
    myTable <- rbind(myTable, 0)
    replaceData(proxy_start, myTable, resetPaging = FALSE, rownames = FALSE)
    rv$break_points[[1]] <- myTable
  })
  
  # output$start_dto_original <-renderDT({
  #   if (is.null(rv$break_points[[1]])) return(NULL)
  #   rv$break_points[[1]]},
  #   rownames = FALSE,
  #   escape   = FALSE
  # )
  
  
  
  proxy_end <- dataTableProxy('end_dto')
  output$end_dto <-renderDT({
    if (is.null(rv$bp_to_show[[2]])) return(NULL)
    cbind(merge(rv$break_points[[2]], rv$bp_to_show[[2]], by.x = c("breakp", "left", "right"), by.y = c("breakpoint", "lower_limit", "upper_limit"), all.x = T, sort = F),
          Remove = unlist(lapply(1:nrow(rv$break_points[[2]]), function(x) getRemoveButton(x, idS = "", lab = "Tabend"))))},
    editable = list(target = "cell", disable = list(columns = c(3,4,5))),
    colnames = c("Breakpoint", "Lower limit", "Upper limit", "# reads at breakpoint (%)", "# of reads in the interval  (%)", ""),
    rownames = FALSE,
    escape   = FALSE
    )
  
  # output$end_dto <-renderDT({rv$break_points[[2]]}, editable = TRUE, rownames = FALSE)
  # Edit row
  observeEvent(input$end_dto_cell_edit, {
    info = input$end_dto_cell_edit
    i = info$row
    j = info$col + 1
    v = info$value
    rv$break_points[[2]][i, j] <<- DT::coerceValue(v, rv$break_points[[2]][i, j])
    replaceData(proxy_end, rv$break_points[[2]], resetPaging = FALSE, rownames = FALSE)
    print("end modification")
  })

  # Remove row
  observeEvent(input$remove_button_Tabend, {
    myTable <- rv$break_points[[2]]
    s <- as.numeric(strsplit(input$remove_button_Tabend, "_")[[1]][2])
    myTable <- myTable[-s,]
    replaceData(proxy_end, myTable, resetPaging = FALSE, rownames = FALSE)
    rv$break_points[[2]] <- myTable
  })
  
  
  # Add row
  observeEvent(input$end_bp_add, {
    myTable <- rv$break_points[[2]]
    myTable <- rbind(myTable, 0)
    replaceData(proxy_end, myTable, resetPaging = FALSE, rownames = FALSE)
    rv$break_points[[2]] <- myTable
  })
  
  
  

  ######################
  # Known break points #
  ######################

  
  observeEvent(ignoreInit=T, c(
    input$known_sites
    ), {
      print("Known break points")

      file <- input$known_sites
      
      print("Known break points before req")
      req(file)
      print("Known break points ater req")
      
      withProgress(message = 'Loading known break points', value = 1, {
        if (input$known_sites_merge == "Replace"){
          rv$break_points <- known_sites_replace(file$datapath, input$selectgene)
        }else if (input$known_sites_merge == "Merge"){
          rv$break_points <- known_sites_merge(file$datapath, rv$break_points, input$selectgene)
        }
        
      })
  })



    
  
  
  #########################
  # Break point asignment #
  #########################
  
  observeEvent(ignoreInit=T, c(
    rv$break_points
  ), {
    print("break point assignment before req(rv$break_points)")
    req(rv$break_points)
    
    
    withProgress(message = 'Break point asignment',value = 1, {
      # Break point asignment
      print("break point assignment")
      rv$bp_asigned_exons <- break_point_asignment(rv$trimmed_data, rv$break_points)
      print("break point assignment (end)")
      
      rv$bp_to_show <- break_point_info_calculator(rv$break_points, rv$num_reads_post_trimming, rv$bp_asigned_exons)
      # print(rv$bp_to_show[[2]])
      print("rv$bp_to_show creation (end)")
    })
  })
    
  
  ###############################
  # exon and isoform definition #
  ###############################

  observeEvent(ignoreInit=T, c(
    rv$break_points_original,
    input$break_point_apply
  ), {
    print("break point assignment before req(rv$break_points)")
    req(rv$bp_asigned_exons)

    
    withProgress(message = 'Exon and isoform definition', {
      # # Break point asignment
      # print("break point assignment")
      # rv$bp_asigned_exons <- break_point_asignment(rv$trimmed_data, rv$break_points)
      # rv$bp_to_show = break_point_info_calculator(rv$break_points, rv$num_reads_post_trimming, rv$bp_asigned_exons)
        
      incProgress(1/3)
      
      # Exon definition
      print("exon definition")
      rv$defined_exons <- exon_definition(rv$bp_asigned_exons)
      rv$exon_information = exon_info_fun(rv$defined_exons[[1]])
      rv$num_reads_final = length(unique(rv$defined_exons[[1]]$id))
      incProgress(1/3)
      
      # Isoform definition
      print("isoform definition")
      rv$isoforms <- isoform_definition(rv$defined_exons[[1]]) # Slow
      incProgress(1/3)
      
    })
  })
 
  
  
  observeEvent(ignoreInit=T, c(
    rv$isoforms,
    input$full_reads
  ), { 
    
    req(rv$isoforms)
    
    withProgress(message = 'Exon and isoform definition', {
      
      # Full length filter and isoform information
      if (input$full_reads){
        rv$isoforms_filtered = isoform_full_length_filter(rv$isoforms, rv$break_points)
        rv$isoforms_information <- isoform_info_fun(rv$isoforms_filtered[[1]], rv$num_reads_initial, rv$num_reads_post_trimming, rv$num_reads_final, length(rv$isoforms_filtered[[3]]))
      }else{
        rv$isoforms_filtered = rv$isoforms
        rv$isoforms_filtered[[3]] = NA
        rv$isoforms_information <- isoform_info_fun(rv$isoforms_filtered[[1]], rv$num_reads_initial, rv$num_reads_post_trimming, rv$num_reads_final)
      }
      incProgress(1/3)
      
      # df to plot
      print("df2plot")
      rv$isoform_df_plot <- isoform_df_plot_fun(rv$isoforms_filtered[[1]])
      incProgress(1/5)
      rv$df_superplot <- break_point_freq_fun(rv$trimmed_data)
      incProgress(1/5)
      
      print( "df2plot finish")
    })
  })

  output$exon_information <-renderDT({
    if (is.null(rv$exon_information)) return(NULL)
    rv$exon_information},
    colnames = c("Exon", "Size (# of bases)", "# of reads", "% of reads (within all mapped reads)"),
    rownames = FALSE,
    editable = FALSE,
    escape   = FALSE
  )

  output$isoform_information <-renderDT({
    if (is.null(rv$isoforms_information)) return(NULL)
    rv$isoforms_information[, -2]},
    colnames = c("Isoform ID", "Size (# of bases)", "# of reads", "Partial % (within clasified reads)", "Total % (within all mapped reads)"),
    rownames = FALSE,
    editable = FALSE,
    escape   = FALSE
  )



  ################
  # Plot funtion #
  ################

  observeEvent(ignoreInit=T, c(
    rv$isoform_df_plot,
    rv$df_superplot,
    input$abundance_apply,
    rv$transcripts_df_to_plot
    
  ), {
    print("plot")
    req(rv$isoform_df_plot)
    req(rv$df_superplot)
    if (is.null(input$abundance)) abundance = 1 else abundance = input$abundance
    if (is.null(input$histo_pixels)) rv$histo_pixels = 250 else rv$histo_pixels = input$histo_pixels
    if (is.null(input$iso_pixels)) rv$iso_pixels = 20 else rv$iso_pixels = input$iso_pixels
    
    withProgress(message = 'Plotting', value = 1, {
      
      # Abundance filter
      rv$isoform_df_plot_filtered = rv$isoform_df_plot[rv$isoform_df_plot$perc >= abundance, c("x_pos_start", "x_pos_end", "y_pos", "star_stop", "grupo")]
      
      
      # Color palete syncronization
      isoform_df_plot_and_gtf = rbind(rv$isoform_df_plot_filtered, rv$transcripts_df_to_plot)
      ordered_start_stop = unique(isoform_df_plot_and_gtf[order(isoform_df_plot_and_gtf$x_pos_start, isoform_df_plot_and_gtf$x_pos_end),"star_stop"])
      isoform_df_plot_and_gtf$star_stop = factor(isoform_df_plot_and_gtf$star_stop, levels = ordered_start_stop)
      
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      
      pal <- col_vector[1:length(unique(isoform_df_plot_and_gtf$star_stop))]
      pal <- setNames(pal, unique(isoform_df_plot_and_gtf$star_stop))
      
      legend_exons_df = data.frame(x_pos_start = rv$isoform_df_plot_filtered$x_pos_start[1],
                                   x_pos_end = rv$isoform_df_plot_filtered$x_pos_start[1],
                                   y_pos = rv$isoform_df_plot_filtered$y_pos[1],
                                   star_stop = as.character(isoform_df_plot_and_gtf$star_stop),
                                   grupo = "legend_colors",
                                   stringsAsFactors = F) 
      rv$isoform_df_plot_filtered = rbind(rv$isoform_df_plot_filtered, legend_exons_df)
      
      
      # Plot
      print("isoform_plot")
      rv$isoform_plot1 <- isoform_plot_fun(rv$isoform_df_plot_filtered, 
                                           "VisoQLR detected isoforms\nIsoform ID | size | % of reads", 
                                           pal, 
                                           ordered_start_stop, 
                                           showlegend = TRUE)
      rv$num_iso1 <- max(length(unique(rv$isoform_df_plot_filtered[,"y_pos"])), 11)
      
      
      
      print("break_point_freq_plot")
      rv$break_point_freq_plot <- break_point_freq_plot_fun(rv$df_superplot, 
                                                            rv$break_points, 
                                                            input$start_dto_rows_selected, 
                                                            input$end_dto_rows_selected)
  
      
      if (is.null(rv$transcripts_df_to_plot)){
        
        rv$multiplot <- multiplot_fun1(rv$isoform_plot1, rv$break_point_freq_plot, rv$num_iso1, rv$iso_pixels, rv$histo_pixels)
        rv$multiplot$height=rv$num_iso1*rv$iso_pixels+rv$histo_pixels

      }else{
        
        rv$isoform_plot2 = isoform_plot_fun(rv$transcripts_df_to_plot, 
                                            gsub(" \\| $", "", paste("Loaded transcripts\nTranscript ID", paste(input$transcriptinfo2display,  sep = "", collapse = " | "), sep = " | ")),
                                            pal,
                                            ordered_start_stop,
                                            showlegend = FALSE)
        rv$num_iso2 <- max(length(unique(rv$transcripts_df_to_plot[,"y_pos"])), 11)
        rv$multiplot <- multiplot_fun2(rv$isoform_plot2, rv$isoform_plot1, rv$break_point_freq_plot, rv$num_iso2, rv$num_iso1, rv$iso_pixels, rv$histo_pixels)
        rv$multiplot$height=rv$num_iso1*rv$iso_pixels+rv$histo_pixels+rv$num_iso2*rv$iso_pixels
        
      }
      
 
  
      print("plot1")
      output$plot1 <- renderPlotly(rv$multiplot)
  
      
      print("ui_plot")
      output$ui_plot <- renderUI({
        plotlyOutput("plot1", height = rv$multiplot$height)
        })
      # cat(file=stderr(), "plot finish\n")
      print("plot finish")
      
    })
    
    })


  
  
  
  ###########################
  # Transcript data loading #
  ###########################

  output$transcripts_data_wait<-renderUI({
    req(input$transcriptinputfile)
    tagList(
      numericInput("num_transcript2display", "Maximum number of transcripts to display", min = 1, max = 100, value = 10),
      selectInput(inputId='transcriptorder2display', label='Sort transcripts by', choices = c("id", rv$transcriptcolumns), selected = "id"),
      selectInput(inputId='decreasingorder2display', label='Sort increasing/decrasing', choices = c("Increasing", "Decreasing"), selected = "Increasing"),
      checkboxGroupInput("transcriptinfo2display", "Information to display", rv$transcriptcolumns, selected="size"),
      actionButton("transcriptinfo2display_apply", "Apply")
    )
  })


  # Data loading
  observeEvent(ignoreInit=T, c(
    input$transcriptinputtype,
    input$transcriptinputfile
  ), {
    file <- input$transcriptinputfile
    ext <- tools::file_ext(file$datapath)

    print("Loading input transcrits (prerequisite)")
    req(file)

    withProgress(message = 'Loading transcripts', value = 1, {
      if (input$transcriptinputtype == "GFF3"){
        validate(need(ext %in% c("gff3", "gff"), "Please upload a GFF3 file with the extension .gff or .gff3"))
        rv$transcriptinput <- transcript_data_loading(file$datapath, input$transcriptinputtype)
      } else {
        validate(need(ext %in% c("gtf"), "Please upload a GTF file with the extension .gtf"))
        rv$transcriptinput <- transcript_data_loading(file$datapath, input$transcriptinputtype)
      }

    })


    print("Loading data")
  })


  # Gene selection
  observeEvent(ignoreInit=T, c(
    input$selectgene,
    rv$transcriptinput
  ),{
    print("Gene selection transcripts (prerequisites)")
    req(rv$transcriptinput)
    req(input$selectgene)
    print("Gene selection transcripts")
    
    withProgress(message = 'Selecting gene', value = 1, {
      
      rv$transcriptcolumns = colnames(rv$transcriptinput[,6:ncol(rv$transcriptinput)])
      rv$gene_filter_transcripts = rv$transcriptinput[rv$transcriptinput$gene == input$selectgene, ]
      
      # Frequency calculation, numeric annotations
      rv$gene_filter_transcripts_unique = rv$gene_filter_transcripts[!duplicated(rv$gene_filter_transcripts$id), ]
      new_col = c()
      for (i in rv$transcriptcolumns){
        if (is.numeric(rv$gene_filter_transcripts[,i])){
          rv$gene_filter_transcripts_unique[,paste0(i, " (%)")] = paste0(round(rv$gene_filter_transcripts_unique[,i]/sum(rv$gene_filter_transcripts_unique[,i])*100, 2), "%")
          new_col = c(new_col, paste0(i, " (%)"))
        }
      }
      rv$gene_filter_transcripts_unique = rv$gene_filter_transcripts_unique[,c("id", new_col), drop = F]
      rv$gene_filter_transcripts = merge(rv$gene_filter_transcripts, rv$gene_filter_transcripts_unique, by = "id")
      rv$transcriptcolumns = c(rv$transcriptcolumns, new_col)
    })
    print("Gene selection")
  })
  
  
  
  observeEvent(ignoreInit=T, c(
    rv$transcriptcolumns,
    input$transcriptinfo2display_apply,
    rv$transcriptinput
  ),{
    print("Transcript data to plot (prerequisites)")
    req(rv$gene_filter_transcripts)
    req(input$num_transcript2display)
    req(input$transcriptorder2display)
    req(input$decreasingorder2display)
    req(input$transcriptinfo2display)
    print("Transcript data to plot")
    
    withProgress(message = 'Data preparation to plot', value = 1, {
      
      rv$gene_filter_transcripts$y_pos = apply( rv$gene_filter_transcripts[ , c("id", input$transcriptinfo2display) , drop = F] , 1 , paste , collapse = " | " )

      rv$transcripts_df_to_plot = data.frame(
        x_pos_start = rv$gene_filter_transcripts$start,
        x_pos_end = rv$gene_filter_transcripts$end,
        y_pos = rv$gene_filter_transcripts$y_pos,
        star_stop = paste(rv$gene_filter_transcripts$start, rv$gene_filter_transcripts$end, sep = "-"),
        grupo = "Loaded transcripts",
        stringsAsFactors = F)

            
      if (input$decreasingorder2display == "Decreasing") decreasing = TRUE else decreasing = FALSE
      factorlevels = unique(rv$gene_filter_transcripts[order(rv$gene_filter_transcripts[,input$transcriptorder2display], decreasing = decreasing), "y_pos"])
      factorlevels = factorlevels[1:min(length(factorlevels),input$num_transcript2display)]
      factorlevels = rev(factorlevels)
      rv$transcripts_df_to_plot = rv$transcripts_df_to_plot[rv$transcripts_df_to_plot$y_pos %in% factorlevels,]
      rv$transcripts_df_to_plot$y_pos = factor(rv$transcripts_df_to_plot$y_pos, levels = factorlevels)
      
    })
    print("Transcript data to plot (end)")
  })




  
  
  
  
  
  
  

  ################
  # Write output #
  ################
  
 
  # observeEvent(ignoreInit=T, c(
  # ),{
  # })
  
  observe({
    if (is.null(input$prefix_output)) rv$prefix_output = gsub("(.gff3|.gff|.bed|.bed6|.bed12)$", "", input$inputfile$name, perl = T) else rv$prefix_output = input$prefix_output

  # observe({})
  # observe({})
  # observe({})
  # observe({})
  if (is.null(input$plot_format)) rv$plot_format = "pdf" else rv$plot_format = input$plot_format
  if (is.null(input$plot_width)) rv$plot_width = 1000 else rv$plot_width = input$plot_width
  if (is.null(input$plot_height)) rv$plot_height = rv$multiplot$height else rv$plot_height = input$plot_height
    
  print(is.null(input$plot_format))
  print(input$plot_format)
  })


  
  
  # Dynamic plot (html)
  output$download_plot_html <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".combined_plot.html")},
    content = function(file){saveWidget(widget = rv$multiplot, file = file, selfcontained = T)}
  )
  # Static plot (pdf)
  output$download_plot_pdf <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".combined_plot.", rv$plot_format)},
    content = function(file){save_image(rv$multiplot, file, height = rv$plot_height, width = rv$plot_width)}
  )
  

  # Break point information
  output$download_bp_info <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".breakpoints_info.tsv")},
    content = function(file){write.table(rv$bp_to_show[[3]], file, sep = "\t", quote = F, row.names = F, col.names = T)}
  )
  # Exon information
  output$download_exon_info <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".exon_info.tsv")},
    content = function(file){write.table(rv$exon_information, file, sep = "\t", quote = F, row.names = F, col.names = T)}
  )
  # Isoform information (TSV)
  output$download_iso_information_tsv <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".isoform_info.tsv")},
    content = function(file){write.table(rv$isoforms_information, file, sep = "\t", quote = F, row.names = F, col.names = T)}
  )
  # Isoform information (GTF)
  output$download_iso_information_gtf <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".isoform_info.gtf")},
    content = function(file){write.table(isoinfo2gtf(rv$isoforms_information, rv$exon_information, input$selectgene), 
                                         file, sep = "\t", quote = F, row.names = F, col.names = F)}
  )
  # Read clasification
  output$download_read_clasification <- downloadHandler(
    filename = function(){paste0(rv$prefix_output, ".read_clasification.tsv")},
    content = function(file){write.table(read_clasification_fun(all_ids = rv$data$id, no_consensus_ids = rv$defined_exons[[2]], classified_ids = rv$isoforms_filtered[[2]], not_full_length = rv$isoforms_filtered[[3]], iso_dict = rv$isoforms_information), 
                                         file, sep = "\t", quote = F, row.names = F, col.names = T)}
  )
  
  
  
  
  
  
  
}



shinyApp(ui, server)





























