library(shiny)
library(ggplot2)
library(vctrs)
library(plotly)
library(DT)
library(htmlwidgets)
library(shinyFiles)
library(reticulate)

# install.packages('reticulate')
# reticulate::install_miniconda()
# reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
# reticulate::use_miniconda('r-reticulate')

rm(list=ls()) 

options(shiny.maxRequestSize=3000*1024^2)



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
  bed6 = bed6[bed6$score == 60,]
  return(bed6)
}


transcript_data_loading <- function(input_path,  inputformat = "GTF"){
  if (inputformat == "GTF") {separation = " "}
  if (inputformat == "GFF3") {separation = "="}
  
  gtf = read.delim(input_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  
  transcripts = gtf[gtf$V3 == "transcript","V9"]
  transcriptinfo = gsub("; ", ";", transcripts)
  transcriptinfo = do.call("rbind", strsplit(x = transcriptinfo, split = ";"))
  cnames = do.call("rbind", strsplit(x = transcriptinfo[1,], split = separation))[,1]
  transcriptinfo = gsub("^.* ", "", transcriptinfo, perl = T)
  colnames(transcriptinfo) = cnames
  
  gtf = gtf[gtf$V3 == "exon", c("V1", "V4","V5","V9", "V6")]
  colnames(gtf) = c("gene", "start", "end", "id", "score")
  gtf$score = as.numeric(gtf$score)
  gtf$id=gsub("^.*transcript_id ", "", gtf$id, perl = TRUE)
  gtf$id=gsub(";.*$", "", gtf$id, perl = TRUE)
  
  gtf = merge(gtf, transcriptinfo, by.x = "id", by.y = "transcript_id")
  
  # Size annotation
  gtf_split = split(gtf, gtf$id)
  gtf_lengths = unlist(lapply(gtf_split, function(x) sum(x$end - x$start)))
  gtf$size = paste0(gtf_lengths[gtf$id], "bp")
  
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
      for (i in 1:(nrow(start_position)-1)){
        
        if (is.element(TRUE, start_position[i, "breakp"] + very_close_bp >= start_position[, "breakp"] & 
                       start_position[i, "breakp"] - very_close_bp <= start_position[, "breakp"] &
                       start_position[i, "freq"] * 2 <= start_position[, "freq"])) {start_pos_to_remove = c(start_pos_to_remove, i)}
      }
      if (!is.null(start_pos_to_remove)) {start_position = start_position[-start_pos_to_remove,]}
      
      end_pos_to_remove=c()
      for (i in 1:(nrow(end_position)-1)){
        if (is.element(TRUE, end_position[i, "breakp"] + very_close_bp >= end_position[, "breakp"] & 
                       end_position[i, "breakp"] - very_close_bp <= end_position[, "breakp"] &
                       end_position[i, "freq"] * 2 <= end_position[, "freq"])) {end_pos_to_remove = c(end_pos_to_remove, i)}
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




############################
# Known break points merge #
############################
known_sites_merge <- function(known_sites_file, break_points_list, gene, raw_exons){
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
  exon_info_df$percent_of_reads = round(exon_info_df$number_of_reads / length(unique(defined_exons$id)) * 100, 2)
  exon_info_df = exon_info_df[order(exon_info_df$exon),]
  
  exon_info_df$Size = unlist(lapply(gsub("_", "+", exon_info_df[,"exon"]), function(x) -eval(parse(text = x))))
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




# Isoform information
isoform_info_fun <- function(isoform_frequencies, n_inicial, n_post_trim, n_final){
  
  n_vector_reads = n_inicial - n_post_trim
  n_no_consensous_breakpoint_reads = n_post_trim - n_final
  
  other_groups_df = data.frame(isoform_id = c("Only vector reads", "No consensous breakpoint reads"),
                               read_isoform = c("-", "-"), 
                               Freq = c(n_vector_reads, n_no_consensous_breakpoint_reads),
                               perc = c("-", "-"), 
                               Size = c(0, 0),
                               stringsAsFactors = FALSE)
  
  isoform_id = paste("Iso", 1:nrow(isoform_frequencies), sep = "")
  isoform_frequencies = cbind(isoform_id, isoform_frequencies)
  isoform_frequencies$Size = unlist(lapply(gsub("_", "+", isoform_frequencies[,2]), function(x) -eval(parse(text = x))))
  
  all_groups = rbind(other_groups_df, isoform_frequencies)
  all_groups$perc_total = round(all_groups$Freq*100/sum(all_groups$Freq), 2)
  colnames(all_groups) =  c("Isoform_id", "Isoform", "N_of_reads", "Prerc_partial", "Size", "Perc_total")
  
  all_groups = all_groups[,c("Isoform_id", "Isoform", "Size", "N_of_reads", "Prerc_partial", "Perc_total")]
  return(all_groups)
}



# Read classification
read_clasification_fun <- function(all_ids, no_consensus_ids, classified_ids){
  classified_df <- data.frame(read_id = names(classified_ids), type = as.character(classified_ids))
  
  only_vector_ids = setdiff(unique(all_ids), c(unique(no_consensus_ids), names(classified_ids)))
  if (length(only_vector_ids) != 0){
    trimmed_df <- data.frame(read_id = only_vector_ids, type = "Only_vector_reads")
    classified_df = rbind(classified_df, trimmed_df)
  }
  
  if (length(unique(no_consensus_ids)) != 0){
    no_consensus_df <- data.frame(read_id = unique(no_consensus_ids), type = "No_consensous_breakpoint_reads")
    classified_df = rbind(classified_df, no_consensus_df)
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
    isosize = paste0(-eval(parse(text = gsub("_", "+", isoform_frequencies[i,1]))), "bp")
    isoform_id = paste("Iso", isoform_num, " | ", isoform_frequencies[i,3], "% | ", isosize, sep = "")
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





























#############################################################################################
###   USER INTERFACE   ######################################################################
#############################################################################################

ui <- fluidPage(
  titlePanel("VIsoQLR"),
  
  sidebarLayout(
    sidebarPanel(
      width=2,
      
      wellPanel(
        h3("Input"),
        radioButtons("inputtype", "Input format", c("GFF3", "BED6"), selected="GFF3"),
        fileInput("inputfile", "Input file"),
        uiOutput("known_sites_wait"),
        uiOutput("transcripts_data_wait")
        # fileInput("reference", "Reference sequence"),
      ),
      uiOutput("input_wait_sidebarpanel"),
    ),
    
    mainPanel(
      width = 10,
      uiOutput("input_wait_mainpanel"),
      ))
)









#############################################################################################
###   SERVER   ##############################################################################
#############################################################################################

server <- function(input, output) {
  
  
  output$input_wait_sidebarpanel<-renderUI({
    req(input$inputfile)
    tagList(
      wellPanel(
        h3("Download"), 
        shinyDirButton('output', 'Download directory', 'Please select a folder', FALSE),
        textOutput("output_dir"),
        textInput(inputId='prefix_output', label='Output prefix', 
                  value = gsub("(.gff3|.gff|.bed|.bed6|.bed12)$", "", input$inputfile$name, perl = T)),
        actionButton("save_apply", "Save")
      ),
      
      wellPanel(
        h3("Analysis window"),
        selectInput(inputId='selectgene', label='Select gene', choices = unique(rv$orig$gene)),
        uiOutput("studied_range")
      ),
      
      
      wellPanel(
        h3("Automatic peak detection"),
        numericInput("peak_threshold", "Read threshold (%)", min = 0, max = 100, value = 3),
        numericInput("padding", "Padding (# of bases)", min = 0, step = 1, value = 5),
        numericInput("very_close_bp", "Merge close splice sites (# of bases)", min = 0, step = 1, value = 3),
        actionButton("peak_detection_apply", "Apply")
      ),
      
      wellPanel(
        h3("Display options"),
        numericInput("abundance", "Isoform abundance threshold to display (%)", min = 0, max = 100, value = 1),
        sliderInput("iso_pixels", "Isoform separation (px)", min = 10, max = 40, value = 20, step = 1, round = TRUE),
        sliderInput("histo_pixels", "Barplot height (px)", min = 100, max = 500, value = 250, step = 10, round = TRUE),
        actionButton("abundance_apply", "Apply")
      ),
    )
  })
  
  
  output$input_wait_mainpanel<-renderUI({
    req(input$inputfile)
    tagList(
      fluidRow(
        column(12, class = "well", uiOutput('ui_plot'))
      ),
      
      
      fluidRow(class = "well",
               column(6, class = "well", h3("Exonic starting points"), DTOutput('start_dto'), actionButton("start_bp_add", "Add row")),
               column(6, class = "well", h3("Exonic ending points"), DTOutput('end_dto'), actionButton("end_bp_add", "Add row")),
               fluidRow(column(12, align = "center",actionButton("break_point_apply", "Apply changes"))),
      ),
      
      fluidRow(class = "well",
               column(7, class = "well", h3("Isoform information"), DTOutput('isoform_information')), #isoform_df_plot
               column(5, class = "well", h3("Exon information"), DTOutput('exon_information'))
      ),
    )
  })
  
  
  
  
  
  
  rv <- reactiveValues(data = NULL, orig=NULL)
  
  
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
      } else {
        validate(need(ext %in% c("bed", "bed6"), "Please upload a 6 column BED file"))
        rv$orig <- bed_data_loading(file$datapath)
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
  output$known_sites_wait<-renderUI({
    req(input$inputfile)
    fileInput("known_sites", "Known break points")
  })
  
  observeEvent(ignoreInit=T, c(
    input$known_sites
    ), {
      print("Known break points")

      file <- input$known_sites
      
      print("Known break points before req")
      req(file)
      print("Known break points ater req")
      
      withProgress(message = 'Loading known break points', value = 1, {
        rv$break_points <- known_sites_merge(file$datapath, rv$break_points, input$selectgene)
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
    req(rv$break_points)

    
    withProgress(message = 'Exon and isoform definition', {
      # # Break point asignment
      # print("break point assignment")
      # rv$bp_asigned_exons <- break_point_asignment(rv$trimmed_data, rv$break_points)
      # rv$bp_to_show = break_point_info_calculator(rv$break_points, rv$num_reads_post_trimming, rv$bp_asigned_exons)
        
      incProgress(1/5)
      
      # Exon definition
      print("exon definition")
      rv$defined_exons <- exon_definition(rv$bp_asigned_exons)
      rv$exon_information = exon_info_fun(rv$defined_exons[[1]])
      rv$num_reads_final = length(unique(rv$defined_exons[[1]]$id))
      incProgress(1/5)
      
      # Isoform definition
      print("isoform definition")
      rv$isoforms <- isoform_definition(rv$defined_exons[[1]]) # Slow
      rv$isoforms_information <- isoform_info_fun(rv$isoforms[[1]], rv$num_reads_initial, rv$num_reads_post_trimming, rv$num_reads_final)
      incProgress(1/5)
      
      # df to plot
      print("df2plot")
      rv$isoform_df_plot <- isoform_df_plot_fun(rv$isoforms[[1]])
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
    if (is.null(input$iso_pixels)) iso_pixels = 20 else iso_pixels = input$iso_pixels
    if (is.null(input$histo_pixels)) histo_pixels = 250 else histo_pixels = input$histo_pixels
    
    withProgress(message = 'Plotting', value = 1, {
      
      # Abundance filter
      rv$isoform_df_plot_filtered = rv$isoform_df_plot[rv$isoform_df_plot$perc >= abundance, c("x_pos_start", "x_pos_end", "y_pos", "star_stop", "grupo")]
      
      
      # Color palete syncronization
      isoform_df_plot_and_gtf = rbind(rv$isoform_df_plot_filtered, rv$transcripts_df_to_plot)
      ordered_start_stop = unique(isoform_df_plot_and_gtf[order(isoform_df_plot_and_gtf$x_pos_start, isoform_df_plot_and_gtf$x_pos_end),"star_stop"])
      isoform_df_plot_and_gtf$star_stop = factor(isoform_df_plot_and_gtf$star_stop, levels = ordered_start_stop)
      
      library(RColorBrewer)
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
                                           "VisoQLR detected isoforms\nIsoform ID | % of reads| size", 
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
        
        rv$multiplot <- multiplot_fun1(rv$isoform_plot1, rv$break_point_freq_plot, rv$num_iso1, iso_pixels, histo_pixels)
        plotheight=rv$num_iso1*iso_pixels+histo_pixels
        rv$multiplot$height = plotheight
      
      }else{
        
        rv$isoform_plot2 = isoform_plot_fun(rv$transcripts_df_to_plot, 
                                            gsub(" \\| $", "", paste("Loaded transcripts\nTranscript ID", paste(input$transcriptinfo2display,  sep = "", collapse = " | "), sep = " | ")),
                                            pal,
                                            ordered_start_stop,
                                            showlegend = FALSE)
        rv$num_iso2 <- max(length(unique(rv$transcripts_df_to_plot[,"y_pos"])), 11)
        rv$multiplot <- multiplot_fun2(rv$isoform_plot1, rv$isoform_plot2, rv$break_point_freq_plot, rv$num_iso1, rv$num_iso2, iso_pixels, histo_pixels)
        plotheight=rv$num_iso1*iso_pixels+histo_pixels+rv$num_iso2*iso_pixels
        rv$multiplot$height = plotheight

      }
      
 
  
      print("plot1")
      output$plot1 <- renderPlotly(rv$multiplot)
  
      
      print("ui_plot")
      output$ui_plot <- renderUI({
        plotlyOutput("plot1", height = plotheight)
        })
      # cat(file=stderr(), "plot finish\n")
      print("plot finish")
      
    })
    
    })


  
  
  
  ###########################
  # Transcript data loading #
  ###########################
  output$transcripts_data_wait<-renderUI({
    req(input$inputfile)
    tagList(
      wellPanel(
        h4("Load transcripts fom file"),
        radioButtons("transcriptinputtype", "Input format", c("GTF", "GFF3"), selected="GTF"),
        fileInput("transcriptinputfile", "Transcript Input file"),
        uiOutput("transcripts_data_wait2")
      )
    )
  })
  
  output$transcripts_data_wait2<-renderUI({
    req(input$transcriptinputfile)
    tagList(
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

    rv$transcriptcolumns = colnames(rv$transcriptinput[,6:ncol(rv$transcriptinput)])

    print("Loading data")
  })


  # Gene selection
  observeEvent(ignoreInit=T, c(
    input$selectgene,
    input$transcriptinfo2display_apply,
    rv$transcriptinput
  ),{
    req(rv$transcriptinput)
    req(input$selectgene)

    withProgress(message = 'Selecting gene', value = 1, {
      rv$gene_filter_transcripts = rv$transcriptinput[rv$transcriptinput$gene == input$selectgene, ]

      rv$transcripts_df_to_plot = data.frame(
        x_pos_start = rv$gene_filter_transcripts$start,
        x_pos_end = rv$gene_filter_transcripts$end,
        y_pos = apply( rv$gene_filter_transcripts[ , c("id", input$transcriptinfo2display) , drop = F] , 1 , paste , collapse = " | " ),
        star_stop = paste(rv$gene_filter_transcripts$start, rv$gene_filter_transcripts$end, sep = "-"),
        grupo = "Loaded transcripts",
        stringsAsFactors = F)
      rv$transcripts_df_to_plot$y_pos = factor(rv$transcripts_df_to_plot$y_pos, levels = sort(unique(rv$transcripts_df_to_plot$y_pos), decreasing = T))
      

    })
    print("Gene selection")
  })




  
  
  
  
  
  
  

  ################
  # Write output #
  ################
  
  observe({
    shinyDirChoose(input, 'output', roots = c('root' = '/'))
    req(input$output)
    if (!"path" %in% names(input$output)) return()
    # print(input$output)
    print(paste(input$output$path, sep = "", collapse = "/"))
    
    rv$output_dir = paste(input$output$path, sep = "", collapse = "/")
    output$output_dir <- renderText({ rv$output_dir })
  })
  
  
  observeEvent(ignoreInit=T, c(
    input$save_apply
  ), {
    print("Saving")
    req(rv$output_dir, rv$multiplot)
    
    if (is.null(input$prefix_output)) prefix_output = gsub("(.gff3|.gff|.bed|.bed6|.bed12)$", "", input$inputfile$name, perl = T) else prefix_output = input$prefix_output
    dir_plus_prefix = paste0(rv$output_dir, "/", prefix_output)
    
    
    withProgress(message = 'Saving data', value = 0, {
      
      if (is.null(input$histo_pixels)) histo_pixels = 250 else histo_pixels = input$histo_pixels
      if (is.null(input$iso_pixels)) iso_pixels = 20 else iso_pixels = input$iso_pixels
      
      
      # Combined plot (selected isoforms)
      saveWidget(rv$multiplot, paste0(dir_plus_prefix, ".combined_plot.html"), selfcontained = T)
      incProgress(1/8)
      save_image(rv$multiplot, paste0(dir_plus_prefix, ".combined_plot.pdf"),format = "pdf", height = rv$multiplot$height)
      incProgress(1/8)
      print("Saving Combined plot (selected isoforms)")
      

      
      # Combined plot (all isoforms)
      ordered_start_stop_all_iso = unique(rv$isoform_df_plot[order(rv$isoform_df_plot$x_pos_start, rv$isoform_df_plot$x_pos_end),"star_stop"])
      
      rv$isoform_plot_all_iso <- isoform_plot_fun(df_pos = rv$isoform_df_plot, 
                                                  ylabel = "VisoQLR detected isoforms\nIsoform ID | % of reads| size", 
                                                  pal = NULL, 
                                                  start_end_levels = ordered_start_stop_all_iso, 
                                                  showlegend = TRUE)
      num_all_iso <- max(length(unique(rv$isoform_df_plot[,"y_pos"])), 11)
      rv$isoform_plot_all_iso$height = num_all_iso*iso_pixels
      saveWidget(rv$isoform_plot_all_iso, paste0(dir_plus_prefix, ".all_isoforms_plot.html"), selfcontained = T)
      incProgress(1/8)
      save_image(rv$isoform_plot_all_iso, paste0(dir_plus_prefix, ".all_isoforms_plot.pdf"),format = "pdf", height = rv$isoform_plot_all_iso$height)
      incProgress(1/8)
      print("Saving Combined plot (all isoforms)")

      
      
      
      
      # Break point information
      write.table(rv$bp_to_show[[3]], paste0(dir_plus_prefix, ".breakpoints_info.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      incProgress(1/8)
      
      
      # Exon information
      write.table(rv$exon_information, paste0(dir_plus_prefix, ".exon_info.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      incProgress(1/8)
      
      
      # Isoform information 
      write.table(rv$isoforms_information, paste0(dir_plus_prefix, ".isoform_info.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      incProgress(1/8)
      
      
      # Read clasification
      read_clasification = read_clasification_fun(all_ids = rv$data$id , no_consensus_ids = rv$defined_exons[[2]], classified_ids = rv$isoforms[[2]])
      write.table(read_clasification, paste0(dir_plus_prefix, ".read_clasification.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
      incProgress(1/8)
    })
    
  })
  
 
  # #############################
  # # Pruebas ................. #
  # #############################
  # # output$dto <-renderDT({rv$data}, editable = TRUE, rownames = FALSE)
  # 
  # proxy <- dataTableProxy('dto')
  # 
  # observeEvent(input$dto_cell_edit, {
  #   info = input$dto_cell_edit
  #   i = info$row
  #   j = info$col + 1
  #   v = info$value
  #   rv$data[i, j] <<- DT::coerceValue(v, rv$data[i, j])
  #   replaceData(proxy, rv$data, resetPaging = FALSE, rownames = FALSE)
  # })
  # 
  # #############################
  # # ................. Pruebas #
  # #############################
  # 
 
  
  # 
  # output$hover <- renderPrint({
  #   d <- event_data("plotly_hover")
  #   if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  # })
  # 
  # output$click <- renderPrint({
  #   d <- event_data("plotly_click")
  #   if (is.null(d)) "Click events appear here (double-click to clear)" else d
  # })
  # 
  # output$brushing <- renderPrint({
  #   d <- event_data("plotly_brushing")
  #   if (is.null(d)) "Brush extents appear here (double-click to clear)" else d
  # })
  # 
  
  
}



shinyApp(ui, server)






























