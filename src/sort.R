
fill_df <- function(res, df, category)
{
  for(i in 1:length(res))
  {
    df_row <- data.frame(t(res[[i]]$df[category]), row.names=NULL)
    colnames(df_row) <-  t(res[[i]]$df["CellID"])
    
    df_row <- cbind(Image = res[[i]]$filename[[1]], df_row)
    df <- merge(df, df_row, all.x = TRUE, all.y= TRUE)
  }
  return(df)
}
# Group individual experiment types, aggregate, and output
group_types <- function(df, groups)
{
  df_new = data.frame(Image="df new placeholder") 
  
  aggregated = vector("list", nrow(df))
  for(group in groups)
  {
    df_row = data.frame(Image = "df row placeholder") 
    rows = which(grepl(group, df[["Image"]], fixed = TRUE))
    print("############################")
    print(group)
    print(rows)
    
    for(row in rows)
    {
      
      df_row <- cbind(df_row,df[row,])
      df[row, "Image"] <- "used"
      }
    df_row <- df_row[, !(colnames(df_row) %in% c("Image"))]
    df_row <- df_row[ , colSums(is.na(df_row)) == 0]
    colnames(df_row)<- 1:ncol(df_row)
    df_row["Image"] = group
    
    df_new <- merge(df_new, df_row, all.x= TRUE, all.y = TRUE)
  }
  df_new <- df_new[-1,]
  print(df)
  return(df_new)
}

# Sort all data from computed images.  Aggregate, output
sort_data <- function(res, groups,title)
{
  df_pm_mpi <- data.frame(Image = character(),stringsAsFactors=FALSE)
  df_vac_mpi <- data.frame(Image = character(),stringsAsFactors=FALSE)
  df_pm_vac_ratio <- data.frame(Image = character(),stringsAsFactors=FALSE)
  
 
  
  df_pm_mpi <- fill_df(res,df_pm_mpi, "cell_mpi")
  df_vac_mpi <- fill_df(res,df_vac_mpi, "vac_mpi")
  df_pm_vac_ratio <- fill_df(res,df_pm_vac_ratio, "PM_vac_ratio")

  write.csv(t(df_pm_mpi), paste0("FinalOutput/Aggregated Spreadsheets/", title, "_pm_mpi_all.csv"), na = "", row.names = FALSE)
  write.csv(t(df_vac_mpi), paste0("FinalOutput/Aggregated Spreadsheets/", title, "_vac_mpi_all.csv"), na = "", row.names = FALSE)
  write.csv(t(df_pm_vac_ratio), paste0("FinalOutput/Aggregated Spreadsheets/", title, "_pm_vac_ratio_all.csv"), na = "", row.names = FALSE)
  
  df_pm_mpi_grouped <- group_types(df_pm_mpi, groups)
  df_vac_mpi_grouped <- group_types(df_vac_mpi, groups)
  df_pm_vac_ratio_grouped <- group_types(df_pm_vac_ratio, groups)
  
  
  write.csv(t(df_pm_mpi_grouped), paste0("FinalOutput/Aggregated Spreadsheets/", title, "_pm_mpi_grouped.csv"), na = "", row.names = FALSE)
  write.csv(t(df_vac_mpi_grouped), paste0("FinalOutput/Aggregated Spreadsheets/", title, "_vac_mpi_grouped.csv"), na = "", row.names = FALSE)
  write.csv(t(df_pm_vac_ratio_grouped), paste0("FinalOutput/Aggregated Spreadsheets/", title, "_pm_vac_ratio_grouped.csv"), na = "", row.names = FALSE)
  

  
  print("Finished sorting data.")
}









