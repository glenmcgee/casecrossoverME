## set correct path to results
mypath="~/Results/"
setwd(mypath)

## function to combine data
combine_dat <- function(ptrn="sim_par0_results_iter*",   ## ptrn should be the common prefix INCLUDING _iter, followed by asterisk 
                        orderby="X.iter"){ ## name of a column by which to order the results (ie iteration number) 
                                           ## set to NA to not reorder
  
  ## get list of files
  txt_files_ls = list.files(path=mypath, pattern=ptrn) 
  
  ## load files
  txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T)})
  
  ## combine
  combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
  
  ## save
  filename=substr(ptrn,1,nchar(ptrn)-6) ## remove the "_iter*"
  if(!is.na(orderby)){ ## reorder by a columnname
    combined_df <- combined_df[order(combined_df[,orderby]),]
  }
  write.table(combined_df,file=paste0("Full/",filename,".txt")) ## Added "Full/" to put files in a separate folder, making things a little tidier
}

# ## Example
# combine_dat("sim_par0_results_iter*")
