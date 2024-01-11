

add_info_repname = function(sub_df){
  
  library("stringr")
  
  # detecting the replictes name vector
  rep_vec = sub_df$Replicate
  
  # split the information in filenames
  y = str_split( rep_vec , "_")
  
  # edit the output to dataframe
  y_df = t(data.frame(y))
  y_df2 = data.frame(y_df)
  colnames(y_df2) = c( "date","Cell.line","Drug","rep","loc")
  
  # add the info to the original table
  sub_df2 = cbind( sub_df, y_df2)
  
  # editing the column types and infos
  sub_df2 = sub_df2[ , -c(2,3,4) ]
  sub_df2$Quantification = str_replace( sub_df2$Quantification, "ng/ml", "" )
  sub_df2$Quantification = as.numeric(sub_df2$Quantification)
  
  return(sub_df2)
}


concat_skyline_tab = function( cell_list, CC_num_str ) {
  
  # create output container
  total_tab = data.frame()
  
  for (cell in cell_list){
   
    # arrange the filenames
    file_name = paste( cell, "_", CC_num_str, ".csv", sep="")
    
    # loading files
    cell_file = read.csv(file_name)
    
    # substracting the important number
    cell_quan = cell_file[ which(cell_file$Quantification != "#N/A") , ]
    cell_sub = cell_quan[ which(cell_quan$Sample.Type == "Unknown") , ]
    
    # adding to total_table
    total_tab = rbind(total_tab, cell_sub)
  }
  
  return(total_tab)
} 



calc_mol_in_cell = function(quant, drug_mol_weight) {
  
  # calculate mass in cell - 2.5 comes from 25 [= dilution factor] * 100uL [extraction volume] to ml [since quant is a concentration in ng/ml]
  mass_in_cell = 2.5 * quant
  
  # convert to gram from nanogram
  mass_in_cell_gram = mass_in_cell / 10^9
  
  # calculate moles [n=m/Mw]
  moles = mass_in_cell_gram / drug_mol_weight
  
  return(moles)
}



calc_mol_conc = function(moles_in_cell, avg_count, avg_diam){
  
  # calculate a single cell volumes
  radius = avg_diam / 2
  vol_per_cell = (4/3) * pi * radius^3
  
  # convert micro-cubic to Litre
  vol_per_cell_litre = vol_per_cell / 10^15
  
  # calculate total volume of cells per well
  total_vol = vol_per_cell_litre * avg_count
  
  # so molar concentration is ...[c=n/v]
  mol_conc = moles_in_cell / total_vol
  
  
  return(mol_conc)
}


