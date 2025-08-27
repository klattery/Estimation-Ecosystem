#########################################################
#     Based on Estimation Ecosystem         
#     Created Nov 2021 for AWS Portal                                 
#     Kevin Lattery                                                       
#     Run Everything to create environments                  
#########################################################

# Create 3 Environments to Store Functions
# Attaching an environment multiple times creates duplicates
if (exists("env_code")) rm(env_code) 
if (exists("env_modfun"))  rm(env_modfun)
if (exists("env_eb")) rm(env_eb)
if (exists("env_stan")) rm(env_stan)
env_code <- new.env(parent = emptyenv())
env_modfun <- new.env(parent = emptyenv())
env_eb <- new.env(parent = emptyenv())
env_stan <- new.env(parent = emptyenv())

###############  Coding Functions Environment ################ 

env_code$read_csv_rds <- function(dir_data, data_file, as_name = NULL){
  if (is.null(data_file)){
    result <- NULL
  } else{
    if (!file.exists(file.path(dir_data, data_file))){
      message(paste0("ERROR: CANNOT FIND YOUR DATA FILE: ", data_file)) 
      result <- NULL
    } else {
      ftype <- toupper(substr(data_file,nchar(data_file)-2, nchar(data_file)))
      if (ftype %in% c("CSV", "RDS")){
        if (ftype == "CSV") {result <- read.csv(file.path(dir_data, data_file), as.is=TRUE, check.names = FALSE)}
        if (ftype == "RDS") {result <- readRDS(file.path(dir_data, data_file))}
        f_inputs <- as.list(match.call(expand.dots = FALSE))
        message(paste0("\nREAD ", eval(f_inputs[[3]]), " INTO R ", as_name)) # [[1]] is name of function
      } else {
        message("ERROR: DATA MUST BE .CSV OR .RDS FILE")
        result <- NULL
      } 
    } # End normal data file 
  } # End null
  return(result)
}

env_code$setup_cores <- function(ncores){
  # R multi-threading
  if (exists("k_multi_core", envir = globalenv())){
    stopCluster(.GlobalEnv$k_multi_core)
    rm("k_multi_core", envir = globalenv())
  }
  if (ncores >= 1){
    .GlobalEnv$k_multi_core <- makeCluster(min(ncores,detectCores()))
    registerDoParallel(.GlobalEnv$k_multi_core)
    message(paste0("\nUsing ", ncores, " cores for multi-threading.  Set lower if RAM problems"))
  } else {
    message("\nRemoving parallel threads in R defined by k_multi_core")
  }
}

env_code$catcode <- function(kdata, vname, codetype = 3, varout = NULL, reflev = NULL, setna = 0, priorcode = c(0, NA), paircon = NULL) {
  #codetype 1 = indicator, 2= dummy, 3 = effects
  #reflev of NULL defaults to most shown level
  colvec <- kdata[[vname]] 
  if (is.null(varout)) varout <- vname
  colvec[colvec == priorcode[1]] <- priorcode[2]
  na_vals <- is.na(colvec)
  kmax <- max(colvec, na.rm = TRUE)
  colvec[na_vals] <- kmax # temp set to max
  labels_in <- sort(unique(colvec))
  newval <- match(colvec,labels_in)
  varnames <- paste(varout,labels_in, sep = "_")
  numlevs <- length(labels_in)
  if (numlevs == 1){
    outcode <- as.matrix(colvec) # no coding for 1 level
    code_matrix <- as.matrix(1)
    colnames(code_matrix) <- varnames
  } 
  if (numlevs >= 2) {
    if (is.null(paircon)){
      code_matrix <- diag(numlevs)
      colnames(code_matrix) <- varnames
      if (is.null(reflev)){ # set to level with highest shows
        shows <- colSums(code_matrix[newval,TRUE,drop = FALSE])
        shows[length(shows)] <- shows[length(shows)] - sum(na_vals) # last value has NAs too
        reflev <- match(max(shows), shows)
      } 
      if (codetype %in% c(2, 3))code_matrix <- as.matrix(code_matrix[, -reflev, drop = FALSE]) # ref lev col dropped
      if (codetype == 3) code_matrix[reflev,] <- -1
    } else { # we have constraints
      paircon[,2] <- match(paircon[,2], labels_in) # recode to seq levels
      paircon[,3] <- match(paircon[,3], labels_in) # recode to seq levels
      paircon <- paircon[(rowSums(is.na(paircon)) == 0),,drop = FALSE]  # drop rows with NA
      cat_con <- code_cat_wcon(paircon,numlevs)
      code_matrix <- cat_con$code_matrix
      con_vec <- cat_con$con
      reflev <- cat_con$reflev
      pairs_add <- cat_con$pairs_add
      colnames(code_matrix) <- varnames[-reflev]
    }  
    outcode <- (code_matrix[newval,TRUE,drop = FALSE])
  } # end numlev >=2
  if (is.null(paircon)){
    con_vec <- rep(0, ncol(code_matrix))
    pairs_add <- NULL
  }
  outcode[na_vals,] <- setna
  # Now get priors
  if (is.null(paircon)){
    if (numlevs <= 2){
      if (codetype == 1) prior <- diag(numlevs)    
      if (codetype > 1) prior <- matrix(1)    
    } 
    if (numlevs >= 3) {
      if (codetype == 1){ # ind
        prior <- diag(numlevs)
      }
      if (codetype == 2){ # dummy
        prior <- matrix(1, nrow = numlevs-1, ncol = numlevs-1) # off-diagonal
        diag(prior) <- 2
      }
      if (codetype == 3){ #effects
        prior <- matrix(-1/numlevs, nrow = numlevs-1, ncol = numlevs-1) # off-diagonal
        diag(prior) <- (numlevs -1)/numlevs
      }
    }
  } else prior <- diag(ncol(code_matrix))
  newval[na_vals] <- setna 
  return(list(outcode = outcode, code_matrix = code_matrix, levels = newval, levels_in = labels_in, con_sign = con_vec, vnames = varnames, reflev = reflev, prior = prior, pairs_add = pairs_add))
}

#### Kees new Code ###
env_code$remove_implicits<-function(constraints){
  if (nrow(constraints) == 0){
    new_constraints <- constraints
  } else {
    result <- NULL
    
    constraints <- data.frame(constraints) # Tibble does not work
    atts<-unique(constraints[,1])
    new_constraints<-c()
    
    for(att in atts){
      attcons<-constraints[constraints[,1]==att,-1]
      orig_prohibs<-create_zero_ones(attcons)
      new_prohibs<-check_implicits(orig_prohibs)
      new_constraints<-rbind(new_constraints,cbind(att,new_prohibs))
    }
    
  }
  colnames(new_constraints)<-colnames(constraints)
  new_constraints<-data.frame(new_constraints)
  return(new_constraints)
}

env_code$create_zero_ones<-function(attcons){
  
  prohibs<-diag(max(attcons))
  for (con in 1:nrow(attcons)){
    prohibs[attcons[con,1],attcons[con,2]]=1
  }
  prevsum=0
  newsum=sum(prohibs)
  while(newsum!=prevsum){
    prevsum=newsum
    prohibs<-prohibs %*% prohibs
    prohibs[prohibs>0]<-1
    newsum=sum(prohibs)
  }  
  diag(prohibs)<-0
  return (prohibs)
}
env_code$zero_one_to_lev<-function(prohibs){
  nlev<-ncol(prohibs)
  result<-c()
  for (rij in 1:nlev){
    for (col in 1:nlev){
      if (prohibs[rij,col]==1){
        result<-rbind(result,c(rij,col))
      }
    }
  }
  return(result)
}

env_code$check_implicits<-function(orig_prohibs){
  
  nlev<-ncol(orig_prohibs)
  new_prohibs<-orig_prohibs
  
  for (rij in 1:nlev){
    for (col in 1:nlev){
      if (orig_prohibs[rij, col] == 1){
        rij2 = col
        vlag <- TRUE
        if (sum(orig_prohibs[rij2,])==0) vlag<-FALSE
        for (col2 in 1:nlev){
          if (orig_prohibs[rij2, col2] == 1 && orig_prohibs[rij, col2] == 0){
            vlag <-FALSE
            break
          }
        }
        if (vlag){
          new_prohibs[rij, orig_prohibs[rij2, ]==1] <- 0
        }
      }
    }
  }
  return(zero_one_to_lev(new_prohibs))
}
##### End Kees new code ############

env_code$ordinal_chain <-function(constraints){
  result <- list()
  nrofcons <- nrow(constraints)
  i <- 0
  nrDone <- 0
  while (nrDone < nrofcons){
    startFound <- FALSE
    for (lvl in unique(constraints[,2])){
      if(sum(constraints[,2]==lvl) > 0 && sum(constraints[,3]==lvl) == 0) {
        i <- i + 1
        A <- lvl
        result[[i]]<-c(A)
        startFound <- TRUE
        break
      }
    }
    if(startFound == FALSE){
      browser(expr = TRUE)
      cat ("loop exist in constraints")
      break
    }
    nextFound <- TRUE
    while (nextFound){
      nextFound <- FALSE
      for (c in 1:nrofcons){
        if(constraints[c, 2] == A){ 
          B <- constraints[c,3]
          result[[i]]<-append(B,result[[i]],after = length(result[[i]]))
          A <- B
          constraints[c,2]<-0
          constraints[c,3]<-0
          nrDone <- nrDone + 1
          nextFound <- TRUE
          break
        }
      }
    }
  }
  #reorder so first element is longest
  result <- result[order(sapply(result,function(x) -length(x)))]
  return(result)
}

# Nominal with constraints
env_code$code_cat_wcon <-function(constraints, numlevs){
  constraint_chain <- ordinal_chain(data.frame(constraints)) # create list of chains
  missing <- NULL
  con <- rep(0,numlevs) # storage for constraint
  X <- matrix(0,numlevs,numlevs) # X will store the coded matrix
  for (constrings in 1:length(constraint_chain)){
    for (c in 1:(length(constraint_chain[[constrings]])-1)){
      A <- constraint_chain[[constrings]][c]
      B <- constraint_chain[[constrings]][c+1]
      diagA <- X[A,A]  # diagonal element A
      sumB <- sum(X[B,]!=0) # Num of coded elements in row B
      if (diagA==0 && sumB==0){
        X[A,A] <- 1
        X[B,A] <- 1
        X[B,B] <- 1
        con[B] <- 1
      }
      if (diagA!=0 && sumB==0){
        X[B,] <- X[A,]
        X[B,B] <- 1
        con[B] <- 1
      }
      if (diagA==0 && sumB!=0){
        X[A,] <- X[B,]
        X[A,A] <- -1
        con[A] <- 1
      }
      if (diagA!=0 && sumB!=0){        
        missing <- rbind(missing,c(B,A))
      }
    }
  }
  for (lev in 1:numlevs){
    if (X[lev,lev]==0) X[lev,lev] <- 1
  }
  kolsum <- colSums(X)
  kolsum[(con != 0)] <- -99 # columns with constraints (!=0) are keepers
  reflev <- match(max(kolsum),kolsum)
  X <- X[,-reflev, drop = FALSE]
  con <- con[-reflev]
  return(list(code_matrix=X,con=con,pairs_add=missing, reflev = reflev))
}


env_code$usercode1 <- function(kdata, vname, varout = NULL, con_sign = 0){
  # Single variable usercode
  if (is.null(varout)) varout <- vname
  con_sign[is.na(con_sign)] <- 0
  outcode <- as.matrix(kdata[[vname]])
  colnames(outcode) <- varout
  code_matrix <- diag(1)
  colnames(code_matrix) <- varout
  prior <- code_matrix
  return(list(outcode = outcode, code_matrix = code_matrix, con_sign = con_sign, vnames = varout, prior = prior))
}

env_code$usercode <- function(kdata, kcol, varout = NULL, con_sign = 0){
  if (is.null(varout)){
    varout <- colnames(kdata)[kcol]
    if (is.null(varout)) varout <- "user_"
  }
  outcode <- as.matrix(kdata[,kcol])
  numcol <- ncol(outcode)
  if ((length(varout) == 1) & (numcol > 1)){
    varout <- paste0(varout, 1:numcol)
  }
  if (!(length(varout) == numcol)){
    varout <- paste0("user_", 1:numcol)
  }
  colnames(outcode) <- varout
  return(list(outcode = outcode, code_matrix = diag(ncol(outcode)), con_sign = rep(con_sign,numcol), vnames = varout, prior = diag(ncol(outcode))))
}

env_code$ordmatrix <- function(num_levels) {
  #num_levels must be >= 2, 2 levels is one variable coded -1,1 or 0,1
  negval <- (num_levels - 1) * -1
  ord_matrix <- matrix(rep(1:(num_levels - 1), each = num_levels), nrow = num_levels)
  ord_matrix[1,] <- (negval:-1)
  if (num_levels > 2){
    for (i in 2:(num_levels - 1)) {
      negval <- negval + 1
      ord_matrix[i,] <- c(1:(i - 1), (negval:-1))
    }
  }
  maxsum <- sum(ord_matrix[num_levels,])
  return(ord_matrix / num_levels)
}

env_code$ordmatrix2 <- function(num_levels) {
  ord_matrix <- matrix(0,num_levels-1, num_levels-1)
  ord_matrix <- diag(num_levels-1)
  ord_matrix[lower.tri(ord_matrix)] <- 1
  ord_matrix <- rbind(rep(0, ncol(ord_matrix)), ord_matrix)
  return(ord_matrix)
}

env_code$ordcode <- function(kdata, vname, cut_pts = NULL, thermcode = TRUE, varout = NULL, setna = 0, con_sign = 0) {
  # xvec must be vector
  # cut_pts must be sequential vector from low to high.  if null then all values except c(0,NA)
  # NA and values outside cut_pts are set to "setna", default = 0
  # varout is prefix for varout_cut
  # uses function ordmatrix
  if (is.null(varout)) varout <- vname
  con_sign[is.na(con_sign)] <- 0
  xvec <- kdata[[vname]]
  if (is.null(cut_pts)){
    na_vals <- is.na(xvec) | xvec == 0
    xvec[na_vals] <- NA
    kmax <- max(xvec, na.rm = TRUE)
    xvec[na_vals] <- kmax # temp set to max
    cut_pts <- sort(unique(xvec)) # Unique values
    newval <- match(xvec,cut_pts)
    if (thermcode){
      code_matrix <- ordmatrix2(length(cut_pts))
    } else code_matrix <- ordmatrix(length(cut_pts))
    outcode <- (code_matrix[newval,TRUE,drop = FALSE])
    outcode[na_vals,] <- setna
    newval[na_vals] <- setna
    vnames <- unlist(lapply(2:length(cut_pts), function(i) paste0(varout, "_",cut_pts[i-1],"_", cut_pts[i])))
    colnames(code_matrix) <- colnames(outcode) <- vnames
    varnames <- paste0(varout, "_", cut_pts)
    return(list(outcode = outcode, code_matrix = code_matrix, levels = newval, con_sign = rep(con_sign, ncol(code_matrix)), vnames = varnames, prior = diag(ncol(code_matrix))))
  } else{
    bad <- xvec < min(cut_pts) | xvec > max(cut_pts) | is.na(xvec)
    xvec[bad] <- min(cut_pts) # temp set to min
    high <- sapply(cut_pts, function(x) xvec <= x)
    high_col <- apply(high, 1, function(x) match(TRUE, x))
    low_col <- pmax(high_col - 1,1)
    low_pt <- cut_pts[low_col] 
    high_pt <- cut_pts[high_col] 
    dist1 <- (high_pt - xvec)/(high_pt - low_pt)
    dist1[is.na(dist1)] <- 0
    dist2 <- 1 - dist1
    mat0 <- matrix(0, length(xvec), length(cut_pts))
    mat_low <- mat0
    mat_high <- mat0
    mat_low[cbind(1:length(xvec),low_col)] <- dist1 
    mat_high[cbind(1:length(xvec),high_col)] <- dist2 
    rowcode <- mat_low + mat_high
    if (thermcode){
      code_matrix <- ordmatrix2(length(cut_pts))
    } else code_matrix <- ordmatrix(length(cut_pts))
    ordcode <- round(rowcode %*% code_matrix, 5)
    ordcode[bad] <- setna # reset initial NA (default 0)
    vnames <- unlist(lapply(2:length(cut_pts), function(i) paste0(varout, "_",cut_pts[i-1],"_", cut_pts[i])))
    colnames(ordcode) <- vnames
    colnames(code_matrix) <- vnames
    varnames <- paste0(varout, "_", cut_pts)
    return(list(outcode = ordcode, code_matrix = code_matrix, con_sign = rep(con_sign, ncol(code_matrix)), vnames = varnames, prior = diag(ncol(code_matrix))))
  }
}

env_code$check_atts_constraints <- function(data_in, att_coding, constraints){
  # Check consistency of specified attributes and constraints with data
  result <- TRUE
  check_att <- att_coding[,1,drop = TRUE] %in% colnames(data_in)
  if (min(check_att) == 0){
    cat("Below are attributes specified that are not in your data.  Please fix.\n")
    print(att_coding[!check_att,1])
    result <- FALSE
  }
  if (nrow(constraints) > 0){
    check_match <- constraints[,1, drop = TRUE] %in% att_coding[,1, drop = TRUE]
    if (min(check_match) == 0){
      cat("You specified constraints that do no match any attribute in attribute file.  Please fix.")
      print(constraints[!check_match,1])
      result <- FALSE
    }
  }
  return(result)
}

env_code$indcode_spec_files <- function(data_in, att_coding, constraints,
                                        add_none = FALSE, dep_col = NULL){
  if (is.null(constraints)) constraints <- data.frame(matrix(ncol = 3, nrow = 0))
  constraints <- remove_implicits(constraints) 
  if (check_atts_constraints(data_in, att_coding,constraints)){
    catcode_types <- c("INDICATOR", "DUMMY","EFFECT","EFFECTS","CATEGORICAL","NOMINAL")
    indcode_spec <- vector("list", length = nrow(att_coding))
    
    for (i in 1:nrow(att_coding)){
      att_name <- att_coding[i,1,drop = TRUE]
      att_type <- toupper(att_coding[i,2,drop = TRUE]) # UPPERCASE
      if (att_type %in% catcode_types){ # CATEGORICAL
        codetype <- min(3, match(att_type, catcode_types))
        if (sum(constraints[,1,drop = TRUE] == att_name) == 0) {
          indcode_spec[[i]] <- catcode(data_in, att_name, codetype) # No constraints
        } else{
          indcode_spec[[i]] <- catcode(data_in, att_name, codetype, paircon = constraints[constraints[,1,drop=TRUE] == att_name,])
        }    
      }
      if (att_type == "ORDINAL"){
        indcode_spec[[i]] <- ordcode(data_in, att_name, thermcode = FALSE, con_sign = att_coding[i,3,drop = TRUE])
      }
      if (att_type == "THERMOMETER"){
        indcode_spec[[i]] <- ordcode(data_in, att_name, thermcode = TRUE, con_sign = att_coding[i,3,drop = TRUE])
      }      
      if (att_type == "USERSPECIFIED"){
        indcode_spec[[i]] <- usercode1(data_in, att_name, con_sign = att_coding[i,3,drop = TRUE])
      }
    }  
    indcode_spec <- setNames(indcode_spec,att_coding[,1,drop = TRUE])
    if(add_none){
      None <- data.frame(None = 0 + apply(data_in[,colnames(data_in) %in% 
                                                          att_coding[,1], drop = FALSE], 1, function(x_row){all(x_row == 0)}))
      if (sum(None$None) >0){
        indcode_spec$None <- usercode1(None, "None")
        indcode_spec$None$message <- c("\nAdded variable None, chosen ",
                                       sum(None$None[data_in[,dep_col]> 0]),
                                       " times out of ", sum(None$None), " shows\n")
      } 
    } 
  } else indcode_spec <- NULL 
  return(indcode_spec)
}

env_code$list_to_matrix <- function(klist){
  # Create one matrix with blocks of matrices from list and 0 fills
  row_sizes <- sapply(klist, nrow)
  col_sizes <- sapply(klist, ncol)
  row_end <- cumsum(row_sizes)
  col_end <- cumsum(col_sizes)
  row_start <- c(1, row_end[-length(row_end)] + 1)
  col_start <- c(1, col_end[-length(col_end)] + 1)
  result <- matrix(0, nrow = sum(row_sizes), ncol = sum(col_sizes))
  for (i in 1:length(klist)){
    result[(row_start[i]:row_end[i]), (col_start[i]:col_end[i])] <- klist[[i]]   
  }
  return(result)  
} 

env_code$make_codefiles <- function(indcode_spec, control_code = .GlobalEnv$control_code){
  # Converts list of codes to matrices:  # code_master, indcode, indprior
  result <- list()
  indcode_list <- indcode_spec[sapply(indcode_spec, length) > 0] # remove NULL elements
  att_names <- names(indcode_list)
  names(indcode_list) <- NULL # to avoid adding to vnames
  
  # Create two files: coded and uncoded levels, and a combined version for R
  ind_coded <- NULL # coded data
  ind_levels <- NULL # uncoded (will be coded in Stan)
  ind_levels_names <- NULL
  ind_coded <- matrix(0, nrow(indcode_list[[1]]$outcode),0)
  ind_levels <- matrix(0, nrow(indcode_list[[1]]$outcode),0)
  col_beg <- 1 # code_matrix col
  row_beg <- 1 # code_matrix row
  coded_beg <- 1 # column of coded data
  code_blocks <- matrix(0, nrow = length(indcode_list), ncol = 5)
  # col_beg, col_end, row_beg, row_end, coded_beg
  for (i in 1:length(indcode_list)){
    x <- indcode_list[[i]]
    col_end <- col_beg + ncol(x$code_matrix) - 1
    row_end <- row_beg + nrow(x$code_matrix) - 1
    if (is.null(x$levels)){
      ind_coded <- cbind(ind_coded, x$outcode)
      code_blocks[i,] <- c(col_beg, col_end, row_beg, row_end, coded_beg)
      coded_beg <- coded_beg + ncol(x$code_matrix) # could also use x$outcode
    } else {
      ind_levels_names <- c(ind_levels_names, att_names[i])
      ind_levels <- cbind(ind_levels, x$levels)
      code_blocks[i,] <- c(col_beg, col_end, row_beg, row_end, 0) # 0 = not in coded file
    }
    col_beg <- col_end + 1 # next column of code_matrix
    row_beg <- row_end + 1 # next row of code_matrix
  }
  colnames(ind_levels) <- ind_levels_names
  
  result$indcode <- do.call(cbind, lapply(indcode_list, function(x) x$outcode)) # coded variables 
  code_master <- list_to_matrix(lapply(indcode_list, function(x) x$code_matrix))
  colnames(code_master) <- colnames(result$indcode) 
  rownames(code_master) <- do.call("c", lapply(indcode_list, function(x) x$vnames))
  
  # Recode cat_var x thermometer if checked
  cat_var <- control_code$recode_therm$cat_var
  therm_var <- control_code$recode_therm$therm_var    
  if (cat_var != "" && therm_var != ""){
    att_num <- match(cat_var, att_names, nomatch = 0)
    goodinput <- TRUE
    if (att_num == 0) {
      cat(paste0("*** Recode Slopes ERROR: Could not find categorical variable ", cat_var, " in specs ***\n"))
      goodinput <- FALSE
    } 
    ind_level_col <- match(toupper(cat_var), toupper(colnames(ind_levels)), nomatch = 0)
    if (ind_level_col == 0) {
      cat(paste0("*** Recode Slopes ERROR: Could not find categorical variable ", cat_var, " in data ***\n"))
      goodinput <- FALSE
    }
    therm_data <- ind_coded[, grepl(therm_var, colnames(ind_coded), ignore.case = T)]
    if (ncol(therm_data) == 0){
      cat(paste0("*** Recode Slopes ERROR: Could not find thermometer variable: ", therm_var, " in data ***\n"))
      goodinput <- FALSE
    }
    therm_cols_change <- match(colnames(therm_data), colnames(code_master))
    if (sum(is.na(therm_cols_change))>0){
      cat(paste0("*** Recode Slopes ERROR: Could not match all data thermometer variables ", colnames(therm_data), "with code master ***\n"))
      goodinput <- FALSE
    }
    if (goodinput){
      cat("Recoding these thermometer variables\n")
      print(colnames(therm_data))
      cat(paste0("Based on: ", cat_var, "\n\n"))
      att_specs <- code_blocks[att_num,]
      cat_values <- ind_levels[, ind_level_col]
      #cat_values_u <- indcode_spec[[att_num]]$levels_in
      cat_values_u <- sort(unique(cat_values))
      beg_row <- att_specs[3]
      row_change <- beg_row
      for (i in cat_values_u){
        if (i >0){
          cat_filter <- (cat_values == i)
          min_cols <- apply(therm_data[cat_filter,], 2 , min)
          max_cols <- apply(therm_data[cat_filter,], 2 , max)
          recodes <- -min_cols * (min_cols == max_cols) * (min_cols != 0) # 0 unless min_cols != 0 and min = max
          # Change code_master
          code_master[row_change,therm_cols_change] <- recodes # Changing code_master
          row_change <- row_change+1    
        } # end if
      } # end for
    } # end if goodinput
  } # end code_master recode
  
  # Create indcode just like we do in Stan
  result$indcode[] <- 0
  n_atts <- nrow(code_blocks)
  lev_col <- 1
  ## Coding of attributes with levels
  for (i in 1:n_atts){
    if (code_blocks[i,5] > 0){ ## already coded, just copy
      col_beg <- code_blocks[i, 1];
      col_end <- code_blocks[i, 2];
      col_coded <- code_blocks[i, 5]; #  beginning of coded column
      result$indcode[,col_beg:col_end] <- result$indcode[,col_beg:col_end] + ind_coded[,col_coded:(col_coded + col_end - col_beg)];
    } else {
      row_beg <- code_blocks[i,3];
      row_end <- code_blocks[i,4];
      #int nrows = row_end - row_beg +1;
      code_att <-  code_master[row_beg:row_end, ]; ## code for this attribute
      for (n in 1:nrow(result$indcode)){
        lev <- ind_levels[n,lev_col];
        if (lev >= 1) result$indcode[n,] = result$indcode[n,] + code_att[lev,];
      } # end 1:N loop coding
      lev_col <- lev_col + 1;
    } # end if  
  } # end for (att coding)
  # End creating indcode
  
  result$ind_coded <- ind_coded
  result$ind_levels <- ind_levels
  result$code_blocks <- code_blocks
  result$con_sign <- do.call("c", lapply(indcode_list, function(x) x$con_sign))
  result$code_master <- code_master
  result$indprior <- list_to_matrix(lapply(indcode_list, function(x) x$prior))
  
  # This is for the pair constraints
  pair_m <- lapply(indcode_list, function(one_list) {
    if (is.null(one_list$pairs_add)){
      result <- t(as.matrix(rep(0, ncol(one_list$code_matrix))))
    } else {
      code <- one_list$code_matrix
      pairs <- one_list$pairs_add
      result <- do.call(rbind, lapply(1:nrow(pairs), function(krow){
        return(code[pairs[krow,1],] - code[pairs[krow,2],])
      } 
      ))
    }
    return(result)
  })
  pair_m2 <- list_to_matrix(pair_m)
  colnames(pair_m2) <- colnames(result$code_master)
  pair_m2 <- pair_m2[(rowSums(pair_m2 != 0) > 0),,drop = FALSE] # keep as matrix
  result$con_matrix <- pair_m2
  
  # Find initial x0 satisfying constraints
  if (is.null(nrow(result$paircon_matrix))){
    result$x0 <- result$con_sign/10
  } else {
    con_matrix <- diag(result$con_sign)
    con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], result$paircon_matrix)
    hinge <- function(xvec, kmatrix){
      kprod <- kmatrix %*% xvec
      result <- 1000 * sum(log(1+ exp(-100 * kprod)))
      return(result)
    }
    x0 <- optim(par = rep(0, ncol(con_matrix)), fn = hinge, method = "L-BFGS-B", lower = -2, upper = 2, kmatrix = con_matrix)$par/5
    check <- con_matrix %*% x0
    if (min(check) <= 0) cat("Could not find initial x0 that satisfies constraints. \nSet indcode_list$x0 manually if using EB")
    result$x0 <- x0
  }
  return(result)
}


env_code$save_codemastercon <- function(indcode_list, dir, out_prefix){
  # Makes a code_master + constraints file and saves it
  # export code matrix, sign, and constraints in 1 file
  df1 <- rbind(indcode_list$code_master, indcode_list$con_sign, indcode_list$con_matrix)
  att_labels <- c(rownames(indcode_list$code_master), "con_sign")
  if (nrow(df1) > length(att_labels)){
    more_names <- rep("", nrow(df1) - length(att_labels))
    more_names[1] <- "con_pair(s)"
    att_labels <- c(att_labels, more_names)
  } 
  write.table(cbind(att_labels, df1), file = file.path(dir, paste0(out_prefix,"_code_master_andcon.csv")), sep = ",", na = ".", row.names = FALSE)
}

env_code$code_covariates <- function(cov_in, cov_coding, resp_id){
  cov_code_spec <- indcode_spec_files(cov_in, cov_coding, NULL)
  cov_code <- do.call(cbind, lapply(cov_code_spec[sapply(cov_code_spec, length) > 0],
                                    function(x) x$outcode)) # coded variables
  code_master <- list_to_matrix(lapply(cov_code_spec, function(x) x$code_matrix))
  colnames(code_master) <- colnames(cov_code) 
  rownames(code_master) <- do.call("c", lapply(cov_code_spec, function(x) x$vnames))
  
  row_match <- match(resp_id, cov_in[,1])
  no_cov <- is.na(row_match) # Respondents in data with no matching covariate
  row_match[is.na(row_match)] <- 1
  result <- as.matrix(cov_code[row_match,])
  if (sum(no_cov) >0){
    message(paste0(sum(no_cov), " respondents had no covariates.  Coded to 0"))
    result[no_cov,] <- 0 # set bad to 0
  } else message(paste0("All respondents matched in covariates file.  ", ncol(result), " coded parameters"))
  return(list(output = result, cov_code_master = code_master))
}

env_code$make_wts <- function(cov_in, resp_id){
  # cov_in has id, weight for first 2 columns
  wts <- cov_in[,2,drop = TRUE]
  result <- wts[match(resp_id, cov_in[,1])]
  if (sum(is.na(result)) > 0){
    cat(paste0(sum(is.na(result)), " respondents in data_stan have no matching weights.  These weights are set to 1"))
  }
  result[is.na(result)] <- 1
  cat("Weights used from column 2 of covariates file:\n")
  print(summary(result))
  return(result)
}

env_code$make_con <- function(con_specs, code_master, x0_try){
  # col_pos & col_neg are VECTORS from COLUMNS of code_master
  # row_rel is LIST of pairs from ROWS of code_master
  diaguse <- diag(ncol(code_master))
  result <- vector("list", length = 3)
  col_pos <- con_specs$col_pos
  col_neg <- con_specs$col_neg
  row_rel <- con_specs$row_rel
  x0 <- x0_try
  if (length(col_pos) > 0){
    result[[1]] <- diaguse[col_pos,]
    message("Positive Variables - Coded")
    print(colnames(code_master)[col_pos])
    x0[col_pos] <- pmax(.001, x0[col_pos]) 
  } 
  if (length(col_neg) > 0){
    result[[2]] <- diaguse[col_neg,] * -1
    message("Negative Variables - Coded")
    print(colnames(code_master)[col_neg])
    x0[col_neg] <- pmin(-.001, x0[col_neg])  
  } 
  if (length(row_rel) > 0){
    vec0 <- rep(0,ncol(code_master))
    message("Relative Non-Coded 1st >= 2nd")
    code_rel <- lapply(row_rel, function(pair){
      vec_new <- code_master[pair[1],] - code_master[pair[2],]
      print(rownames(code_master)[pair])
      return(vec_new)
    })
    result[[3]] <- do.call(rbind,code_rel)  
  }
  if (is.null(col_pos) & is.null(col_neg) & is.null(row_rel)){
    result_mat <- con_trivial(ncol(code_master))  
  } else {result_mat <- cbind(0, do.call(rbind, result))}
  colnames(result_mat) <- c("C", colnames(code_master))
  if (nrow(result_mat) == 1) result_mat <- rbind(result_mat, result_mat)
  return(list(constrain = result_mat, x0 = x0))  
}

env_code$check_beta <- function(beta, constrain){
  #checks whether vector beta is within constraints
  good <- FALSE
  Diff <- (constrain[,-1] %*% beta) - constrain[,1]
  result <- as.data.frame(Diff)
  colnames(result) <- "RowVal"
  result$NoViol <- (Diff >= 0)
  result$Inside <- (Diff > 0)
  result$RowVal <- round(Diff, 5)
  print(result)
  if (sum(!result$Inside) == 0){
    message("x0 is within constraints")
    message("No need to set manually below")
    good <- TRUE
  } else {
    if (sum(!result$NoViol) == 0){
      message("x0 does not violate constraints")
      message("But it is on Constraint Boundary (will not work for intial value)")
      message("Must set x0 manually below")
    } else {
      message("x0 violates constraints")
      message("Must set x0 manually below")
    } 
  } 
  return(list(result = result, good = good))
}

###############  Modeling Functions Environment ###############

env_modfun$PredMNL <- function(x, data_list, model_env) {
  U <- exp(data_list$ind %*% x)
  tasksum <- rowsum(U, data_list$idtask_r) # Summarize by task
  esum <- tasksum[data_list$idtask_r,]
  predprob = U / esum
}

env_modfun$LL_Neg <- function(x, data_list, model_env) {
  predprob <- model_env$func$pred(x = x, data_list = data_list, model_env = model_env)
  LLTot <- -1 * sum((log(predprob) * data_list$dep * data_list$wts))
  return(LLTot)
}

env_modfun$grad_MNL <- function(x, data_list, model_env) {
  # Gradient for LogLIke of MNL
  # To be accutae requires f.pred is MNL
  predprob <- model_env$func$pred(x = x, data_list = data_list, model_env = model_env)
  diff <- (predprob - data_list$dep) * data_list$wts
  grad <- t(data_list$ind) %*% diff
  return(grad)
}

env_modfun$LL_wPriorPDF <- function(x, data_list, model_env) {
  predprob <- model_env$func$pred(x = x, data_list = data_list, model_env = model_env)
  LLTot <- -1 * (sum((log(predprob) * data_list$dep * data_list$wts)) +
                   model_env$func$logpdf(x[model_env$prior$upper_model], model_env$prior$alpha, model_env$prior$cov_inv * model_env$prior$scale))
  return(LLTot)
}
# x is n x 1 column vector

env_modfun$grad_MNLwMVN <- function(x, data_list, model_env) {
  predprob <- model_env$func$pred(x = x, data_list = data_list)
  diff <- (predprob - data_list$dep) * data_list$wts # MNL Gradient
  diff2 <- (model_env$prior$cov_inv * model_env$prior$scale) %*% (x[model_env$prior$upper_model] - model_env$prior$alpha)
  diff2_all <- rep(0,length(x))
  diff2_all[model_env$prior$upper_model] <- diff2
  grad <- (t(data_list$ind) %*% (diff)) + diff2_all #MNL + MVN
  #grad[grad < -100] <- -100
  #grad[grad > 100] <- 100
  return(grad)
}

env_modfun$grad_num_make <- function(model_env){
  # model_env is unquoted name of list
  kstr <- bquote(
    function(x, data_list, model_env) {
      result <- grad(func = .(model_env$func$min), x = x, method = "simple", data_list = data_list, model_env = model_env)
      return(result)
    }  
  )
  eval(kstr)
}

env_modfun$PredProb_wPI <- function(x, data_list, model_env) {
  V <- as.matrix(data_list$ind[, c(-1, -2)]) %*% as.matrix(x[c(-1, -2)]) # No scale, none
  U <- exp(V)
  tasksum <- rowsum(U, data_list$idtask_r) # Summarize by task
  esum <- as.matrix(tasksum[data_list$idtask_r,])
  predprob1 <- U / esum # regular conjoint
  U <- exp((V * x[1]) + (x[2] * data_list$ind[, 2])) # (Beta x * scale) + (None * None Ind)
  tasksum <- rowsum(U, data_list$idtask_r) # Summarize by task
  esum <- as.matrix(tasksum[data_list$idtask_r,])
  predprob2 <- U / esum # PI tasks
  predprob <- (predprob1 * (data_list$ind[, 1] == 0)) +
    (predprob2 * (data_list$ind[, 1] == 1))
  return(predprob)
}

# Helper functions
env_modfun$logpdf_mvnorm <- function(beta, alpha, cov_inv) {
  result <- -.5 * ((t(beta - alpha) %*% (cov_inv)) %*% (beta - alpha))
}

env_modfun$logpdf_mvt <- function(beta, alpha, cov_inv, df = 100) {
  #multivariate t dist
  kexp <- (0 + df)/-2
  result <- kexp * log(1 + ((t(beta - alpha) %*% cov_inv) %*% (beta - alpha))/df)
  return(result)
}

###############  Stan Functions Environment ###############
env_stan$prep_file_stan <- function(idtaskdep, indcode_list, train = TRUE,
                                    data_cov, specs_cov_coding,
                                    check_collinearity = FALSE, other_data = NULL,
                                    force_sumtask_dep1 = TRUE,
                                    holdouts = 0, seed = NULL) {
  result <- list(tag = 0, N = 0, P = 0, I = 0, T = 0, dep = 0, ind = 0, sizes = 0,
                 code_master = indcode_list$code_master, n_atts = nrow(indcode_list$code_blocks), code_blocks = indcode_list$code_blocks)
  sort_order <- order(idtaskdep[, 1], idtaskdep[, 2], -idtaskdep[,3])
  sort_order[!train] <- 0 # Non-training gets order = 0, which removes
  result$row_in <- (1:nrow(idtaskdep))[sort_order] # Initial order of data with non-training removed
  result$ind <- as.matrix(indcode_list$indcode[sort_order,])
  result$ind_coded <- as.matrix(indcode_list$ind_coded[sort_order,])
  result$ind_levels <- as.matrix(indcode_list$ind_levels[sort_order,])
  colnames(result$ind_levels) <- colnames(indcode_list$ind_levels)
  result$sizes <- c(ncol(result$ind_coded), ncol(result$ind_levels), nrow(indcode_list$code_master))
  
  dep <- as.vector(as.matrix(idtaskdep[sort_order, 3]))
  idtask <- data.frame(idtaskdep[sort_order, 1:2])
  idtask_u <- as.matrix(unique(idtask))
  idtask_r <- (match(data.frame(t(idtask)), data.frame(t(idtask_u)))) # unique tasks
  resp_id <- as.vector(unique(idtask_u[, 1]))
  match_id <- match(idtask[, 1], as.matrix(resp_id))
  # Next 3 lines recodes dep to sum to 1
  if (force_sumtask_dep1){
    depsum <- rowsum(dep, idtask_r) # Sum dep each task
    depsum_match <- (depsum[idtask_r,]) # Map Sum to rows
    dep <- dep / depsum_match # sum of dep will add to 1
  }
  # Recode NAs to 0
  result$ind[is.na(result$ind)] <- 0
  dep[is.na(dep)] <- 0
  result$dep <- dep
  result$N <- nrow(result$ind); result$P <- ncol(result$ind);
  result$I <- length(resp_id); result$T <- max(idtask_r);
  result$task_scale_group <- rep(1, result$T)
  
  # Add Stan ragged array stuff
  result$end <- c(which(diff(idtask_r)!=0), length(idtask_r))
  result$start <- c(1, result$end[-length(result$end)]+1)
  result$task_individual <- match_id[result$start]
  id_task_end <- c(which(diff(result$task_individual, lag = 1) != 0),
                   length(result$task_individual))
  id_task_beg <- c(1, id_task_end[-length(id_task_end)]+1)
  id_row_beg <- result$start[id_task_beg]
  id_row_end <- result$end[id_task_end]
  result$id_ranges <- cbind(task_beg = id_task_beg, task_end = id_task_end,
                            row_beg = id_row_beg, row_end = id_row_end, nrows = id_row_end - id_row_beg + 1)
  result$npos_m1 <- aggregate(dep>0, list(idtask_r), sum)[,2]- 1
  result$npos_m1[result$npos_m1 < 0] <- 0 # In case all deps for task are 0
  if (!is.null(other_data)) {
    result$other_data <- as.matrix(other_data)[sort_order,]
  } else result$other_data <- 0
  
  # Friendly output
  cat("Prepared data_stan with coded data and constraints\n")
  cat(paste0("    ",length(resp_id)), " Respondents\n")
  cat(paste0("    ",sprintf("%.1f",max(idtask_r)/length(resp_id)),
             " Tasks per Respondent\n"))
  khead <- data.frame("   # Concepts in Tasks:  " = "   ", check.names=FALSE)
  t_table <- as.data.frame(table(result$end - result$start + 1, dnn = "N_Concepts"),
                           responseName = "N_Tasks")
  print(cbind(khead, t_table), row.names = F)
  cat(paste0("There are ", ncol(indcode_list$code_master)), " coded parameters\n")  
  cat(paste0("The final utilities will have ", nrow(indcode_list$code_master), " parameters:\n"))
  print(rownames(indcode_list$code_master))
  
  # Covariates, weights
  result$P_cov <- 0
  result$i_cov <- matrix(0, length(resp_id), 0)
  result$wts <- rep(1, length(result$task_individual))
  if (!is.null(data_cov)){ # Code respondent covariates and weights
    if (nrow(specs_cov_coding) >0){
      covariates_out <- code_covariates(data_cov, specs_cov_coding, resp_id)
      result$i_cov <- covariates_out$output 
      result$P_cov <- ncol(result$i_cov) # Num of coded parameters
      result$covariates_code <- covariates_out$cov_code_master
    }
    if (toupper(colnames(data_cov)[2]) %in% c("WTS","WT","WEIGHTS","WEIGHT")){
      result$wts <- make_wts(data_cov, resp_id)[result$task_individual]
    } else cat("Optional weights variable not found in covariates data.  All data weighted at 1\n")
  }
  
  # Select # holdouts corresponding to option > 0
  result$task_idtask <- idtask[result$start,]
  result$holdout <- as.integer(result$task_idtask[,2] < 0) # negative task numbers are holdouts
  if (holdouts >0){
    if (!is.null(seed)) set.seed(seed)
    for (i in 1:nrow(result$id_ranges)){
      task_set <- c(result$id_ranges[i,1]:result$id_ranges[i,2]) # tasks for resp
      holds_left <- holdouts - sum(result$holdout[task_set] == 1) # num holdouts still to pick
      if (holds_left > 0){
        picks <- sample(task_set, holds_left, prob = 1 - result$holdout[task_set])
        result$holdout[picks] <- 1
      }
    }
  } 
  
  # Check for collinearity
  if (check_collinearity){
    if (nrow(result$ind) > 1100000){
      result$collinear <- check_collinear(result$ind[c(1:500000, (nrow(result$ind) - 499999):nrow(result$ind)),])
      cat("Large data file, so checking collinearity for subset of 1 million rows\n")
    } else result$collinear <- check_collinear(result$ind)
  }  
  
  # Get initial x0 that satisfies constraints
  con_matrix <- diag(indcode_list$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], indcode_list$con_matrix)
  fcon <- function(x, con_matrix){
    result <- sum(log(1+exp(-100 * (con_matrix %*% x))))
    return(result)
  }
  x_initial <- optim(par = rep(0, ncol(indcode_list$code_master)), fn = fcon,
                     lower = -.5, upper = .5, method = "L-BFGS-B", control = list(maxit = 20),
                     con_matrix = con_matrix)
  if (nrow(con_matrix)>0){
    check_con <- min(con_matrix %*% x_initial$par)
    if (check_con <= 0) message("Please check constraints. Could not find initial value that satisfied your constraints. HB estimation in Stan can run, but constrained optimization models in R like Empirical Bayes cannot.\n")    
  }
  
  
  return(modifyList(result, list(
    con_sign = indcode_list$con_sign,
    paircon_rows = nrow(indcode_list$con_matrix),
    paircon_matrix = indcode_list$con_matrix,
    df = 2,
    prior_alpha = rep(0, ncol(result$ind)),
    a_sig = 10,
    prior_cov = indcode_list$indprior,
    cov_block = matrix(1, ncol(result$ind), ncol(result$ind)),
    prior_cov_scale = 1,
    x0 = x_initial$par,
    threads_rec = min(max(1,(detectCores() - 2)/2), round(.5 + result$T/500), round(.5 + result$I/16)),
    idtask = idtask, idtask_r = idtask_r,
    resp_id = resp_id, match_id = match_id))) 
}  

  

env_stan$check_collinear <- function(x, add_int = TRUE, vnames = NULL){
  #add_int adds an intercept which we usually want
  if(!is.matrix(x)) x <- as.matrix(x)
  if (is.null(vnames)){
    if (!is.null(colnames(x))){
      vnames <- colnames(x)
    } else vnames <- 1:ncol(x)
  }
  if (add_int){
    x <- cbind(x,1)
    vnames <- c(vnames, "intercept")
  } 
  mu <- colMeans(x)
  qr_obj <- qr(x) #QR decomp
  R <- qr.R(qr_obj) 
  cov_qr <- crossprod(R)/(nrow(x)-1) - mu %*%t(mu) # covariance of data from QR
  if (add_int){
    P <- ncol(cov_qr)
    cov_qr <- cov_qr[-P,-P] # remove last row and column because they are intercept
  } 
  result <- list()
  rank <- qr_obj$rank
  if (rank < ncol(x)){
    # This code for bad_combos from R package findLinearCombos
    pivot <- qr_obj$pivot              
    p1 <- 1:rank
    X <- R[p1, p1]                    # extract the independent columns
    Y <- R[p1, -p1, drop = FALSE]     # extract the dependent columns
    b <- qr(X)                        # factor the independent columns
    b <- qr.coef(b, Y)                # get regression coefficients of dependent columns
    b[abs(b) < 1e-6] <- 0             # zap small values
    # generate a list with one element for each dependent column
    bad_combos <- lapply(1:ncol(Y),
                         function(i) vnames[sort(c(pivot[rank + i], pivot[which(b[,i] != 0)]))]
    )
    
    message(paste0("\n############################################",
                   "\nWARNING!!!! YOUR DESIGN IS DEFICIENT",
                   "\nThe following variables are perfectly collinear:",
                   paste("\n", bad_combos, collapse = ""),
                   "\n############################################\n"))
    # result$bad_combos <- bad_combos
    result$cor_eigen <- 0
  } else {
    result$cor_eigen <- eigen(cov2cor(cov_qr))$value
    if (min(result$cor_eigen) < 1e-10){
      cat(paste0("\n############################################",
                 "\nWARNING!!!! YOUR DESIGN MAY BE DEFICIENT.",
                 "\nSmallest eigenvalue of correlation matrix is:", min(result$cor_eigen),
                 "\n############################################\n"))
    } else{
      cat(paste0("\nYour design is not deficient. Smallest eigenvalue is: ", min(result$cor_eigen))) 
    }
  }  
  return(result)
}


env_stan$create_tempdir <- function(dir, out_folder, out_prefix, stan_outname, save_specs = FALSE, code_master = NULL){
  dir.create(my_temp <- file.path(dir$work, out_folder))
  # dir.create(stan_out <- file.path(my_temp, "stan_out"))
  if (save_specs) saveRDS(object = list(specs_att_coding = specs_att_coding,
                                        specs_pair_constraints = specs_pair_constraints,
                                        specs_cov_coding  = specs_cov_coding),
                          file.path(my_temp, paste0(out_prefix,"_specs.rds")))
  if (!is.null(code_master)) write.csv(code_master, file.path(my_temp, paste0(out_prefix,"_code_master.csv")))
  cat("\nWhile Stan runs, you may check convergence with Stan csv output.\n")
  cat("To create subset 'stan_part.csv' of first 300 columns using Linux terminal:\n")
  cat(paste0("cd ", dir$stanmodel, "   # Stan output directory and then:\n",
             "tail -n +45 '",stan_outname,"-1.csv'  | cut -d, -f 1-300 > stan_part.csv"))
  return(my_temp)
}


env_stan$plot_draws_df <- function(draws, vnames = NULL, ylab = "Draw", pdf_path = NULL,
                                   chain_colors = rep(c("red","blue","green","black"),2),makepdf=NULL){
  if (is.null(makepdf)) makepdf <- TRUE
  # Assume draws is a list of data frames
  # Each data frame is draws for one chain of results with rows for iterations 
  # We plot each column
  if (is.null(vnames)) vnames <- colnames(draws[[1]]) # take names from first data frame
  fit_stats <- data.frame(
    variable = vnames,
    mean = NA,
    sd =  NA,
    rhat = NA,
    ESS = NA
  )
  if (!is.null(pdf_path) && makepdf){
    pdf(file = pdf_path,   # The directory you want to save the file in
        width = 7, # The width of the plot in inches
        height = 5) # The height of the plot in inches
  }
  for (i in 1:ncol(draws[[1]])){
    x <- sapply(1:length(draws), function(chain){
      draws[[chain]][,i]     
    })
    fit_stats$mean[i] <- round(mean(x), 5)
    fit_stats$sd[i] <- round(sd(x),5)
    fit_stats$rhat[i] <- round(rhat(x),2)
    fit_stats$ESS[i] <- round(ess_basic(x),1)
    if (makepdf==TRUE){
      plot(x[,1], type = "l", col = chain_colors[1], ylim = c(min(x), max(x)),
           xlab = "Sample Iteration", ylab = ylab,
           main = paste(vnames[i],
                        "| rhat = ", round(rhat(x),2),
                        "| ESS = ", round(ess_basic(x),1)
           ))
      if (ncol(x) > 1){
        for (chain in 2:ncol(x)){
          lines(x[,chain], type = "l", col = chain_colors[chain])
        }
      }
    }
  } # end for
  if (!is.null(pdf_path) && makepdf==TRUE) dev.off()
  return(fit_stats)
}

env_stan$check_draws_vector <- function(output_files, stan_vector, vector_labels = NULL, out_dir, out_prefix,makepdf=NULL){
  # stan_vector is vector/1D array or scalar
  # vector labels are names for each item in vector, NULL auto assigns
  draws_stan_var <- read_cmdstan_csv(output_files, variables = stan_vector, format = "draws_list")
  draws_stan_df <- lapply(draws_stan_var$post_warmup_draws, function(x) do.call(cbind, x))
  fit_stats_var <- plot_draws_df(draws_stan_df, vector_labels, toupper(stan_vector),
                                 pdf_path = file.path(out_dir,paste0(out_prefix,"_trace_",stan_vector,".pdf")),makepdf=makepdf)
  write.table(fit_stats_var, file = file.path(out_dir, paste0(out_prefix,"_fit_stats_", stan_vector,".csv")), sep = ",", na = ".", row.names = FALSE)
} 

env_stan$get_mean_beta <- function(draws_beta){
  nchains <- length(draws_beta$post_warmup_draws)
  result <- matrix(
    Reduce("+",lapply(draws_beta$post_warmup_draws, function(x) sapply(x, mean)))/nchains,
    draws_beta$metadata$stan_variable_sizes$beta_ind[2], # Respondents
    draws_beta$metadata$stan_variable_sizes$beta_ind[1], # Parameters
    byrow = TRUE) # First P entries are resp 1
  return(result)
}

env_stan$get_covariates <- function(output_files, data_stan, dir_out, out_prefix){
  draws_i_cov_load <- read_cmdstan_csv(output_files, variables = "i_cov_load", format = "draws_list")
  nchains <- length(draws_i_cov_load$post_warmup_draws)
  result <- matrix(
    Reduce("+",lapply(draws_i_cov_load$post_warmup_draws, function(x) sapply(x, mean)))/nchains,
    draws_i_cov_load$metadata$stan_variable_sizes$i_cov_load[2], # Covariates
    draws_i_cov_load$metadata$stan_variable_sizes$i_cov_load[1], # Parameters
    byrow = TRUE) # First P entries are Covariate level 1
  result <- (data_stan$code_master %*% t(result)) %*% t(data_stan$covariates_code) # Back code parameters
  result <-  cbind(parameters = rownames(result), data.frame(result))
  write.table(result, file = file.path(dir_out, paste0(out_prefix,"_covariates.csv")), sep = ",", na = ".", row.names = FALSE)
  return(result)
}

env_stan$checkconverge_export <- function(draws_beta, vnames, out_prefix, dir_run,
                                          export_draws_betas, export_draws_means,
                                          resp_id = data_stan$resp_id, code_master = data_stan$code_master,
                                          makepdf = TRUE){
  nchains <- length(draws_beta$post_warmup_draws)
  nresp <- draws_beta$metadata$stan_variable_sizes$beta_ind[2] # num respondents
  npar <- draws_beta$metadata$stan_variable_sizes$beta_ind[1]
  
  ### Save output files and check convergence ###
  draws_name <- paste0(out_prefix,"_draws_beta.rds")
  pdf_name <- paste0(out_prefix,"_trace_plots.pdf")
  fit_name <-  paste0(out_prefix,"_fit_stats.csv")
  cat("\nSaving Results in folder:", dir_run, "\n")
  message(paste0(
    "\nSaving post warm-up files for:\n",
    " convergence stats of mean:     ", fit_name, "\n",
    " PDF of detailed traceplots:    ", pdf_name,"\n"
  ))
  
  try(hist(as.vector(sapply(draws_beta$post_warmup_sampler_diagnostics, function(x) x$accept_stat__)), breaks = 30, main = "Acceptance Rate - Sampling", xlab = "", xlim = c(0,1)))
  if (export_draws_betas){
    stacked_draws <- do.call(rbind, lapply(1:nchains, function(i){
      result <- cbind(chain = i,
                      convert_draws(draws_beta$post_warmup_draws[[i]],
                                    nresp = nresp,
                                    P = npar,
                                    resp_id = resp_id, # Only if stacking
                                    code_master = code_master, # Only if back coding
                                    by_row = TRUE, 
                                    back_code = TRUE, stack = TRUE))
    }))    
    saveRDS(stacked_draws, file.path(dir_run, paste0(out_prefix,"_draws_beta.rds")))
  } 
  
  # Get mean of respondent betas for each iteration (like alpha in constrained space)
  ndraws <- draws_beta$metadata$iter_sampling
  draws_beta_mu <- list()
  for (chain_i in (1:nchains)){
    x <- do.call(rbind, draws_beta$post_warmup_draws[[chain_i]]) # bind list  PxResp rows, iter in cols 
    draws_beta_mu[[chain_i]] <- t(aggregate.data.frame(data.frame(x), by = list(rep(1:npar, nresp)), FUN = "mean")[,-1])
  }
  if (export_draws_means){
    draws_beta_mu_export <- do.call(rbind, lapply(1:length(draws_beta_mu), function(i){
      result <- cbind(chain = i, iter = 1:nrow(draws_beta_mu[[i]]), draws_beta_mu[[i]])
    }))
    colnames(draws_beta_mu_export)[3:ncol(draws_beta_mu_export)] <- vnames
    write.table(draws_beta_mu_export, file = file.path(dir_run, paste0(out_prefix,"_draws_beta_means.csv")), sep = ",", na = ".", row.names = FALSE)
  }
  
  # Create pdf of acceptance rate + traceplots
  if (makepdf){
    pdf(file = file.path(dir_run, pdf_name),   # The directory you want to save the file in
        width = 7, # The width of the plot in inches
        height = 5) # The height of the plot in inches
    hist(as.vector(sapply(draws_beta$post_warmup_sampler_diagnostics, function(x) x$accept_stat__)), breaks = 30, main = "Acceptance Rate - Sampling", xlab = "", xlim = c(0,1))
    for (chain_i in (1:nchains)){
      matplot(1:nrow(draws_beta_mu[[chain_i]]), draws_beta_mu[[chain_i]],
              type = "l" , lty = 1, lwd = 1, main = paste0("Chain ", chain_i), xlab = "Iteration", ylab = "Mean Beta")   
    }
    fit_stats <- plot_draws_df(draws= draws_beta_mu, vnames = vnames, ylab = "Mean Beta") # Plots each column
    dev.off()
    write.table(fit_stats, file = file.path(dir_run, paste0(out_prefix,"_fit_stats.csv")), sep = ",", na = ".", row.names = FALSE)
  }
}


env_stan$obs_vs_pred <- function(wts_obs_pred, cat_vars, predvar_out = "pred_per"){
  ind_levels <- as.data.frame(cat_vars) # in case cat_vars is a vector
  obs_vs_pred <- do.call(rbind, lapply(1:ncol(ind_levels), function(i){
    ksplit <- split(wts_obs_pred, ind_levels[,i], drop = TRUE)
    result <- do.call(rbind,lapply(ksplit, function(x) {
      base <- sum(x[,1])
      return(c(base, sum(x[,2]), sum(x[,1] * x[,2])/base, sum(x[,1] * x[,3])/base))
    }))
    result_df <- data.frame(cat_var= colnames(ind_levels)[i],
                            level = (rownames(result)),
                            shows = result[,1],
                            chosen = result[,2],
                            obs_per = result[,3],
                            pred_per = result[,4])
    return(result_df)
  }))
  
  return(obs_vs_pred)
}

env_stan$predict_fromutil <- function(data_stan, utilities, task_type = NULL, inv_logit_thresh = NULL){
  if (is.null(task_type)){ # Standard MNL
    pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                      function(t){
                                        U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                   utilities[data_stan$task_individual[t],])
                                        pred <- U/sum(U)
                                        return(pred)
                                      }
    )) # lapply, rbind
  } else {
    pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                      function(t){
                                        if (task_type[t] == 1){ # MNL
                                          U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                     utilities[data_stan$task_individual[t],])
                                          pred <- U/sum(U)
                                        }
                                        if (task_type[t] == 2){ # Inv Logit
                                          U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                     utilities[data_stan$task_individual[t],])
                                          if (is.null(inv_logit_thresh)){
                                            pred <- U/(1+U)
                                          } else pred <- U/(exp(inv_logit_thresh[t]) + U)
                                        }
                                        return(pred)
                                      }
    )) # lapply, rbind                                   
  }
  return(pred_all)
}

env_stan$process_utilities <- function(data_stan, utilities, out_prefix, dir_run,
                                       task_type = NULL, inv_logit_thresh = NULL){
  # Compute predictions
  row_weights <- data_stan$wts[data_stan$idtask_r] # April 2025 task weights * holdout -> row weights
  row_holdout <- data_stan$holdout[data_stan$idtask_r]
  n_holdouts <- sum(data_stan$holdout)
  pred_all <- predict_fromutil(data_stan, utilities, task_type, inv_logit_thresh)   
  utilities_r <- utilities %*% t(data_stan$code_master)
  util_name <- paste0(out_prefix,"_utilities_r.csv")
  pred_name <- paste0(out_prefix,"_preds_meanpt.csv")
  fit_name <- paste0(out_prefix,"_Hit_LL.csv")
  failcon_name <- paste0(out_prefix,"_utilities_failcon.csv")
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred.csv")
  
  LL_id <- rowsum(log(pred_all) * data_stan$dep * row_weights, data_stan$match_id)
  sum_wts <- rowsum(data_stan$dep * row_weights, data_stan$match_id)
  sum_wts[sum_wts == 0] <- 1
  header <- data.frame(id = data_stan$resp_id, rlh = exp(LL_id/sum_wts))
  message(paste0(
    "Saving: \n",
    " Respondent mean utilities: ", util_name, "\n",
    " Predictions for data     : ", pred_name, "\n",
    " Observed vs Predicted    : ", obs_vs_pred_name, "\n",
    " Hit Rate and LogLike     : ", fit_name
  ))
  
  write.table(cbind(header, utilities_r), file = file.path(dir_run, util_name), sep = ",", na = ".", row.names = FALSE)
  
  pred_all_export <- cbind(data_stan$idtask, holdout = row_holdout, wt = row_weights, dep = data_stan$dep, pred = pred_all)
  if (ncol(data_stan$ind_levels) >0){
    obs_vs_pred <- obs_vs_pred(pred_all_export[,4:6], data_stan$ind_levels)
    write.table(obs_vs_pred, file = file.path(dir_run, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  } else message("No Categorical Variables found for observed vs predicted")
  if ("row_in" %in% names(data_stan)){
    pred_all_export <- cbind(row_in = data_stan$row_in, pred_all_export)
    pred_all_export <- pred_all_export[order(data_stan$row_in),]
  } else pred_all_export <- cbind(row_in = NA, pred_all_export)
  write.table(pred_all_export, file = file.path(dir_run, pred_name), sep = ",", na = ".", row.names = FALSE)
  
  # Check if utilities meet constraints
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  bad_ids <- rowSums(((utilities %*% t(con_matrix)) < 0)) > 0
  if (sum(bad_ids) > 0){
    message(paste0(sum(bad_ids), " Respondents had reversals from constraints.\n",
                   "Reversals saved to: ", failcon_name))
    write.table(cbind(header, utilities_r)[bad_ids,], file = file.path(dir_run, failcon_name), sep = ",", na = ".", row.names = FALSE)
  } else message(" All respondent mean utilities obey constraints")
  fit_hit_LL <- hit_rate_LL(pred_all_export ,col_id = 2,col_task = 3, col_dep = 6,
                            col_pred = 7, col_split = 4, col_wts = 5)
  write.table(fit_hit_LL, file = file.path(dir_run, fit_name), sep = ",", na = ".", row.names = FALSE)      
}

env_stan$hit_rate_LL <- function(data_wpreds, col_id, col_task, col_dep, col_pred, col_split = NULL, col_wts = NULL){
  # Note task_wts are the weights of each 
  ksplit <- split(data_wpreds, data_wpreds[,c(col_task, col_id)], drop = TRUE) # No empty elements
  fit_tasks <- data.frame(do.call(rbind, lapply(ksplit, function(x){ # wts for id x task
    if (is.null(col_wts)){
      wt <- 1
    } else {
      wt <- mean(x[,col_wts])
    }
    if (is.null(col_split)){
      split <- 1
    } else {
      split1 <- min(x[,col_split])
      split2 <- max(x[,col_split])
      split <- NA
      if (split1 == split2) split <- split1
    }
    preds <- x[,col_pred, drop = TRUE] # Force vector
    dep <- x[,col_dep, drop = TRUE]
    row_max <- which(preds == max(preds, na.rm = TRUE)) # Allows for ties
    hit <- mean(dep[row_max]) 
    LL <- sum(log(preds) * dep)
    result <- c(split, wt, hit, LL)
    return(result)
  })))
  knames <- colnames(data_wpreds)
  split_name <- "Group"
  if (!is.null(col_split)) split_name <- knames[col_split]
  colnames(fit_tasks) <- c(split_name, "wt", "hit", "LL")  
  
  split_group <- split(fit_tasks, fit_tasks[,1], drop = TRUE)
  fit_agg <- data.frame(do.call(rbind, lapply(split_group, function(x){
    sum_wts <- sum(x[,2])
    hit_rate <- sum(x[,3] * x[,2])/sum_wts
    LL_mean <- sum(x[,4] * x[,2])/sum_wts
    rlh <- exp(LL_mean)
    result <- c(x[1,1],nrow(x),sum_wts,hit_rate, LL_mean, rlh)
    return(result)
  })))
  colnames(fit_agg) <- c(split_name,"num_tasks", "sum_wts","hit_rate","LL_mean", "RLH")
  
  # cat("Total of ", length(ksplit), " tasks\n")
  # cat("Used the following for\n")
  # cat(paste("  defining tasks  :", knames[col_id], knames[col_task],"\n"))
  # cat(paste("  observed vs pred:", knames[col_dep], knames[col_pred],"\n"))
  # print.data.frame(fit_agg, row.names = FALSE)
  return(fit_agg)
}

env_stan$process_dualtypes <- function(data_stan, utilities, util_thresh, scale, lb, range, out_prefix, dir_run, ind_thresh = NULL){
  # Compute predictions
  row_weights <- data_stan$wts[data_stan$idtask_r] # convert task weights to row weights
  task_type <- data_stan$task_type
  row_task_type <- task_type[data_stan$idtask_r] # Carry task_type to row
  choice_tasks <- row_task_type == 1
  if (is.null(ind_thresh)) ind_thresh <- as.matrix(rep(1,data_stan$N))
  pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                    function(t){
                                      if (task_type[t] == 1){ # MNL
                                        U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                   utilities[data_stan$task_individual[t],])
                                        pred <- U/sum(U)
                                      }
                                      if (task_type[t] == 2){ # Rating = lb + range * Inv Logit
                                        util <- scale[data_stan$task_individual[t]] * data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                          utilities[data_stan$task_individual[t],]
                                        thresh <- ind_thresh[data_stan$start[t]:data_stan$end[t],] %*%
                                          (util_thresh[data_stan$task_individual[t],])
                                        pred <- lb + range/(1+ exp(thresh - util)) # = lb + range * e^util/(e^util + e^thresh)
                                        
                                      }
                                      return(pred)
                                    }
  )) # lapply, rbind                                   
  
  utilities_r <- cbind(utilities %*% t(data_stan$code_master), beta_scale, utilities_thresh)
  util_name <- paste0(out_prefix,"_utilities_r.csv")
  pred_name <- paste0(out_prefix,"_preds_meanpt.csv")
  failcon_name <- paste0(out_prefix,"_utilities_failcon.csv")
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred.csv")
  
  LL_id <- rowsum(log(pred_all) * data_stan$dep * row_weights * choice_tasks, data_stan$match_id)
  sum_wts <- rowsum(data_stan$dep * row_weights * choice_tasks, data_stan$match_id)
  sum_wts[sum_wts == 0] <- 1
  SUM_AE_id <- rowsum(abs(pred_all - data_stan$dep) * row_weights * (!choice_tasks), data_stan$match_id)
  sum_wts_2 <- rowsum(row_weights * (!choice_tasks), data_stan$match_id)
  sum_wts_2[sum_wts_2 == 0] <- 1
  
  header <- data.frame(id = data_stan$resp_id, rlh = exp(LL_id/sum_wts), MAE_rate = SUM_AE_id/sum_wts_2)
  message(paste0(
    "Saving: \n",
    " Respondent mean utilities: ", util_name, "\n",
    " Predictions for data     : ", pred_name, "\n",
    " Observed vs Predicted    : ", obs_vs_pred_name 
  ))
  
  write.table(cbind(header, utilities_r), file = file.path(dir_run, util_name), sep = ",", na = ".", row.names = FALSE)
  
  pred_all_export <- cbind(data_stan$idtask, wts = row_weights, dep = data_stan$dep, pred = pred_all)
  if (ncol(data_stan$ind_levels) >0){
    obs_vs_pred <- obs_vs_pred(pred_all_export[choice_tasks,3:5], data_stan$ind_levels[choice_tasks,])
    write.table(obs_vs_pred, file = file.path(dir_run, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  } else message("No Categorical Variables found for observed vs predicted")
  if ("row_in" %in% names(data_stan)){
    pred_all_export <- cbind(row_in = data_stan$row_in, pred_all_export)
    pred_all_export <- pred_all_export[order(data_stan$row_in),]
  } else pred_all_export <- cbind(row_in = NA, pred_all_export)
  write.table(pred_all_export, file = file.path(dir_run, pred_name), sep = ",", na = ".", row.names = FALSE)
  
  # Check if utilities meet constraints
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  bad_ids <- rowSums(((utilities %*% t(con_matrix)) < 0)) > 0
  if (sum(bad_ids) > 0){
    message(paste0(sum(bad_ids), " Respondents had reversals from constraints.\n",
                   "Reversals saved to: ", failcon_name))
    write.table(cbind(header, utilities_r)[bad_ids,], file = file.path(dir_run, failcon_name), sep = ",", na = ".", row.names = FALSE)
  } else message(" All respondent mean utilities obey constraints")
}

env_stan$process_utilities_utility_add <- function(data_stan, utilities, out_prefix, dir_run,
                                                   task_type = NULL, inv_logit_thresh = NULL){
  # Compute predictions
  row_weights <- data_stan$wts[data_stan$idtask_r] # convert task weights to row weights
  if (is.null(task_type)){ # Standard MNL
    pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                      function(t){
                                        U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                   utilities[data_stan$task_individual[t],] +
                                                   data_stan$utility_add[data_stan$start[t]:data_stan$end[t]])
                                        pred <- U/sum(U)
                                        return(pred)
                                      }
    )) # lapply, rbind
  } else {
    pred_all <- do.call(rbind, lapply(1:data_stan$T,
                                      function(t){
                                        if (task_type[t] == 1){ # MNL
                                          U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                     utilities[data_stan$task_individual[t],] +
                                                     data_stan$utility_add[data_stan$start[t]:data_stan$end[t]])
                                          pred <- U/sum(U)
                                        }
                                        if (task_type[t] == 2){ # Inv Logit
                                          U <- exp(data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                                                     utilities[data_stan$task_individual[t],] +
                                                     data_stan$utility_add[data_stan$start[t]:data_stan$end[t]])
                                          if (is.null(inv_logit_thresh)){
                                            pred <- U/(1+U)
                                          } else pred <- U/(exp(inv_logit_thresh[t]) + U)
                                        }
                                        return(pred)
                                      }
    )) # lapply, rbind                                   
  }
  utilities_r <- utilities %*% t(data_stan$code_master)
  util_name <- paste0(out_prefix,"_utilities_r.csv")
  pred_name <- paste0(out_prefix,"_preds_meanpt.csv")
  failcon_name <- paste0(out_prefix,"_utilities_failcon.csv")
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred.csv")
  
  LL_id <- rowsum(log(pred_all) * data_stan$dep * row_weights, data_stan$match_id)
  sum_wts <- rowsum(data_stan$dep * row_weights, data_stan$match_id)
  sum_wts[sum_wts == 0] <- 1
  header <- data.frame(id = data_stan$resp_id, rlh = exp(LL_id/sum_wts))
  message(paste0(
    "Saving: \n",
    " Respondent mean utilities: ", util_name, "\n",
    " Predictions for data     : ", pred_name, "\n",
    " Observed vs Predicted    : ", obs_vs_pred_name 
  ))
  
  write.table(cbind(header, utilities_r), file = file.path(dir_run, util_name), sep = ",", na = ".", row.names = FALSE)
  
  pred_all_export <- cbind(data_stan$idtask, wts = row_weights, dep = data_stan$dep, pred = pred_all)
  if (ncol(data_stan$ind_levels) >0){
    obs_vs_pred <- obs_vs_pred(pred_all_export[,3:5], data_stan$ind_levels)
    write.table(obs_vs_pred, file = file.path(dir_run, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  } else message("No Categorical Variables found for observed vs predicted")
  if ("row_in" %in% names(data_stan)){
    pred_all_export <- cbind(row_in = data_stan$row_in, pred_all_export)
    pred_all_export <- pred_all_export[order(data_stan$row_in),]
  } else pred_all_export <- cbind(row_in = NA, pred_all_export)
  write.table(pred_all_export, file = file.path(dir_run, pred_name), sep = ",", na = ".", row.names = FALSE)
  
  # Check if utilities meet constraints
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  bad_ids <- rowSums(((utilities %*% t(con_matrix)) < 0)) > 0
  if (sum(bad_ids) > 0){
    message(paste0(sum(bad_ids), " Respondents had reversals from constraints.\n",
                   "Reversals saved to: ", failcon_name))
    write.table(cbind(header, utilities_r)[bad_ids,], file = file.path(dir_run, failcon_name), sep = ",", na = ".", row.names = FALSE)
  } else message(" All respondent mean utilities obey constraints")
}


env_stan$convert_draws <- function(draws_chain, # draws for one chain
                            nresp = data_stan$I,
                            P = data_stan$P,
                            resp_id = data_stan$resp_id, # Only if stacking
                            code_master = data_stan$code_master, # Only if back coding
                            by_row = TRUE,
                            back_code = FALSE, 
                            stack = FALSE) { # Convert list of matrices to one stacked matrix
    # converts draws for one chain to list of matrices
    niter <- length(draws_chain[[1]])
    result <- lapply(1:niter, function(i){
      util <- matrix(sapply(draws_chain, "[[", i), nresp, P, byrow = by_row)
      if (back_code) util <- util %*% t(code_master)
      return (util)           
    })
    # result is list of matrices.  Option to stack that list
    if(stack){
      result <- do.call(rbind, lapply(1:length(result), function(i){
        addid <- cbind(iter = i, id = resp_id, result[[i]])
      }))
    }
    return(result)
}

env_stan$zip_allout <- function(files_dir, out_dir, out_name){
  # Zip all files in in_dir to file.path(out_dir, out_names)
  # Default: Remove directory where initial files were stored
  options(warn = -1)
  zip::zip(zipfile = file.path(out_dir, out_name), files = file.path(files_dir, list.files(files_dir)),
           include_directories = FALSE, recurse = FALSE, mode = "cherry-pick") 
  options(warn = 0)
  message(paste0("Files zipped into folder", out_dir, " : ", out_name, "\n"))
}

env_stan$est_agg_model <- function(data_list, maxit = 100, reltol = 1e-5, con_use = 0){
  data_list$wts <- data_list$wts[data_list$idtask_r] # Convert task weights to row weights
  model_agg <- list(
    func = list(pred = PredMNL, min = LL_Neg, gr = grad_MNL),
    x0 =  rep(0, data_list$P) # start at 0
  )
  cat("Estimating aggregate MNL model for checking\n")
  if (con_use == 0){
    agg_beta <- optim(par = model_agg$x0, fn = model_agg$func$min, gr = model_agg$func$gr, method ="BFGS",
                      data_list = data_list, model_env = model_agg, control = list(maxit = maxit, reltol = reltol, trace = 1, REPORT = 1))
  } else {
    con_matrix <- diag(data_list$con_sign)
    con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_list$paircon_matrix)
    check_con <- min(con_matrix %*% data_list$x0)
    message(check_con)
    if (check_con <= 0){
      message("Constrained model not estimated. Starting value in data$x0 does not satisfy constraints")
      message("Estimating unconstrained aggregate model")
      agg_beta <- optim(par = model_agg$x0, fn = model_agg$func$min, gr = model_agg$func$gr, method ="BFGS",
                        data_list = data_list, model_env = model_agg, control = list(maxit = maxit, reltol = reltol, trace = 1, REPORT = 1))
    } else {
      agg_beta <- constrOptim(theta = data_list$x0, f = model_agg$func$min, grad = model_agg$func$gr,
                              ui = con_matrix, ci = rep(0,nrow(con_matrix)), mu = 1e-02, 
                              method = "BFGS", outer.iterations = 100, outer.eps = 1e-05, 
                              data_list = data_list,
                              model_env = model_agg,
                              control = list(maxit = maxit, reltol = reltol, trace = 1, REPORT = 1))   
    }
  }
  return(agg_beta)
}


env_stan$eb_betas_est <- function(data_stan, draws_beta, x0, r_cores, out_prefix, dir_work, cov_scale, linux = TRUE, nids_core = 5){
  cat("\n")
  cat("Computing Empirical Bayes point estimates using respondent draws and constraints:")
  
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  model_eb <- list(
    func = list(
      pred = env_modfun$PredMNL,
      min = env_modfun$LL_wPriorPDF,
      gr = env_modfun$grad_MNLwMVN,
      logpdf = env_modfun$logpdf_mvnorm
    ),
    prior = list(
      alpha = NULL,
      cov_inv = NULL, # Cov Inverse 
      upper_model = rep(TRUE, length(x0)),
      scale = 1 # scale factor that we will solve for
    ),
    con = as.matrix(con_matrix), # must be matrix, 
    x0 = x0 # starting point inside constraints - use overall mean
  )
  
  id_eb <- function(idseq){
    end <- idseq * data_stan$P
    start <- end - data_stan$P + 1
    resp_draws <- do.call(rbind, lapply(draws_beta$post_warmup_draws, function(x) do.call(cbind, x[start:end])))
    model_id <- model_eb
    resp_mu <- colMeans(resp_draws)
    model_id$prior$alpha <- resp_mu
    # model_id$prior$cov_inv <- solve(cov(resp_draws) * cov_scale)
    # Replaced solve with SVD March 21, 2023
    ksvd <- svd(cov(resp_draws) * cov_scale)
    ksvd$d[ksvd$d < 1e-12] <- 1e-12
    model_id$prior$cov_inv <- ksvd$v %*% diag(1/ksvd$d) %*% t(ksvd$u)
    
    id_list <- list()
    id_filter <- (data_stan$match_id == idseq)
    task_id <- data_stan$idtask_r[id_filter]
    task_id_u <- unique(task_id)
    id_list$idtask_r <- (match(data.frame(t(task_id)), data.frame(t(task_id_u)))) # unique tasks
    id_list$ind = data_stan$ind[id_filter,]
    id_list$dep <- as.matrix(data_stan$dep[id_filter]) # Need matrix
    id_list$wts <- data_stan$wts[task_id]
    
    eb_solve <- constrOptim(theta = model_id$x0, f = model_id$func$min, grad = model_id$func$gr,
                            ui = model_id$con, ci = rep(0,nrow(model_id$con)), mu = 1e-02, control = list(),
                            method = "BFGS", outer.iterations = 100, outer.eps = 1e-05, 
                            data_list = id_list,
                            model_env = model_id)
    
    predprob <- model_id$func$pred(x = eb_solve$par, data_list = id_list)
    predprob_mu <- model_id$func$pred(x = resp_mu, data_list = id_list)
    sum_wt <- sum(id_list$dep *id_list$wts)
    rlh <- exp(sum(log(predprob) * id_list$dep * id_list$wts)/sum_wt)
    betas <- c(idseq, rlh, round(eb_solve$par, 6))
    names(betas) <- c("idseq", "rlh_eb", colnames(data_stan$ind))
    preds <- data.frame(data_stan$idtask[id_filter,], id_list$wts, id_list$dep, predprob)
    result <- list(betas, preds)
    return(result)
  }
  
  if (linux){
    idseq_all <- seq(1,data_stan$I, by = r_cores * nids_core)
    if (idseq_all[length(idseq_all)] < data_stan$I){
      idseq_all <- c(idseq_all, data_stan$I)
    } 
    nruns <- length(idseq_all) - 1
    betas_eb <- vector(mode = "list", length = nruns)
    preds <- vector(mode = "list", length = nruns)
    for (i in 1:nruns){
      if (i == nruns){
        sub_result <- mclapply(idseq_all[i]:(idseq_all[i+1]), id_eb, mc.cores = r_cores)
      } else {
        sub_result <- mclapply(idseq_all[i]:(idseq_all[i+1]-1), id_eb, mc.cores = r_cores)      
      } 
      betas_eb[[i]] <- do.call(rbind, lapply(sub_result, function(x) x[[1]]))
      preds[[i]] <- do.call(rbind, lapply(sub_result, function(x) x[[2]]))
      cat(paste0("\n", sprintf("%.1f", 100 * i/nruns), "%"))
    }
    betas_eb <- do.call(rbind, betas_eb)
    preds <- do.call(rbind, preds)
  } else {
    setup_cores(r_cores)
    result <- list()
    result <- foreach(idseq = 1:length(data_stan$resp_id)) %dopar% {
      result[[idseq]] <- id_eb(idseq)
    }
    setup_cores(0) # Close cores
    betas_eb <- do.call(rbind, lapply(result, function(x) x[[1]]))
    preds <- do.call(rbind, lapply(result, function(x) x[[2]]))
  }
  if (nrow(betas_eb) < data_stan$I) cat("Some respondents not estimated for EB")
  header <- data.frame(betas_eb[,1:2])
  header[,1] <- data_stan$resp_id[header[,1]] # convert idseq to id
  utilities_r_eb <- as.matrix(betas_eb[,-1:-2])  %*% t(data_stan$code_master) # id, rlh_eb
  util_eb_name <- paste0(out_prefix,"_utilities_r_eb.csv")
  write.table(cbind(header, utilities_r_eb), file = file.path(dir_work, util_eb_name), sep = ",", na = ".", row.names = FALSE)
  message(paste0("\nEB point estimates       : ",util_eb_name))
  colnames(preds) <- c("id","task","wts", "dep","pred_eb")
  preds_name <- paste0(out_prefix,"_preds_EB.csv")

  obs_vs_pred <- obs_vs_pred(preds[,3:5], data_stan$ind_levels)
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred_EB.csv")
  write.table(obs_vs_pred, file = file.path(dir_work, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  message(paste0("EB Observed vs Predicted: ", obs_vs_pred_name))  
  
  if ("row_in" %in% names(data_stan)){
    preds <- cbind(row_in = data_stan$row_in, preds)
    preds <- preds[order(data_stan$row_in),]
  } else preds <- cbind(row_in = NA, preds)
  write.table(preds, file = file.path(dir_work, preds_name), sep = ",", na = ".", row.names = FALSE)  
  message(paste0("EB predictions for data : ", preds_name))
}

# Version 6.0 Improvements
env_stan$LC_est <- function(data_stan, data_model, out_folder, stanout, out_prefix, dir_run){
  if (min(data_model$LC_segs) > 0){
    LC_outname <- paste0(out_folder, "_LCOut")
    LC_betas <- c()
    LC_probs0 <- data.frame(id = data_stan$resp_id)
    LC_probsno0 <- data.frame(id = data_stan$resp_id)
    LC_LL <- c()
    LC_pred <- cbind(row_in = data_stan$row_in, dep = data_stan$dep)
    for (nseg in data_model$LC_segs){
      data_model$Nseg <- nseg #fancy formula as standard for the number of segments?
      LC_optim <- LC_model$optimize(modifyList(data_stan, data_model),
                                    seed = NULL,
                                    refresh = 10,
                                    init = 0.1,
                                    save_latent_dynamics = TRUE,
                                    output_dir = stanout,
                                    output_basename = LC_outname,
                                    sig_figs = NULL,
                                    threads =  (detectCores() ),
                                    opencl_ids = NULL,
                                    algorithm = "lbfgs",
                                    init_alpha = NULL,
                                    iter = 20000)
      LC_lp <- unlist(read_cmdstan_csv(LC_optim$output_files(), variables = "lp__", format = "draws_list")$point_estimates)
      resp_seg_prob_posterior <- read_cmdstan_csv(LC_optim$output_files(), variables = "resp_seg_prob_posterior", format = "draws_list")
      resp_seg_prob_posterior <- matrix((unlist(resp_seg_prob_posterior$point_estimates)),nrow=data_stan$I)
      # resp_seg_prob_prior <- read_cmdstan_csv(LC_optim$output_files(), variables = "resp_seg_prob_prior", format = "draws_list")
      # resp_seg_prob_prior <- matrix((unlist(resp_seg_prob_prior$point_estimates)),nrow=data_stan$I)
      colnames(resp_seg_prob_posterior) <- paste0("prob_seg",nseg,"_",0:nseg)
      resp_seg_x0 <- resp_seg_prob_posterior[,-1]
      resp_seg_x0 <- sweep(resp_seg_x0, 1, rowSums(resp_seg_x0),"/")
      
      beta_seg <- read_cmdstan_csv(LC_optim$output_files(), variables = "beta_seg", format = "draws_list")
      beta_seg <- matrix((unlist(beta_seg$point_estimates)),nrow=data_stan$P)
      
      # Get respondent predictions and overall fit
      pred <- data_stan$dep
      for (t in (1:data_stan$T)){
        data_rows <- data_stan$start[t]:data_stan$end[t]
        util_byseg <- exp(data_stan$ind[data_rows,] %*% beta_seg)
        pred_t <- util_byseg[,-1] %*% as.matrix(resp_seg_x0[data_stan$task_individual[t],]/colSums(util_byseg[,-1]))
        pred[data_rows] <- pred_t
      }
      
      LC_utilities <- data_stan$code_master %*% beta_seg
      colnames(LC_utilities)<-paste0("seg",nseg,"_",0:nseg)
      LC_utilities_no0 <- LC_utilities[,-1] # remove LC segment with all 0 utilities
      # Combine outputs
      LC_LL <- c(LC_LL, LC_lp)
      LC_betas <- cbind(LC_betas,LC_utilities_no0)
      LC_probs0 <- cbind(LC_probs0,resp_seg_prob_posterior[,1])
      LC_probsno0 <- cbind(LC_probsno0,resp_seg_x0)
      LC_pred <- cbind(LC_pred, pred)
    }
    names(LC_LL) <- paste0(c("LL_seg_"),data_model$LC_segs)
    LC_pred[data_stan$row_in,] <- LC_pred
    message(paste0(
      "LC results in LC_fit \n",
      "Total Loglikelihood of each LC solution: \n"
    ))
    print(as.data.frame(LC_LL))
    write.table(data.frame(LC_LL), file = file.path(dir_run, paste0(out_prefix,"_LC_LL.csv")), sep = ",", na = ".", row.names = TRUE)
    write.table(LC_probsno0, file = file.path(dir_run, paste0(out_prefix,"_LC_probs.csv")), sep = ",", na = ".", row.names = FALSE)
    LC_betas <- cbind(data.frame(att = row.names(data_stan$code_master)), SawPrize_LC_betas)
    write.table(LC_betas, file = file.path(dir_run, paste0(out_prefix,"_LC_betas.csv")), sep = ",", na = ".", row.names = FALSE)
    result <- list(LC_LL = LC_LL, LC_betas = LC_betas, LC_probs = LC_probsno0, LC_probs0 = LC_probs0, LC_pred = LC_pred, LC_prob_good = 1 - apply(LC_probs0[,-1],1,min))
    saveRDS(result, file.path(dir_run, paste0(out_prefix, "_LC_fit.rds")))
  } else result <- NULL
  return(result)
}


env_stan$predict_task_base <- function(t, ind, utilities, scale_factor, task_type = NULL, inv_logit_thresh = NULL, draws = FALSE, pred_list_more = NULL){
  # utilities is        Px1 vector or P x # draws 
  # inv_logit_thresh is Tx1 vector or t x # draws
  # scale_factor is scale_factor for that task (single pt or vector of draws)
  if (!draws){ # U is nalt x 1 vector
    U <- exp(ind %*% (utilities * scale_factor)) # utilities is vector, scale_factor is scalar, U is vector
    if (is.null(task_type)){
      pred <- U/sum(U)
    } else {
      if (task_type[t] == 1){ # MNL
        pred <- U/sum(U)
      }
      if (task_type[t] == 2){ # Inv Logit
        if (is.null(inv_logit_thresh)){
          pred <- U/(1 + U)
        } else {
          pred <- U/(exp(inv_logit_thresh[t]) + U)
        }
      }
    }
  } # End points
  
  if (draws) {
    # DRAWS U is nalt x #draws
    U <- exp(ind %*% sweep(utilities,2,scale_factor,"*")) # (nalt x p) * (p x ndraws utilities)
    if (is.null(task_type)){
      pred <- sweep(U, 2, colSums(U), "/")
    } else {
      if (task_type[t] == 1){ # MNL
        pred <- sweep(U, 2, colSums(U), "/")
      }
      if (task_type[t] == 2){ # Inv Logit
        if (is.null(inv_logit_thresh)){
          pred <- U/(1 + U)
        } else {
          UminusT <- sweep(log(U), 2, inv_logit_thresh[t,], "-")
          pred <- 1/(1 + exp(-UminusT))
        }
      }
    }  
  } # End draws 
  return(pred) # For draws this returns n_alt x n_draws
}

env_stan$process_HB_old <- function(data_stan, draws_beta, out_prefix, dir_run, prob_bad = NULL,
                       scale_factor, task_type = NULL, inv_logit_thresh = NULL,
                       predict_task = predict_task_base, pred_list_more = NULL){
  P <- data_stan$P
  row_weights <- data_stan$wts[data_stan$idtask_r] 
  row_holdout <- data_stan$holdout[data_stan$idtask_r]
  n_holdouts <- sum(data_stan$holdout)
  
  utilities <- matrix(NA, data_stan$I, data_stan$P) # Hold utilities
  pred_all <- matrix(NA, data_stan$N, 2) # Hold predictions points, draws
  fit_id <- matrix(NA, data_stan$I, 4) # RlH (pts, draws), Prob_Good (HB, LC)
  if(is.null(prob_bad)){
    prob_bad_use <- rep(.0001, draws_beta$metadata$iter_sampling * length(draws_beta$post_warmup_draws))
  } else prob_bad_use <- prob_bad
  scale_factor <- as.matrix(scale_factor)
  scale_factor_mean <- exp(colMeans(log(scale_factor)))
  for (id in (1:data_stan$I)){
    draw_start <- 1 + (id-1) * P
    draw_end <- draw_start + P - 1
    resp_draws <- do.call(cbind, lapply(draws_beta$post_warmup_draws, function(x) do.call(rbind, x[draw_start:draw_end])))
    resp_pt_mean <- rowMeans(resp_draws)
    
    LL_0 <- log(prob_bad_use) # Vector(ndraws)
    LL_U <- log(1 - prob_bad_use)
    LL_preds <- c(0,0) # LL of point estimate, mean of draw preds
    sum_dep <- 0
    id_specs <- data_stan$id_ranges[id,] # task_beg, task_end, row_beg, row_end, nrows
    for (t in (id_specs[1]:id_specs[2])){
      row_beg <- data_stan$start[t]
      row_end <- data_stan$end[t]
      ind <- data_stan$ind[row_beg:row_end,]
      pred_draws_t <- predict_task(t, ind, resp_draws, scale_factor[,data_stan$task_scale_group[t]], task_type,
                                   inv_logit_thresh, draws = TRUE, pred_list_more) # n_alt x n_draws for task 
      pred_pt_t <-    predict_task(t, ind, resp_pt_mean, scale_factor_mean[data_stan$task_scale_group[t]], task_type,
                                   inv_logit_thresh, draws = FALSE, pred_list_more)
      preds <- cbind(pred_pt_t, rowMeans(pred_draws_t)) # For draws, take mean of n_alt x n_draws
      pred_all[row_beg:row_end,] <- preds
      dep <- data_stan$dep[row_beg:row_end] * data_stan$wts[t] * (1 - data_stan$holdout[t]) # Don't count holdouts
      LL_0 <- LL_0 + (sum(dep) * -log(length(dep))) # LL of 0 vector
      LL_U <- LL_U + (t(log(pred_draws_t)) %*% dep) # LL of each draw
      LL_preds <- LL_preds +  dep %*% log(preds)
      sum_dep <- sum_dep + sum(dep)
    }
    
    id_prob_good <- mean(1/(1 + exp(LL_0 - LL_U))) # = e^LL_U/(e^LL_U + e^LL_0)
    if(is.null(prob_bad)) id_prob_good <- NA
    rlh <- exp(LL_preds/sum_dep)
    fit_id[id,] <- c(rlh, id_prob_good, NA)
    utilities[id,] <- resp_pt_mean
    colnames(fit_id) <- c("rlh_pt","rlh_draws","prob_good_hb","prob_good_LC")
    colnames(pred_all) <- c("pred_pt","pred_draws")
    
  }
  # Export and save
  if (!is.null(LC_fit)) fit_id[,4] <- LC_fit$LC_prob_good
  util_name <- paste0(out_prefix,"_utilities_r.csv")
  pred_name <- paste0(out_prefix,"_preds.csv")
  fit_name <- paste0(out_prefix,"_Hit_LL.csv")
  failcon_name <- paste0(out_prefix,"_utilities_failcon.csv")
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred.csv")
  message(paste0(
    "Saving: \n",
    " Respondent mean utilities: ", util_name, "\n",
    " Predictions for data     : ", pred_name, "\n",
    " Observed vs Predicted    : ", obs_vs_pred_name, "\n",
    " Hit Rate and LogLike     : ", fit_name
  ))
  utilities_r <- utilities %*% t(data_stan$code_master)
  header <- data.frame(id = data_stan$resp_id, fit_id)
  write.table(cbind(header, utilities_r), file = file.path(dir_run, util_name), sep = ",", na = ".", row.names = FALSE)
  
  pred_all_export <- cbind(data_stan$idtask, holdout = row_holdout, wt = row_weights, dep = data_stan$dep, pred_pt = pred_all[,1], pred_draws = pred_all[,2])
  if (ncol(data_stan$ind_levels) >0){
    obs_vs_pred <- obs_vs_pred(pred_all_export[,4:6], data_stan$ind_levels)
    obs_vs_pred2 <- obs_vs_pred(pred_all_export[,c(4,5,7)], data_stan$ind_levels)
    obs_vs_pred$pred_per_draws <- obs_vs_pred2$pred_per
    write.table(obs_vs_pred, file = file.path(dir_run, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  } else message("No Categorical Variables found for observed vs predicted")
  if ("row_in" %in% names(data_stan)){
    pred_all_export <- cbind(row_in = data_stan$row_in, pred_all_export)
    pred_all_export <- pred_all_export[order(data_stan$row_in),]
  } else pred_all_export <- cbind(row_in = NA, pred_all_export)
  write.table(pred_all_export, file = file.path(dir_run, pred_name), sep = ",", na = ".", row.names = FALSE)
  
  fit_hit_LL <- hit_rate_LL(pred_all_export ,col_id = 2,col_task = 3, col_dep = 6,
                            col_pred = 7, col_split = 4, col_wts = 5)
  fit_hit_LL <- cbind(type = "MeanPts", data.frame(fit_hit_LL))
  fit_hit_LL2 <- hit_rate_LL(pred_all_export ,col_id = 2,col_task = 3, col_dep = 6,
                             col_pred = 8, col_split = 4, col_wts = 5)
  fit_hit_LL2 <- cbind(type = "Draws", data.frame(fit_hit_LL2))
  fit_hit_LL <- rbind(fit_hit_LL,fit_hit_LL2)
  write.table(fit_hit_LL, file = file.path(dir_run, fit_name), sep = ",", na = ".", row.names = FALSE) 
  write.table(scale_factor_mean, file = file.path(dir_run, paste0(out_prefix,"_scale_factor_mean.csv")), sep = ",", na = ".", row.names = FALSE) 
  
  # Check if utilities meet constraints
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  bad_ids <- rowSums(((utilities %*% t(con_matrix)) < 0)) > 0
  if (sum(bad_ids) > 0){
    message(paste0(sum(bad_ids), " Respondents had reversals from constraints.\n",
                   "Reversals saved to: ", failcon_name))
    write.table(cbind(header, utilities_r)[bad_ids,], file = file.path(dir_run, failcon_name), sep = ",", na = ".", row.names = FALSE)
  } else message(" All respondent mean utilities obey constraints")
  
  return(result = list(utilities = utilities, fit_id = fit_id, pred_all = pred_all, scale_factor_pts = scale_factor_mean))
}

env_stan$process_HB <- function(data_stan, data_model, control_code, meta_data,
                                predict_task = predict_task_base, pred_list_more = NULL,
                                task_type = NULL, inv_logit_thresh = NULL){
  if (min(meta_data$return_codes) == 0){
    cat("Reading draws from Stan csv output into R (large files take time)...")
    draws_beta <- read_cmdstan_csv(meta_data$output_files, variables = "beta_ind", format = "draws_list")
    draws_beta$warmup_draws <- NULL # to save space
    checkconverge_export(draws_beta, vnames = colnames(data_stan$code_master), control_code$out_prefix,
                         control_code$dir_run, control_code$export_draws_betas, control_code$export_draws_means, makepdf = ("draws_beta" %in% control_code$makepdf))
    if (data_model$check_bad == 1){
      prob_bad <- read_cmdstan_csv(meta_data$output_files, variables = "prob_bad", format = "draws_list")
      prob_bad <- do.call(rbind, lapply(prob_bad$post_warmup_draws, function(x) do.call(cbind, x)))  # draws x n_scales
      check_draws_vector(meta_data$output_files, "prob_bad", "prob_bad_agg", control_code$dir_run, control_code$out_prefix,makepdf = ("prob_bad" %in% control_code$makepdf))
    } else prob_bad <- NULL
    scale_factor_list <- read_cmdstan_csv(meta_data$output_files, variables = "scale_factor_final", format = "draws_list")
    scale_factor <- do.call(rbind, lapply(scale_factor_list$post_warmup_draws, function(x) do.call(cbind, x)))  # draws x n_scales
    if (max(data_stan$task_scale_group > 1)) check_draws_vector(meta_data$output_files, "z_log_scale_xref",
                                                                data_stan$resp_id, control_code$dir_run, control_code$out_prefix,makepdf = ("scale_factor" %in% control_code$makepdf))
    out_prefix <- control_code$out_prefix
    dir_run <- control_code$dir_run
    P <- data_stan$P
    row_weights <- data_stan$wts[data_stan$idtask_r] 
    row_holdout <- data_stan$holdout[data_stan$idtask_r]
    n_holdouts <- sum(data_stan$holdout)
  
    utilities <- matrix(NA, data_stan$I, data_stan$P) # Hold utilities
    pred_all <- matrix(NA, data_stan$N, 2) # Hold predictions points, draws
    fit_id <- matrix(NA, data_stan$I, 4) # RlH (pts, draws), Prob_Good (HB, LC)
    if(is.null(prob_bad)){
     prob_bad_use <- rep(.0001, draws_beta$metadata$iter_sampling * length(draws_beta$post_warmup_draws))
    } else prob_bad_use <- prob_bad
    scale_factor <- as.matrix(scale_factor)
    scale_factor_mean <- exp(colMeans(log(scale_factor)))
    for (id in (1:data_stan$I)){
     draw_start <- 1 + (id-1) * P
     draw_end <- draw_start + P - 1
     resp_draws <- do.call(cbind, lapply(draws_beta$post_warmup_draws, function(x) do.call(rbind, x[draw_start:draw_end])))
     resp_pt_mean <- rowMeans(resp_draws)
    
     LL_0 <- log(prob_bad_use) # Vector(ndraws)
     LL_U <- log(1 - prob_bad_use)
     LL_preds <- c(0,0) # LL of point estimate, mean of draw preds
     sum_dep <- 0
     id_specs <- data_stan$id_ranges[id,] # task_beg, task_end, row_beg, row_end, nrows
     for (t in (id_specs[1]:id_specs[2])){
      row_beg <- data_stan$start[t]
      row_end <- data_stan$end[t]
      ind <- data_stan$ind[row_beg:row_end,]
      pred_draws_t <- predict_task(t, ind, resp_draws, scale_factor[,data_stan$task_scale_group[t]], task_type,
                                   inv_logit_thresh, draws = TRUE, pred_list_more) # n_alt x n_draws for task 
      pred_pt_t <-    predict_task(t, ind, resp_pt_mean, scale_factor_mean[data_stan$task_scale_group[t]], task_type,
                                   inv_logit_thresh, draws = FALSE, pred_list_more)
      preds <- cbind(pred_pt_t, rowMeans(pred_draws_t)) # For draws, take mean of n_alt x n_draws
      pred_all[row_beg:row_end,] <- preds
      dep <- data_stan$dep[row_beg:row_end] * data_stan$wts[t] * (1 - data_stan$holdout[t]) # Don't count holdouts
      LL_0 <- LL_0 + (sum(dep) * -log(length(dep))) # LL of 0 vector
      LL_U <- LL_U + (t(log(pred_draws_t)) %*% dep) # LL of each draw
      LL_preds <- LL_preds +  dep %*% log(preds)
      sum_dep <- sum_dep + sum(dep)
    }
    
    id_prob_good <- mean(1/(1 + exp(LL_0 - LL_U))) # = e^LL_U/(e^LL_U + e^LL_0)
    if(is.null(prob_bad)) id_prob_good <- NA
    rlh <- exp(LL_preds/sum_dep)
    fit_id[id,] <- c(rlh, id_prob_good, NA)
    utilities[id,] <- resp_pt_mean
    colnames(fit_id) <- c("rlh_pt","rlh_draws","prob_good_hb","prob_good_LC")
    colnames(pred_all) <- c("pred_pt","pred_draws")
    
  }
  # Export and save
  if (!is.null(LC_fit)) fit_id[,4] <- LC_fit$LC_prob_good
  util_name <- paste0(out_prefix,"_utilities_r.csv")
  pred_name <- paste0(out_prefix,"_preds.csv")
  fit_name <- paste0(out_prefix,"_Hit_LL.csv")
  failcon_name <- paste0(out_prefix,"_utilities_failcon.csv")
  obs_vs_pred_name <- paste0(out_prefix,"_obs_vs_pred.csv")
  message(paste0(
    "Saving: \n",
    " Respondent mean utilities: ", util_name, "\n",
    " Predictions for data     : ", pred_name, "\n",
    " Observed vs Predicted    : ", obs_vs_pred_name, "\n",
    " Hit Rate and LogLike     : ", fit_name
  ))
  utilities_r <- utilities %*% t(data_stan$code_master)
  header <- data.frame(id = data_stan$resp_id, fit_id)
  write.table(cbind(header, utilities_r), file = file.path(dir_run, util_name), sep = ",", na = ".", row.names = FALSE)
  
  pred_all_export <- cbind(data_stan$idtask, holdout = row_holdout, wt = row_weights, dep = data_stan$dep, pred_pt = pred_all[,1], pred_draws = pred_all[,2])
  if (ncol(data_stan$ind_levels) >0){
    obs_vs_pred <- obs_vs_pred(pred_all_export[,4:6], data_stan$ind_levels)
    obs_vs_pred2 <- obs_vs_pred(pred_all_export[,c(4,5,7)], data_stan$ind_levels)
    obs_vs_pred$pred_per_draws <- obs_vs_pred2$pred_per
    write.table(obs_vs_pred, file = file.path(dir_run, obs_vs_pred_name), sep = ",", na = ".", row.names = FALSE)
  } else message("No Categorical Variables found for observed vs predicted")
  if ("row_in" %in% names(data_stan)){
    pred_all_export <- cbind(row_in = data_stan$row_in, pred_all_export)
    pred_all_export <- pred_all_export[order(data_stan$row_in),]
  } else pred_all_export <- cbind(row_in = NA, pred_all_export)
  write.table(pred_all_export, file = file.path(dir_run, pred_name), sep = ",", na = ".", row.names = FALSE)
  
  fit_hit_LL <- hit_rate_LL(pred_all_export ,col_id = 2,col_task = 3, col_dep = 6,
                            col_pred = 7, col_split = 4, col_wts = 5)
  fit_hit_LL <- cbind(type = "MeanPts", data.frame(fit_hit_LL))
  fit_hit_LL2 <- hit_rate_LL(pred_all_export ,col_id = 2,col_task = 3, col_dep = 6,
                             col_pred = 8, col_split = 4, col_wts = 5)
  fit_hit_LL2 <- cbind(type = "Draws", data.frame(fit_hit_LL2))
  fit_hit_LL <- rbind(fit_hit_LL,fit_hit_LL2)
  write.table(fit_hit_LL, file = file.path(dir_run, fit_name), sep = ",", na = ".", row.names = FALSE) 
  write.table(scale_factor_mean, file = file.path(dir_run, paste0(out_prefix,"_scale_factor_mean.csv")), sep = ",", na = ".", row.names = FALSE) 
  
  # Check if utilities meet constraints
  con_matrix <- diag(data_stan$con_sign)
  con_matrix <- rbind(con_matrix[rowSums(con_matrix !=0) > 0,,drop = FALSE], data_stan$paircon_matrix)
  bad_ids <- rowSums(((utilities %*% t(con_matrix)) < 0)) > 0
  if (sum(bad_ids) > 0){
    message(paste0(sum(bad_ids), " Respondents had reversals from constraints.\n",
                   "Reversals saved to: ", failcon_name))
    write.table(cbind(header, utilities_r)[bad_ids,], file = file.path(dir_run, failcon_name), sep = ",", na = ".", row.names = FALSE)
  } else message(" All respondent mean utilities obey constraints")
  
  return(result = list(utilities = utilities, fit_id = fit_id, pred_all = pred_all, scale_factor_draws = scale_factor, scale_factor_pts = scale_factor_mean, draws_beta = draws_beta))
  } else message("Stan Estimation Did not Finish")
}

env_eb$numder_2 <- function(x, pos, delta = .01){
  # 2nd derivative
  xup <- x
  xup[pos] <- x[pos] + delta
  xdown <- x
  xdown[pos] <- x[pos] - delta
  up <- model_list$func$min(xup, data_list, model_env)
  down <- model_list$func$min(xdown, data_list, model_env)
  result <- (up + down - (2* model_list$func$min(x, data_list, model_env)))/(delta^2)
  return(result)
}

env_eb$agg_solve <- function(data_list, model_list, fish_inf = FALSE) {
  model_env <- list2env(model_list, parent = emptyenv())
  ksolve <- constrOptim(theta = model_env$x0, f = model_env$func$min, grad = model_env$func$gr, ui = model_env$con[, -1], ci = model_env$con[, 1], mu = 1e-02, control = list(trace = 1, REPORT = 1),
                        method = "BFGS", outer.iterations = 100, outer.eps = 1e-05, hessian = fish_inf,
                        data_list = data_list,
                        model_env = model_env)
  ksolve$nresp <- length(data_list$resp_id)
  if (fish_inf) {
    ksolve$fish_inf <- (ksolve$hessian / ksolve$nresp) # unit fisher
  }
  return(ksolve)
}

env_eb$SolveID <- function(idseq, data_list, model_env, sample_est = TRUE) {
  id_list <- list()
  id_filter <- (data_list$match_id == idseq)
  task_id <- data_list$idtask_r[id_filter]
  task_id_u <- unique(task_id)
  id_list$idtask_r <- (match(data.frame(t(task_id)), data.frame(t(task_id_u)))) # unique tasks
  id_list$ind = data_list$ind[id_filter,]
  id_list$dep <- as.matrix(data_list$dep[id_filter]) # Need matrix
  sample_est <- id_filter & sample_est
  wts_use <- data_list$wts
  wts_use[!sample_est] <- 0 # set dep_use to 0 for cases we want to exclude    
  id_list$wts <- wts_use[id_filter]
  #id_list$util_mult <- data_list$util_mult[id_filter] #ADDED
  sample_est <- sample_est[id_filter]
  eb_solve <- constrOptim(theta = model_env$x0, f = model_env$func$min, grad = model_env$func$gr, ui = model_env$con[, -1], ci = model_env$con[, 1], mu = 1e-02, control = list(),
                          method = "BFGS", outer.iterations = 100, outer.eps = 1e-05, 
                          data_list = id_list,
                          model_env = model_env)
  rlh <- exp(-1 * (eb_solve$value/sum(id_list$dep *id_list$wts)))
  betas <- c(data_list$resp_id[idseq], rlh, round(eb_solve$par, 8))
  names(betas) <- c("id", "rlh", colnames(data_list$ind))
  predprob <- model_env$func$pred(x = eb_solve$par, data_list = id_list)
  predprob <- cbind(data_list$idtask[id_filter,], sample_est, predprob, dep = id_list$dep, wts = data_list$wts[id_filter])
  result <- list(betas = betas, predprob = predprob)
  return(result)
}  

env_eb$SolveID_Nest <- function(idseq, emp_bayes_list, scale_fish = 1, sample_est = TRUE, fn = LLn_EB_GetBeta) {
  
  ##### New Code ####
  create_nest_link(data_nest_in[id_filter,], nest_struc, as.matrix(data_list$idtask_r[id_filter]))
  #######################  Also changed fn = LLn_GetBeta
  return(result)
}

env_eb$mleb <- function(data_list, model_list, mleb_control){
  .GlobalEnv$model_env <- list2env(model_list, parent = emptyenv())
  model_env$SolveID <- env_eb$SolveID
  model_env$mleb_result <- list()
  model_env$mleb_result$timebeg <- Sys.time()
  Call <- match.call(expand.dots = FALSE)
  model_env$mleb_result$name_data_list <- as.list(Call)[[2]] # data_list used
  model_env$mleb_result$name_model_list <- as.list(Call)[[3]] # model_list used
  model_env$mleb_result$model_list <- model_list
  
  if (mleb_control$solveonly) {
    eb_betas <- list(length(data_list$resp_id))
    eb_betas <- foreach(i = 1:length(data_list$resp_id)) %dopar% {
      eb_betas[[i]] <- model_env$SolveID(idseq = i, data_list=data_list, model_env = model_env)
    } 
    model_env$mleb_result$eb_betas <- do.call(rbind, lapply(eb_betas, function(x) x$betas))
    model_env$mleb_result$predprob <- do.call(rbind, lapply(eb_betas, function(x) x$predprob))
    colnames( model_env$mleb_result$predprob) <- c("id", "task", "est", "pred", "dep", "wt")
  } else {
    kscale_hist <- matrix(NA, mleb_control$maxiter, 5)
    colnames(kscale_hist) <- c("scale", "rlh_cal", "rlh_inthold", "rlh_exthold", "cov_size")
    model_env$mleb_result$kscale_hist <- kscale_hist
    model_env$mleb_result$eb_betas_hist <- list(); model_env$mleb_result$prior_cov_hist <- list(); model_env$mleb_result$prior_alpha_hist <- list(); model_env$mleb_result$scalefind_detail <- list()
    internal_hold <- jack_knife(idtask = data_list$idtask, num_round = mleb_control$hold_tasks, num_resp = mleb_control$hold_resp, hold_poss = mleb_control$hold_poss) # select holdout tasks across iterations  
    rlh_fac <- 1 # Initial factor for how well cov fits
    
    plot_setup()     
    
    cat("\014")
    message( "   Created Environment model_env for computations")
    cat("-----------------------------------------------------------\n")
    print(kscale_hist[NULL,], right = TRUE, row.names = FALSE)
    cat("------------------------------------------------------------\n")
    iter <- 1
    converge <- FALSE
    while (iter <= mleb_control$maxiter & !converge) {
      time_beg_iter <- Sys.time()
      message(paste0("   Iteration ", iter))
      message("   Estimating Scale Factor")
      screen(3)
      plot(1,.5, xlab = paste0("Scale Iter ", iter), ylab = "RLH", type = "n") # blank plot
      text(x =1, y=.5, labels = "Testing Initial Values", col = "blue")
      cal_filter <- jack_knife(idtask = data_list$idtask, num_round = mleb_control$cal_tasks, num_resp = mleb_control$cal_resp, hold_poss = mleb_control$hold_poss) # select calibration tasks per respondent
      resp_est <- cal_filter$resp_pick #
      cal_filter <- cal_filter$jk_filter # Boolean for all rows
      ub_k <- 2
      if (iter == 1) {ub_k <- 10}
      scale_find <- line_search_gr(f = function(x) - 1 * ll_out(x, cal_filter, resp_est,
                                                                data_list=data_list, model_env = model_env),
                                   lb = .1, ub = ub_k, tolerance = mleb_control$tolerance, iter = iter
      ) # 
      
      message("   Computing Internal Holdouts")
      kscale <- scale_find$optim_x
      model_env$prior$scale <- kscale # Set scale of cov_inv    
      rlh_hold <- exp(ll_out(kscale = kscale, cal_filter = internal_hold$jk_filter, resp_est = internal_hold$resp_pick,
                             data_list=data_list, model_env = model_env))
      if (iter == 1) rlh_base <- rlh_hold
      rlh_fac <- rlh_hold
      
      message("   Estimating Respondent Level Betas")
      #Solution of Betas
      eb_betas <- foreach(i = 1:length(data_list$resp_id), .combine='rbind') %dopar% {
        result <- model_env$SolveID(idseq = i, data_list=data_list, model_env = model_env)
        result$betas
      }
      
      # store results
      #extfit <- test_util(eb_betas[, 1], eb_betas[, -1:-2], data_conj_hold)$rlh_mean # optional ext holdout
      extfit <- NA
      prior_cov <- solve(model_env$prior$scale * model_env$prior$cov_inv)
      
      cov_size <- mean(diag(prior_cov))
      kscale_hist[iter,] <- c(kscale, scale_find$optim_y, rlh_hold, extfit, cov_size)
      model_env$mleb_result$kscale_hist <- kscale_hist
      model_env$mleb_result$eb_betas_hist[[iter]] <- eb_betas    
      model_env$mleb_result$prior_cov_hist[[iter]] <- prior_cov
      model_env$mleb_result$prior_alpha_hist[[iter]] <- model_env$prior$alpha
      model_env$mleb_result$scalefind_detail[[iter]] <- scale_find$pts
      
      #update priors
      betas_upper <- (eb_betas[, -1:-2])[, model_env$prior$upper_model]
      alpha_r <- next_alpha(betas_upper, model_env$prior$alpha)
      cov_r <- next_cov(prior_alpha = model_env$prior$alpha, betas = betas_upper, prior_cov = prior_cov, v0 = (ncol(betas_upper) + 1 + (sqrt(rlh_hold) * nrow(betas_upper)))) # 
      plot_iter(iter, model_list, setup = TRUE)
      dev.copy(device=pdf, file.path(mleb_control$dir_pdf, paste0("MLEB_Plot_iter_", iter, ".pdf")), width = 10, height = 6)
      dev.off()
      
      # model_env$prior$cov_inv <- solve(cov_r)
      ksvd <- svd(cov_r)
      ksvd$d[ksvd$d < 0] <- 1e-15
      model_env$prior$cov_inv <- ksvd$v %*% diag(1/ksvd$d) %*% t(ksvd$u)
      
      model_env$prior$alpha <- alpha_r
      model_env$prior$scale <- 1 # set to 1 just for cleanliness
      model_env$x0 <- alpha_r # Update x0 to match alpha (Nov 2020)
      save(mleb_result, file = file.path(mleb_control$dir_pdf, "mleb_result.RData"), envir = model_env)
      
      # print
      #cat(rep("\n", 3))
      cat("\014")
      cat("-----------------------------------------------------------\n")
      print(kscale_hist[1:iter,], right = TRUE, row.names = FALSE)
      cat("------------------------------------------------------------\n")
      tpi <- (Sys.time() - time_beg_iter) 
      units(tpi) <- "mins"
      message(paste0("  Time for iteration ", iter, " was: ", format(tpi, digits = 3)))
      
      # Check for convergence
      fit <- kscale_hist[,3]
      fit <- fit[!is.na(fit)]
      best_fit <- max(fit)
      best_iter <- match(best_fit, fit) # Picks iteration with best rlh
      if ((length(fit) - best_iter) >= mleb_control$conv_n){
        message(" Completed: Holdout Fit no longer improved")
        converge <- TRUE
      } 
      message("")
      iter <- iter + 1
    }}
}

env_eb$plot_setup <- function(x) {
  # setup graphics
  if (!is.null(dev.list())) dev.off()    
  m <- rbind(c(0, 0.5, 0.55, 1), c(0.5, 1, 0.55, 1),
             c(0, 0.5, 0, 0.55), c(0.5, 1, 0, 0.55), c(.5,.55,.57,.66))
  split.screen(m)
  screen(1)
  par(mar = c(4, 4, 1, 1))
  screen(2)
  par(mar = c(4, 4, 1, 1))
  screen(3)
  par(mar = c(4, 4, 2, 1))
  screen(4)
  par(mar = c(4, 4, 2, 1))
}


env_eb$plot_iter <- function(iter, model_list, setup = FALSE){
  if (setup) plot_setup()
  scalepts <- model_env$mleb_result$scalefind_detail[[iter]]
  rlh <- model_env$mleb_result$kscale_hist[iter,3]
  a <- formatC(rlh, digits = 4, format = "f")
  a <- substr(a,2,nchar(a))
  
  screen(3)
  plot(scalepts,xlab=paste0("Scale Iter ", iter), ylab= "RLH", col = "blue", fg = "blue")
  
  if (iter > 1){
    screen(1)
    plot(model_list$prior$alpha, model_env$mleb_result$prior_alpha_hist[[iter]],xlab="Alpha Initial",
         ylab=paste0("Alpha Iter ", iter))
    
    screen(2)
    plot(model_env$mleb_result$prior_alpha_hist[[iter - 1]], 
         model_env$mleb_result$prior_alpha_hist[[iter]],
         xlab=paste0("Alpha Iter ", iter - 1),
         ylab=paste0("Alpha Iter ", iter))
    
    screen(4)
    cov_keep <- lower.tri(model_env$mleb_result$prior_cov_hist[[iter]], diag = TRUE)
    plot(model_env$mleb_result$prior_cov_hist[[iter-1]][cov_keep], 
         model_env$mleb_result$prior_cov_hist[[iter]][cov_keep], 
         xlab=paste0("Cov Iter ", iter-1), ylab=paste0("Cov Iter ", iter))
  }
  screen(5)
  mtext(paste0("Iter ", iter, "\n", "Hold RLH ","\n", a), side = 3, col = "coral3")
  
}
###############  Updating Upper Level ###############
env_eb$next_alpha <- function(betas, prior_alpha, k0 = nrow(betas)) {
  result <- ((k0 * prior_alpha) + colSums(betas))/(k0 + nrow(betas))
  return(result)
}

env_eb$chol_qr <- function(kmatrix) {
  ksvd <- svd(kmatrix)
  k_qr <- qr(t((ksvd$u %*% sqrt(diag(ksvd$d)))))
  Q <- qr.Q(k_qr)
  R <- qr.R(k_qr)
  L <- (Q %*% R)
  return(L) # equivalent of upper traingle cholesky
}

env_eb$next_cov <- function(prior_alpha, betas,
                            prior_cov = diag(length(alpha)),
                            v0 = nrow(betas),
                            k0 = nrow(betas)) {
  # prior mean and cov are alpha and prior_cov
  # betas are respondent level betas
  # v0 is total degrees of freedom (must be > ncol(betas)), k0 is n size basis of alpha
  
  n <- nrow(betas)
  p <- ncol(betas)
  add_z <- matrix(rnorm(n * p)/1000, nrow = n, ncol = p) # Add small random error
  betas_r <- betas + add_z
  
  xbar <- colMeans(betas_r)
  xbar_diff <- (xbar - prior_alpha) 
  alpha_covn <- ((k0*n)/(k0+n))* (xbar_diff %*% t(xbar_diff)) # comp 1
  
  ab_diff <- t(apply(betas_r, 1, function(x) x - xbar))
  beta_covn <- t(ab_diff) %*% ab_diff # comp 2
  
  cov_new_mean <- (prior_cov * v0 + beta_covn + alpha_covn)/(v0 + n - p - 1) #denom defaults to 2n
  return(cov_new_mean) # returns new cov
}

###############  Scaling Factor ###############
env_eb$ll_out <- function(kscale, cal_filter, resp_est, 
                          data_list=data_list, model_env = model_env){
  # jack knife 
  # Global: cal_filter, resp_est, emp_bayes
  result <- list()
  model_env$prior$scale <- kscale
  #data_list <- data_list
  for (i in 1:ncol(cal_filter)){
    pred_s <- foreach(j = 1:length(resp_est), .combine='rbind') %dopar% {
      model_env$SolveID(idseq = resp_est[j], data_list=data_list, model_env = model_env, sample_est = cal_filter[, i])$predprob
    } 
    pred_s_hold <- pred_s[!pred_s$sample_est,]
    mean_ll <- sum((log(pred_s_hold$predprob) * pred_s_hold$dep * pred_s_hold$wts)) / sum(pred_s_hold$dep * pred_s_hold$wts)
    result[[i]] <- mean_ll
  }
  return(mean(unlist(result)))
}

env_eb$line_search_gr <- function(f, lb, ub, tolerance, iter) {
  # ASSUMES MINIMIZATION
  #golden ratio line search
  gr <- 2/(sqrt(5) + 1)
  
  ### Use the golden ratio to set the initial test points
  x1 <- ub - gr*(ub - lb)
  x2 <- lb + gr*(ub - lb)
  
  ### Evaluate the function at the test points
  f1 <- f(x = x1)
  f2 <- f(x = x2)
  
  iteration <- 0
  hist <- c(lb, ub, x1, x2, exp(f1*-1), exp(f2*-1))
  
  while (abs(ub - lb) > tolerance) {
    #iter <- iteration + 1
    if (f2 > f1) {
      ### Set the new upper bound
      ub <- x2
      x2 <- x1
      f2 <- f1
      
      ### Set the new lower test point
      x1 = ub - gr*(ub - lb)
      f1 <- f(x1)
    } else {
      ### Set the new lower bound
      lb <- x1
      x1 <- x2
      f1 <- f2
      
      ### Set the new upper test point
      x2 <- lb + gr*(ub - lb)
      f2 <- f(x2)
    }
    newpts <- c(lb, ub, x1, x2, exp(f1*-1), exp(f2*-1))
    hist <- rbind(hist, newpts)
    pts <- rbind(hist[,c(3,5)],hist[,c(4,6)])
    colnames(pts) <- c("kscale", "rlh")
    if (nrow(pts) > 2){
      screen(3)
      plot(pts,xlab=paste0("Scale Iter ", iter), ylab= "RLH", col = "blue", fg = "blue")
    }
  }
  optim_x <- (lb + ub)/2
  optim_y <- exp(f(optim_x) * -1)
  colnames(hist) <- c("lb","ub","x1","x2", "rlh1", "rlh2")
  pts <- rbind(hist[,c(3,5)],hist[,c(4,6)],c(optim_x,optim_y))
  result <- list(optim_x = optim_x, optim_y = optim_y, pts = pts, hist = hist)
}

env_eb$jack_knife <- function(idtask, num_round = 1, num_resp = NULL, hold_poss = TRUE) {
  # num_resp = NULL gives all respondents
  # hold_poss means task is possible holdout
  idtask_uhold <- unique(idtask[hold_poss,]) # Only use those 
  ksplit <- split(idtask_uhold, idtask_uhold[, 1])
  if (is.null(num_resp)) {
    resp_pick <- 1:length(ksplit)
  } else {
    num_resp <- min(num_resp,length(ksplit))
    resp_pick <- sample(length(ksplit),num_resp)      
  }
  pick_hold <- lapply(ksplit[resp_pick], function(x) x[sample(1:nrow(x), num_round),])
  jk_filter <- sapply(1:num_round, function(i){
    this_round <- do.call(rbind,lapply(pick_hold, function(x) x[i,]))
    kmatch <- match(data.frame(t(idtask)), data.frame(t(this_round))) # match holdout tasks
    result <- is.na(kmatch)
    return(result)
  })
  return(list(jk_filter = jk_filter, resp_pick = resp_pick))
}

env_eb$best_x_fromy <- function(x, y) {
  #x,y vectors - find wt avg of x for best y values +- 1
  best_y <- match(max(y), y)
  result <- x[best_y] # initial
  if ((best_y > 1) & (best_y < length(y))) {
    seq_use <- (best_y - 1):(best_y + 1)
    kcoef <- cbind((x[seq_use] ^ 2), x[seq_use], 1)
    betas <- solve(kcoef, y[seq_use])
    der_0 <- -1 * betas[2] / (2 * betas[1]) # -b/2a is optimal x for quadratic
    if (abs(der_0 - result) < .5) result <- der_0
  }
  return(result)
}

###############  General Prep ###############
env_eb$prep_file <- function(idtaskdep, indcode_list, train = TRUE, other_data = NULL) {
  sort_order <- order(idtaskdep[, 1], idtaskdep[, 2])
  sort_order[!train] <- 0 # Non-training gets order = 0, which removes
  ind <- as.matrix(indcode_list$indcode[sort_order,])
  
  dep <- as.vector(as.matrix(idtaskdep[sort_order, 3]))
  idtask <- data.frame(idtaskdep[sort_order, 1:2])
  idtask_u <- as.matrix(unique(idtask))
  idtask_r <- (match(data.frame(t(idtask)), data.frame(t(idtask_u)))) # unique tasks
  resp_id <- as.vector(unique(idtask_u[, 1]))
  match_id <- match(idtask[, 1], as.matrix(resp_id))
  # Next 3 lines recodes dep to sum to 1
  depsum <- rowsum(dep, idtask_r) # Sum dep each task
  depsum_match <- (depsum[idtask_r,]) # Map Sum to rows
  dep <- dep / depsum_match # sum of dep will add to 1
  # Recode NAs to 0
  ind[is.na(ind)] <- 0
  dep[is.na(dep)] <- 0
  wts <- rep(1, length(dep)) # initial weights are 1
  
  # Add Stan stuff
  end <- c(which(diff(idtask_r)!=0), length(idtask_r))
  start <- c(1, end[-length(end)]+1)
  # Return list of data (depends on whether other data was also chosen)
  if (!is.null(other_data)){
    return(list(tag = 0, N = nrow(ind), P = ncol(ind), T = max(idtask_r), I = length(resp_id),
                dep = dep, ind = ind, idtask = idtask, idtask_r = idtask_r, resp_id = resp_id, match_id = match_id,
                task_individual = match_id[start],
                start = start,
                end = end,
                prior_cov = indcode_list$indprior,
                code_master = indcode_list$code_master,
                wts = wts,
                other_data = as.matrix(other_data)[sort_order,])) # Only Added item in list vs below
  } else {
    return(list(tag = 0, N = nrow(ind), P = ncol(ind), T = max(idtask_r), I = length(resp_id),
                dep = dep, ind = ind, idtask = idtask, idtask_r = idtask_r, resp_id = resp_id, match_id = match_id,
                task_individual = match_id[start],
                start = start,
                end = end,
                prior_cov = indcode_list$indprior,
                code_master = indcode_list$code_master,
                wts = wts))
  }  
}

env_eb$con_trivial <- function(num_par, umin = -999){
  # Make trivial constraint: first util > -999
  con_null <- t(as.matrix(rep(0, 1 + num_par)))
  con_null2 <- con_null
  con_null[,1:2] <- c(-999, 1)
  con_null2[,1:2] <- c(-999, -1)
  con_null <- rbind(con_null, con_null2)
  return(con_null)
}


# Custom function for Unspoken.  Test weights in seq_test for each respondent
env_eb$SolveID_TestWts <- function(idseq, data_list, model_env, seq_test, hold_poss = TRUE) {
  # Scale for time by respondent based on LOO
  # Create data_list for specific respondent
  id_list <- list()
  id_filter <- (data_list$match_id == idseq)
  task_id <- data_list$idtask_r[id_filter] # Tasks for this id (all rows)
  task_id_u <- unique(task_id)
  id_list$idtask_r <- (match(data.frame(t(task_id)), data.frame(t(task_id_u)))) # unique tasks
  id_list$ind = data_list$ind[id_filter,]
  id_list$dep <- as.matrix(data_list$dep[id_filter]) # Need matrix
  
  task_hold_u <- unique(data_list$idtask_r[id_filter & hold_poss]) 
  cal_n <- min(floor(length(task_hold_u) * 1), 10) # 10 Validations is plenty 
  cal_tasks <- sample(task_hold_u, cal_n)
  sample_est <- as.matrix(!sapply(cal_tasks, function(x)(task_id == x)))
  # sample_est is matrix with >= 1 column specifying training and holdout
  wts_use <- data_list$wts
  id_list$wts <- wts_use[id_filter]
  
  # Test different weights
  ll_out <- function(time_scale) {
    time_r <- pnorm(time_scale * time_std_id[id_filter]) # Cumulative norm
    id_list_timewt <- id_list
    id_list_timewt$wts <- time_r / mean(time_r) #wts for ind model
    
    ksolve <- function(ksample, id_list_timewt_hold = id_list_timewt) {
      id_list_timewt_hold$dep <- id_list_timewt$dep * ksample 
      eb_solve <- constrOptim(theta = model_env$x0, f = model_env$func$min, grad = model_env$func$gr, ui = model_env$con[, -1], ci = model_env$con[, 1], mu = 1e-02, control = list(),
                              method = "BFGS", outer.iterations = 100, outer.eps = 1e-05, 
                              data_list = id_list_timewt_hold,
                              model_env = model_env)
      
      predprob <- model_env$func$pred(x = eb_solve$par, id_list_timewt_hold, model_env)
      betas <- c(data_list$resp_id[idseq], eb_solve$convergence, round(eb_solve$par, 8))
      result <- list(betas = betas, predprob = predprob)
      return(result)
    }
    solve_list <- apply(sample_est, 2, ksolve) #list for each column of sample_est
    
    # Compute fit to holdouts
    pred_all <- sapply(solve_list, function(x) x$predprob)
    dep_all <- apply(pred_all, 2, function(x) id_list$dep) # use all raw dep
    n_est <- sum(dep_all * sample_est)
    n_holdout <- sum(dep_all * !sample_est)
    ll_hold <- sum(log(pred_all) * dep_all * !sample_est) / n_holdout
    ll_est <- sum(log(pred_all) * dep_all * sample_est) / n_est
    return(exp(ll_hold))
  }
  
  result <- sapply(seq_test, ll_out) # seq_test is global
  return(result)
}

env_eb$JeffPrior <- function(xvec, data_list, model_list, ID_Var = FALSE, ObsFish = TRUE, ExpFish = FALSE, score_n = 300){
  result <- list()
  model_env <- list2env(model_list, parent = emptyenv())
  nresp <- length(data_list$resp_id)
  if (ID_Var) {
    score_id <- do.call(rbind, lapply(1:nresp, function(i){
      id_list <- list()
      id_filter <- (data_list$match_id == i)
      task_id <- data_list$idtask_r[id_filter]
      task_id_u <- unique(task_id)
      id_list$idtask_r <- (match(data.frame(t(task_id)), data.frame(t(task_id_u)))) # unique tasks
      id_list$ind = data_list$ind[id_filter,]
      id_list$dep <- as.matrix(data_list$dep[id_filter]) # Need matrix
      id_list$wts <- data_list$wts[id_filter]
      gr_sub <- as.vector(model_env$func$gr(xvec, id_list, model_env))
      return(gr_sub)
    }))
    result$ID_Var <- sqrt(apply(score_id, 2, var))
  }
  if (ObsFish) {
    message("Computing 2nd Derivatives (ObsFish) for Observed Fisher Information")
    derive_2 <- function(x, pos, delta = .01){
      xup <- x
      xup[pos] <- x[pos] + delta
      xdown <- x
      xdown[pos] <- x[pos] - delta
      up <- model_list$func$min(xup, data_list, model_env)
      down <- model_list$func$min(xdown, data_list, model_env)
      result <- (up + down - (2* model_env$func$min(x, data_list, model_env)))/(delta^2)
      return(result)
    }
    result$ObsFish <- sqrt(sapply(1:length(xvec), function(i) derive_2(xvec, i, .01))/nresp)
  }
  
  if (ExpFish){
    k_div <- 2 # 2 default means use half the tasks
    message(paste0("Computing Expected Fisher Information (ExpFish) from ", score_n, " scores" ))
    wts_in <- data_list$wts
    tot_n <- max(data_list$idtask_r)
    npicks <- floor(tot_n/2)
    score_samp <- do.call(rbind,lapply(1:score_n, function(x){
      task_pick <- sample(1:tot_n, npicks)
      data_list$wts <- (!is.na(match(data_list$idtask_r, task_pick))) * wts_in
      gr_sub <- as.vector(model_env$func$gr(xvec, data_list, model_env))
      return(gr_sub)
    }
    )) # end score
    result$ExpFish <- sqrt(apply(score_samp * k_div,2,var)/nresp)
  }
  
  JeffPrior <- as.data.frame(do.call(cbind, result))
  minvals <- round(apply(JeffPrior, 2, min), 4)
  message(paste0("Min values are: ", paste(minvals, collapse = " ")))
  return(JeffPrior)
}


# Do not attach multiple times -- detach first
attach(env_code)
attach(env_modfun)
attach(env_eb)
attach(env_stan)

