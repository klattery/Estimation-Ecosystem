#########################################################
#     Based on Estimation Ecosystem         
#     Created Nov 2021 for AWS Portal                                 
#     Kevin Lattery                                                       
#     Run Everything to create environments                  
#########################################################

# Create 3 Environments to Store Functions
# Attaching an environment multiple times creates duplicates
if (exists("env_code")) detach(env_code)
if (exists("env_modfun")) detach(env_modfun)
if (exists("env_eb")) detach(env_eb)
env_code <- new.env(parent = emptyenv())
env_modfun <- new.env(parent = emptyenv())
env_eb <- new.env(parent = emptyenv())

###############  Coding Functions Environment ################ 
env_code$setup_cores <- function(ncores){
  if (file.exists(".GlobalEnv$k_multi_core")){
    stopCluster(.GlobalEnv$k_multi_core)
    remove(.GlobalEnv$k_multi_core)
  }
  if (ncores >= 1){
    .GlobalEnv$k_multi_core <- makeCluster(min(ncores,detectCores()))
    registerDoParallel(.GlobalEnv$k_multi_core)
  } else {
    message(" No parallel threads. Setting ncores < 1 removes parallel threads ")
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
        reflev <- match(max(shows), shows)
      } 
      if (codetype %in% c(2, 3))code_matrix <- as.matrix(code_matrix[, -reflev, drop = FALSE]) # ref lev col dropped
      if (codetype == 3) code_matrix[reflev,] <- -1
    } else { # we have constraints
      cat_con <- code_cat_wcon(paircon,numlevs)
      code_matrix <- cat_con$code_matrix
      con_vec <- cat_con$con
      reflev <- cat_con$reflev
      pairs_add <- cat_con$pairs_add
      colnames(code_matrix) <- varnames[-reflev]
    }  
    outcode <- (code_matrix[newval,TRUE,drop = FALSE])
    if (is.null(paircon)){
      con_vec <- rep(0, ncol(code_matrix))
      pairs_add <- NULL
    } 
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
  return(list(outcode = outcode, code_matrix = code_matrix, con_sign = con_vec, vnames = varnames, reflev = reflev, prior = prior, pairs_add = pairs_add))
}

env_code$remove_implicits<-function(constraints){
  result <- c()
  constraints <- as.matrix(constraints)
  for (c in 1:nrow(constraints)){
    att = constraints[c,1]
    A = constraints[c,2]
    B = constraints[c,3]
    if(!(isImplicit(constraints, c, att, A, B))){
      result<-rbind(result,constraints[c,])
    }
  }
  return(result)
}
env_code$isImplicit <-function(constraints,conrow, att, A, B){
  nrofcons = nrow(constraints)
  vlag <- FALSE
  for (c in 1:nrofcons){
    if(constraints[c, 1] == att && constraints[c, 3] == B && c != conrow){ 
      bnew <- constraints[c, 2]
      if(A == bnew){
        vlag <- TRUE
        return(vlag)
      } else {
        vlag <- isImplicit(constraints,c, att, A, bnew)
      } 
    }
  }
  return(vlag)
}

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
          B <-constraints[c,3]
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
  result<-result[order(sapply(result,function(x) -length(x) ))]
  return(result)
}

# Nominal with constraints
env_code$code_cat_wcon <-function(constraints, numlevs){
  constraint_chain <- ordinal_chain(as.matrix(constraints)) # create list of chains
  missing <- NULL
  con<-rep(0,numlevs) # storage for constraint
  X<-matrix(0,numlevs,numlevs) # X will store the coded matrix
  for (constrings in 1:length(constraint_chain)){
    for (c in 1:(length(constraint_chain[[constrings]])-1)){
      A<-constraint_chain[[constrings]][c]
      B<-constraint_chain[[constrings]][c+1]
      diagA<-X[A,A]  # diagonoal element A
      sumB<-sum(X[B,]!=0) # Num of coded elements in row B
      if (diagA==0 && sumB==0){
        X[A,A]<-1
        X[B,A]<-1
        X[B,B]<-1
        con[B]<-1
      }
      if (diagA!=0 && sumB==0){
        X[B,]<-X[A,]
        X[B,B]<-1
        con[B]<-1
      }
      if (diagA==0 && sumB!=0){
        X[A,]<-X[B,]
        X[A,A]<--1
        con[A]<-1
      }
      if (diagA!=0 && sumB!=0){        
        missing<-rbind(missing,c(B,A))
      }
    }
  }
  for (lev in 1:numlevs){
    if (X[lev,lev]==0) X[lev,lev]<-1
  }
  kolsum<-colSums(X)
  reflev=match(max(kolsum),kolsum)
  X<-X[,-reflev]
  con<-con[-reflev]
  return(list(code_matrix=X,con=con,pairs_add=missing, reflev = reflev))
}


env_code$usercode1 <- function(kdata, vname, varout = NULL){
  if (is.null(varout)) varout <- vname
  outcode <- as.matrix(kdata[[vname]])
  colnames(outcode) <- varout
  code_matrix <- diag(1)
  colnames(code_matrix) <- varout
  prior <- code_matrix
  return(list(outcode = outcode, code_matrix = code_matrix, vnames = varout, prior = prior))
}


env_code$usercode <- function(kdata, kcol, varout = NULL){
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
  return(list(outcode = outcode, code_matrix = diag(ncol(outcode)), vnames = varout, prior = diag(ncol(outcode))))
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

env_code$ordcode <- function(kdata, vname, cut_pts = NULL, thermcode = TRUE, varout = NULL, setna = 0) {
  # xvec must be vector
  # cut_pts must be sequential vector from low to high.  if null then all values except c(0,NA)
  # NA and values outside cut_pts are set to "setna", default = 0
  # varout is prefix for varout_cut
  # uses function ordmatrix
  if (is.null(varout)) varout <- vname
  xvec <- kdata[[vname]]
  if (is.null(cut_pts)){
    cut_pts <- sort(unique(xvec)) # Unique values
    cut_pts <- cut_pts[!(cut_pts %in% c(0, NA))]  # But exclude values in vector of setna (default 0,NA)
  }
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
  return(list(outcode = ordcode, code_matrix = code_matrix, vnames = varnames, prior = diag(ncol(code_matrix))))
}

env_code$check_atts_constraints <- function(data_in, att_coding, constraints){
  # Check consistency of specified attributes and constraints with data
  result <- TRUE
  check_att <- att_coding[,1,drop = TRUE] %in% colnames(data_in)
  if (min(check_att) == 0){
    message("You specified attributes that are not in your data.  Please fix.")
    att_coding[!check_att,1]
    result <- FALSE
  }
  if (!is.null(constraints)){
    check_match <- constraints[,1, drop = TRUE] %in% att_coding[,1, drop = TRUE]
    if (min(check_match) == 0){
      message("You specified constraints that do no match any attribute in attribute file.  Please fix.")
      constraints[!check_match,1]
      result <- FALSE
    }
  }
  return(result)
}

env_code$indcode_spec_get <- function(data_in, att_coding,constraints){
  if (check_atts_constraints(data_in, att_coding,constraints)){
    catcode_types <- c("INDICATOR", "DUMMY","EFFECT","EFFECTS","NOMINAL")
    indcode_spec <- list(nrow(att_coding))
    for (i in 1:nrow(att_coding)){
      att_name <- att_coding[i,1,drop = TRUE]
      att_type <- toupper(att_coding[i,2,drop = TRUE]) # UPPERCASE
      if (att_type %in% catcode_types){
        codetype <- match(att_type, catcode_types)
        codetype <- min(codetype,3) # anything listed after effect will default to effect
        if (max(constraints[,1,drop = TRUE] %in% att_name) == 0) {
          indcode_spec[[i]] <- catcode(data_in, att_name, codetype) # No constraints
        } else{
          indcode_spec[[i]] <- catcode(data_in, att_name, codetype, paircon = constraints[constraints[,1] == att_name,])
        }    
      }
      if (att_type == "ORDINAL"){
        indcode_spec[[i]] <- ordcode(data_in, att_name)
        indcode_spec[[i]]$con_sign <- rep(att_coding[i,3], ncol(indcode_spec[[i]]$code_matrix))  
      }
      if (att_type == "USERSPECIFIED"){
        indcode_spec[[i]] <- usercode1(data_in, att_name)
        indcode_spec[[i]]$con_sign <- att_coding[i,3]
      }
    }  
    indcode_spec <- setNames(indcode_spec,att_coding[,1,drop = TRUE])
  } else indcode_spec <- "Error: Data file, attributes, constraints do not match" 
  return(indcode_spec)
}

env_code$list_to_matrix <- function(klist){
  # converts list of matrices to single matrix with blocks and 0's filling
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

env_code$make_codefiles <- function(indcode_list){
  # Converts list of codes to matrices:  # code_master, indcode, indprior
  result <- list()
  indcode_list <- indcode_list[sapply(indcode_list, length) > 0] # remove NULL elements
  names(indcode_list) <- NULL
  result$indcode <- do.call(cbind, lapply(indcode_list, function(x) x$outcode)) # coded variables 
  result$con_sign <- do.call(c, lapply(indcode_list, function(x) x$con_sign))
  result$code_master <- list_to_matrix(lapply(indcode_list, function(x) x$code_matrix))
  result$indprior <- list_to_matrix(lapply(indcode_list, function(x) x$prior))
  colnames(result$code_master) <- colnames(result$indcode) 
  rownames(result$code_master) <- do.call(c, lapply(indcode_list, function(x) x$vnames))
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

###############  EB Functions Environment ###############
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
  # score_n 300 default means compute 300 score stats
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

# 1) Scale as mean(previous prior,  scale * current) 
# 2) Future test Likelihood of beta is average of likelihoods from multiple Priors
# 3) Use predicted values of 1st solution and boost
# 4) Boost during iterations w1 * (prior + (w2 * newprior)

