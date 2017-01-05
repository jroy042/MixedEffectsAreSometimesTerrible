
compare_lme4_stan_imbalanced_sims <- function(B = 10000, threads = 16, dir = getwd(),
gau_file = "imbalanced_gaussian.stan", bin_file = "imbalanced_logistic.stan"){
  # B is the number of models to run for each of the 2 model types.
  # threads is the number of threads to split this task by.
  # set threads to 1 to not parallelize the iterations,
  # in which case the individual Stan models will be parallelized

  single_thread <- function(pardat){
    require(lme4)
    require(rstan)
    require(stringr)
    require(gtools)
    require(clusterGeneration)
    options(contrasts=c("contr.sum","contr.poly"))

    ### Parameters of interest and lme4 formula
    keep <- c("coef","s0","s1","s2","s3",
      "r01","r02","r03","r12","r13","r23")
    formula <- y ~ x1 + x2 + (1 + x1 + x2 | subject)

    ### Generate dataset and model
    generate_model <- function(isgau){
      balance <- function(dat){
        nobs <- nrow(dat)
        counts <- as.data.frame(xtabs(~.,dat))
        cells <- nrow(counts)
        non_zero_counts <- subset(counts, Freq > 0)
        non_empty_cells <- nrow(non_zero_counts)
        pct_non_empty_cells <- non_empty_cells / cells
        expected <- nobs / non_empty_cells
        freq <- non_zero_counts$Freq
        ratio <- freq / expected
        ratio[ratio < 1] <- 1 / ratio[ratio < 1]
        bal <- mean(ratio) / pct_non_empty_cells
        bal <- 2 * (1 - bal / (1 + bal))
        return(bal)
      }

      mod <- list(fam = ifelse(isgau,"Gaussian","Logistic"))

      # S is the total number of subjects
      # L is the mean observations per subject
      # NS is a vector of the number of observations for each subject
      mod$S <- S <- sample(30:60,1)
      mod$L <- L <- runif(1,20,30)
      NS <- rpois(S,L)
      
      # Ensure all subjects have observations
      # This loop should rarely execute
      while(any(NS==0)){
        NS <- rpois(S,L)
      }

      # Generate data
      subject <- NULL
      for(s in 1:S) subject <- c(subject,rep(s,NS[s]))
      x1 = sample(c("a","b"), sum(NS),
        prob=rdirichlet(1,c(1,1)), replace=T)
      x2 = sample(c("a","b","c"), sum(NS),
        prob=rdirichlet(1,c(1,1,1)), replace=T)
      
      # Ensure all levels have observations, which should
      # almost never be required
      x1[sample(1:sum(NS),2,replace=F)] <- c("a","b")
      x2[sample(1:sum(NS),3,replace=F)] <- c("a","b","c")
      
      # Create data frame and compute balance
      mod$frame <- data.frame(x1,x2,subject)
      mod$N <- N <- nrow(mod$frame)
      mod$balance <- balance(mod$frame)
      mod$frame$y <- 0
      
      # Create model matrix
      x <- model.matrix(nobars(formula),mod$frame)
      z <- mkReTrms(findbars(formula),mod$frame)
      mod$mat <- cbind(x,t(z$Ztlist[[1]]))
      sparse <- extract_sparse_parts(mod$mat)
      mod$stan <- list(N = N, S = S, P = 4, QS = 4,
        nz = length(sparse$w), x_w = sparse$w, x_v = sparse$v,
        x_u = sparse$u)
      
      # Generate model
      mod$b <- b <- c(runif(1,-2,2),runif(3,-1,1))
      mod$ss <- ss <- c(runif(1,0,1),runif(3,0,0.5))
      mod$os <- os <- rcorrmatrix(4)
      cholfail <- tryCatch(Ls <- t(chol(os)), error = function(e) e)
      while(inherits(cholfail, "error")){
        mod$os <- os <- rcorrmatrix(4)
        cholfail <- tryCatch(Ls <- t(chol(os)), error = function(e) e)
      }
      gs <- Ls %*% matrix(rnorm(S*4),4,S)
      for(i in 1:4) gs[i,] <- gs[i,] * ss[i]

      y <- as.vector(mod$mat %*% c(b,as.vector(gs)))

      if(isgau){
        mod$se <- se <- runif(1,0,1)
        y <- y + rnorm(N,0,se)
        mod$pext <- NA
        mod$frame$y <- mod$stan$y <- scale(y)[,1]
        mod$b[1] <- mod$b[1] - mean(y)
        mod$b <- mod$b / sd(y)
        mod$ss <- mod$ss / sd(y)
        mod$se <- mod$se / sd(y)
      } else {
        mod$se <- NA
        mod$pext <- mean(abs(y)>5)
        mod$frame$y <- mod$stan$y <- rbinom(N,1,plogis(y))
      }

      return(mod)
    }

    ### function to perform rank-deficiency test from Bates et al. 2015
    my.repca <- function(theta, cutoff = 1){
      m <- matrix(0,4,4)
      m[lower.tri(m,diag=T)] <- theta
      vv <- svd(m,nv=0L)
      names(vv) <- c("sdev","rotation")
      vv$center <- FALSE
      vv$scale <- FALSE
      class(vv) <- "prcomp"
      return(1 * (which(summary(vv)$importance[3,] >= cutoff)[1] < 4))
    }

    ### Fit single model with both lme4 and stan
    single_model <- function(linear){
      # generate_model shouldn't throw errors, but this is included in case
      genmod_fails <- 0
      possfail_genmod <- tryCatch(mod <- generate_model(linear), error = function(e) e)
      while(inherits(possfail_genmod, "error")){
        genmod_fails <- genmod_fails + 1
        possfail_genmod <- tryCatch(mod <- generate_model(linear), error = function(e) e)
      }

      if(linear){
        possfail_lme4 <- tryCatch(m_lme4 <- summary(lmer(formula, mod$frame)), error = function(e) e)
        possfail_stan <- tryCatch(m_stan <- sampling(object = gaustan, data = mod$stan, chains = 3,
          pars = c(keep,"res"), control = list(adapt_delta = 0.99)), error = function(e) e)
      } else {
        possfail_lme4 <- tryCatch(m_lme4 <- summary(glmer(formula, mod$frame, family = binomial)), error = function(e) e)
        possfail_stan <- tryCatch(m_stan <- sampling(object = binstan, data = mod$stan, chains = 3,
          pars = keep, control = list(adapt_delta = 0.99)), error = function(e) e)
      }
      
      if(inherits(possfail_lme4, "error") | inherits(possfail_stan, "error")){
        save(mod, possfail_lme4, possfail_stan, file = paste(getwd(),"/Imbalanced_",
          threadchar,"_",bchar,"_",mod$fam,"_Error.rda",sep=""))
        return(NULL)
        
      } else {
        ## assess lme4 convergence
        mess_lme4 <- m_lme4$optinfo$conv$lme4
        if(!is.null(names(mess_lme4))){
          mess_lme4 <- paste(mess_lme4$messages,collapse="; ")
          mess_lme4 <- str_replace_all(mess_lme4,",","")
          mess_lme4 <- str_replace_all(mess_lme4,"\n","")
          conv_lme4 <- 1 * (nrow(str_locate_all(mess_lme4,"negative")[[1]])==0
            & nrow(str_locate_all(mess_lme4,"ratio")[[1]])==0
            & nrow(str_locate_all(mess_lme4,"unable")[[1]])==0)
          if(nrow(str_locate_all(mess_lme4,"converge with max")[[1]])>0){
            j <- str_split(mess_lme4,"converge with max")[[1]][2]
            j <- substr(j,10,nchar(j))
            grad_lme4 <- as.numeric(str_split(j," ")[[1]][1])
            conv_lme4 <- conv_lme4 * (grad_lme4 < .01)
          } else grad_lme4 <- .002
        } else {
          mess_lme4 <- ""
          conv_lme4 <- 1
          grad_lme4 <- .002
        }
        
        ## get ranef cov and chol from lme4
        vc <- data.frame(m_lme4$varcor)
        theta <- m_lme4$optinfo$val

        ## assess stan convergence and get parameter estimates
        samp <- data.frame(do.call(rbind,
          args=get_sampler_params(m_stan,inc_warmup=FALSE)))
        ndiv <- sum(samp$divergent__)
        nmtd <- sum(samp$treedepth__ > 10)
        m_stan <- summary(m_stan,probs=c(.025,.975))$summary
        rhat <- quantile(m_stan[,"Rhat"],c(0,.25,.5,.75,1))
        neff <- quantile(m_stan[,"n_eff"],c(0,.25,.5,.75,1))
        m_stan <- m_stan[,1]
        conv_stan <- 1 * (ndiv == 0 & as.numeric(rhat[5]) < 1.1)
        
        if(conv_stan==0 | conv_lme4==0){
          save(mod, file = paste(getwd(),"/Imbalanced_",
            threadchar,"_",bchar,"_",mod$fam,"_Unconverged.rda",sep=""))
        }

        ## assemble results
        resblme4 <- data.frame(
          thread = thread, iteration = b, fam = mod$fam, reg = "lme4", minvar = min(mod$ss),
          pext = mod$pext, N = mod$N, S = mod$S, L = mod$L, balance = mod$balance,
          b0 = mod$b[1], b0_pred = m_lme4$coef[1,1],
          b1 = mod$b[2], b1_pred = m_lme4$coef[2,1],
          b2 = mod$b[3], b2_pred = m_lme4$coef[3,1],
          b3 = mod$b[4], b3_pred = m_lme4$coef[4,1],
          s0 = mod$ss[1], s0_pred = vc[1,5],
          s1 = mod$ss[2], s1_pred = vc[2,5],
          s2 = mod$ss[3], s2_pred = vc[3,5],
          s3 = mod$ss[4], s3_pred = vc[4,5],
          r01 = mod$os[2,1], r01_pred = vc[5,5],
          r02 = mod$os[3,1], r02_pred = vc[6,5],
          r03 = mod$os[4,1], r03_pred = vc[7,5],
          r12 = mod$os[3,2], r12_pred = vc[8,5],
          r13 = mod$os[4,2], r13_pred = vc[9,5],
          r23 = mod$os[4,3], r23_pred = vc[10,5],
          sres = mod$se, sres_pred = ifelse(linear,vc[11,5],NA),
          conv = conv_lme4, grad = grad_lme4, mess = mess_lme4, ndiv = NA, nmtd = NA,
          rhat_min = NA, rhat_q25 = NA, rhat_med = NA, rhat_q75 = NA, rhat_max = NA,
          neff_min = NA, neff_q25 = NA, neff_med = NA, neff_q75 = NA, neff_max = NA,
          genmod_fails, stringsAsFactors = FALSE, row.names = NULL)
        resblme4$repca <- my.repca(theta)
        resblme4[,paste("theta",str_pad(1:10,width=2,pad="0",side="left"),sep="_")] <- theta

        resbstan <- data.frame(
          thread = thread, iteration = b, fam = mod$fam, reg = "Stan", minvar = min(mod$ss),
          pext = mod$pext, N = mod$N, S = mod$S, L = mod$L, balance = mod$balance,
          b0 = mod$b[1], b0_pred = m_stan["coef[1]"],
          b1 = mod$b[2], b1_pred = m_stan["coef[2]"],
          b2 = mod$b[3], b2_pred = m_stan["coef[3]"],
          b3 = mod$b[4], b3_pred = m_stan["coef[4]"],
          s0 = mod$ss[1], s0_pred = m_stan["s0"],
          s1 = mod$ss[2], s1_pred = m_stan["s1"],
          s2 = mod$ss[3], s2_pred = m_stan["s2"],
          s3 = mod$ss[4], s3_pred = m_stan["s3"],
          r01 = mod$os[2,1], r01_pred = m_stan["r01"],
          r02 = mod$os[3,1], r02_pred = m_stan["r02"],
          r03 = mod$os[4,1], r03_pred = m_stan["r03"],
          r12 = mod$os[3,2], r12_pred = m_stan["r12"],
          r13 = mod$os[4,2], r13_pred = m_stan["r13"],
          r23 = mod$os[4,3], r23_pred = m_stan["r23"],
          sres = mod$se, sres_pred = ifelse(linear,m_stan["res"],NA),
          conv = conv_stan, grad = NA, mess = NA, ndiv, nmtd,
          rhat_min = rhat[1], rhat_q25 = rhat[2], rhat_med = rhat[3], rhat_q75 = rhat[4], rhat_max = rhat[5],
          neff_min = neff[1], neff_q25 = neff[2], neff_med = neff[3], neff_q75 = neff[4], neff_max = neff[5],
          genmod_fails, stringsAsFactors = FALSE, row.names = NULL)
        resbstan$repca <- NA
        resbstan[,paste("theta",str_pad(1:10,width=2,pad="0",side="left"),sep="_")] <- NA

        return(rbind(resblme4,resbstan))
      }
    }

    ## unpack arguments
    B <- pardat[[1]]
    thread <- pardat[[2]]
    gaucode <- pardat[[3]]
    bincode <- pardat[[4]]
    threadchar <- str_pad(paste(thread),width=nchar(paste(pardat[[5]])),side="left",pad="0")

    gaustan <- stan_model(model_name = "linear_sim", model_code = gaucode, save_dso = TRUE)
    binstan <- stan_model(model_name = "logistic_sim", model_code = bincode, save_dso = TRUE)

    ### run B iterations of each of the 2 models and output progress to thread-specific file
    results <- NULL
    progfile <- paste(getwd(),"/","Imbalanced_",threadchar,".progress",sep="")
    csvfile <- paste(getwd(),"/","Imbalanced_",threadchar,".csv",sep="")

    cat(paste(Sys.time(),"Started Thread",threadchar),file=progfile,sep="\n",append=FALSE)
    for(b in 1:B){  # Each iteration has 2 calls to single_model, one for each model family
      bchar <- str_pad(paste(b),width=nchar(paste(B)),side="left",pad="0")
      cat("\n",paste(Sys.time(),"Iteration",bchar,"/",B),file=progfile,sep="",append=TRUE)
      
      temp <- single_model(TRUE)
      if(is.null(temp)){
        results[nrow(results)+1,1:2] <- c(thread,b)
        results[nrow(results),3] <- "Gaussian"
        results[nrow(results),4] <- "lme4"
        results[nrow(results)+1,1:2] <- c(thread,b)
        results[nrow(results),3] <- "Gaussian"
        results[nrow(results),4] <- "Stan"
        cat(" ERROR ",file=progfile,sep="",append=TRUE)
      } else {
        results <- rbind(results, temp)
        cat(" . ",file=progfile,sep="",append=TRUE)
      }
      
      temp <- single_model(FALSE)
      if(is.null(temp)){
        results[nrow(results)+1,1:2] <- c(thread,b)
        results[nrow(results),3] <- "Logistic"
        results[nrow(results),4] <- "lme4"
        results[nrow(results)+1,1:2] <- c(thread,b)
        results[nrow(results),3] <- "Logistic"
        results[nrow(results),4] <- "Stan"
        cat(" ERROR ",file=progfile,sep="",append=TRUE)
      } else {
        results <- rbind(results, temp)
        cat(" . ",file=progfile,sep="",append=TRUE)
      }
      
      write.csv(results,row.names=F,quote=F,file=csvfile)
    }
    cat("\n\n",paste(Sys.time(),"Finished Thread",threadchar),file=progfile,sep="",append=TRUE)
    
    ### return results
    return(results)
  }
  
  setwd(dir)
  gaucode <- paste(readLines(gau_file),collapse="\n")
  bincode <- paste(readLines(bin_file),collapse="\n")

  if(threads==1){
    options(mc.cores = parallel::detectCores())
    results <- single_thread(list(B,threads,gaucode,bincode,threads))
  } else {
    require(parallel)
    cl = makeCluster(rep("localhost", threads))
    
    # simulations per thread and thread number
    left <- B %% threads
    B <- rep(floor(B/threads),threads)
    if(left>0) B[1:left] <- B[1:left] + 1
    pardat <- list()
    for(i in 1:threads) pardat[[i]] <- list(B[i],i,gaucode,bincode,threads)
    
    results <- parLapply(cl,pardat,single_thread)
    
    stopCluster(cl)
    
    results <- do.call(rbind,args=results)
  }
  
  return(results)
}
