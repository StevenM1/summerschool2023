# data generation 

# Used in make_data and make_samplers
add_trials <- function(dat) 
  # Add trials column, 1:n for each subject   
{
  n <- table(dat$subjects)
  if (!any(names(dat)=="trials")) dat <- cbind.data.frame(dat,trials=NA)
  for (i in names(n)) dat$trials[dat$subjects==i] <- 1:n[i]
  dat
}


# model=NULL;trials=NULL;data=NULL;expand=1;
# mapped_p=FALSE;LT=NULL;UT=NULL;LC=NULL;UC=NULL
# trials=10; design=design_B_MT
# 
# p_vector=pars[[i]];design=design[[j]];model=model[[j]];data=data[[j]]

make_data <- function(p_vector,design,model=NULL,trials=NULL,data=NULL,expand=1,
                      mapped_p=FALSE,LT=NULL,UT=NULL,LC=NULL,UC=NULL,Fcovariates=NULL)
  # Simulates data using rfun from model specified by a formula list (Flist) 
  # a factor contrast list (Clist, if null data frame creation defaults used) 
  # using model (list specifying p_types, transform, Mtransform and rfun).
  # If data is supplied that determines the design, with data sorted by subjects
  #   and a trials = 1:n trials/subject factor added (or overwriting any existing)
  #   NB: trials differs from generated where numbering is cell specific. 
  # Otherwise list cfactors and vector rfactor used to create a design fully 
  # crossed in cfactors with response (R) factor specified by rfactor, 
  #  e.g for simple  design cfactors = list(subjects=1:2,trials=1:2,S=1:2) and 
  #  rfactor= c(R=1:2).
  #  NB Factor order level is always the same as cfactors and rfactor
  # expand replicates the design expand times
  # matchfun scores correct response accumulator lM, e.g., for scoring latent 
  #   response = stimulus  matchfun = function(d)d$S==d$lR
# If mapped_p is true instead returns cbind(da,pars), augmented data and 
# mapped pars.

{
  
  missingFilter <- function(data,LT,UT,LC,UC,Rmissing) 
  {
    
    makeExclude <- function(data,L=NULL,U=NULL) 
    {
      
      makeCrit <- function(data,bound) {
        exclude <- factor(data$subjects)
        levels(exclude) <- bound[levels(exclude)]
        as.numeric(as.character(exclude))
      }
      
      if (!is.null(L)) {
        if (length(L)==1) exclude <- data$rt < L else {
          if (!all(sort(levels(data$subjects))==names(L)))
            stop("Missing subjects")
          exclude <- data$rt < makeCrit(exclude)
        }
      } else {
        if (length(U)==1) exclude <- data$rt > U else {
          if (!all(sort(levels(data$subjects))==names(U)))
            stop("Missing subjects")
          exclude <- data$rt > makeCrit(exclude)
        }
      }
      exclude
    }
     
    if (!is.null(LT)) data <- data[!makeExclude(data,L=LT),] 
    if (!is.null(UT)) data <- data[!makeExclude(data,U=UT),]
    if (!is.null(LT)) {
      exclude <- makeExclude(data,L=LC)
      data$rt[exclude] <- -Inf
      if (Rmissing) data$R[exclude] <- NA 
    }
    if (!is.null(UT)) {
      exclude <- makeExclude(data,U=UC)
      data$rt[exclude] <- -Inf
      if (Rmissing) data$R[exclude] <- NA 
    }
    data
  }
  
  if (is.null(model)) if (is.null(design$model)) 
    stop("Must specify model as not in design") else model <- design$model

  if ( is.null(data) ) {
    if (mapped_p) trials <- 1
    if ( is.null(trials) )
      stop("If data is not provided need to specify number of trials")
    
    Ffactors=c(design$Ffactors,list(trials=1:trials))
    data <- as.data.frame.table(array(dim=unlist(lapply(Ffactors,length)),
                                    dimnames=Ffactors))
    for (i in names(design$Ffactors)) 
      data[[i]] <- factor(data[[i]],levels=design$Ffactors[[i]])
    names(data)[dim(data)[2]] <- "R"
    data$R <- factor(data$R,levels=design$Rlevels)
    data$trials <- as.numeric(as.character(data$trials))
    if (!is.null(design$Ffunctions)) data <- 
      cbind.data.frame(data,data.frame(lapply(design$Ffunctions,function(f){f(data)})))
    LT <- UT <- LC <- UC <- Rmissing <- NULL
    if (!is.null(Fcovariates)) {
      if (is.null(design$Fcovariates)) stop("Design does not contain Fcovariates")
      if (!(all(sort(names(Fcovariates))==sort(design$Fcovariates))))
        stop("Fcovariates and design$Fcovariates have different names")
      if (!is.data.frame(Fcovariates)) {
        if (!all(unlist(lapply(Fcovariates,is.function))))
          stop("Fcovariates must be either a data frame or list of functions")
        nams <- names(Fcovariates)
        Fcovariates <- do.call(cbind.data.frame,lapply(Fcovariates,function(x){x(data)}))
        names(Fcovariates) <- nams 
      }
      n <- dim(Fcovariates)[1]
      if (!(n==dim(data)[1])) stop("Fcovariates must specify ",dim(data)[1]," values per covariate")
      data <- cbind.data.frame(data,Fcovariates)
    }
  } else {
    LT <- attr(data,"LT"); UT <- attr(data,"UT")
    LC <- attr(data,"LC"); UC <- attr(data,"UC")
    Rmissing <- any(is.na(data$R))
    data <- add_trials(data[order(data$subjects),])
  }
  if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
  if ( is.null(model$p_types) ) stop("model$p_types must be specified")
  if ( is.null(model$transform) ) model$transform <- identity
  if ( is.null(model$Mtransform) ) mdoel$Mtransform <- identity
  data <- design_model(
    add_accumulators(data,design$matchfun,simulate=TRUE,type=model$type),
    design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
    rt_check=FALSE)
  pars <- model$Mtransform(map_p(
    model$transform(add_constants(p_vector,design$constants)),data
  ))
  if (mapped_p) 
    return(cbind(data[,!(names(data) %in% c("R","rt"))],pars))
  if (expand==1) 
    Rrt <- model$rfun(data$lR,pars) else
      Rrt <- model$rfun(rep(data$lR,expand),apply(pars,2,rep,times=expand))
  data <- data[data$lR==levels(data$lR)[1],!(names(data) %in% c("lR","lM"))]
  if (expand>1) data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                              data.frame(lapply(data,rep,times=expand))) 
  data[,c("R","rt")] <- Rrt
  missingFilter(data[,names(data)!="winner"],LT,UT,LC,UC,Rmissing)
}

# hyper=FALSE;n_post=100;expand=1; filter="sample";subfilter=1;thin=1;
# use_par=c("random","mean","median")[1]; n_cores=1
# samples=rdm_B_MT; n_post=1;subfilter=1000
post_predict <- function(samples,hyper=FALSE,n_post=100,expand=1,
                         filter="sample",subfilter=0,thin=1,n_cores=1,
                         use_par=c("random","mean","median")[1]) 
  # Post predictions for samples object, based on random samples or some
  # central tendency statistic. 
  # n_post is number of parameter vectors used
  # expand=1 gives exact data design, larger values replicate whole design
  # filter/subfilter/thin as for as_mcmc.list
  # hyper=FALSE draws from alphas (participant level)
  # hyper=TRUE draws from hyper 
{

  data <- attr(samples,"data_list")
  design <- attr(samples,"design_list")
  model <- attr(samples,"model_list")
  if(length(data) > 1){
    jointModel <- TRUE
    all_samples <- samples
  } else{
    jointModel <- FALSE
  }
  post_out <- vector("list", length = length(data))
  for(j in 1:length(data)){
    if(jointModel) samples <- single_out_joint(all_samples, j)
    subjects <- levels(data[[j]]$subjects)
    
    if (hyper) {
      pars <- vector(mode="list",length=n_post) 
      for (i in 1:n_post) {
        pars[[i]] <- get_prior_samples(samples,selection="alpha",
                                       filter=filter,thin=thin,subfilter=subfilter,n_prior=length(subjects))
        row.names(pars[[i]]) <- subjects   
      }
    } else {
      samps <- lapply(as_mcmc.list(samples,selection="alpha",
                                   filter=filter,subfilter=subfilter,thin=thin),function(x){do.call(rbind,x)})
      if (use_par != "random") {
        p <- do.call(rbind,lapply(samps,function(x){apply(x,2,use_par)}))
        row.names(p) <- subjects
      }
      pars <- vector(mode="list",length=n_post) 
      for (i in 1:n_post) {
        if (use_par != "random") pars[[i]] <- p else {
          pars[[i]] <- do.call(rbind,lapply(samps,function(x){x[sample(1:dim(x)[1],1),]})) 
          row.names(pars[[i]]) <- subjects
        }
      }
    }
    if (n_cores==1) {
      simDat <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        cat(".")
        simDat[[i]] <- make_data(pars[[i]],design=design[[j]],model=model[[j]],data=data[[j]],expand=expand)
      }
      cat("\n")
    } else {
      simDat <- mclapply(1:n_post,function(i){
        make_data(pars[[i]],design=design[[j]],model=model[[j]],data=data[[j]],expand=expand)
      },mc.cores=n_cores)
    }
    out <- cbind(postn=rep(1:n_post,each=dim(simDat[[1]])[1]),do.call(rbind,simDat))
    if (n_post==1) pars <- pars[[1]]
    attr(out,"pars") <- pars 
    post_out[[j]] <- out
  }
  if(!jointModel) post_out <- post_out[[1]]
  return(post_out)
}


