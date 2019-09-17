
#' BCF predictions after propensity score estimation with internal parallelization
#' 
#' \code{mc.wbart_bcf} returns the (continuous) predicted values for each datapoint with treatment = 0 and treatment = 1 with BART-based prognostic score and BART-based treatment heterogenity. It just calls wbart_bcf.
#' 
#' @param mc.cores Integer, how many processes to spawn.
#' @param nice Integer, what nice level to apply (to prevent it locking up the system)
#' @param seed Integer, seed for special parallel RNG
#' 
#' @inheritParams wbart_bcf
#' @return Identical to wbart_bcf
#' @export

mc.wbart_bcf <- function(
    x.train, y.train, x.test=matrix(0.0,0,0),
    numtreated =1L , extravar=TRUE, use_cauchy=FALSE,
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=FALSE,
    xinfo2=matrix(0.0,0,0),
    cont=FALSE, rm.const=TRUE,
    sigest=NA, sigdf=3, sigquant=0.90,
    k=2.0, power=2.0, base=.95,
    sigmaf=NA, lambda=NA,
    fmean=mean(y.train),
    w=rep(1,length(y.train)),
    ntree=200L, numcut=100L,
    ntree_treated = NULL , score_index = ncol(x.train ) ,
    ndpost=1000L, nskip=100L,
    keepevery=1L, printevery=100L,
    keeptrainfits=TRUE, transposed=FALSE,
    ##treesaslists=FALSE,
    mc.cores = 2L, nice = 19L,
    seed = 99L, ci_frac = 0.99
    )
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()
    n = length(y.train)
    score_index <- score_index*2L/2L ## defeat lazy evaluation

    if( is.null(ntree_treated)) {
      ntree_treated = ceiling(sqrt(numtreated) / sqrt(n) *ntree)
       ntree_treated <- max(ntree_treated, 1)
    }

    if(!transposed) {
        temp2 = bartModelMatrix(x.train[1:numtreated,], numcut, usequants=usequants, cont=cont,  rm.const=FALSE)
        
        xinfo2 = temp2$xinfo
        numcut2 =temp2$numcut
        rm(temp2)

        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               cont=cont, xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        
        if(length(x.test)>0)
            x.test = t(bartModelMatrix(x.test, rm.const=FALSE)[ , temp$rm.const])
        rm.const <- temp$rm.const
        numcut_treated <- numcut_treated[rm.const]
        xinfo2 <- xinfo2[rm.const, ]
        rm(temp)
    }

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected
        ## warning(paste0('The number of cores requested, mc.cores=', mc.cores,
        ##                ',\n exceeds the number of cores detected via detectCores() ',
        ##                'which yields ', mc.cores.detected, ' .'))

    mc.ndpost <- ceiling(ndpost/mc.cores)

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
                   wbart_bcf(x.train=x.train, y.train=y.train, x.test=x.test,
                         numtreated=numtreated, extravar=extravar, use_cauchy=use_cauchy,
                         sparse=sparse, theta=theta, omega=omega,
                         a=a, b=b, augment=augment, rho=rho,
                         xinfo=xinfo,
                         xinfo2=xinfo2,
                         sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                         k=k, power=power, base=base,
                         sigmaf=sigmaf, lambda=lambda, fmean=fmean, w=w,
                         ntree=ntree, numcut=numcut, numcut_treated=numcut2,
                         ntree_treated=ntree_treated , score_index=score_index,
                         ndpost=mc.ndpost, nskip=nskip, keepevery=keepevery,
                         printevery=printevery, transposed=TRUE)},
                         ##treesaslists=treesaslists)},
                   silent=(i!=1))
                   ## to avoid duplication of output
                   ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    ## sigma.len <- length(post$sigma)
    ## if(sigma.len>mc.ndpost) {
    ##     sigma.beg <- 1+sigma.len-mc.ndpost
    ##     post$sigma <- post$sigma[sigma.beg:sigma.len]
    ## }

    if(mc.cores==1 ) return(post)
    else {
        if(class(rm.const)!='logical') post$rm.const <- rm.const
    
        p <- nrow(x.train[post$rm.const, ])
        ##p <- nrow(x.train[ , post$rm.const])

        ## if(length(rm.const)==0) rm.const <- 1:p
        ## post$rm.const <- rm.const

        old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree), ' ', as.character(p))
        old.text2 <- paste0(as.character(mc.ndpost), ' ', as.character(ntree_treated), ' ', as.character(p))
        old.stop <- nchar(old.text)
        old.stop2 <- nchar(old.text2)

        post$treedraws$trees <- sub(old.text,
                                    paste0(as.character(mc.cores*mc.ndpost), ' ', as.character(ntree), ' ',
                                           as.character(p)),
                                    post$treedraws$trees)
        post$treedraws_treated$trees <- sub(old.text2,
                                    paste0(as.character(mc.cores*mc.ndpost), ' ', as.character(ntree_treated), ' ',
                                           as.character(p)),
                                    post$treedraws_treated$trees)

        for(i in 2:mc.cores) {
            post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)
            post$yhat.train.treated <- rbind(post$yhat.train.treated, post.list[[i]]$yhat.train.treated)

            if(length(post$yhat.test)>0) {
                post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)
                post$yhat.test.treated <- rbind(post$yhat.test.treated, post.list[[i]]$yhat.test.treated)                
            }

            ## if(sigma.len>0)
            ##     post$sigma <- c(post$sigma, post.list[[i]]$sigma[sigma.beg:sigma.len])

            post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)

            post$treedraws$trees <- paste0(post$treedraws$trees,
                                           substr(post.list[[i]]$treedraws$trees, old.stop+2,
                                                  nchar(post.list[[i]]$treedraws$trees)))
            post$treedraws_treated$trees <- paste0(post$treedraws_treated$trees,
                                           substr(post.list[[i]]$treedraws_treated$trees, old.stop2+2,
                                                  nchar(post.list[[i]]$treedraws_treated$trees)))

            ## if(treesaslists) post$treedraws$lists <-
            ##                      c(post$treedraws$lists, post.list[[i]]$treedraws$lists)

            if(length(post$varcount)>0) {
                post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
                post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
                post$varcount.treated <- rbind(post$varcount.treated, post.list[[i]]$varcount.treated)
                post$varprob.treated <- rbind(post$varprob.treated, post.list[[i]]$varprob.treated)
                
            }
        }

        if(length(post$yhat.train.mean)>0) {
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
            post$yhat.train.mean.treated <- apply(post$yhat.train.treated, 2, mean)
        }
            

        if(length(post$yhat.test.mean)>0)
            post$yhat.test.mean.treated <- apply(post$yhat.test.treated, 2, mean)
            post$yhat.test.mean.treated <- apply(post$yhat.test.treated, 2, mean)

        if(length(post$varcount)>0) {
            post$varcount.mean.treated <- apply(post$varcount.treated, 2, mean)
            post$varprob.mean.treated <- apply(post$varprob.treated, 2, mean)
            post$varcount.mean.treated <- apply(post$varcount.treated, 2, mean)
            post$varprob.mean.treated <- apply(post$varprob.treated, 2, mean)
        }
        
        ate_vector <- rowMeans(post$yhat.train.treated - post$yhat.train)
        post$ate_est <- mean(ate_vector)
        post$ate_sd <- sd(ate_vector)
        post$ate_ci <- HDInterval::hdi(ate_vector, credMass = ci_frac)


        att_vector <- rowMeans(post$yhat.train.treated[,seq.int(numtreated)] - post$yhat.train[,seq.int(numtreated)])
        post$att_est <- mean(att_vector)
        post$att_sd <- sd(att_vector)
        post$att_ci <- HDInterval::hdi(att_vector, credMass = ci_frac)

        att_modified_vector <- mean(y.train[seq.int(numtreated)]) - rowMeans( post$yhat.train[,seq.int(numtreated)])
        post$att_modified_est <- mean(att_modified_vector)
        post$att_modified_sd <- sd(att_modified_vector)
        post$att_modified_ci <- HDInterval::hdi(att_modified_vector, credMass = ci_frac)

        
        late_set <- which( (x.train[ score_index, ] < 0.6) && (x.train[ score_index, ] > 0.4))
        if(length(late_set) > 20) {
        late_vector <- rowMeans(post$yhat.train.treated[ , late_set] - post$yhat.train[ , late_set])
        post$late_est <- mean(late_vector)
        post$late_sd <- sd(late_vector)
        post$late_ci <- HDInterval::hdi(late_vector, credMass = ci_frac)
        } else {
        post$late_est <- post$late_sd <- post$late_ci <- NA
        }
        
        


#         attr(post, 'class') <- 'wbart'

        return(post)
    }
}


