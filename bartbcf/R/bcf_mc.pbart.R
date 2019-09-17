
#' BCF predictions after propensity score estimation with internal parallelization
#' 
#' \code{mc.pbart_bcf} returns the (latent) predicted values for each datapoint with treatment = 0 and treatment = 1 with BART-based prognostic score and BART-based treatment heterogenity for binary outcomes. It just calls pbart_bcf.
#' 
#' @param mc.cores Integer, how many processes to spawn.
#' @param nice Integer, what nice level to apply (to prevent it locking up the system)
#' @param seed Integer, seed for special parallel RNG
#' 
#' @inheritParams pbart_bcf
#' @return Identical to pbart_bcf
#' @export
#' @importFrom parallel detectCores
#' @importFrom tools psnice

mc.pbart_bcf <- function(
    x.train, y.train, x.test = matrix(0.0, 0L, 0L),
    numtreated =1L , extravar=TRUE, use_cauchy=FALSE,
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=FALSE,
    xinfo2=matrix(0.0,0,0),
    cont=FALSE, rm.const=TRUE,
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = 0,
    ntree=50L, numcut=100L,
    ntree_treated = NULL , score_index = ncol(x.train ) ,
    ndpost=1000L, nskip=100L,
    keepevery=1L, printevery=100L,
    keeptrainfits=TRUE, transposed=FALSE,
##    treesaslists=FALSE,
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
    score_index <- score_index*1L ## defeat lazy evaluation
    if( is.null(ntree_treated)) {
      temp <- c(sum(y.train[seq.int(numtreated)]) , sum(y.train) )
      ntree_treated = sqrt( min(temp[1], numtreated-temp[1]) ) / sqrt( min(temp[2], n-temp[2]) )
      ntree_treated <- max(ceiling(ntree_treated*ntree) ,1)
    }

    if(!transposed) {
    temp2 = bartModelMatrix(x.train[1:numtreated,], numcut_treated, usequants=usequants, cont=cont,  rm.const=FALSE)
        
        xinfo2 = temp2$xinfo
        numcut_treated =temp2$numcut
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
    ## mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

#     while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+1

    ## mc.nkeep <- mc.ndpost %/% keepevery

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
                  pbart_bcf(x.train=x.train, y.train=y.train, x.test=x.test,
                         numtreated=numtreated, extravar=extravar, use_cauchy=use_cauchy,
                        sparse=sparse, theta=theta, omega=omega,
                        a=a, b=b, augment=augment, rho=rho,
                        xinfo=xinfo, usequants=usequants,
                         xinfo2=xinfo2,
                        cont=FALSE, rm.const=rm.const,
                        k=k, power=power, base=base,
                        binaryOffset=binaryOffset,
                        ntree=ntree, numcut=numcut, numcut_treated=numcut_treated,
                         ntree_treated=ntree_treated ,  score_index=score_index,
                        ndpost=mc.ndpost, nskip=nskip, keepevery=keepevery,
                        ## nkeeptrain=mc.nkeep, nkeeptest=mc.nkeep,
                        ## nkeeptestmean=mc.nkeep, nkeeptreedraws=mc.nkeep,
                        printevery=printevery, transposed=TRUE)},
                        ##treesaslists=treesaslists)},
                  silent=(i!=1))
                  ## to avoid duplication of output
                  ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    if(mc.cores==1 ) return(post)
    else {
        if(class(rm.const)!='logical') post$rm.const <- rm.const

        p <- nrow(x.train[post$rm.const, ])
        ##p <- nrow(x.train[ , post$rm.const])

        ## if(length(rm.const)==0) rm.const <- 1:p
        ## post$rm.const <- rm.const

        old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree), ' ', as.character(p))
        old.text2 <- paste0(as.character(mc.ndpost), ' ', as.character(ntree_treated), ' ', as.character(p))
        
        ##old.text <- paste0(as.character(mc.nkeep), ' ', as.character(ntree), ' ', as.character(p))
        old.stop <- nchar(old.text)
        old.stop2 <- nchar(old.text2)


        post$treedraws$trees <- sub(old.text,
                                    paste0(as.character(mc.cores*mc.ndpost), ' ', as.character(ntree), ' ',
                                    ##paste0(as.character(mc.cores*mc.nkeep), ' ', as.character(ntree), ' ',
                                           as.character(p)),
                                    post$treedraws$trees)
        post$treedraws_treated$trees <- sub(old.text2,
                                    paste0(as.character(mc.cores*mc.ndpost), ' ', as.character(ntree_treated), ' ',
                                           as.character(p)),
                                    post$treedraws_treated$trees)

        keeptestfits <- length(x.test)>0

        for(i in 2:mc.cores) {
            if(keeptrainfits) {
                post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)
                post$prob.train <- rbind(post$prob.train, post.list[[i]]$prob.train)
                post$yhat.train.treated <- rbind(post$yhat.train.treated, post.list[[i]]$yhat.train.treated)
                post$prob.train.treated <- rbind(post$prob.train.treated, post.list[[i]]$prob.train.treated)
                
            }

            if(keeptestfits) {
                post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)
                post$prob.test <- rbind(post$prob.test, post.list[[i]]$prob.test)
            }

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
            post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
            post$varcount.treated <- rbind(post$varcount.treated, post.list[[i]]$varcount.treated)
                post$varprob.treated <- rbind(post$varprob.treated, post.list[[i]]$varprob.treated)

            post$treedraws$trees <- paste0(post$treedraws$trees,
                                           substr(post.list[[i]]$treedraws$trees, old.stop+2,
                                                  nchar(post.list[[i]]$treedraws$trees)))
            post$treedraws_treated$trees <- paste0(post$treedraws_treated$trees,
                                           substr(post.list[[i]]$treedraws_treated$trees, old.stop2+2,
                                                  nchar(post.list[[i]]$treedraws_treated$trees)))

            ## if(treesaslists) post$treedraws$lists <-
            ##                      c(post$treedraws$lists, post.list[[i]]$treedraws$lists)
        }

        ## if(length(post$yhat.train.mean)>0)
        ##     post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

        ## if(length(post$yhat.test.mean)>0)
        ##     post$yhat.test.mean <- apply(post$yhat.test, 2, mean)





        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)
        post$varcount.mean.treated <- apply(post$varcount.treated, 2, mean)
        post$varprob.mean.treated <- apply(post$varprob.treated, 2, mean)
        
        ate_vector <- rowMeans(pnorm(post$yhat.train.treated) - pnorm(post$yhat.train))
        att_vector <- rowMeans(pnorm(post$yhat.train.treated[,seq.int(numtreated)]) - pnorm(post$yhat.train[,seq.int(numtreated)]) )
        att_modified_vector <- mean(y.train[seq.int(numtreated)]) - rowMeans( pnorm(post$yhat.train[,seq.int(numtreated)]) )

        post$ate_est <- mean(ate_vector)
        post$ate_sd <- sd(ate_vector)
        post$ate_ci <- HDInterval::hdi(ate_vector, credMass = ci_frac)


        post$att_est <- mean(att_vector)
        post$att_sd <- sd(att_vector)
        post$att_ci <- HDInterval::hdi(att_vector, credMass = ci_frac)

        post$att_modified_est <- mean(att_modified_vector)
        post$att_modified_sd <- sd(att_modified_vector)
        post$att_modified_ci <- HDInterval::hdi(att_modified_vector, credMass = ci_frac)

        
        late_set <- which( (x.train[ score_index, ] < 0.6) && (x.train[ score_index, ] > 0.4))
        if(length(late_set) > 20) {
        late_vector <- rowMeans(pnorm(post$yhat.train.treated[ , late_set]) - pnorm(post$yhat.train[ , late_set]))

        post$late_est <- mean(late_vector)
        post$late_sd <- sd(late_vector)
        post$late_ci <- HDInterval::hdi(late_vector, credMass = ci_frac)
        } else {
        post$late_est <- post$late_sd <- post$late_ci <- NA
        }

#         attr(post, 'class') <- 'pbart'

        return(post)
    }
}
