setwd("/work2/gpeters/ozonevge/Programming/RPackages/rnome/R")
library(Rcpp)

#Sys.setenv("PKG_CXXFLAGS"="")
#sourceCpp("../src/rcpp_hello_world.cpp")

binding_models <- list(list("PROTECT_PROB" = rep(0.99,147),
                            "PRIOR" = 1e-3,
                            "NAME" = "Nucleosome"),
                       list("PROTECT_PROB" = rep(0.99,30),
                            "PRIOR" = 1e-4,
                            "NAME" = "TF"))

## Load test NOMEseq data

laura.data <- readRDS("/work2/gpeters/gasplaur/NOMe-seq/EXPERIMENTS/EXP2_Spermatids_1506R/Rdata/04_sm-unique-allamplicons_allsamples.rds")


cell.stages <- c("1506F1" = "Rounds1",
                 "1506F2" = "Rounds2",
                 "1506F3" = "Elongating_10",
                 "1506F4" = "Elogated_13",
                 "1506F5" = "ESC_JM8",
                 "1506F6" = "gRounds")


ESC.data <- laura.data[["83_nonTSS_meCGI_H_S_11_P0R"]][["1506F5"]]
#ESC.data <- laura.data[["2_Cntrl_5"]][["1506F5"]]
Meth.mat.GpC <- ESC.data[["Meth.mat"]]
Meth.mat.GpC[,ESC.data$CpG.positions] <- NA
Meth.mat.GpC[,ESC.data$GpC.positions] <- 1 - Meth.mat.GpC[,ESC.data$GpC.positions]
row.names(Meth.mat.GpC) <- paste0("Read",1:nrow(Meth.mat.GpC))
create_data_list <- function(matr){
  stopifnot(is.matrix(matr))
  if(is.null(row.names(matr))){
    row.names(matr) <- 1:nrow(matr)
  }
  
  dat.out <- lapply(row.names(matr),function(nm){
    seq <- matr[nm,]
    seq[is.na(seq)] <- 2
    seq <- paste(seq,collapse = "")
    list("DATA" = seq,
         "NAME" = nm)
  })
  dat.out
}
nome_data <- create_data_list(Meth.mat.GpC)


sourceCpp("../src/rcpp_hello_world.cpp")

test.out <- rcpp_hello_world(data=nome_data,
                 binding_models = binding_models[1],
                 bgcoverprob = 0.01,
                 bgprior = 0.5,
                 bound_fit_tol = 1e-5,
                 bound_min_fderiv_val = 1e-5,bound_max_steps = 1000,
                 run_priorEM = FALSE,
                 priorEM_fit_tol = 1e-5,
                 priorEM_max_steps = 1000)


a <- as.data.frame(test.out[[2]],stringsAsFactors = F)
library(reshape2)
dat <- melt(Meth.mat.GpC)
library(ggplot2)
a$dat <- dat[match(paste(a$seq,a$pos),
                   paste(as.character(dat$Var1),dat$Var2)),"value"]
ggplot(subset(a,seq %in% paste0("Read",1:20))) + 
  facet_wrap(~seq)+
  geom_point(aes(x=pos,y=Nucleosome),color="red")+
  #geom_point(aes(x=pos,y=TF),color="blue")+
  geom_point(aes(x=pos,y=dat),color="black")
  

cover.prob <- as.data.frame(test.out[[2]],stringsAsFactors = F)
cover.prob$dat <- dat[match(paste(cover.prob$seq,cover.prob$pos),
                            paste(as.character(dat$Var1),dat$Var2)),"value"]

start.prob <-  as.data.frame(test.out[[1]],stringsAsFactors = F)
nucl.cover.prob.matr <- dcast(data=cover.prob,seq ~ pos,value.var = "Nucleosome")

tf.cover.prob.matr <- dcast(data=cover.prob,seq ~ pos,value.var = "TF")
library(ggplot2)

cov.plot <- ggplot(subset(cover.prob,seq %in% paste0("Read",1))) + 
  facet_wrap(~seq)+
  geom_point(aes(x=pos,y=Nucleosome),color="red")+
  geom_point(aes(x=pos,y=background),color="green")+
  geom_line(aes(x=pos,y=Nucleosome),color="red")+
  geom_line(aes(x=pos,y=background),color="green")+
  
  geom_vline(xintercept = c(-76,-76 + 147),color="grey",linetype="dashed")+
  #geom_point(aes(x=pos,y=TF),color="blue")+
  geom_point(aes(x=pos,y=dat),color="black")+
  theme_bw()
start.plot <- ggplot(subset(start.prob,seq %in% paste0("Read",1))) + 
  facet_wrap(~seq)+
  geom_point(aes(x=pos,y=Nucleosome),color="red")+
  geom_point(aes(x=pos,y=background),color="green")+
  geom_line(aes(x=pos,y=Nucleosome),color="red")+
  geom_line(aes(x=pos,y=background),color="green")+
  
  geom_vline(xintercept = c(-76,-76 + 147),color="grey",linetype="dashed")+
  #geom_point(aes(x=pos,y=TF),color="blue")+
  #geom_point(aes(x=pos,y=dat),color="black")+
  theme_bw()

library(gridExtra)
grid.arrange(cov.plot,start.plot,ncol=1)

plot(colMeans(nucl.cover.prob.matr[,-1]),pch=20,col="grey",ylim=c(0,1))
points(colMeans(tf.cover.prob.matr[,-1]),pch=20,col="blue")
points(colMeans(Meth.mat.GpC,na.rm = T),col="red",pch=20)



#### Most likely configuration? ####

maxLKHconf <- function(cov.prof,model.names = c("Nucleosome","background")){
  #browser()
  colnm <- colnames(cov.prof)[colnames(cov.prof) %in% model.names]
  maxlk <- apply(cov.prof[,colnm],1,which.max)
  cov.prof$Max.cov <- colnm[maxlk]
  cov.prof
}

max.conf <- maxLKHconf(subset(cover.prob,pos >= -147 & pos <= nrow(Meth.mat.GpC) + 147))

nucl.len <- do.call(rbind,by(max.conf,max.conf$seq,function(dat){
  #browser()
  rle.d <- rle(dat$Max.cov)
  len <- rle.d$lengths
  val <- rle.d$values
  n.pos <- which(val == "Nucleosome")
  pos.in.read <- dat$pos[cumsum(len) - len + 1]
  if(length(n.pos) == 0)
    return(NULL)
  data.frame(seq=dat$seq[1],
             start.pos = pos.in.read[n.pos],
             nucl.len = len[n.pos],
             stringsAsFactors = F)
  
}))


tail(subset(max.conf,pos >= -77 & pos <= -76+143+1 & seq=="Read1"))


nucl.len$end.pos <- nucl.len$start.pos + 147
nucl.len$outside <- ifelse(nucl.len$start.pos <=0 | nucl.len$end.pos > nrow(Meth.mat.GpC),"out","in" )


st.change <- diff(ifelse(max.conf$Max.cov == "background",0,1))





