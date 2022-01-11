#####################################################################################################
# this code file contains all the functions and the main program used to implement the MCCI method; #              #
# users don't have to modify this code file                                                         #
#####################################################################################################

library(MplusAutomation)
library(doSNOW)
library(foreach)
curpath <<- getwd()
simulate_datapath <- 'random_samples_generation/'
datalist <- 'datalist.dat'

## function for creating new folders to store Mplus output files
createDirForOutput <- function(c_dirname) {
  c_outpath <- c_dirname
  if (!file.exists(c_outpath)){ 
    dir.create(c_outpath, recursive=TRUE, showWarnings=FALSE) 
  }
  simulate_outpath <<- paste0(c_outpath,'/simulateddata_output')
  if (!file.exists(simulate_outpath)){ 
    dir.create(simulate_outpath, recursive=TRUE, showWarnings=FALSE) 
  }
  empirical_outpath <<- paste0(c_outpath,'/empiricaldata_output')
  if (!file.exists(empirical_outpath)){ 
    dir.create(empirical_outpath, recursive=TRUE, showWarnings=FALSE) 
  }
}

## function for analyzing the combined sample with the CFA model
combined_sample_est <- function(c_inputdata, out_file_name, cur_savepath, model_var, model_model, cur_title, cur_output, cur_savedata) {
  setwd(cur_savepath)
  retfile1 <- paste0(out_file_name,".inp")
  pathmodel <- mplusObject(
    TITLE = cur_title,
    VARIABLE = model_var,
    MODEL = model_model,
    OUTPUT = cur_output,
    SAVEDATA = cur_savedata,
    rdata = c_inputdata)
  fit <- mplusModeler(pathmodel, modelout = retfile1, run = 1L)
  unlink(fit$results$input$data$file)
  unlink(retfile1)
  return (fit)
}

## function for generating random samples 
random_samples_generation <- function(out_file_name, cur_savepath, cur_title, cur_montecalo, cur_population) {
  setwd(cur_savepath)
  retfile1 <- paste0(out_file_name,".inp")
  pathmodel <- mplusObject(
    TITLE = cur_title,
    MONTECARLO = cur_montecalo,
    MODELPOPULATION = cur_population)
  fit <- mplusModeler(pathmodel, modelout = retfile1, run = 1L)
  unlink(fit$results$input$data$file)
  unlink(retfile1)
}

## function for analyzing the empirical data and the random samples with the two-group CFA model (e.g., the configural model for the loading tests)
multigroup_cfa_est <<- function(c_inputdata, out_file_name, cur_savepath, model_var, model_model, cur_title, cur_output) {
  library(MplusAutomation)
  setwd(cur_savepath)
  retfile1 <- paste0(out_file_name,".inp")
  pathmodel <- mplusObject(
    TITLE = cur_title,
    VARIABLE = model_var,
    MODEL = model_model,
    OUTPUT = cur_output,
    rdata = c_inputdata)
  fit <- mplusModeler(pathmodel, modelout = retfile1, run = 1L)
  unlink(fit$results$input$data$file)
  unlink(retfile1)
  return (fit)
}

## functions (getdiff, dodiffcore, and dodiff) for obtaining the between-group differences in measurement parameters
getdiff <- function(fit1) {
  para_est <- fit1$parameters$unstandardized$est
  para_num2 <- length(para_est)
  para_num1 <- para_num2 / 2
  est1 <- para_est[1:para_num1]
  est2 <- para_est[(para_num1+1):(para_num2)]
  diff1 <- est1 - est2
  return (diff1)
}

dodiffcore <- function(){
  curname <- outfile_names[i]
  cur1 <- unlist(strsplit(curname, split="\\."))
  dataname <- cur1[1]
  curoutfile <- paste0(curpath,'/',fileout,'/',curname)
  cur_result <- readModels(curoutfile,what="all")
  curdiff <- getdiff(cur_result)
  c(dataname,curdiff);
}

dodiff <- function(pinpath,numpara) {
  simu_diff <- data.frame(matrix(NA,numpara,0))
  fileout <- paste0(pinpath)
  outfile_names <- list.files(fileout, pattern='*.out')
  funclst <- c("getdiff","curpath","readModels","simu_diff","fileout","dodiffcore","outfile_names")
  simu_diff <- foreach(i=1:length(outfile_names), .verbose=T, .export=funclst, .combine="cbind") %dopar% dodiffcore()
  colnames(simu_diff) <- simu_diff[c(1),]
  simu_diff <- simu_diff[c(2:nrow(simu_diff)),]
  write.csv(simu_diff, file=paste0(pinpath,"/diff.csv"))
}

## function for obtaining the population values for random samples generation 
main_flow_combine_sample_est <- function() {
  setwd(curpath)
  if (!file.exists(simulate_datapath)){ 
    dir.create(simulate_datapath, recursive=TRUE, showWarnings=FALSE) 
  }
  empirical_data <- read.csv(empirical_datapath, header = FALSE)
  empiricalfit <- combined_sample_est(empirical_data,'combined_sample',simulate_datapath, model_variable, model_model, model_title, model_output, model_savedata)
  
  combined_estimates <- read.csv("combined_sample_estimates.dat",header = FALSE)
  copy_estimates <- combined_estimates
  population_values<-rbind(combined_estimates,copy_estimates)
  write.table(population_values,"combined_sample_estimates.dat",row.names = FALSE, col.names = FALSE,quote = FALSE)
}

## function for performing random samples generation
main_flow_random_samples <- function() {
  setwd(curpath)
  random_samples_generation("random_samples_generation", simulate_datapath, model_title, model_montecarlo, model_population)
}

#######################################################################################
#           main program for the implementation of the MCCI method                    #
#######################################################################################
meta_flow_core <- function() {
  setDefaultClusterOptions(master="localhost")
  clusterobj <- makeSOCKcluster(njobs)
  registerDoSNOW(clusterobj)
  setwd(curpath)
  para_num<<-3*item_num+choose(factor_num,2)+2*factor_num

############## Step 1:analyze the observed data and obtain the observed differences in the parameters ##############
  empirical_data <- read.csv(empirical_datapath, header = FALSE)
  empiricalfit <- multigroup_cfa_est(empirical_data,'empirical',empirical_outpath, model_variable, model_model, model_title, model_output)
  setwd(curpath)
  empirical_results <- readModels(paste0(empirical_outpath,"/empirical.out"))
  paraheader <- paste(empirical_results$parameters$unstandardized$paramHeader,empirical_results$parameters$unstandardized$param)
  empirical_diff <- getdiff(empirical_results)
  
############## Step 2:use the Monte Carlo method to approximate the null distribution of the differences in the parameters ##############
  ## analyze each of the random samples
  setwd(curpath)
  datalist_path_tmp <<- paste0(simulate_datapath,datalist)
  simu_datalist <<- as.matrix(read.table(datalist_path_tmp))
  funclst <- c("simu_datalist", "multigroup_cfa_est", "model_variable", "model_model", "model_title", "model_output", "curpath", "simulate_datapath", "simulate_outpath")
  foreach(r=1:nrow(simu_datalist), .verbose=T, .export=funclst, .combine="c") %dopar% {
    setwd(curpath)
    simu_dataname <- simu_datalist[r]
    fileone <- paste0(simulate_datapath,simu_dataname)
    cur1 <- read.table(fileone)
    curdatacc <- data.frame(cur1)
    filetmpgroup <- unlist(strsplit(simu_dataname, split="\\."))
    curfname <- filetmpgroup[1]
    multigroup_cfa_est(curdatacc, curfname, simulate_outpath, model_variable, model_model, model_title, model_output)
  }
  
  # calculate the between-group differences in parameters for each of the random samples
  setwd(curpath)
  dodiff(simulate_outpath,para_num)
  
############## Step 3:establish confidence intervals for the differences in the parameters under the null hypothesis ##############
  diff<-read.csv(paste0(simulate_outpath,"/diff.csv"))
  diff2<-diff[,2:1001]
  interval<- matrix(NA,para_num,3)
  for (i in 1:para_num) {
    interval[i,1] <- quantile(diff2[i,],probs=0.025,na.rm = TRUE)
    interval[i,2] <- quantile(diff2[i,],probs=0.975,na.rm = TRUE)

############## Step 4: evaluate the invariance of each item ##############    
    if(empirical_diff[i] >= interval[i,1] && empirical_diff[i] <= interval[i,2]){
      interval[i,3] <- 1
    }
    else{
      interval[i,3] <- 0
    }
  }
  ## save the testing results
  empirical_and_CI<- cbind(paraheader[1:para_num],empirical_diff,interval)
  colnames(empirical_and_CI)<-c("Parameter","Observed between-group difference","2.5% quantile","97.5% quantile","0-noninvariant,1-invariant")
  result_by_tag <- NULL
  if(process_tag == 'loading') {
    result_by_tag <- empirical_and_CI[1:item_num,]
  }else if(process_tag == 'intercept') {
    result_by_tag <- empirical_and_CI[(item_num+choose(factor_num,2)+factor_num+1):(2*item_num+choose(factor_num,2)+factor_num),]
  }else if(process_tag == 'residualvariance') {
    result_by_tag <- empirical_and_CI[(para_num-item_num+1):para_num,]
  }
  write.csv(result_by_tag, file=paste0(process_tag,"_testresults.csv"))
  stopCluster(clusterobj)
}


#################### testing for metric invariance (factor loading) ################
main_flow_loading_1 <- function() {
  setwd(curpath)
  process_tag <<- 'loading'
  loading_outpath <- 'loading_outputfile/'
  createDirForOutput(loading_outpath)
  meta_flow_core()
}

#################### testing for scalar invariance (intercept) #####################
main_flow_intercept_2 <- function() {
  setwd(curpath)
  process_tag <<- 'intercept'
  loading_outpath <- 'intercept_outputfile/'
  createDirForOutput(loading_outpath)
  meta_flow_core()
}

#################### testing for strict invariance (residual variance) #############
main_flow_residualvariance_3 <- function() {
  setwd(curpath)
  process_tag <<- 'residualvariance'
  loading_outpath <- 'residualvariance_outputfile/'
  createDirForOutput(loading_outpath)
  meta_flow_core()
}