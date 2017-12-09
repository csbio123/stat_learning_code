#P value analysis conrad
pval_anal <- function(path) {
  library(ggplot2)
  
  filenames <-
    list.files(path,
               pattern = "*acc*",
               full.names = TRUE,
               recursive = TRUE)
  files = sapply(
    filenames,
    read.csv,
    sep = " ",
    simplify = FALSE,
    header = FALSE
  )
  print(length(filenames))
  
  
  my_files = sapply(c(1:length(files)),  function(f) {
    name = basename(dirname(filenames[f]))
    my_f = files[[f]]
    my_f[, 8] = rep(name, 1100)
    my_f
  }, simplify = FALSE)
  
  results  = do.call(rbind, my_files)
  groups = unique(results[, 8])
  
  groups_ann = sapply(1:length(groups), function(x) {
    samples = length(which(results[, 8] == groups[x])) / 1100
    ann = unlist(sapply(1:samples,  function(f) {
      rep(f, 1100)
    }, simplify = FALSE))
    
    unlist(ann)
  }, simplify = FALSE)
  
  results[9] = unlist(groups_ann)
  colnames(results)  =  c(
    "lambda_min",
    "error_min",
    "lambda_dev",
    "error_dev",
    "alpha",
    "iteration",
    "parameter_combi",
    'run_type',
    'dataset'
  )
  
  
  mean_se <- function(x, mult = 100) {
    x <- na.omit(x)
    se <- mult * sqrt(var(x) / length(x))
    mean <- mean(x)
    data.frame(y = mean,
               ymin = mean - se,
               ymax = mean + se)
  }
  
  mean_t = mean(results[, 2])
  sd_t = sd(results[, 2])
  #line_plot = ggplot(data = results, aes(dataset, (error_dev), group = run_type, color =
                                         #  run_type)) + stat_summary(alpha = 0.1) +  geom_smooth() +
    #theme(
    #  panel.grid.major = element_blank(),
   #   panel.grid.minor = element_blank(),
    #  panel.background = element_blank(),
   #   axis.line = element_line(colour = "black")
   # ) + scale_y_continuous(limits = c(0, 100))
  #ggsave(paste0(path, "/gene_runs_lineplot.pdf"), line_plot)
  
  
  #scale_y_continuous(limits = quantile(log10(results$error_dev), c(0.00, 0.95))) removes the 5% outliers
  my_plot = ggplot(data = results, aes(run_type, (error_dev)), outlier.shape = NA) + geom_boxplot(aes(color =
                                                                                 run_type)) + scale_y_continuous(limits = quantile((results$error_dev), c(0.00, 0.90)))
  
  
  # compute lower and upper whiskers
 # ylim1 = boxplot.stats(results$error_dev)$stats[c(1, 5)]
  
  # scale y limits based on ylim1
  #p1 = my_plot + coord_cartesian(ylim = ylim1*1.05)
  
  ggsave(paste0(path, "/gene_runs_boxplot.pdf"), my_plot)
  
}
#arg function below is in-built. It allows you to interact with the script from the command line.

#####USE AS FOLLOWS#####:
#Rscript p_value_analysis.R /users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/all_transcripts 476####

args <-commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc
print(paste0('Input directory for p value analysis:', args[1]))#This line will tell you the input directory

args[1] = "C:/Users/spjtcoi/Google Drive/bmi/categorical/classification_bmi/plot/plot"
#args[1] = "/Users/ti1/Documents/conrad/results/snps2"
if (length(args) == 1) {
  pval_anal(args[1])#if a single argument is entered then that paramater should be the directory path
} 

