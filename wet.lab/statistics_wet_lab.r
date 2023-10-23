t_test <- function(data
                  ,...){
    t.test(value~group
           ,data
           ,...
          )
}

welch_test <- function(data
                  ,...){
    t.test(value~group
           ,data
           ,var.equal = FALSE
           ,...
          )
}

mann_whitney_u_test <- function(data
                               ,...){
    wilcox.test(value~group
                ,data
                ,exact = NULL
                ,correct = TRUE
                ,conf.int = FALSE
                ,...
               )
}

moods_median_test <- function(data
                             ,...){
    mood.test(value~group
              ,data
              ,...)
}

run_statistics <- function(experiment
                         ,groups_to_test = NA 
                          # If there are more than two groups in the dataset, by default all vs all groups will be tested. 
                          # Provide character vector in form of e.g. c("het vs WT", "KO vs WT") to perform only selected comparisons
                         ,adjustment_method = "bonferroni"
                          ){
    # read in data
    data <- read.csv(file = paste0("./input/"
                                  ,experiment
                                  ,".tsv"
                                  )
                    ,sep = "\t"
                     ,dec = ","
                    ,header = TRUE)
    print(str(data)) # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cat('\n')
    
    # how many groups we have?
    groups <- unique(data$group)
    print("we have following groups in the data:")
    print(groups) # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cat('\n')
    
    all_vs_all <- (length(groups >2)) & is.na(groups_to_test) # TRUE or FALSE
    
    # do we have biological replicates?
    ### if yes, paired version will be used
    paired <- sum(!is.na(data$bio_rep)) # TRUE or FALSE
    if(paired) {print("data are paired observations")}
    cat('\n')

    # pick the correct test
    ### check the normality of distribution assumption in all groups
    print("run Shapiro test to check the assumption of the distribution normality...")
    not_norm <- sapply(groups
                   ,function(group){
                       idx_group <- data$group == group
                       pvalue <- shapiro.test(data$value[idx_group])$p.value
                       
                       print(group) # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       ifelse(pvalue <= 0.05,print("non-normal"),print("normal")) # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       pvalue <= 0.05
                   })
    norm_dist <- sum(not_norm)== 0 # TRUE or FALSE
    print("overall distribution is:") # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ifelse(norm_dist,print("normal"),print("non-normal")) # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cat('\n')
    
    ### check the homogeneity of variance assumption in all groups
    not_homogeneic <- if(norm_dist){
        print("run Bartlett test to check the assumption of the variance homogeneity...")
        pvalue <- bartlett.test(value~group, data = data)$p.value 
        pvalue <= 0.05
    } else {
        print("run Fligner-Killeen test to check the assumption of the variance homogeneity...")
        pvalue <- fligner.test(value~group, data = data)$p.value 
        pvalue <= 0.05
    }
    
    homo_var <- sum(not_homogeneic)== 0 # TRUE or FALSE
    print("overall variance is:") # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ifelse(homo_var,print("stable"),print("non-stable")) # REMOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ### pick the test
    which_test <- if(norm_dist){ 
        ifelse(homo_var
              ,"t-test"  # distribution is normal, variance is stable
              ,"Welch test"  # distribution is normal, variance is not stable
              )
    }else {
        ifelse(homo_var
              ,"Mann-Whitney U-test"  # distribution is not normal, variance is stable
              ,"Mood's Median-test" # distribution is not normal, variance is not stable
              )
    }
    
    ### run the test
    run_one_comparison <- function(data){
        if(which_test == "t-test"){
            t_test(data, paired = paired)$p.value
    } else if(which_test == "Welch test"){
        welch_test(data, paired = paired)$p.value
    } else if(which_test == "Mann-Whitney U-test"){
        mann_whitney_u_test(data, paired = paired)$p.value
    } else if(which_test == "Mood's Median-test"){
        moods_median_test(data)$p.value
    } else error("ERROR: undefined test")
        }
    
    p.values <- sapply(groups_to_test
                      ,function(comparison){
                          group1 <- sub(" vs .*","",comparison)
                          group2 <- sub(".* vs ","",comparison)
                          
                          idx <- data$group %in% c(group1, group2)
                          run_one_comparison(data[idx,])
                      })
    
    # adjust for multiple testing if needed
    p.adjs <- p.adjust(p.values
                      ,method = adjustment_method)
        
    
    output <- data.frame(comparison = groups_to_test
                        ,p.value = p.values)
        
    if(length(groups_to_test)> 1){output$p.adj = p.adjs}
    
    # print mesages
    
    cat('\n')
    print(paste0("Since the normality of distribution assumption is "
                ,ifelse(norm_dist
                       ,"met"
                       ,"not met"
                       )
                 ," and the homogeneity of variace assumption is "
                 ,ifelse(homo_var
                        ,"met"
                        ,"not met"
                        )
                 ,", we use the "
                 ,which_test
                 ,"."
                )
         )
    
    cat('\n')
    if(paired){print(paste0("Since we have several independent biological replicates as paired observations across groups, we use the paired version of the "
                           ,which_test
                           ," test.")
                    )}
    
    
    cat('\n')
    print("We are doing the comparison of the following group(s):")
    print(ifelse(all_vs_all
          ,"all vs all"
          ,groups_to_test))
    
    
    cat('\n')    
    print(output)

    # export results
        write.table(output
                 ,file = paste0("./output/"
                                ,experiment
                                ,"_stats.tsv"
                                )
                 ,sep = "\t"
                 ,quote = FALSE
                 ,col.names = TRUE
                 ,row.names = FALSE)
    
}

run_statistics("Fig1C_DGEP"
              ,groups_to_test = "MUT vs WT")

run_statistics("Fig1D_QMI"
              ,groups_to_test = "MUT vs WT")

run_statistics("Fig2B_FACS_OCI-Ly1"
              ,groups_to_test = c("het vs WT"
                                 ,"KO vs WT"
                                 )
              )

run_statistics("Fig2B_FACS_OCI-Ly8"
              ,groups_to_test = c("het vs WT"
                                 ,"KO vs WT"
                                 ,"het+ARID1A vs WT")
              )

run_statistics("Fig2E_qPCR_OCI-Ly1_raw"
              ,groups_to_test = c("het vs WT"
                                 ,"KO vs WT"
                                 )
              )

run_statistics("Fig2E_qPCR_OCI-Ly8_raw"
              ,groups_to_test = c("het vs WT"
                                 ,"KO vs WT"
                                 ,"het+ARID1A vs WT")
              )

run_statistics("Fig3E_qPCR_OCI-Ly1"
              ,groups_to_test = c("het vs WT"
                                 ,"KO vs WT"
                                 )
              )

run_statistics("Fig3E_qPCR_OCI-Ly8"
              ,groups_to_test = c("het vs WT"
                                 ,"KO vs WT"
                                 ,"het+ARID1A vs WT")
              )

run_statistics("Fig4C_luc_promoter1"
              ,groups_to_test = c("50 ng vs 0 ng"
                                 ,"200 ng vs 0 ng"
                                 ,"500 ng vs 0 ng")
              )

run_statistics("Fig4C_luc_promoter2"
              ,groups_to_test = c("50 ng vs 0 ng"
                                 ,"200 ng vs 0 ng"
                                 ,"500 ng vs 0 ng")
              )

run_statistics("Fig4F_qPCR_OCI-Ly8"
              ,groups_to_test = c("het+RUNX3 vs het"
                                 ,"het+RUNX3 vs WT"
                                 )
              )

run_statistics("Fig4G_FACS_OCI-Ly8"
              ,groups_to_test = c("het+RUNX3 vs het"
                                 ,"het+RUNX3 vs WT"
                                 )
              )

run_statistics("Fig5B_FACS_OCI-Ly1"
              ,groups_to_test = c("het 0 ng vs WT 0 ng"
                                 ,"het+RUNX3 0 ng vs WT 0 ng"
                                 ,"KO 0 ng vs WT 0 ng"
                                  ,"het 3 ng vs WT 3 ng"
                                 ,"het+RUNX3 3 ng vs WT 3 ng"
                                 ,"KO 3 ng vs WT 3 ng"
                                  ,"het 30 ng vs WT 30 ng"
                                 ,"het+RUNX3 30 ng vs WT 30 ng"
                                 ,"KO 30 ng vs WT 30 ng"
                                  ,"het 300 ng vs WT 300 ng"
                                 ,"het+RUNX3 300 ng vs WT 300 ng"
                                 ,"KO 300 ng vs WT 300 ng"
                                 )
              )

run_statistics("Fig5B_FACS_OCI-Ly8"
              ,groups_to_test = c("het 0 ng vs WT 0 ng"
                                 ,"het+RUNX3 0 ng vs WT 0 ng"
                                 ,"KO 0 ng vs WT 0 ng"
                                  ,"het 3 ng vs WT 3 ng"
                                 ,"het+RUNX3 3 ng vs WT 3 ng"
                                 ,"KO 3 ng vs WT 3 ng"
                                  ,"het 30 ng vs WT 30 ng"
                                 ,"het+RUNX3 30 ng vs WT 30 ng"
                                 ,"KO 30 ng vs WT 30 ng"
                                  ,"het 300 ng vs WT 300 ng"
                                 ,"het+RUNX3 300 ng vs WT 300 ng"
                                 ,"KO 300 ng vs WT 300 ng"
                                 )
              )

run_statistics("Fig5E_FACS_OCI-Ly8"
              ,groups_to_test = c("het+RUNX3 vs het"
                                 ,"het+RUNX3 vs WT"
                                 )
              )

run_statistics("FigS1B_FACS"
              ,groups_to_test = c("MUT vs WT")
              )

run_statistics("FigS2D_qPCR_OCI-Ly1"
              ,groups_to_test = c("het vs WT"
                                 ,"het+RUNX3 vs WT"
                                 )
              )

run_statistics("FigS2E_FACS_OCI-Ly1"
              ,groups_to_test = c("het+RUNX3 vs het"
                                 ,"het+RUNX3 vs WT"
                                 )
              )

sessionInfo()


