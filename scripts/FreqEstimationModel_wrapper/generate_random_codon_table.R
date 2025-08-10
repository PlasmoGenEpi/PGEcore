setwd("C:/Users/zoeor/OneDrive/Documents/Life/Collaboration with Dr. Wesolowski/")
amino_acids <- read.csv("codon-table-grouped.csv")
rownames(amino_acids) <- amino_acids$codon
create_random_dataset <- function(){
  targets <- c("csp", "gama", "protein_x")
  target_positions <- list(c(72, 81), c(300, 27, 309), c(13))
  normal_AA <- list(c("ATC", "GCG"), c("ATG", "CTC", "AGC"), c("CCC"))
  
  samples <- c("Sample 1", "Sample 2", "Sample 3")
  random_dataset <- data.frame(sample_id=c("A"), 
                               target=c("A"), 
                               position=c("A"), 
                               codon=c("A"), 
                               AA=c("A"))
  for(sample in rep(samples,2)){
    for(target_ix in 1:length(target_positions)){
      target_position <- target_positions[[target_ix]]
      normal_AA_target <- normal_AA[[target_ix]]
      for(target_position_ix in 1:length(target_position)){
        randomize <- FALSE
        if(runif(1) < 0.1){
          randomize <- TRUE
        }
        skip <- FALSE
        if(runif(1)<0.){
          skip <- TRUE
        }
        if(randomize){
          codon_ix <- sample(x=1:64, size=1)
          codon <- amino_acids[codon_ix, 2]
        }
        else{
          codon <- normal_AA_target[target_position_ix]
        }
        print(codon)
        letter <- amino_acids[codon, 1] 
        print(target_position[[target_position_ix]])
        if(!skip){
          random_dataset[nrow(random_dataset)+1,] <- c(sample, targets[target_ix], 
                                                       target_position[target_position_ix],codon, letter)
          
        }
      }
    }
  }
  random_dataset <- random_dataset[random_dataset$sample_id != "A",]
  return(random_dataset)
}

test_dataset <- create_random_dataset()
write.csv(test_dataset, "test_dataset.csv")
