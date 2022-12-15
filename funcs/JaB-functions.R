########################
### Needed Libraries
########################
library(igraph)
library(snowboot)
library(dplyr)

########################
### Internal Functions 
########################

find_central_stat <- function(g, central_stat){
  ## Calculate given centrality statistics of the network 
  
  ### input
  # g = igraph object 
  #central_stat = centrality statistic we want to consider, options are listed in null_dist_algorithm function
  #   options: "betweenness", "closeness" "eigen", "authority", "hub"
  
  ### returns centrality statistic of graph g
  
  if(central_stat == "betweenness"){
    return(betweenness(g))
  }else if(central_stat == "closeness"){
    return(closeness(g))
  }else if(central_stat == "eigen"){
    return(eigen_centrality(g)$vector)
  }else if(central_stat =="authority"){
    return(authority_score(g)$vector)
  }else if(central_stat == "hub"){
    return(hub_score(g)$vector)
  }else if(central_stat == "degree"){
    return(unname(degree(g)))
  }
}


get_bootstrap <- function(g, B, U, num_seed = NA, num_wave = NA, central_stat, orig_stat, mult = TRUE){
  ##Bootstrapping step of JaB algorithm with LSMI method
  
  ### Inputs 
  # g = igraph object of original network
  # U = number of unique networks to sample
  # B = number of bootstrap samples within each u for i in 1:U, only needed if type = "snowboot" 
  # num_seed = number of seeds to use in snowboot, only needed if type = "snowboot" 
  # num_waves = number of waves to use in snowboot, only needed if type = "snowboot" 
  # centrality_stat = type of centrality statistic we will consider. 
  # orig_stat = vector of centrality statistics from g, output of `find_central_stat`
  # mult = if multiplicity is to be used, TRUE if using multiplicity, FALSE otherwise 
  
  ### Outputs a list of 
  # 1. list of nodes in each bootstrap sample 
  # 2. centrality statisitc of all nodes in bootstrap sample 
  # 3. number of times each node was in each bootstrap sample
  
  
  ### Convert igraph object to snowboot object 
  net <- igraph_to_network(g)
  
  ### Initialize storage 
  boot_node <-  vector(mode = "list", length = B*U)
  boot_central_stat <- vector(mode = "list", length = B*U)
  boot_num_times <- vector(mode = "list", length = B*U)
  
  num_v <- gorder(g) #number of nodes 

  for(u in 1:U){
    sample_network <- lsmi(net, n.seed=num_seed, n.wave = num_wave)
    
    #get list of seeds and non-seed nodes 
    wave_nodes <- c()
    seed_nodes <- rep(NA, num_seed)
    for(i in 1:num_seed){
      patch_i <- unlist(sample_network[[i]])
      wave_nodes <- c(wave_nodes,  patch_i[-1]) #remove the seed 
      seed_nodes[i] <- patch_i[1] #just the seed 
    }
    
    #make numeric
    wave_nodes <- as.numeric(wave_nodes)
    seed_nodes <- as.numeric(seed_nodes)
    
    #size of the wave set
    wave_size <- length(wave_nodes)
    
    #bootstrap
    for(b in 1:B){
      sample_seed <- sample(seed_nodes, num_seed, replace = TRUE)
      sample_wave <- sample(wave_nodes, wave_size, replace = TRUE)
      
      current_index <- (u-1)*B + b
      
      if(mult){ #Has multiplicity in it
        boot_node[[current_index]] <- c(sample_seed, sample_wave)
      } else{ #does not have multiplicity in it
        boot_node[[current_index]] <- unique(c(sample_seed, sample_wave))
      }
      
      #if the seed has no neighbors the waves will just say NA
      # we remove the NA's 
      # this is equivalent to never include that seed's waves in the sampling process 
      boot_node[[current_index]] <- boot_node[[current_index]][ !is.na(boot_node[[current_index]])]
      
      boot_central_stat[[current_index]] <- unname(orig_stat[boot_node[[current_index]]] )
      boot_num_times[[current_index]] <- table(boot_node[[current_index]])
      
    }#end B for
  }#end U for
  
  return(list(boot_node, boot_central_stat, boot_num_times))
  
}


get_jackknife <- function(boot_result, n, B, U, quants){
  ## Jackknife-after step in JaB
  
  ### Inputs 
  # boot_result = output from get_bootstrap function
  # n = number of nodes in original graph 
  # B = # of bootstrap samples 
  # U = number of samples
  # quants = quantile of null distribution
  
  ### Outputs 
  # 1. Lower and Upper quantiles 
  # 2. Null Distribution for each node 
  
  node_list <- boot_result[[1]]
  stat_list <- boot_result[[2]]
  numtime_list <- boot_result[[3]]
  
  JaB <- matrix(NA, nrow = n, ncol = 2)
  
  #centrality_storage <- list(n) #used for exploration purposes, not actually needed in algorithm
  
  for(i in 1:n){
    #get index of all bootstrap samples that do not have node i
    keep_index <- rep(NA,length(node_list))
    for(b in 1:length(node_list)){
      keep_index[b] <- !( i %in% node_list[[b]])
    }
    
    #Get centrality statistic for all networks in keep_index
    #This repeats each vertex each number of times it was samples 
    # i.e. if in bootstrap sample b, node i was sampled 3 times, there will be 3 repetitions of the i^th statistics in the null distribution 
    centrality_vector <- c(unlist(stat_list[keep_index]))
    
    #sum(is.na(centrality_vector))
    
    #just for exploring null distribution shape, not actually needed to run the algorithm
    #centrality_storage[[i]] <- centrality_vector
    
    #Get Quantiles 
    JaB[i, ] <- quantile(centrality_vector, quants, na.rm = T) 
    
  }
  
  return(JaB)
}



#########################
### JaB Algorithm
### Generate null distribution 
###########################

JaB <- function(graph,  
                B=1000 , U = 10,
                num_seed=NA, num_wave=NA, 
                centrality_stat,
                quants =  c(0, 95)/100,
                mult = TRUE,
                print_null_dist = T){
  
  ### Inputs 
  # graph = igraph object of original network
  # B = number of bootstrap samples within each u
  # U = number of unique networks to sample
  # num_seed = number of seeds to use in snowboot, only needed if boot_type = "snowboot" 
  # num_waves = number of waves to use in snowboot, only needed if boot_type = "snowboot" 
  # centrality_stat = type of centrality statistic we will consider. 
  #   options: "betweenness", "closeness" "eigen", "authority.score", "hub.score"
  # quants = c(lower quantile, upper quantile) for jackknife 
  # mult = if multiplicity is to be used, TRUE if using multiplicity, FALSE otherwise 
  # print_null_dist = T if you want to include the simulated null distribution for each node, 
  #                   F if you do not want to include this in functino output.
  #                   This storage can be very large, especially for large networks 
  
  ### Output a list of 
  # 1. "$Results" : Table of node number, original statistic, JaB quantiles ordered from most influential to least influential 
  # 2. "$Null_Distribution": Entire empirical null distribution for each node 
  #     2. is not actually needed to complete the analysis, it is used for exploratory purposes 
  
  
  #number of nodes in graph
  num_v <- vcount(graph)
  orig_nodes <- V(graph)$name
  V(graph)$name <- 1:num_v
  
  ### Get original centrality statistic 
  orig_stat <- find_central_stat(graph, centrality_stat)
  
  ### Bootstrap step 
  boot_result <- get_bootstrap(graph, B , U, num_seed, num_wave, centrality_stat, orig_stat, mult)
  
  ## Jackknife step
  jack_result <- get_jackknife(boot_result, num_v, B,U, quants)
  
  ### Return something useful
  ret <- data.frame(1:num_v, orig_nodes, orig_stat, jack_result)
  colnames(ret) <- c("Node_Number", "Node_Name", "Orig_Stat", "Lower", "Upper")
  
  #order from most outside on CI to least
  ret <- ret %>%
    mutate(diff = Orig_Stat - Upper) %>%
    arrange(desc(diff)) %>%
    mutate(Influential = diff > 0) %>%
    mutate(diff = abs(diff)) %>%
    mutate(rank = row_number()) %>%
    select(-diff)
  
  #for exploratory purposes - this part is too much storage for large networks 
  #ret1 <- list(ret, jack_result[[2]])
  #names(ret1) <- c("Results", "Null_Distribution")
  
  #for exploratory purposes - 
  if(print_null_dist){
    ret1 <- list("Results" = ret, 
                 "Bootstrap_Nodes" = boot_result[[1]])
  }else{
    ret1 <- list("Results" = ret)
  }
  
  
  return(ret1)
  
}



######################################
### Function to plot null distribution 
#######################################

null_dist_large_plot <- function(large_output, nn, orig_line=TRUE, upper_line= TRUE, breaks=30){
  ### Inputs 
  # large_output = output of null_dist_algorithm2_large
  # nn = node number you want to plot 
  # orig_line = if TRUE, will plot a red vertical line on the histogram at the loaction of the original statistic of nn,
  #             if FALSE, does not plot line
  # upper_line = if TRUE, will plot a green vertical line on the histogram at the loaction of the upper quantile from the JaB of of nn,
  #             if FALSE, does not plot line  
  # breaks = number of breaks in the histogram 
  
  ### Output
  # a histogram of the estimated null distribuion from the JaB of nn 
  
  #get correct row 
  r <- which(large_output$Results$Node_Number ==nn)
  #get correct data from Results 
  values <- large_output$Results[r, ]
  
  #get which bootstrap samples do not have nn 
  B <- length(large_output$Bootstrap_Nodes)
  keep_index <- rep(NA, B)
  for(b in 1:B){
    keep_index[b] <- !( nn %in% large_output$Bootstrap_Nodes[[b]])
  }
  
  orig_stats <- large_output$Results[, 1:2]%>%
    arrange(Node_Number)
  
  stats <- orig_stats$Orig_Stat[unlist(large_output$Bootstrap_Nodes[keep_index])]
  
  hist(stats, freq = F,
       main = paste0("Estimated Null Distribution of Node ", nn),
       breaks = breaks, 
       xlab = "Boostrap Statisitcs")
  if(orig_line){
    abline(v = large_output$Results$Orig_Stat[r], col = 2)
  }
  if(upper_line){
    abline(v = large_output$Results$Upper[r], col = 3)
  }
  
  
  
}






