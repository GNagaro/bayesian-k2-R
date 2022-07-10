library(tidyverse)

log.f <- function(node, parents, dataset) {

    node.nunique <- dataset[node] |>
                        unique() |>
                        nrow()

    contingency_table <- dataset[c(node, parents)] |>
                            table() |>
                            as.data.frame() |>
                            as_tibble() |>
                            filter(Freq > 0)

    sum_over_x <- 1:(node.nunique - 1) |>
                    log() |>
                    sum()

    sum_over_y <- contingency_table |>
                    rowwise() |>
                    mutate(Freq = ((1:Freq) |> log() |> sum())) |>
                    ungroup() |>
                    group_by_at(parents) |>
                    summarise(Freq = sum(Freq), .groups='drop_last') |>
                    ungroup() |>
                    select(Freq) |>
                    deframe()

    sum_over_z <- contingency_table |>
                    group_by_at(parents) |>
                    summarise(Freq = sum(Freq), .groups='drop_last') |>
                    ungroup() |>
                    rowwise() |>
                    mutate(Freq = ((1:(Freq + node.nunique - 1)) |> log() |> sum())) |>
                    ungroup() |>
                    select(Freq) |>
                    deframe()

    dataset.log.prob <- sum(sum_over_x + sum_over_y - sum_over_z)
    
    return(dataset.log.prob)

}


k2 <- function(dataset, parents.nmax, f = log.f) {
        
        
        #####################
        ##  PREPROCESSING  ##
        #####################
        # Since bnstruct uses custom objects, we do some stuff to use
        #  our function in this environment.
        nodes <- names(dataset)
        n <- length(nodes)
        
        # creating a named vector to select the correct index in dag matrix
        candidate.idx <- 1:n
        names(candidate.idx) <- nodes
        
        # allocating an empty adjacency matrix of the network  (see AllClasses.R)
        adjm <- matrix(0, n, n)
        
        
        
        ####################
        ##  K2 ALGORITHM  ##
        ####################
        
        tic <- Sys.time()

        #dag <- empty.graph(nodes=nodes)
        for (i in 2:length(nodes)) {
            
            # select current node & all the previous
            node <- nodes[i]
            node.previous <- nodes[1:i-1]
            
            parents <- c()
            P_old <- f(node, parents, dataset)

            proceed <- ifelse(length(node.previous) > 0, T, F)
            
            while (proceed & (length(parents) < parents.nmax)) {

                candidates <- setdiff(node.previous, parents)

                P_new <- P_old

                for (candidate in candidates) {
                    score <- f(node, c(parents, candidate), dataset)
                    if (score > P_new) {
                        candidates.best <- candidate
                        P_new <- score
                    }
                }

                if (P_new > P_old) {
                    P_old <- P_new
                    parents <- c(parents, candidates.best)
                    #dag <- set.arc(dag, from=candidates.best, to=node)
                    adjm[ candidate.idx[candidates.best], i] <- 1 # set the arc in matrix
                } else {
                    proceed = F
                }

            }

            #cat('Node:', node, '| Parents:', parents, '\n')

        }

        toc <- Sys.time()

        cat('Execution time:', toc - tic, 's')

        return(adjm)
}
