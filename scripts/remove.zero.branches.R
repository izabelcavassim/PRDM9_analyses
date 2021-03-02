library(ape)

setwd("/Users/Zbaker/Downloads/")

### function for counting number of losses
dollo.num.losses = function(species.tree, tip.values) {
  # in the tree, edges describe the branches. 
  # Each node is listed in the first column twice 
  # each value in the second column corresponds to the more recent tips or nodes
  node.labels = unique(species.tree$edge[,1])
  node.labels[order(node.labels)]
  
  values = c(tip.values, rep(0,length(node.labels)))
  ancestors = lapply(which(tip.values == 1), FUN = function(x) {
    anc = species.tree$edge[which(species.tree$edge[,2] == x),1]
    done = FALSE
    while (!done) {
      new.anc = which(species.tree$edge[,2] == anc[length(anc)])
      if (length(new.anc) > 0) {
        anc = c(anc, species.tree$edge[which(species.tree$edge[,2] == anc[length(anc)]),1])
      } else {
        done = TRUE
      }
    }
    anc
  })
  values[unique(unlist(ancestors))] = 1
  
  ## identify common ancestors of all tips with a value of 1
  common.ancestors = names(table(unlist(ancestors))[which(table(unlist(ancestors)) == sum(tip.values))])
  anc.of.ancs = lapply(common.ancestors, FUN = function(x) {
    anc = species.tree$edge[which(species.tree$edge[,2] == x),1]
    if (length(anc) > 0) {
      done = FALSE
      while (!done) {
        new.anc = which(species.tree$edge[,2] == anc[length(anc)])
        if (length(new.anc) > 0) {
          anc = c(anc, species.tree$edge[which(species.tree$edge[,2] == anc[length(anc)]),1])
        } else {
          done = TRUE
        }
      }
    }
    anc
  })
  num.ancs = sapply(anc.of.ancs,length)
  mrca = common.ancestors[which(num.ancs == max(num.ancs))]
  values[as.numeric(common.ancestors)] = 0
  values[as.numeric(mrca)] = 1
  
  new.edges = NULL
  for (x in 1:length(species.tree$edge.length)) {
    new.edges = rbind(new.edges, 
                      values[species.tree$edge[x,]])
  }
  transitions = new.edges[,2] - new.edges[,1]
  
  sum(transitions == -1)
}

### function for looking at species contained within each loss
characterize.losses = function(tree, tip.values, mrca.sub = TRUE) {
  
  values = c(tip.values, rep(0,length(unique(tree$edge[,1]))))
  ancestors = lapply(which(tip.values == 1), FUN = function(x) {
    anc = tree$edge[which(tree$edge[,2] == x),1]
    done = FALSE
    while (!done) {
      new.anc = which(tree$edge[,2] == anc[length(anc)])
      if (length(new.anc) > 0) {
        anc = c(anc, tree$edge[which(tree$edge[,2] == anc[length(anc)]),1])
      } else {
        done = TRUE
      }
    }
    anc
  })
  values[unique(unlist(ancestors))] = 1
  if (mrca.sub) {
    ## identify common ancestors of all tips with a value of 1
    common.ancestors = names(table(unlist(ancestors))[which(table(unlist(ancestors)) == sum(tip.values))])
    anc.of.ancs = lapply(common.ancestors, FUN = function(x) {
      anc = tree$edge[which(tree$edge[,2] == x),1]
      if (length(anc) > 0) {
        done = FALSE
        while (!done) {
          new.anc = which(tree$edge[,2] == anc[length(anc)])
          if (length(new.anc) > 0) {
            anc = c(anc, tree$edge[which(tree$edge[,2] == anc[length(anc)]),1])
          } else {
            done = TRUE
          }
        }
      }
      anc
    })
    num.ancs = sapply(anc.of.ancs,length)
    mrca = common.ancestors[which(num.ancs == max(num.ancs))]
    values[as.numeric(common.ancestors)] = 0
    values[as.numeric(mrca)] = 1
  }
  
  new.edges = NULL
  for (x in 1:length(tree$edge.length)) {
    new.edges = rbind(new.edges, 
                      values[tree$edge[x,]])
  }
  transitions = new.edges[,2] - new.edges[,1]
  
  result = list()
  result[[1]] = which(transitions == -1)
  result[[2]] = sapply(which(transitions == -1), FUN = function(x) {
    descendents = tree$edge[x,2]
    incomplete = TRUE
    des.length = length(descendents)
    while (incomplete) {
      descendents = unique(c(descendents,tree$edge[which(tree$edge[,1] %in% descendents),2]))
      if (length(descendents) == des.length) {
        incomplete = FALSE
      } else {
        des.length = length(descendents)
      }
    }
    tree$tip.label[descendents[descendents %in% (1:length(tree$tip.label))]]
  })
  result[[3]] = sapply(which(transitions == -1), FUN = function(x) {
    descendents = tree$edge[x,2]
    incomplete = TRUE
    des.length = length(descendents)
    while (incomplete) {
      descendents = unique(c(descendents,tree$edge[which(tree$edge[,1] %in% descendents),2]))
      if (length(descendents) == des.length) {
        incomplete = FALSE
      } else {
        des.length = length(descendents)
      }
    }
    descendents
  })
  
  result
}

### reading in tree and PRDM9 calls
tree = read.tree("/Users/PM/Dropbox/PRDM9_evolution_paper/Simulations/data/species_tree_339_changed_replacements_to_old_names.nwk")
species.table = read.csv("/Users/PM/Dropbox/PRDM9_submission/Final_calls/PRDM9.species.table.csv",stringsAsFactors = FALSE)

# need to adjust Cebus capucinus to Cebus imitator
tree$tip.label[which(tree$tip.label == "Cebus_capucinus")] = "Cebus_imitator"

# removing 32 species due to busco 
busco.species = read.csv("/Users/zbaker/Dropbox/PRDM9_evolution_paper/Simulations/data/species_to_exclude_busco.txt", header = TRUE)
busco.species = sapply(busco.species[,2], FUN = function(x) {
  paste(strsplit(x,"_")[[1]][1:2], collapse = "_")
})
tree = drop.tip(tree, which(tree$tip.label %in% busco.species))

tree.species = as.character(sapply(tree$tip.label, FUN = function(x) {
  paste(strsplit(x,"_")[[1]],collapse=" ")
}))
PRDM9.calls = as.logical(sapply(tree.species, FUN = function(x) {
  species.table$Plan.A[which(species.table$species == x)]
}))

# checking the number of losses of PRDM9 observed in the tree initially
# to do so I temporarily discount the NA calls
PRDM9.calls[which(is.na(PRDM9.calls))] = TRUE
num.losses = dollo.num.losses(tree, PRDM9.calls)
# there are 8 losses (the 9 from before minus scleropages)

# setting what branch lengths to remove
rm.length = 0  ## this can be changed if we want to exclude very short branches that have non-zero values
zero.branches = which(tree$edge.length <= rm.length)

# determining which species to remove
descendants = sapply(zero.branches, FUN = function(x) {
  from = tree$edge[x,2] # the descendent node of the branch with length zero
  from = tree$edge[which(tree$edge[,1] == from),2] # the 2 descendents of that node
  from1 = from[1]
  from2 = from[2]
  
  complete = FALSE
  while(!complete) {
    len = length(from1)
    from1 = unique(c(from1, tree$edge[which(tree$edge[,1] %in% from1),2]))
    if (length(from1) == len) {
      complete = TRUE
    }
  }
  
  complete = FALSE
  while(!complete) {
    len = length(from2)
    from2 = unique(c(from2, tree$edge[which(tree$edge[,1] %in% from2),2]))
    if (length(from2) == len) {
      complete = TRUE
    }
  }
  
  from1 = from1[which(from1 %in% (1:length(tree$tip.label)))]
  from2 = from2[which(from2 %in% (1:length(tree$tip.label)))]
  
  ### the sister node of the zero branch
  
  anc = tree$edge[x,1]
  decs = tree$edge[which(tree$edge[,1] == anc),2]
  from3 = decs[decs != tree$edge[x,2]]
  
  complete = FALSE
  while(!complete) {
    len = length(from3)
    from3 = unique(c(from3, tree$edge[which(tree$edge[,1] %in% from3),2]))
    if (length(from3) == len) {
      complete = TRUE
    }
  }
  from3 = from3[which(from3 %in% (1:length(tree$tip.label)))]
  
  ### identify how many species are lost for each of three groups, and order them by how many species are lost (increasing)
  options = c(length(from1),length(from2),length(from3))
  option.order = order(options)
  
  ### figuring out how many transitions are removed per each of 3 possible groups of species to remove
  temp.tree = drop.tip(tree, from1)
  lost1 = num.losses - dollo.num.losses(temp.tree, PRDM9.calls[!(1:length(PRDM9.calls)) %in% from1])
  
  temp.tree = drop.tip(tree, from2)
  lost2 = num.losses - dollo.num.losses(temp.tree, PRDM9.calls[!(1:length(PRDM9.calls)) %in% from2])
  
  temp.tree = drop.tip(tree, from3)
  lost3 = num.losses - dollo.num.losses(temp.tree, PRDM9.calls[!(1:length(PRDM9.calls)) %in% from3])
  
  ### which options result in the minimum number of lost transitions
  best.options = which(c(lost1,lost2,lost3) == min(c(lost1,lost2,lost3)))
  
  ### choosing the group that results in the minimum number of lost transitions, and among those, with the fewest species
  if (option.order[1] %in% best.options) {
    option = option.order[1]
  } else {
    if (option.order[2] %in% best.options) {
      option = option.order[2]
    } else {
      option = option.order[3]
    }
  }
  
  if (option == 1) {
    result = from1
  } else {
    if (option == 2) {
      result = from2
    } else {
      result = from3
    }
  }
  
  result
})
rm.species = unique(unlist(descendants))

# looking at which species get removed
new.tree = drop.tip(tree, rm.species)
new.calls = PRDM9.calls[!(1:length(PRDM9.calls)) %in% rm.species]
num.losses = dollo.num.losses(new.tree, new.calls)

########### there is still 1 zero branch length (this issue might be explaind by 2 adjacent branches with length of zero?)

zero.branches = which(new.tree$edge.length <= rm.length)
descendants = sapply(zero.branches, FUN = function(x) {
  from = new.tree$edge[x,2] # the descendent node of the branch with length zero
  from = new.tree$edge[which(new.tree$edge[,1] == from),2] # the 2 descendents of that node
  from1 = from[1]
  from2 = from[2]
  
  complete = FALSE
  while(!complete) {
    len = length(from1)
    from1 = unique(c(from1, new.tree$edge[which(new.tree$edge[,1] %in% from1),2]))
    if (length(from1) == len) {
      complete = TRUE
    }
  }
  
  complete = FALSE
  while(!complete) {
    len = length(from2)
    from2 = unique(c(from2, new.tree$edge[which(new.tree$edge[,1] %in% from2),2]))
    if (length(from2) == len) {
      complete = TRUE
    }
  }
  
  from1 = from1[which(from1 %in% (1:length(new.tree$tip.label)))]
  from2 = from2[which(from2 %in% (1:length(new.tree$tip.label)))]
  
  ### the sister node of the zero branch
  
  anc = new.tree$edge[x,1]
  decs = new.tree$edge[which(new.tree$edge[,1] == anc),2]
  from3 = decs[decs != new.tree$edge[x,2]]
  
  complete = FALSE
  while(!complete) {
    len = length(from3)
    from3 = unique(c(from3, new.tree$edge[which(new.tree$edge[,1] %in% from3),2]))
    if (length(from3) == len) {
      complete = TRUE
    }
  }
  from3 = from3[which(from3 %in% (1:length(new.tree$tip.label)))]
  
  ### identify how many species are lost for each of three groups, and order them by how many species are lost (increasing)
  options = c(length(from1),length(from2),length(from3))
  option.order = order(options)
  
  ### figuring out how many transitions are removed per each of 3 possible groups of species to remove
  temp.new.tree = drop.tip(new.tree, from1)
  lost1 = num.losses - dollo.num.losses(temp.new.tree, new.calls[!(1:length(new.calls)) %in% from1])
  
  temp.new.tree = drop.tip(new.tree, from2)
  lost2 = num.losses - dollo.num.losses(temp.new.tree, new.calls[!(1:length(new.calls)) %in% from2])
  
  temp.new.tree = drop.tip(new.tree, from3)
  lost3 = num.losses - dollo.num.losses(temp.new.tree, new.calls[!(1:length(new.calls)) %in% from3])
  
  ### which options result in the minimum number of lost transitions
  best.options = which(c(lost1,lost2,lost3) == min(c(lost1,lost2,lost3)))
  
  ### choosing the group that results in the minimum number of lost transitions, and among those, with the fewest species
  if (option.order[1] %in% best.options) {
    option = option.order[1]
  } else {
    if (option.order[2] %in% best.options) {
      option = option.order[2]
    } else {
      option = option.order[3]
    }
  }
  
  if (option == 1) {
    result = from1
  } else {
    if (option == 2) {
      result = from2
    } else {
      result = from3
    }
  }
  
  result
})

rm.species = as.numeric(descendants)

new.new.tree = drop.tip(new.tree, rm.species)
new.new.calls = new.calls[!(1:length(new.calls)) %in% rm.species]

new.num.losses = dollo.num.losses(new.new.tree,new.new.calls)
PRDM9.losses = characterize.losses(new.new.tree,new.new.calls)

## this will show you the species lost per loss of PRDM9
PRDM9.losses[[2]]

## Define which species to keep from those PRDM9 losses (3 species per loss when possible)
## We will pick species based on their genome assembly quality


contigs_stat = read.csv( "/Users/PM/Desktop/339_species/stats/stats_assemblies.csv", sep=",")

picking_species_based_genome_assembly <- function(PRDM9_loss){
  rownames(contigs_stat) <- contigs_stat$Species
  sub = contigs_stat[PRDM9_loss,]
  sub %>%
  group_by(contig_N50, contigs) %>%
    arrange(contigs,desc(contig_N50)) -> ordered
  return(ordered)
}

for(i in 1:8){
  loss = as.vector(PRDM9.losses[[2]][i])
  results = picking_species_based_genome_assembly(loss[[1]])
  print(head(results))
}

## in the end, removed 29 species to remove zero branch lengths
write.tree(new.new.tree,"modified.tree")

tree.species = as.character(sapply(new.new.tree$tip.label, FUN = function(x) {
  paste(strsplit(x,"_")[[1]],collapse=" ")
}))
final.PRDM9.calls.A = sapply(tree.species, FUN = function(x) {
  species.table$Plan.A[which(species.table$species == x)]
})
final.PRDM9.calls.B = sapply(tree.species, FUN = function(x) {
  species.table$Plan.B[which(species.table$species == x)]
})

write.table(final.PRDM9.calls.A, "modified.calls.plan.A", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(final.PRDM9.calls.B, "modified.calls.plan.B", quote = FALSE, row.names = FALSE, col.names = FALSE)

