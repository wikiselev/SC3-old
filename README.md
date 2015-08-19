### Installation

```{R}
install.packages("devtools")
devtools::install_github("wikiselev/SC3")
```

### Input format

The input expression matrix should be of the following format:

gene1 1 2 3 4 5  
gene2 1 2 3 4 5  
gene3 1 2 3 4 5  

It should not have a header and separators should be either spaces or tabs.

### Running CS3

```{R}
library(CS3)
run_cs3(filename, k.min:k.max)
```

where filename is the path to your expression matrix, k.min is the minimum number of clusters, k.max is the maximum number of clusters.
