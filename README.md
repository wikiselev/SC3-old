### Installation

```{R}
install.packages("devtools")
devtools::install_github("hemberg-lab/SC3")
```

### Input format

The input expression matrix should be of the following format:

gene1 1 2 3 4 5  
gene2 1 2 3 4 5  
gene3 1 2 3 4 5  

It should not have a header and separators should be either spaces or tabs. If separators are commas (,) and the file is in the proper csv format, then the extension of the file should be .csv

### Running SC3

```{R}
library(SC3)
run_sc3(filename, k.min:k.max)
```

where __filename__ is the path to your expression matrix, __k.min__ is the minimum number of clusters, __k.max__ is the maximum number of clusters. For example, if you would like to check clustering of the cells in "expression-matrix.txt" file, where __k__ is in the region from 2 to 5, then you need to run the following:

```{R}
run_sc3("expression-matrix.txt", 2:5)
```
