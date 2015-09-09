### 1. Installation

Start R and then type:

```{R}
install.packages("devtools")
devtools::install_github("hemberg-lab/SC3")
```

### 2. Test run

To test that the package has been installed successfully please run SC3 on a [published dataset](http://www.nature.com/nature/journal/v509/n7500/full/nature13173.html):

```{R}
library(SC3)
run_sc3("quake_all_fpkm", 3:7)
```

It should open SC3 in a browser window without providing any error. If there is any error please send it to [Vladimir Kiselev](mailto:vk6@sanger.ac.uk).

### 3. Input format

To run SC3 on your own data one need to prepare an input file with an expression matrix. The expression matrix should be of the following format:

cell1 cell2 cell3 cell4 cell5  
gene1 1 2 3 4 5  
gene2 1 2 3 4 5  
gene3 1 2 3 4 5  

The first row of the expression matrix should contain one fewer field than all other rows and separators should be either spaces or tabs. If separators are commas (,) and the file is in the proper csv format, then the extension of the file should be .csv.

### 4. Running SC3

```{R}
library(SC3)
run_sc3(filename, k.min:k.max)
```

where __filename__ is the path to your input file, __k.min__ is the minimum number of clusters, __k.max__ is the maximum number of clusters. For example, if you would like to check clustering of the cells in "expression-matrix.txt" file, where __k__ is in the region from 2 to 5, then you need to run the following:

```{R}
run_sc3("expression-matrix.txt", 2:5)
```
