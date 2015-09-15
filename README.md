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
run_sc3(quake_all_fpkm, 3:7)
```

It should open SC3 in a browser window without providing any error. If there is any error please send it to [Vladimir Kiselev](mailto:vk6@sanger.ac.uk).

### 3. Running SC3

```{R}
library(SC3)
run_sc3(dataset, k.min:k.max)
```

where __dataset__ is either an R matrix / data.frame / data.table object OR a path to your input file containing an expression matrix, __k.min__ is the minimum number of clusters, __k.max__ is the maximum number of clusters. For example, if you would like to check clustering of your __dataset__ for __k__ from 2 to 5, then you need to run the following:

```{R}
run_sc3(dataset, 2:5)
```

### 4. Input file format

To run SC3 on an input file containing an expression matrix one need to preprocess the input file so that it looks as follows:

cell1 cell2 cell3 cell4 cell5  
gene1 1 2 3 4 5  
gene2 1 2 3 4 5  
gene3 1 2 3 4 5  

The first row of the expression matrix (with cell labels, e.g. cell1, cell2, etc.) should contain one fewer field than all other rows. Separators should be either spaces or tabs. If separators are commas (,) then the extension of the file must be .csv. If a path to your input file is "/path/to/input/file/expression-matrix.txt", to run it:

```{R}
run_sc3("/path/to/input/file/expression-matrix.txt", 2:5)
```
