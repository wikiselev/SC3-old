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
sc3(quake_all_fpkm, ks = 3:7, cell.filter = TRUE)
```

It should open SC3 in a browser window without providing any error. If there is any error please send it to [Vladimir Kiselev](mailto:vk6@sanger.ac.uk).

### 3. Running SC3

```{R}
sc3(dataset, ks = k.min:k.max, cell.filter = FALSE, interactivity = TRUE, svm.num.cells = 1000, cell.filter.genes = 2000)
```

* __dataset__ is either an R matrix / data.frame / data.table object OR a path to your input file containing an expression matrix
* __ks__ is a range of the number of clusters that needs to be tested. __k.min__ is the minimum number of clusters (default is 3). __k.max__ is the maximum number of clusters (default is 7).
* (optional) __cell.filter__ is used to filter cells that express less than __cell.filter.genes__ genes (_lowly expressed cells_). By default it is OFF. To switch it ON please use __TRUE__ value as in the __Test run__ above. Should be used if it is not possible to properly cluster original cells - filtering of _lowly expressed cells_ usually makes clustering better.
* (optional) __interactivity__ defines whether a browser interactive window should be open after all computation is done. By default it is ON. To switch it OFF please use __FALSE__ value. This option can be used to separate clustering calculations from visualisation, e.g. long and time-consuming clustering of really big datasets can be run on a farm cluster and visualisations can be done a personal laptop afterwards. If __interactivity__ is OFF then all clustering results will be saved to __dataset__.rds file. To run interactive visulisation with the precomputed clustering results please use `sc3_interactive(readRDS("dataset.rds"))`.
* (optional) __svm.num.cells__ - if number of cells in your dataset is more than this parameter, then an SVM prediction will be used. The default is 1000.
* (optional) __cell.filter.genes__ - if __cell.filter__ is used that this parameter defines the minimum number of genes that have to be expressed in each cell. If there are less, the cell will be removed from the analysis. The default is 2000.

Usage example: if you would like to check clustering of your __dataset__ for __ks__ from 2 to 5, then you need to run the following:

```{R}
sc3(dataset, ks = 2:5)                        # without filtering of lowly expressed cells
sc3(dataset, ks = 2:5, cell.filter = TRUE)    # with filtering of lowly expressed cells
sc3(dataset, ks = 2:5, interactivity = FALSE) # without interactive visualisation
```

### 4. Input file format

To run SC3 on an input file containing an expression matrix one need to preprocess the input file so that it looks as follows:


|  | cell1 | cell2 | cell3 | cell4 | cell5 
--- | --- | --- | --- | --- | ---
| __gene1__ | 1 | 2 | 3 | 4 | 5 
| __gene2__ | 1 | 2 | 3 | 4 | 5 
| __gene3__ | 1 | 2 | 3 | 4 | 5 


The first row of the expression matrix (with cell labels, e.g. __cell1__, __cell2__, etc.) should contain one fewer field than all other rows. Separators should be either spaces or tabs. If separators are commas (,) then the extension of the file must be .csv. If a path to your input file is "/path/to/input/file/expression-matrix.txt", to run it:

```{R}
sc3("/path/to/input/file/expression-matrix.txt", ks = 2:5)
```
