# comboSC
**comboSC - Personalized tumor combination therapy prediction based on single cell RNA-seq**

![Fig_workflow](https://user-images.githubusercontent.com/37855187/221109476-43db2314-85f1-474e-af96-719f61978a38.jpg)
.jpg)
#### Description

1.comboSC is a scalable toolkit for personalized combination therapy recommendation based on the single-cell sequencing data of cancer patient.

2.comboSC uses a large collection of molecule/drugs, which come from CMAP and GDSC databases.

3.comboSC designs an efficient immune score for personalized immunity profile evaluation based on the single-cell sequencing data of cancer patient.
- For high immune score samples, comboSC recommends to use routine immunotherapy like immune checkpoint inhibitors.
- For middle immune score samples, comboSC recommends to use combination therapy by combining immunotherapy with certaixn small molecule/drugs, which are predicted to regulate the immune microenvironment and boost the immunotherapy effect.
- For low immune score samples, comboSC recommends to use combination therapy by combining small molecule/drugs to eliminate malignant cells directly.

4.The current version of comboSC have been applied to the following 10 cancer types,including: basal cell carcinoma (BCC), breast invasive carcinoma (BRCA), colorectal cancer (CRC), head and neck cancer (HNSC), non-small-cell lung cancer (NSCLC), pancreatic adenocarcinoma (PAAD), skin cutaneous melanoma (SKCM), liver hepatocellular carcinoma (LIHC), and adult acute myeloid leukemia (AML),Endometrioid Carcinoma（UCEC）.

|Cancer |Number of samples|PMID of reference|
|:---:|:---:|:---:|
|BCC|6 |31359002|
|BRCA|1 |bioRxiv|
|BRCA |1 |34493872|
|CRC|7 |31359002|
|HNSC|11 |29198524|
|NSCLC|4 |31033233|
|NSCLC |3 |32103181|
|NSCLC |5 |31979687|
|NSCLC |3 |bioRxiv|
|PAAD|23 |31273297|
|SKCM|1 |29198524|
|SKCM |5 |3038455|
|SKCM |1 |30250229|
|SKCM |6 |27124452|
|LIHC|1 |32103181|
|LIHC |4 |31588021|
|AML|1 |34493872|
|UCEC|2 |32103181|

## Accessible
### Web
The comboSC tools can be accessed from the webserver (www.comboSC.top).
### Github
Local comboSC can be installed from https://github.com/bm2-lab/comboSC. It required R(v3.6.0 or updated),Python(v3.6 or updated), and [Tres](https://github.com/data2intelligence/Tres).The source code is available on this page,  the full and executable comboSC with dependency data can be downloaded from this [link](http://www.combosc.top/combsc/csv/example?name=comboSC.zip). The required R packages are listed in `resources/Allpackage.R`.
## Usage
### Input

The model required two files:

1. **expression_matrix**, a expression profile of the samples where rows are genes, columns are cells, and the value is counts or TPM. 
2. **cell_metadata**, a data frame, where rows are cells, and columns are cell attributes (such as cell type, detail cluster, tissues, samples, batch, etc.)

![input.png](http://www.combosc.top/combsc/static/images/metedata.png)

### Output
The result of the model will be a dataframe that contains the highest-scoring drug combination and their detail information.  

|Drug combination|Personalized score level|Score value|
|:---:|:---:|:---:|
|Sepantronium bromide & AZD6738|Low|31.86558|
|Tretinoin & rTRAIL|Low|23.7361276|
|Sepantronium bromide & BIX02189|Low|23.5998664|


## Running comboSC
Please CD into the comboSC folder in the terminal, such as `cd /home/user/comboSC/`.
### Command in terminal

```
Rscript comboSC.R [exp] [meta] [sim] [patid]
```
- exp: Expression_matrix
- meta:Cell metadate in patients
- sim :Threshold of similarity between query cell and reference expression profile.
- patid:If the input contains more than one patient, we need to specify the id of the patient to be calculated.

### Example
```
Rscript ./comboSC.R "bcc006.exp.gz" "bcc006.meta.gz" 0.5 282
```
The output of the model is in `./comboSC/Auxiliary`. The test data is availabe in [example data](http://www.combosc.top/combsc/csv/example?name=example.zip).

### Metadata specifications
The cell metadata is a file in csv or csv.gz format, where the `"cell.id"` must be in the first column, corresponding to the column name of the gene expression profile. There are other columns to describe the cell information. The columns that must be included are: 

1, `"patient"`, the patient id where the cell is located, such as `"Su006", "Lung01"`.if the sample has only one patient, the same value is sufficient. 

2, `"cancertype"`, the cancer type of the sample, must be one of the ten cancer types mentioned in the description, including, `"BCC", " BRCA", "CRC", "HNSC", "SKCM", "NSCLC", "PAAD", "UCEC", "LIHC", "AML"`. 

3, `"treatment"`, the time period of sampling, including, `"Pre ", "Post", "Duration"`

Lastly, due to the limitations of existing automated annotation tools in identifying cell types with similarity features, we recommend to use manual cell annotations in the metadata. The annotated cell types need to be stated in the `"cluster"` column. comboSC can recognize the following cell type names:
|Input cell type names|Full name|
|:---:|:---:|
|B|B cells|
|CD4_T_cells|CD4+ T cells|
|CD8_mem_T_cells|CD8+ memory T cells|
|CD8_ex_T_cells|CD8+ exhausted T cells|
|CD8_act_T_cells|CD8+ activated T cells|
|Endothelial|endothelial cells|
|Fibroblasts|fibroblasts|
|CAFs|cancer-associated fibroblasts|Malignant cells|
|Plasma| plasma cells|
|Tprolif| proliferating T cells|
|DC|dendritic cells|
|Treg|Regulatory T Cells|
|M1|M1 macrophages|
|M2|M2 macrophages|
|Mast|Mast Cells|
|Myofibroblasts|Myofibroblasts|
|NK_cells|Natural Killer Cells|

## Citation  
Personalized tumor combination therapy prediction based on single cell RNA-seq(submitted)
## Contributing
Bug reports, pull requests and other contributions are welcomed!
