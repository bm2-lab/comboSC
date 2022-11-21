# comboSC
**comboSC - Optimized therapy combination for pan-cancer based scRNA-Seq**

![项目总流程图.jpg](https://github.com/EularTang/Figure_bed/blob/master/Project_workflow5.jpeg)
#### Description

1.comboSC is a scalable toolkit for personalized combination therapy recommendation based on the single-cell sequencing data of cancer patient.

2.comboSC collects more than 30,000 small molecule/drugs, which come from cMAP and GDSC databases.

3.comboSC designed an efficient immune score for personalized immunity profile evaluation based on the single-cell sequencing data of cancer patient.
- For high immune score samples, comboSC recommends to use routine immunotherapy like immune checkpoint inhibitors.
- For middle immune score samples, comboSC recommends to use combination therapy by combining immunotherapy with certaixn small molecule/drugs, which are predicted to regulate the immune microenvironment and boost the immunotherapy effect.
- For low immune score samples, comboSC recommends to use combination therapy by combining small molecule/drugs to eliminate malignant cells directly.

4.The current version of comboSC is tested for the following 10 cancer types:basal cell carcinoma (BCC), breast invasive carcinoma (BRCA), colorectal cancer (CRC), head and neck cancer (HNSC), melanoma, non-small-cell lung cancer (NSCLC), pancreatic adenocarcinoma (PAAD), skin cutaneous melanoma (SKCM), liver hepatocellular carcinoma (LIHC), and adult acute myeloid leukemia (AML).
## Accessiable
### Web
The comboSC tools can be access from the webserver (www.comboSC.top).
### Github
Newly comboSC can be installed from GitHub. It required R, version 3.6.0 or greater. Rstudio is also recommended.The required R packages are listed in `resources/Allpackage.R`.
### Docker
If you are having trouble installing the package on your system, you can build a docker instance that can be used on a wide range of systems and cloud environments.

`docker pull eularg/combsc:lates`

## Test Data
|Cancer type|Number of samples|Reference|
|:---:|:---:|:---:|
|Melanoma|1|Nat Commun . 2018|
|NSCLC|10|Cancer Med . 2019|
|NSCLC|5|Nature . 2020|
|NSCLC|3|bioRxiv . 2019|
|SKCM|10|Cell . 2017|
|SKCM|5|Nat Med . 2019|
|LIHC|2|Cancer Cell . 2019|
|HNSC|13 |Cell . 2017|
|BRCA|1|Nat Genet . 2021|
|BRCA|1|biorxiv . 2019|
|CRC|9|Nat Med . 2019|
|BCC|9 |Nat Med . 2019|
|PAAD|17 |Cell Res . 2019|
|AML|11 |Nat Genet . 2021|


## Usage
### Input

The model required three input files:

1. **expression_matrix**, a expression profile of the samples where rows are genes, columns are cells, and the value is counts or TPM. 
2. **cell_metadata**, a data frame, where rows are cells, and columns are cell attributes (such as cell type, detail cluster, tissues, samples, batch, etc.)
3. **gene_metadata**, an data frame, where rows are gene id (e.g. genes), and columns are gene attributes( such as biotype, gc content, whether key gene, etc.)

![input.png](http://www.combosc.top/combsc/static/images/metedata.png)

### Output
The result of the model will be a dataframe that contains the highest-scoring drug combination and their detail information.  

|Drug combination|Personalized score level|Score value|
|:---:|:---:|:---:|
|Sepantronium bromide & AZD6738|Low|31.86558|
|Tretinoin & rTRAIL|Low|23.7361276|
|Sepantronium bromide & BIX02189|Low|23.5998664|


## Running comboSC
### Command in terminal
```
Rscript new_sctc.R [exp] [meta] [sim] [Patid]
```
- exp: expression_matrix
- meta:cell metadate in patients
- sim :Threshold of similarity between query cell and reference expression profile.
- patid:If the input contains more than one patient, we need to specify the id of the patient to be calculated.
### Metadata specifications
The cell metadata is a file in csv or csv.gz format, where the `"cell.id"` must be in the first column, corresponding to the column name of the expression profile. There are other columns to describe the cell information. The columns that must be included ar: 1, `"patient"`, the patient number where the cell is located, if the sample has only one patient, the same value is sufficient. 2, `"cancertype"`, the cancer type of the sample, must be one of the ten cancer types mentioned in the description, including, `"BCC", " BRCA", "CRC", "HNSC", "SKCM", "NSCLC", "PAAD", "SKCM", "LIHC", "AML"`. 3, `"treatment"`, the time period of sampling, including, `"Pre ", "Post", "Duration"`

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
|CAFs|cancer-associated fibroblastsMalignant|Malignant cells|
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
doi:
## Contributing
Bug reports, pull requests and other contributions are welcomed!
