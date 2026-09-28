# Big Five Intelligence Genetic Overlap

## Analysis code

The [`scripts`](scripts/) folder contains all code used in our main and supplementary analyses. The scripts are separated into folders by analysis, and each one contains a `README` file which explains briefly what each script does.

## Data

[`data`](data/) contains gene sets, gene coordinates, and the excluded genes lists, which are necessary to run any script which uses a gene-level block jackknife.

Most data necessary to run our analyses are not included in this file, specifically the GWAS summary statistics for the Big Five and intelligence, and all of the gene sets we used that we did not create ourselves.

### GWAS summary statistics

- **Big Five summary statistics** can be obtained from <https://osf.io/hgnsm/files/osfstorage>.
- **Intelligence summary statistics** can be obtained from <https://vu.data.surf.nl/index.php/s/9tgwxmO5yosQkmb>.

### Gene sets

Gene set sources are described in the MAGMA part of our methods section.

All but one of the collections we did not generate can be downloaded from <https://www.gsea-msigdb.org/gsea/msigdb>.

The other comes from <https://github.com/WonuAkingbuwa/Gene_expression_SCZ/blob/main/data/final_genesets.xlsx>.
