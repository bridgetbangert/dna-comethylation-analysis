# DNA Co-Methylation Analysis in Breast and Endometrial Cancer

## Overview
This project investigates DNA co-methylation patterns on the X chromosome in breast cancer (BRCA) and endometrial cancer (UCEC) samples. DNA methylation (DNAm) is an epigenetic mechanism involving the addition of a methyl group to cytosine at CpG sites, which can impact gene expression. Previous studies have highlighted unique methylation patterns on the X chromosome, suggesting a potential role in cancer progression. Our study builds upon this research by analyzing within-sample (WS) co-methylation patterns in matched normal and tumor samples from both alive and deceased individuals.

## Motivation
Epigenetic changes, including DNAm, play a critical role in cancer development and progression. Prior research has identified unique co-methylation characteristics on the X chromosome, particularly in relation to breast and endometrial cancers. By further investigating these patterns at the CpG site level, we aim to uncover deeper insights into cancer-specific methylation disruptions and their potential implications for detection, prevention, and treatment.

## Data
- The dataset was provided by Dr. Shuying Sun and is the same dataset used in previous research, filtered to focus exclusively on the X chromosome.
- Data includes matched normal and tumor samples for breast cancer (BRCA) and endometrial cancer (UCEC), stratified by alive and deceased patient status.
- The dataset originates from Illumina 450K methylation array data, which provides high-resolution DNAm profiling.

## Methods
1. **Co-Methylation Analysis:**
   - Identification of highly correlated CpG sites (|r| ≥ 0.8) within each sample.
   - Distinction between highly positively correlated (r ≥ 0.8) and highly negatively correlated (r ≤ -0.8) sites.
2. **Distance Analysis:**
   - Measurement of genomic distances between highly correlated CpG sites.
   - Comparison of distance distributions between normal and tumor samples.
3. **Statistical Evaluation:**
   - Quantification and comparison of co-methylation trends across different sample groups.
   - Assessment of co-methylation differences between BRCA and UCEC.

## Results
- **Distinct Co-Methylation Patterns:** BRCA and UCEC exhibit different co-methylation behaviors.
- **Higher Co-Methylation in UCEC Tumor Samples:** UCEC tumors had more highly correlated CpG sites than normal samples, whereas this pattern was not observed in BRCA.
- **Prevalence of Positive Correlations:** The majority of high correlations were positive, except in UCEC deceased samples, where negative correlations were also present.
- **Distance Reduction in Tumor Samples:** Highly correlated CpG sites were closer together in tumor samples compared to normal samples. This effect was more pronounced in BRCA than in UCEC.

## Key Findings
- Breast and endometrial cancer have distinct co-methylation patterns, emphasizing the potential for cancer-specific epigenetic markers.
- Tumor progression appears to be associated with a disruption of normal co-methylation patterns, particularly in CpG site clustering.
- The X chromosome displays unique co-methylation trends that warrant further exploration in cancer research.

## Tools & Technologies
- **Programming Language:** R

## Future Directions
- Extend analyses to autosomal chromosomes, particularly chromosome 22, to compare patterns of negative correlation.
- Investigate potential biomarkers for cancer detection based on co-methylation differences.
- Apply machine learning techniques to classify cancer subtypes based on co-methylation profiles.

## References
- Sun, S., Dammann, J., Lai, P., Tian, C. (2022). Thorough statistical analyses of breast cancer co-methylation patterns. [DOI: 10.1186/s12863-022-01046-w](https://doi.org/10.1186/s12863-022-01046-w)
- Zhang, Y., Huang, S. (2017). Co-methylation and its association with cancer.

## Contributors
- Dr. Shuying Sun
- Bridget Bangert
- Madison Glenwinkel
