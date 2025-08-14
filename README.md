# Neural Organoid Hormonal Atlas (NOHA)

Repository for the [Neural Organoid Hormonal Atlas](biorxiv doi).

> A molecular cell atlas of endocrine signalling in human neural organoids 

![ExperimentalDesign](FigSchemes/ExperimentalDesign.jpg)

![Hallmarks](FigSchemes/Hallmarks.jpg)


## Abstract

The development of the human brain is profoundly shaped by hormonal signalling, and endocrine disruption is implicated in various neuropsychiatric conditions. However, a comprehensive understanding of how hormonal pathways orchestrate human neurodevelopment remains elusive. Here, we generated a multi-scale atlas of endocrine signalling in human neural organoids through systematic perturbations with agonists and inhibitors of seven key hormonal pathways: androgen (AND), estrogen (EST), glucocorticoid (GC), thyroid (THY), retinoic acid (RA), liver X (LX), and aryl hydrocarbon (AH). By integrating bulk and single-cell transcriptomics, high-throughput imaging, and targeted steroidomics, we mapped the molecular and cellular consequences of these perturbations. Our analyses revealed that retinoic acid exerted the most profound effects, promoting neuronal differentiation and maturation, consistent with its established role as a patterning factor. We confirmed previously reported effects, such as induction of mTOR signalling by AND activation, alteration of disease relevant genes by GC, and enhanced differentiation by TH. Furthermore, we observed that LX activation upregulated genes involved in cholesterol metabolism, while AH inhibition promoted neuronal differentiation. We uncovered extensive crosstalk between these endocrine pathways. For instance, we identified significant convergence in gene expression changes induced by AND agonist and inhibitors of GC, TH, and LX, particularly affecting genes related to protein folding, and metabolic regulation, as also highlighted by weighted gene co-expression network analysis. Single-cell analyses pinpointed cell-type-specific responses to hormonal challenges, such as the caudalisation of progenitors and neurons upon RA activation, and the depletion of a developmental neuronal state upon AH activation. Finally, we dissected the cytoarchitectural and morphometric impact of hormonal perturbations and demonstrated that neural organoids possess active steroidogenic pathways that are functionally modulated by the tested compounds. This foundational resource provides a systematic quantification of hormonal impact on human neurodevelopment, representing a framework for investigating unknown aspects of the developmental origins of neuropsychiatric traits. This will be key to better understand how environmental factors, especially endocrine-disrupting chemicals, and genetic endocrine disorders, contribute to neurodevelopmental pathogenesis, as well as for training advanced generative models to improve their overall predictive power.


## Report Cards and Mining

Acknowledging the wealth and multi-scale nature of the data produced in this study, which integrates bulk and single-cell transcriptomics, high-throughput imaging, and targeted steroidomics, we have developed a public resource to ensure maximum accessibility and utility for the scientific community. 

For each of the seven hormonal pathways investigated, we provide a dedicated "report card" containing a schematic summary of the main observations drawn by our analyses. This allows fellow researchers to gain a rapid, high-level understanding of the key findings.

A more comprehensive asset complements this initial overview: a mining of the data organized in this public repository serving as a portal for deep data exploration, offering detailed descriptions of our results and, crucially, contextualizing them within the current body of scientific knowledge. This transforms our atlas from a static publication into a dynamic resource designed to be actively used, explored, and updated. While we have meticulously curated this initial release, we envision it as a living atlas that will evolve and improve in the future. This was elaborated as a collective effort distributed across scientists with different backgrounds and seniority. Consequently, even if the same structure with thematic sub-chapters was elaborated for each hormonal pathway, the level of details included for each pathway is still heterogeneous, mainly because the amount of available public literature is very different across pathways.  We plan to harmonise better the info included in the near future, and we hope the scientific community will contribute to it. 
Our vision is that this resource will foster further discoveries and provide a framework for investigating the impact of endocrine signalling in human neurodevelopment, and more broadly, the developmental origins of neuropsychiatric traits.

If you notice any typos or glitches, feel welcome to contact us to contribute and improve this resource.

The work is accessible [here](1_bulkRNASeq/7_DataMining/0.Index.md).



## Containers

- for the bulkRNAseq and steroidomics analyses it can be retrieved via `docker pull testalab/downstream:EndPoints-1.1.5`.
- for the scRNAseq analyses it can be retrieved via `docker pull alessiavalenti/sc:sc-noha-0.0.1`.
- for the imaging analyses it can be retrieved via `docker pull alessiavalenti/imaging:TMA-0.0.1`.


## html notebooks

An html version of the notebooks is accessible [here](link to the github.io page).


## Code folders

* `1_bulkRNASeq`: codes and notebooks related to the bulk transcriptomic analysis.

* `2_scRNASeq`: codes and notebooks related to the single-cell transcriptomic analysis.

* `3_ImageAnalysis`: codes and notebooks related to the image analysis.

* `4_steroidomics`: codes and notebooks related to the targeted steroidomics analysis.


