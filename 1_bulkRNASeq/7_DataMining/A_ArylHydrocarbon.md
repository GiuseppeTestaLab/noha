## Table of Content

* [Compounds tested](#Compounds-tested)
* [Results from reference compounds](#Results-from-reference-compounds)
* [General pathway information](#General-pathway-information)
* [Physiological effect of the hormone](#Physiological-effects-of-AHR)
* [Impact on neurodevelopment](#Impact-on-neurodevelopment)
* [In vitro experiment with iPSC and neural models](#In-vitro-experiments-with-iPSC-and-neural-models)
* [Epidemiological data](#Epidemiological-data)
* [Evidence of EDC disruption of the pathway](#Evidence-of-EDC-disruption-of-the-pathway)
* [Crosstalk with other pathways](#Crosstalk-with-other-pathways)

-----


## Compounds Tested
- **Agonist:** benzo-a-pyrene  
- **Antagonist:** 3-methoxy-4-nitroflavone  

-----
  
## Results from reference compounds

### Summary scheme

[Hormonal Card - Aryl Hydrocarbon](Schemes/Cards_Hormones_Aryl%20Hydrocarbon.jpg)


### Result Description

* [Gene Signatures](1.GeneSignatures.md#Aryl-Hydrocarbon)
* [Differential Expression Analysis - Top Genes](2.DEATopGenes.md#Aryl-Hydrocarbon)
* [Differential Expression Analysis - Functional Enrichment](3.DEAFunctional.md#Aryl-Hydrocarbon)
* [WGCNA](4.WGCNA.md#Aryl-Hydrocarbon)




---

## General pathway information
* The **Aryl Hydrocarbon Receptor (AHR)** is a ligand-activated transcription factor in a superfamily involved in sensing endogenous and environmental signals, mediating cellular adaptation to the microenvironment.
* :books: [Nature Review on AHR](https://www.nature.com/articles/s41577-019-0125-8)
* AHR detects environmental chemicals like polyaromatic hydrocarbons and toxins. Initial ligands include exogenous aromatic pollutants like **TCDD** and **benzo-a-pyrene**. It also binds **tryptophan catabolites** (kynurenine, indoxyl sulfate, indole-3-carbinol), expanding its role to metabolic sensing.
* AHR functions by:
  + Binding DNA at **xenobiotic responsive elements (XREs)** as a **heterodimer with ARNT (HIF1β)**
  + Competing with **HIF1α** for ARNT
* **Species differences** in ligand specificity limit rodent models for studying human AHR, though **TCDD** is a full agonist in mice and commonly used in research.

---

## Physiological effects of AHR
AHR is involved in:
- **Detoxification**
- **Tissue homeostasis**
- **Immune response**
- **Carcinogenesis**

Despite early discovery as a receptor for dioxin-like compounds, AHR now is recognized for binding endogenous substrates and regulating various physiological processes.

:books: [Further reading](https://www.nature.com/articles/nrc3846)

---

## Impact on neurodevelopment 
- AHR activation **disrupts neurodevelopment**, especially in early stages due to an immature blood-brain barrier (BBB)
- AHR is **highly expressed** in **endothelial cells** of the BBB
- Activation induces:
  - **CYP1A1**, **CYP1B1**
  - **P-glycoprotein**
  - Detoxification and efflux mechanisms

However:
- **Rodent fetal brains** show **100x higher TCDD transfer** due to immature BBB
- **Wnt/β-catenin** pathway, critical for CNS vascularization, is affected

Overactivation of AHR can impair:
- **Neuronal migration**, **maturation**, **connectivity**  
- Inhibit **CREB-dependent transcription**, reducing **Bdnf** and **Npas4** expression
- :books: [Study 1](https://www.nature.com/articles/srep26386); [Study 2](https://www.ahajournals.org/doi/10.1161/circulationaha.114.011394)

Additional roles:
- Mediates **gut-brain axis** communication - :books: [Nature Article](https://www.nature.com/articles/s41423-020-00585-5)
- Involved in **brain tumor proliferation**  :books: [MDPI Study](https://www.mdpi.com/1422-0067/21/8/2863)

>  There is limited evidence on AHR **inhibition effects**, especially in human-relevant models.

---

## In vitro experiments with iPSC and neural models
#### TCDD exposure
  - Inhibits proliferation in **human (SK-N-SH)** and **murine (C17.2)** neuronal precursors
  - Impairs **NPC proliferation** via **G1/S checkpoint blockade**

#### iPSC-derived motor neurons
  - Susceptible to TCDD-activated AHR
  - Show upregulation of **CYP1A1**, **CYP1B1**  :books: [Study](https://pubmed.ncbi.nlm.nih.gov/35890127/)

> These align with **animal and epidemiological data** showing impaired neurogenesis due to AHR ligands.

---

## Epidemiological data
#### Japanese coastal study 
  - Measured dioxins in breast milk
  - **2,3,7,8-TCDD** inversely correlated with **newborn head circumference**
  - :books: [Study Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC2684781/)
    
>  In our data (CTL04 bulk expo), **TOPGO enrichment term: “head development”**

#### Seveso Women's Health Study (SWHS)
  - :books: [Study Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC3891592/)
  - No significant link between TCDD and pregnancy outcomes
  - Peak lifetime TCDD dose correlated **non-significantly** with lower birthweight
  - **Neuropsychological changes** due to prenatal exposure are observed in animals and **suggested** in human Seveso data. :books: [Study Link](https://www.sciencedirect.com/science/article/abs/pii/S1438463918306400)

#### Autism Spectrum Disorder and AHR
  - Children with ASD have **increased AhR-mediated inflammatory cytokines** (IL-6, STAT3)
  - **STAT3** activates **AhR promoter**
  - **Vitamin D deficiency** common in ASD; metabolized by **CYP1B1**
  - **CYP1B1 expression** reduced in ASD (possibly due to epigenetic silencing)
  - :books: [Study Link](https://www.mdpi.com/1422-0067/22/17/9258)

---

## Evidence of EDC disruption of the pathway

---

## Crosstalk with Other Pathways
- **AHR-ARNT** can indirectly regulate genes **without XREs** via:
  - **Estrogen receptors**
  - **Retinoic acid receptors**

- **Testosterone (T)** vs. **Dihydrotestosterone (DHT)**:
  - **T** induces AHR expression and forms AHR/AR complex
  - Promotes **LRH-1** expression (not mimicked by DHT)
  - Requires **ERK activation** for AHR/AR complex formation  
     :books: [Study Link](https://www.tandfonline.com/doi/full/10.1128/MCB.00011-13)

- **Glucocorticoid (GC)**-induced neuroinflammation:
  - Activates **let7b**, **AHR/ARNT**, **HMGB1/RAGE**, **NRF2/Keap1**
  - Elevates **GFAP**, **TNF-α**, **TBARS**
  - Reduces **Keap1**, **GSH**  
     :books: [Study Link](https://www.sciencedirect.com/science/article/pii/S0006295224006932)
