# Crosstalk among pathways
[Hormones do not act in a strictly specific way on downstream pathways](https://www.sciencedirect.com/science/article/pii/S0160412025001515), and each hormone can regulate different signalling pathways at the same time. 
Therefore, studying hormonal crosstalk is essential for deciphering the complex mechanisms that guide human neurodevelopment. 

## Table of Content

* [1. ArylHydrocarbon](#ArylHydrocarbon)
* [2. Androgen](#Androgen)  
* [3. Estrogen](#Estrogen)
* [4. Glucocorticoid](#Glucocorticoid)
* [5. LiverX](#LiverX)
* [6. RetinoicAcid](#RetinoicAcid)
* [7. Thyroid](#Thyroid)

-----------------

## ArylHydrocarbon
* **ArHy, ESTR and RA**:
  * **AHR-ARNT** can indirectly regulate genes **without XREs** via:
    * **Estrogen receptors**
    * **Retinoic acid receptors**

* **ArHy and ANDR**:
  * **Testosterone (T)** vs. **Dihydrotestosterone (DHT)**:
    * **T** induces AHR expression and forms AHR/AR complex
    * Promotes **LRH-1** expression (not mimicked by DHT)
    * Requires **ERK activation** for AHR/AR complex formation  
      🔗 [Study Link](https://www.tandfonline.com/doi/full/10.1128/MCB.00011-13)

* **ArHy and GC**:
  * **Glucocorticoid (GC)**-induced neuroinflammation:
    * Activates **let7b**, **AHR/ARNT**, **HMGB1/RAGE**, **NRF2/Keap1**
    * Elevates **GFAP**, **TNF-α**, **TBARS**
    * Reduces **Keap1**, **GSH**  
      🔗 [Study Link](https://www.sciencedirect.com/science/article/pii/S0006295224006932)

## Androgen
Androgens do not operate in isolation; rather, they engage in complex crosstalk with several key hormonal pathways. 
* **ANDR and ESTR**:
   * Androgens are not just ligands for the androgen receptor, crucially, they also serve as precursors for estrogen synthesis via the enzyme aromatase, enabling a dynamic hormonal interplay. Estrogens are far more potent than their androgen precursors, and this local conversion supports dynamic hormone interconversion and regulatory balance [Blakemore & Naftolin, 2016](https://journals.physiology.org/doi/full/10.1152/physiol.00054.2015).
* **ANDR and GC**:
   * However, the interplay deepens: androgens are produced not only in the gonads but also in the adrenal cortex, where they share critical biosynthetic and regulatory pathways with glucocorticoids and mineralcorticoids. A classic example is Congenital Adrenal Hyperplasia, most caused by genetic mutations in the enzyme 21-hydroxylase. When this enzyme is deficient, cortisol synthesis is blocked, resulting in persistently low cortisol levels. This deficiency disrupts negative feedback on ACTH. As cortisol precursors accumulate, such as 17 hydroxyprogesterone, they are rerouted into the androgen synthesis pathway, causing a dramatic increase in adrenal androgens [Turcu et al., 2014](https://onlinelibrary.wiley.com/doi/10.1002/cphy.c140006). The androgen receptor (AR) sits at the heart of a complex steroid receptor network, influencing gene expression through extensive chromatin-level interactions with receptors such as ER, PR, GR. The crosstalk between signal-dependent TFs is an important step in the reprogramming of chromatin sites. a signal-activated TF can expand or restrict the chromatin binding of another TF. This crosstalk can rewire gene programs and thus alter biological processes and influence the progression of diseases like breast or prostate cancer. For example, in prostate cancer, AR-targeted treatments enable GR to increase its expression and occupy AR-binding sites on pre-accessible chromatin, a switch that sustains androgen-response transcription programs and drives therapy resistance [Paakinaho & Palvimo, 2021](https://academic.oup.com/nar/article/49/4/1951/6125668?login=false).
* **ANDR and THYR**:
   * Layered on top of these interactions is another hormonal dimension: thyroid hormones (TH). T3 and T4 bind to the AR gene promoter and elevate its expression. They also increase levels of 5α-reductase, promoting conversion of testosterone into more potent dihydrotestosterone (DHT), which further strengthens AR signaling. AR, in turn, binds to androgen response elements (AREs) in promoters of thyroid-related genes, such as deiodinases and thyroid receptor isoforms, reinforcing a feedback loop between the androgen and thyroid systems [Torabinejad et al., 2023](https://linkinghub.elsevier.com/retrieve/pii/S0099239923002145).
   
Together, these intertwined pathways illustrate a biological equilibrium: androgens, estrogens, glucocorticoids, and thyroid hormones constantly intersect through shared precursors, receptor-level competition and cooperation, and chromatin-level interactions.

## Estrogen
* **ESTR and GC**:
  * [Stress, Sex, and Sugar: Glucocorticoids and Sex-Steroid Crosstalk in the Sex-Specific Misprogramming of Metabolism](https://pmc.ncbi.nlm.nih.gov/articles/PMC7382384/#s5)
    * ER and GR have been shown to affect each other’s action both by altering receptor and ligand availability as well as by modulating each other’s genomic binding and transcriptional end points.

  * [Estradiol Rapidly Rescues Synaptic Transmission from Corticosterone-induced Suppression via Synaptic/Extranuclear Steroid Receptors in the Hippocampus](https://academic.oup.com/cercor/article/22/4/926/424336?login=true)

* **ESTR and THYR**:
  * [Thyroid Hormone Causes Mitogen-Activated Protein Kinase-Dependent Phosphorylation of the Nuclear Estrogen Receptor](https://academic.oup.com/endo/article/145/7/3265/2500000)

  * [Thyroid hormone- and estrogen receptor interactions with natural ligands and endocrine disruptors in the cerebellum](https://www.sciencedirect.com/science/article/pii/S0091302217300602)

* **ESTR and ArHy**:
  * [AhR and ARNT modulate ER signaling](https://www.sciencedirect.com/science/article/pii/S0300483X09004697)

  * [Chemical activation of estrogen and aryl hydrocarbon receptor signaling pathways and their interaction in toxicology and metabolism](https://www.tandfonline.com/doi/10.1080/17425255.2019.1569627?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)

* **ESTR and RA**:
  * [Limiting Effects of RIP140 in Estrogen Signaling: POTENTIAL MEDIATION OF ANTI-ESTROGENIC EFFECTS OF RETINOIC ACID](https://www.jbc.org/article/S0021-9258(19)30550-2/fulltext#FN1)

  * [An R package for generic modular response analysis and its application to estrogen and retinoic acid receptor crosstalk](https://www.nature.com/articles/s41598-021-86544-0#Sec9)

* **ESTR and LivX**:
  * [LXRβ/estrogen receptor-α signaling in lipid rafts preserves endothelial integrity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3726156/)



## Glucocorticoid
* **GC and THYR**:
  * [Maternal Prenatal Stress, Thyroid Function and Neurodevelopment of the Offspring: A Mini Review of the Literature](https://www.frontiersin.org/articles/10.3389/fnins.2021.692446/full)
     * Both the Hypothalamic-Pituitary-Adrenal (HPA) and the Hypothalamic-Pituitary-Thyroid (HPT) axes are involved in stress responses, whereas, their final effectors, the Glucocorticoids (GCs) and the Thyroid Hormones (THs), mediate several fundamental processes involved in neurodevelopment.
     * The effects of these hormones on brain development are time and dose-dependent. Inadequate or excess concentrations of both GCs and THs can cause abnormalities in the neuronal and glial structures and functions.
     * Maternal stress and GC excess can impact on growth and neurodevelopment of the offspring. Chronic stress and alterations of the HPA axis interacts and influences HPT axis and TH production.
     * Animal studies: increased GC concentrations related to maternal stress, most likely reduce maternal and thus fetal circulating THs, either directly or through modifications in the expression of placental enzymes responsible for regulating hormone levels in fetal microenvironment.

  * [Maternal hormonal milieu influence on fetal brain development](https://pubmed.ncbi.nlm.nih.gov/29484271/)
     * Subtle changes in fetal brain development have been observed even for maternal hormone levels within the currently accepted physiologic ranges.
     * Thyroid hormones are required for normal brain development. Despite serum TSH appearing to be the most accurate indicator of thyroid function in pregnancy, maternal serum free T4 levels in the first trimester of pregnancy are the major determinant of postnatal psychomotor development.
     * Even a transient period of maternal hypothyroxinemia at the beginning of neurogenesis can confer a higher risk of expressive language and nonverbal cognitive delays in offspring.
     * Corticosteroids are determinant in suppressing cell proliferation and stimulating terminal differentiation, a fundamental switch for the maturation of fetal organs.
     * Intrauterine exposure to stress or high levels of glucocorticoids, endogenous or synthetic, has a molecular and structural impact on brain development and appears to impair cognition and increase anxiety and reactivity to stress.
     * Limbic regions, such as hippocampus and amygdala, are particularly sensitive.
     * Repeated doses of prenatal corticosteroids seem to have short-term benefits of less respiratory distress and fewer serious health problems in offspring. Nevertheless, the impact on neurodevelopmental outcomes needs further clarification.
     * In adults, glucocorticoids generally inhibit thyroid functions. In fetal development however the interplay is much more complex.

**GC and ANDR**:
  * [Restricted effects of androgens on glucocorticoid signaling in the mouse prefrontal cortex and midbrain](https://pmc.ncbi.nlm.nih.gov/articles/PMC10830692/)
      * chronic treatment with corticosterone, dihydrotestosterone, a combination of both, and corticosterone in combination with the AR antagonist enzalutamide, we compared the expression of glucocorticoid receptor target genes in brain regions where AR and GR are co-expressed, namely: prefrontal cortex, hypothalamus, hippocampus, ventral tegmental area and substantia nigra. _Androgen affected glucocorticoid signaling only in the prefrontal cortex and the substantia nigra_. 
      * Dihydrotestosterone and corticosterone independently and inversely regulated expression of Sgk1 and Tsc22d3 in prefrontal cortex. 
      * _AR antagonism_ with enzalutamide _attenuated corticosterone-induced expression of Fkbp5 in the prefrontal cortex_ and of Fkbp5 and Sgk1 in the substantia nigra
      * _AR antagonism increased expression of Th and Slc18a1_, two genes coding for key components of the dopaminergic system

**GC and RA**:
   * [Cross-talk between glucocorticoid and retinoic acid signals involving glucocorticoid receptor interaction with the homoeodomain protein Pbx1](https://pmc.ncbi.nlm.nih.gov/articles/PMC1223238/)
        * in P19 embryonal carcinoma cells, GCs potentiate RA-induced expression of the murine Hoxb -1 gene through an autoregulatory element, b1-ARE, recognized by the Pbx1 and HOXB1 homoeodomain proteins
        * A physical interaction between the GR and Pbx1 proteins was demonstrated in vivo by co-immunoprecipitation experiments. These results are compatible with a model in which the GC/RA synergy is mediated by a direct interaction between the GR and Pbx1
        * On the basis of the ubiquitous expression of both GR and Pbx1, a number of genes regulated by Pbx are likely to be important targets for GC-mediated 'cross-talk'
        
        
## LiverX
* **LivX and RA**:
   * LXRs form obligate heterodimers with retinoid X receptors (RXRs) to regulate target genes. In differentiating stem cells, RXR partners orchestrate neuronal gene programs: [for example](https://doi.org/10.1016/j.mce.2017.07.033), chromatin profiling shows LXR:RXR (and RAR:RXR) complexes drive distinct gene networks during ESC→neuron differentiation . Thus LXR activation often overlaps with retinoic acid (RAR/RXR) signaling in neurogenesis.

* **LivX and THYR**:
   * Beyond RXR, LXRs intersect with other nuclear receptors. In mouse cortex, LXR $\beta$ and thyroid hormone receptor $\alpha$ (TR $\alpha$) share DNA response elements and [co-regulate genes critical for cortical layering](https://doi.org/10.1073/pnas.1006162107). Loss of LXR $\beta$ delays neuronal migration, but rising thyroid hormone postnatally can compensate via TR$\alpha$ , suggesting dynamic LXR–TR crosstalk in cortical development. Several types of crosstalk between TRs and LXRs have been identified and [crosstalk has also been observed in other physiological systems](https://doi.org/10.1507/endocrj.ej11-0114) such as central nervous system rather than lipid metabolism.

* **LivX and GC**:
   * Similarly, in metabolic tissues LXR activation has been shown to directly modulate glucocorticoid receptor (GR) signaling. LXRs bind to glucocorticoid response elements (GREs) and [can inhibit GR-driven transcription](https://doi.org/10.1371/journal.pone.0026751). For example, GW3965 (our LXR agonist) attenuated dexamethasone-induced gene expression of gluconeogenic enzymes by recruiting LXR/RXR to GREs and displacing GR. Although this work was in liver cells, it shows that LXR activation can repress glucocorticoid pathways via competitive DNA binding, a mechanism that could extend to neural progenitors exposed to stress hormones.

* **LivX and ESTR**:
   * Estrogen signaling is also modulated by LXR. In rat sensory neurons, [LXR agonism counteracted ER$\alpha$ -driven inflammatory gene expression](https://doi.org/10.1016/j.bbih.2024.100757). For instance, 24S,25-epoxycholesterol (an LXR ligand) suppressed ER$\alpha$ -mediated IL-6 and $Cav\alpha2\delta$ induction in DRG neurons, reducing estrogen’s pro-nociceptive effects. This implies that LXR activation can dampen estrogen receptor pathways in neurons, which may influence neuronal maturation or responses to estrogenic endocrine disruptors.


Importantly, these hormone cross-regulations influence neural stem cell behavior. LXR agonists promote progenitor proliferation: [in neural cultures](https://doi.org/10.1016/j.bbrc.2016.12.163), GW3965 and LXR623 significantly stimulate NPC expansion (via MEK/ERK signaling) in wild-type but not LXR-knockout cells. 

In human models, the recent comprehensive screen using [fetal human neurospheres](https://doi.org/10.1016/j.envint.2025.109400) show LXR (and RAR, RXR, GR) agonists changed NPC proliferation and differentiation, and RNA-seq revealed activation of conserved pathways (Notch, Wnt) and extensive receptor crosstalk. 

Together these data suggest that LXR signaling intermingles with other hormonal programs to shape gene expression, progenitor proliferation, and neuronal differentiation. This crosstalk also affects vulnerability to chemicals. The human neurosphere study highlights that perturbing LXR (or GR, RXR, etc.) at physiological concentrations disrupts key neurodevelopmental processes. In other words, endocrine disruptors that target one receptor may indirectly influence LXR-dependent pathways (and vice versa), altering cortical neurogenesis.


## Retinoic Acid
* **RA and LivX**:
  * **LXRs** form obligate heterodimers with retinoid X receptors (RXRs) to regulate target genes. In differentiating stem cells, RXR partners orchestrate neuronal gene programs: [for example](https://doi.org/10.1016/j.mce.2017.07.033), chromatin profiling shows LXR:RXR (and RAR:RXR) complexes drive distinct gene networks during ESC→neuron differentiation. Thus LXR activation often overlaps with retinoic acid (RAR/RXR) signaling in neurogenesis.

## Thyroid
* **THYR and RA**:
   * Retinoic acids share carrier proteins with THs and can increase MCT8 expression to increase their import.
     
* **THYR and LivX**:
  * Nuclear liver X receptor beta interacts with TH signaling in regulating cortical layering by acting on the receptor ApoER2.

