## Table of Content

* [Compounds tested](#Compounds-tested)
* [Results from reference compounds](#Results-from-reference-compounds)
* [General pathway information](#General-pathway-information)
* [Physiological effect of the hormone](#Physiological-effect-of-the-hormone)
* [Impact on neurodevelopment](#Impact-on-neurodevelopment)
* [In vitro experiment with iPSC and neural models](#In-vitro-experiment-with-iPSC-and-neural-models)
* [Epidemiological data](#Epidemiological-and-clinical-data)
* [Evidence of EDC disruption of the pathway](#Evidence-of-EDC-disruption-of-the-pathway)
* [Crosstalk with other pathways](#Crosstalk-with-other-pathways)

-----


## Compounds tested
- GR agonist: Dexamethasone (1.0 $\mu M$)
- GR antagonist: RU38486 (10.0 $\mu M$)

-----


## Results from reference compounds

### Summary scheme

[Hormonal Card - Glucocorticoid](Schemes/Cards_Hormones_Glucocorticoids.jpg)

### Result Description
* [Gene Signatures](1.GeneSignatures.md#glucocorticoid)
* [Differential Expression Analysis - Top Genes](2.DEATopGenes.md#glucocorticoid)
* [Differential Expression Analysis - Functional Enrichment](3.DEAFunctional.md#glucocorticoid)
* [WGCNA](4.WGCNA.md#glucocorticoid)

-----

## General pathway information


####  :books: [A General Introduction to Glucocorticoid Biology](https://www.frontiersin.org/articles/10.3389/fimmu.2019.01545/full)

* __GCs are mainly synthesized in the cortex of the adrenal gland__. Adrenal GC production is __regulated by__ the hypothalamic-pituitary-adrenal __(HPA) axis__
* Under basal conditions GCs are released from the adrenal glands in the bloodstream in a circadian rhythm. When the HPA-axis is stimulated, corticotropin-releasing hormone (CRH), and arginine vasopressin (AVP) are released from the hypothalamus. CRH and AVP bind their receptor in the anterior pituitary inducing the release of adrenocorticotrophic hormone (ACTH) in the circulation. ACTH will in turn stimulate the adrenal gland to synthetize and secrete GC hormones (cortisol) in the circulation.
* The HPA axis is subject to a negative feedback inhibition by GCs.

* Once secreted in the bloodstream GCs are bound to and transported by plasma proteins which keep the GCs inactive. __Corticosteroid-binding globulin (CBG) is the main GC-binding protein in the plasma__, with about 80–90% of the GCs bound to it.

* The bioavailability of GCs in the cytoplasm is regulated by the balance between active and inactive forms of GCs. __Two enzymes are responsible for the conversion between inactive cortisone and active cortisol__: (I)  11β-hydroxysteroid dehydrogenase 1 (11β-HSD1) catalyzes the conversion of cortisone to cortisol: (II) 11β-HSD2 carries out the opposite reaction.

* __GCs acts mainly through the glucocorticoid receptor GR (gene NR3C1)__.  It belongs to the nuclear receptor superfamily of transcription factors (TFs) and is a 97 kDa protein that is constitutively and  expressed throughout the body. GCs exert cellular and tissue-specific effects due to the existence of different GR isoforms on the one hand and cell- and context-specific allosteric signals influencing GR function on the other hand

* In the absence of intracellular bioactive GCs, the GR finds itself as a monomer in the cytoplasm where it resides in a multiprotein complex with HSP. GC binding induces a change in the chaperone complex bound to GR, after which it translocates to the nucleus to transactivate (+) or transrepress (-) gene transcription as a monomer or a dimer. In particular, the chaperone FK506-binding protein 51 (FKBP51, encoded by the FKBP5 gene), binds to heat-shock protein 90 (Hsp90) complexes and through them induces conformational changes of GR, leading to its *reduced affinity for cortisol*. Single nucleotide polymorphisms (SNPs) within the FKBP5 gene locus have been reported to impact stress responsivity, and risk or resilience to psychiatric disorders. :books: [FKBP5 polymorphisms induce differential glucocorticoid responsiveness in primary CNS cells – First insights from novel humanized mice](https://doi.org/10.1111/ejn.14999)

* GCs can also act thorough non-genomic mechanisms by binding membrane receptor (probably a G-coupled receptor). The main non-genomic effects described in :books: [Non-genomic Effects of Glucocorticoids: An Updated View](10.1016/j.tips.2018.11.002) are: (I) impairment of calcium homeostasis, (II) altered reactive oxygen and nitrogen species levels, (III) interference of the inflammatory pathways and (IV) the impact on mitochondrial metabolism and apoptosis.

--------

## Physiological effect of the hormone


#### :books: [Glucocorticoid and Mineralocorticoid Receptors in the Brain: A Transcriptional Perspective](https://pmc.ncbi.nlm.nih.gov/articles/PMC6777400/)

* Cortisol activity on the brain is mediated by the _high-affinity mineralocorticoid receptor (MR, encoded by the gene NR3C2)_ and the _lower affinity glucocorticoid receptor (GR encoded by the gene NR3C1)_. MRs and GRs can mediate distinct, sometimes opposite, effects of glucocorticoids.
* Given its high affinity, MR is occupied at basal hormone levels, whereas GR is activated at the circadian peak of glucocorticoid secretion and during stress. Selective brain regions contain MR, which binds aldosterone to control physiology and behavior in relation to salt balance.
* GR is present in most brain regions and cell types, whereas MR is mainly expressed in limbic areas such as the hippocampus, amygdala, and prefrontal cortex.
* GR and MR control a wide range of processes, ranging from _neuronal differentiation and excitability to behavioral reactivity, mood, and cognition_. In addition to complementary actions of GR and MR, they can also exert opposing effects, even within the same cell type.
* Both receptor types can mediate nongenomic steroid effects, but they mainly act as ligand-activated transcription factors. MR and GR protein structure is similar; the receptors can form heterodimers on the DNA at glucocorticoid response elements (GREs), and they share a number of target genes. The transcriptional basis for opposite effects on cellular physiology remains largely unknown.
* Potential mechanisms of transcriptional specificity for MRs and GRs: unique GR binding to “negative GREs,” direct binding to other transcription factors, and binding to specific DNA sequences in conjunction with other transcription factors (e.g MRs and NeuroD proteins). MR- and GR-specific effects may also depend on specific interactions with transcriptional coregulators, downstream mediators of transcriptional receptor activity.


#### :books: [Genomic and epigenomic mechanisms of glucocorticoids in the brain](https://www.nature.com/articles/nrendo.2017.97)
* DNA methylation patterns in brain have been studied to find association with psychological resilience or vulnerability (depending on the methylation pattern) and with familiarity of psychiatric disorders. Furthermore, treatment with antidepressants [showed a change in the methylation pattern](10.3389/fncel.2012.00018)
* high level of mental resilience to early life stress (ELS), expressed as lower incidence of psychiatric outcomes, was associated with [attenuated NR3C1 DNAm levels](https://www.tandfonline.com/doi/full/10.1080/08039488.2024.2436987#abstract)

:books: [Other interesting papers](https://www.nature.com/articles/s41586-025-09083-y#Sec2)


-----------------------

## Impact on neurodevelopment

### Epidemiological and clinical data

#### :books: Review: [The Link Between Perinatal Glucocorticoids Exposure and Psychiatric Disorders](https://www.nature.com/articles/pr9201189). 
Altered levels of glucocorticoids or the hypothalamic-pituitary-adrenal axis activity in either the mother or offspring in the perinatal period cam impoct on physiological, endocrinological and behavioural functions. 

#### :books: [Effects of antenatal glucocorticoids on the developing brain](https://www.sciencedirect.com/science/article/abs/pii/S0039128X16300563?via%3Dihub)

* Throughout development, the placenta tightly regulates fetal GC levels, while the HPA axis regulates maternal GC
* Endogenous GC levels in fetal circulation increase rapidly late in gestation to promote the development of fetal organs.
* This acts through the receptor for cortisol, GR, which is expressed throughout gestation.
* In preterm birth, the developing fetus does not receive sufficient exposure to endogenous GCs _in utero_ for proper organ development predisposing the neonate to complications (e.g intraventricular hemorrhage, respiratory distress syndrome and necrotizing enterocolitis).
* Synthetic GCs (sGCs) have proven useful in the prevention of these complications since they are able to promote the rapid maturation of underdeveloped organs.
* However, they may also trigger adverse developmental side effects. 

#### :books: [Effects of postnatal glucocorticoids on brain structure in preterm infants, a scoping review](https://www.sciencedirect.com/science/article/abs/pii/S0149763423000039?via%3Dihub)
* Glucocorticoids are used in neonatal intensive care units to prevent or reduce the severity of chronic lung disease in preterm infants.
* Searched scientific literature for original research on human preterm infants, postnatal GCs, and brain structure.
* 11 studies assessed the effects of GCs on structural brain outcomes. 56 studies reported brain injury, but not structure.
* Dexamethasone was consistently associated with decreased total and regional brain volumes, including cerebellar volumes.
* Hydrocortisone was often, but not always associated with absence of brain volume differences. 


#### :books: [Glucocorticoids as Mediators of Adverse Outcomes of Prenatal Stress](https://www.cell.com/trends/neurosciences/fulltext/S0166-2236(20)30068-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0166223620300680%3Fshowall%3Dtrue)
* Prenatal stress is associated with neurodevelopmental alterations that persist after birth and manifest at the behavioral level and the physiological one.
* Elevated glucocorticoid signaling in utero may be one of the key mediators of prenatal stress effects on the offspring.
* Prenatal glucocorticoids may impact the activity of the fetal hypothalamic-pituitary-adrenal (HPA) axis, disrupt neurodevelopmental processes and alter the epigenetic landscape of the fetus.

#### :books: [Human cortical neurogenesis is altered via glucocorticoid-mediated regulation of ZBTB16 expression](https://www.cell.com/neuron/fulltext/S0896-6273(24)00089-8)
* GCs increase a specific type of basal progenitors (co-expressing PAX6 and EOMES) that has been shown to contribute to cortical expansion in gyrified species in human cortical brain organoids
* ZBTB16 mediates the increase of neurons in the cortex in response to GCs
* rs648044 is a SNP present within the ZBTB16 locus and it correlates with neurobehavioral and brain structural outcomes (especially in the cortex, where ZBTB16 plays a pivotal role)
* rs648044 SNP showed low methylation levels suggesting that, the enhancer containing it, modulates dex-mediated transcriptional effects on ZBTB16

-----------

## In vitro experiment with iPSC and neural models

#### :books: [In vitro modeling of the neurobiological effects of glucocorticoids: A review](https://doi.org/10.1016/j.ynstr.2023.100530)
Extensive table reporting reference to different types of experimental models in which GR signalling has been studied. 

#### :books: Review: [Impact of glucocorticoid on neurogenesis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5558474/)
* Both embryonic and adult neurogenesis are influenced by glucocorticoids
* Proliferation/differentiation of NSC is influenced by glucocorticoids through action on PI3K/Akt, SHH, WNT. 
* MR and GR differentially affect the proliferation and differentiation of NSPCs. 
* Catalytic conversion of glucocorticoid to cortisone by placental 11 beta-hydroxysteroid dehydrogenase type 2 prevents maternal glucocorticoid transfer to the fetus
* Maternal prolonged elevation of glucocorticoid levels (chronic severe stress) may exceeds catalytic conversion capacity, with GCs reaching the fetus (Reynolds, 2013).
* Fetal exposure to GCs can also result from therapeutic administration of synthetic glucocorticoids to promote fetal lung maturation.

#### :books: [Prenatal programming of stress responsiveness and behaviours: Progress and perspectives](https://onlinelibrary.wiley.com/doi/10.1111/jne.12674)

* Parental exposure to stress or glucocorticoids either before or during pregnancy can have profound influences on neurodevelopment, neuroendocrine function and behaviours in offspring.
* Specific outcomes are dependent on factors such as nature, intensity and timing of the exposure, sex and age of the subject.
* There is strong emerging evidence that epigenetic processes are involved in mediation long-term outcomes of these early exposures. 

#### :books: [Emerging roles of the unfolded protein response (UPR) in the nervous system: A link with adaptive behavior to environmental stress?](https://www.sciencedirect.com/science/article/abs/pii/S1937644820300046?via%3Dihub)

* Chronic exposition to stressors operates as a risk factor for psychiatric diseases.
* Alterations to organelle function at the level of mitochondria and endoplasmic reticulum (ER) are emerging as possible factors contributing to neuronal dysfunction.
* ER proteostasis alterations elicit the unfolded protein response (UPR). UPR has been associated to neurodevelopment, synaptic plasticity and neuronal connectivity.
* Recent studies suggest a role of the UPR in the adaptive behavior to stress, suggesting a mechanistic link between environmental and cellular stress. 

#### :books: [Circadian glucocorticoids throughout development](https://www.frontiersin.org/articles/10.3389/fnins.2023.1165230/full) 

* The developing circadian clock is shaped by maternal GCs. GC deficits, excess, or exposure at the wrong time of day leads to persisting effects later in life.
* The master pacemaker of the circadian system is located in the hypothalamic suprachiasmatic nucleus (SCN)
* During the first two thirds of pregnancy, the embryo/fetus is protected by the placenta from excessive levels of GCs by expressing (among others) the enzyme 11β-Hydroxysteroid dehydrogenase (11β-HSD2) which inactivates GCs.
* Even though only 10% of maternal GCs reach the fetus, a circadian rhythmicity of GCs can still be detected in fetal blood.
* During the last third of pregnancy, there is a gradual increase of GC levels in fetal blood that helps the fetus to continue developing and growing.
* Evidence shows that maternal stress or circadian rhythm disruption leads to a higher risk of developing metabolic, behavioral, and sleep disorders later in life
* During adulthood, GCs are one of the main hormonal outputs of the circadian system, peaking at the beginning of the active phase (i.e., the morning in humans and the evening in nocturnal rodents) and contributing to the coordination of complex functions such as energy metabolism and behavior

#### :books: [Chronic exposure to glucocorticoids amplifies inhibitory neuron cell fate during human neurodevelopment in organoids](https://www.science.org/doi/10.1126/sciadv.adn8631)
* Chronic GC exposure influenced lineage specification primarily by priming the inhibitory neuron lineage through transcription factors like PBX3

-------------------

## Evidence of EDC disruption of the pathway


#### :books: [Prenatal drug exposure and neurodevelopmental programming of glucocorticoid signalling](https://onlinelibrary.wiley.com/doi/10.1111/jne.12786)
Prenatal neurodevelopment is dependent on precise functioning of multiple signalling pathways in the brain, including those mobilised by glucocorticoids (GC) and endocannabinoids (eCBs).
* Prenatal exposure to drugs of abuse has been shown to not only impact GC signalling, but also alter functioning of the hypothalamic-pituitary-adrenal (HPA) axis.
* Such exposures can have long-lasting neurobehavioural consequences.
* Prenatal drug exposure can impact on HPA axis and stress response; interactions between GC and EC signalling in the developing brain has potential for neurodevelopmental consequences.

#### :books: [Glucocorticoid and mineralocorticoid receptors and corticosteroid homeostasis are potential targets for endocrine-disrupting chemicals](https://www.sciencedirect.com/science/article/pii/S016041201931390X?via%3Dihub)
Very comprehensive review on scientific evidence of disturbances of corticosteroid homeostasis and action by EDCs of different classes: metals, metalloids, pesticides, bisphenol analogues, flame retardant, other chemicals.


-------------

## Crosstalk with other pathways

#### :books: [Maternal Prenatal Stress, Thyroid Function and Neurodevelopment of the Offspring: A Mini Review of the Literature](https://www.frontiersin.org/articles/10.3389/fnins.2021.692446/full)
* Both the Hypothalamic-Pituitary-Adrenal (HPA) and the Hypothalamic-Pituitary-Thyroid (HPT) axes are involved in stress responses, whereas, their final effectors, the Glucocorticoids (GCs) and the Thyroid Hormones (THs), mediate several fundamental processes involved in neurodevelopment.
* The effects of these hormones on brain development are time and dose-dependent. Inadequate or excess concentrations of both GCs and THs can cause abnormalities in the neuronal and glial structures and functions.
* Maternal stress and GC excess can impact on growth and neurodevelopment of the offspring. Chronic stress and alterations of the HPA axis interacts and influences HPT axis and TH production.
* Animal studies: increased GC concentrations related to maternal stress, most likely reduce maternal and thus fetal circulating THs, either directly or through modifications in the expression of placental enzymes responsible for regulating hormone levels in fetal microenvironment.

#### :books: Review: [Maternal hormonal milieu influence on fetal brain development](https://pubmed.ncbi.nlm.nih.gov/29484271/)
* Subtle changes in fetal brain development have been observed even for maternal hormone levels within the currently accepted physiologic ranges. 
* Thyroid hormones are required for normal brain development. Despite serum TSH appearing to be the most accurate indicator of thyroid function in pregnancy, maternal serum free T4 levels in the first trimester of pregnancy are the major determinant of postnatal psychomotor development.
* Even a transient period of maternal hypothyroxinemia at the beginning of neurogenesis can confer a higher risk of expressive language and nonverbal cognitive delays in offspring. 
* Corticosteroids are determinant in suppressing cell proliferation and stimulating terminal differentiation, a fundamental switch for the maturation of fetal organs.
* Intrauterine exposure to stress or high levels of glucocorticoids, endogenous or synthetic, has a molecular and structural impact on brain development and appears to impair cognition and increase anxiety and reactivity to stress.
* Limbic regions, such as hippocampus and amygdala, are particularly sensitive.
* Repeated doses of prenatal corticosteroids seem to have short-term benefits of less respiratory distress and fewer serious health problems in offspring. Nevertheless, the impact on neurodevelopmental outcomes needs further clarification. 
* In adults, glucocorticoids generally inhibit thyroid functions. In fetal development however the interplay is much more complex.


#### :books: [Restricted effects of androgens on glucocorticoid signaling in the mouse prefrontal cortex and midbrain](https://pmc.ncbi.nlm.nih.gov/articles/PMC10830692/)
* chronic treatment with corticosterone, dihydrotestosterone, a combination of both, and corticosterone in combination with the AR antagonist enzalutamide, we compared the expression of glucocorticoid receptor target genes in brain regions where AR and GR are co-expressed, namely: prefrontal cortex, hypothalamus, hippocampus, ventral tegmental area and substantia nigra. _Androgen affected glucocorticoid signaling only in the prefrontal cortex and the substantia nigra_. 
* Dihydrotestosterone and corticosterone independently and inversely regulated expression of Sgk1 and Tsc22d3 in prefrontal cortex. 
* _AR antagonism_ with enzalutamide _attenuated corticosterone-induced expression of Fkbp5 in the prefrontal cortex_ and of Fkbp5 and Sgk1 in the substantia nigra
* _AR antagonism increased expression of Th and Slc18a1_, two genes coding for key components of the dopaminergic system

#### :books: [Cross-talk between glucocorticoid and retinoic acid signals involving glucocorticoid receptor interaction with the homoeodomain protein Pbx1](https://pmc.ncbi.nlm.nih.gov/articles/PMC1223238/)
* in P19 embryonal carcinoma cells, GCs potentiate RA-induced expression of the murine Hoxb -1 gene through an autoregulatory element, b1-ARE, recognized by the Pbx1 and HOXB1 homoeodomain proteins
* A physical interaction between the GR and Pbx1 proteins was demonstrated in vivo by co-immunoprecipitation experiments. These results are compatible with a model in which the GC/RA synergy is mediated by a direct interaction between the GR and Pbx1
* On the basis of the ubiquitous expression of both GR and Pbx1, a number of genes regulated by Pbx are likely to be important targets for GC-mediated 'cross-talk'


