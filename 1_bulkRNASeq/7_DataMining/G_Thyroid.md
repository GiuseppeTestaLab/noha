## Table of Content

* [Compounds tested](#Compounds-tested)
* [Results from reference compounds](#Results-from-reference-compounds)
* [General pathway information](#General-pathway-information)
* [Physiological effect of the hormone](#Physiological-effect-of-the-hormone)
* [Impact on neurodevelopment](#Impact-on-neurodevelopment)
* [In vitro experiment with iPSC and neural models](#In-vitro-experiment-with-iPSC-and-neural-models)
* [Epidemiological data](#Epidemiological-data)
* [Evidence of EDC disruption of the pathway](#Evidence-of-EDC-disruption-of-the-pathway)
* [Crosstalk with other pathways](#Crosstalk-with-other-pathways)

-----

## Results from reference compounds

### Summary scheme

[Hormonal Card - Thyroid](Schemes/Cards_Hormones_Thyroid.jpg)

### Result Description
* [Gene Signatures](1.GeneSignatures.md#thyroid)
* [Differential Expression Analysis - Top Genes](2.DEATopGenes.md#thyroid)
* [Differential Expression Analysis - Functional Enrichment](3.DEAFunctional.md#thyroid)
* [WGCNA](4.WGCNA.md#thyroid)

-----

## Compounds tested

- __Agonist__: T3 (0.003 $\mu M$)
- __Antagonist__: NH3 (0.3 $\mu M$)


## General pathway information

* The thyroid gland synthesises and secretes two primary hormones, triiodothyronine (T3) and its inactive form tetraiodothyronine (T4), collectively referred to as thyroid hormones (THs). THs are lipophilic molecules, derived from tyrosine and containing iodine. :books: [Thyroid hormone biosynthesis and release](https://doi.org/10.1016/j.mce.2017.01.038).

* Conversion from T4 to T3 is predominantly carried out in target tissues by selenoprotein enzymes deiodinase I (DIO1) and II (DIO2). However, a modest amount of T3 is secreted by the thyroid gland as well. THs modifications (decarboxylation, deamination, ether-link cleavage, sulfation and glucoronidation) affect their bioactivity and metabolism. Most lead to inactivation, for example, sulfation and glucuronidation are considered detoxification reactions as they increase the solubility of the product. Both T3 and T4 inactivation is catalysed by the enzyme deiodinase III (DIO3) that removes an iodine atom. Inactivation of THs is key to maintaining a homeostatic balance, preventing the phenomena of thyrotoxicosis caused by an excess of THs, which may lead to cardiovascular dysfunctions and osteoporosis. :books: [Transport, Metabolism, and Function of Thyroid Hormones in the Developing Mammalian Brain](https://doi.org/10.3389/fendo.2019.00209);  [Thyrotoxicosis](https://pubmed.ncbi.nlm.nih.gov/29489233/)
  
* THs production and delivery to target tissues throughout the body is meticulously controlled by the hypothalamic-pituitary-thyroid (HPT) axis. Indeed, the release of THs from the thyroid gland is prompted by a signal transduction cascade originating from the hypothalamus. :books: [Genetic Determination of the Hypothalamic-Pituitary-Thyroid Axis: Where Do We Stand? ](https://doi.org/10.1210/er.2014-1081)
  
* Once released into the circulation, THs are bound to carrier proteins such as the thyroxine binding globulin (TBG), transthyretin (TTR) or albumin, which contribute to buffer the amount of free hormones available for tissue uptake. :books: [T3 levels and thyroid hormone signaling](https://doi.org/10.3389/fendo.2022.1044691).
  
* Being lipophilic molecules, THs may enter the target cells via passive diffusion, nevertheless their uptake can be mediated by ATP-dependent mecanisms involving monocarboxylate transporters MCT8 and MCT10, the organic anion transporter proteins (OATPs), NTCP, LAT1 (SLC7A5) and LAT2 (SLC7A8). Other secondary carriers from the SLC family are SLCO1C1 and SLC16A2.
  
* At the cellular level, THs modulate gene expression by binding with the nuclear thyroid receptors (THRA or THRB). Thyroid receptors (TRs) are encoded by two genes, TRα and TRβ, that can undergo alternative processing, generating (among others) TRα1, TRβ1 and TRβ2. Their expression patterns are tissue and cell type-specific. For example, TRα1, is the major isoform expressed in neurons. TRβ is abundant in specific cell types, and in particular, TRβ1 is also expressed in the germinal zones of the cerebral cortex. The balance between TRα and TRβ receptors has been proposed to control the proliferation/differentiation balance.
  
* TRs are structured in  the following regions: 
  + N-terminal region, containing a consitutive activation function
  + a DNA binding domain (DBD) composed of two zinc fingers responsible for DNA binding
  + a region containing the ligand binding domain (LBD)
  + a regions that connect DBD with LBD
  + a region responsible for dimerisation 

* After dimerisation with the retinoid X receptors (RXRs), TRs act as transcription factors and recognise sequences on the DNA, termed thyroid hormone response elements (TRE) :books: [Thyroid hormone receptor localization in target tissues](https://doi.org/10.1530/JOE-17-0708)

* TRs are constitutively bound to corepressors such as NCoR and SMRT, forming complexes containing histone deacetylases, which keeps the chromatin closed. Detrimental effects of hypothyroidism are thought to occur mostly due to this repressive activity. When the ligand binds the receptor, the conformational change releases the corepressors, allowing the TRs to:
  + recruit coactivator complexes, amongst which there may also be proteins that catalyse histone acetylation or arginine methylation. Some coactivators directly recruit the RNA polymerase.
  + repress gene expression by binding to negative TREs or by negative interference on the activity of other transcription factors through protein-to-protein interactions (ex., AP-1, CREB, NF-kB mediated transcription).

* Thyroid hormones may also have non-genomic regulatory effect including regulation of actin polymerisation and ion transport; activation of DIO2 and the Akt/PKB and mTOR pathways; increase of cell proliferation and survival through the binding of integrin α-v β-3 which signalling has been implicated in neocortical development as it upregulates progenitor proliferation.


## Physiological effect of the hormone

* THs are known to have pleiotropic effects by regulating physiological processes associated with growth, metabolic activity and differentiation across nearly all human tissues, including heart, central and autonomic nervous systems, bone and gastrointestinal. :books:  [Thyroid hormone receptor localization in target tissues](https://doi.org/10.1530/JOE-17-0708).

* Typically, the activation of THR implicates an increase in expression of genes that raise the basal metabolic rate and thermogenesis, such as the Na+/K+ ATPase resulting in a ubiquitous higher consumption of oxygen and energy. In the lungs, for example, this results in the stimulation of the respiratory centres to increase oxygenation and perfusion. :books: [Physiology, Thyroid Hormone](https://www.ncbi.nlm.nih.gov/books/NBK500006/).

* Additionally, THs can induce both the anabolism and catabolism of lipids and proteins, depending on the metabolic state, and influence the metabolism of glucose. 

* During development, thyroid hormones act synergistically with growth hormones to stimulate bone formation and growth in children while during foetal development THs role is crucial for the proper formation and maturation of several organs, including the brain. :books: [Transport, Metabolism, and Function of Thyroid Hormones in the Developing Mammalian Brain](https://doi.org/10.3389/fendo.2019.00209)

## Impact on neurodevelopment

* THs are actively transported through tissue barriers, including placenta and BBB. In the developing brain, the most important transporters involved are MCT8 and OATP1C1. The former is predominantly expressed by neurons and preferentially takes up T3, while the latter is expressed in astrocytes and has a higher affinity for T4, which is converted to T3 intracellularly. Low OATP1C1 and high MCT8 mRNA levels were observed in intermediate progenitors.

* The hormones levels in these cells are regulated by the balance of DIO2 and DIO3 in astrocytes and neurons: indeed, while the first expresses DIO2, that is not found in neurons where instead DIO3 helps regulate the levels of active T3 to protect against an excess of it.

* During prenatal development, the embryo relies on the maternal source of THs until mid-gestation, when the thyroid gland is able to supply the fetus with significant Tha levels. In case of fetal THa production deficiencies, maternal THs are able to substitute for them. For this reason, early events such as proliferation of neural progenitors, neuronal migration in the neocortex, hippocampus and MGE depend on the maternal THs, while later events such as neurogenesis, migration and axon growth, dendritic arborization, synaptogenesis and early myelination are under control of both fetal and maternal THs. Even later events (cortex pyramidal cell, hippocampal granule cell and cerebellar granule and Purkinje cell migration, gliogenesis, and myelination) are controlled exclusively by the fetus. 

* Maternal hypothyroxinemia retards fetal glutamatergic neuron migration, resulting in reduced neocortical thickness and decreased neuron number, especially in upper cortical layers. Hypothyroidism causes cell cycle disruption and increases apoptosis, while reducing the progenitor pool and causing defects in neuronal differentiation. Maternal hyperthyroidism has also been associated with neurological problems, with an increased risk of developing epilepsy in the child.
  
* The agonism of THs instead upregulates genes involved in the cell cycle and sustains cell proliferation in the cortex, but this may depend on the specific THs and pathway affected (T4 binding to alpha-v beta-3 upregulates progenitor proliferation, T3 promotes neuronal differentiation). 

* THs can regulate oligodendrocyte lineage development and maturation, inducing more oligodendrocytes from neural stem cells, enhancing the expression of oligo precursors markers and their maturation. In particular, THs, together with the Platelet-Derived Growth Factor (PDGF), induce oligodendrocyte precursor cells to exit the cell cycle (downregulation of cyclin D1 and c-Myc expression, upregulation of CKI) and differentiate. In general, THs are fundamental for myelination, as deficiencies in THs during development result in altered myelin gene expression, impaired oligodendrocyte cell-cycle and consequently in defective myelination. For example, a lack of the transporter MCT8 leads to the Allan-Herndon-Dudley syndrome (AHDS), of which one of the main symptoms is an abnormal white matter content and lower signal of PLP, MBP and MOG.

* On neural stem cells, TH inhibit the proliferation and T3 decreases phosphorilation and DNA binding activity of STAT3. Proliferation is also reduced in SVZ, with a shift towards the oligo lineage. 

* SHH signalling leads to an increase in DIO3 expression, while decreasing DIO2. Through a negative feedback loop, SHH is upregulated by T3 in both adult and fetal brain. Moreover, other factors that are impacted by THs are Reelin, and Ephrins, which are well-characterised morphogens and regulators of migration. Other TH-regulated targets involved in neuronal specification are EMX1, TBR1 and NRGN. 

* Several lines of evidence suggest a direct interplay between TH and the GABAergic system. Developmental hypothyroidism is known to reduce central GABA levels, e.g., by affecting the GABA producing enzyme GAD. Reduced cortical GABA levels were observed in adult-onset hypothyroid patients. Interneuron neurogenesis critically depends on the presence of thyroid hormone, and, interestingly, the PV+ subtype appears to be particularly sensitive towards TH (mostly studies in mice). A decreased number of PV+ interneurons is, however, frequently paralleled by an increase in the SST+ and CR+ subtypes, arguing for more complex changes, such as altered fate decisions.  This puts T3 at the top of a signaling cascade that eventually governs the generation of various interneuronal subtypes. TH has further been suggested to govern migration events, making it especially important for interneuron progenitors that have to traverse long distances before they eventually integrate at their predestined position. Consequently, developmental hypothyroidism interferes with the proper migration of immature interneurons. Critical factors involved in these processes that may be impaired in hypothyroidism are LHX6 and BDNF.

__:books: References:__
- [Transport, Metabolism, and Function of Thyroid Hormones in the Developing Mammalian Brain](https://doi.org/10.3389/fendo.2019.00209)
- [Deficient thyroid hormone transport to the brain leads to impairments in axonal caliber and oligodendroglial development](https://doi.org/10.1016/j.nbd.2021.105567)
- [Thyroid hormone influences brain gene expression programs and behaviors in later generations by altering germ line epigenetic information](https://doi.org/10.1038/s41380-018-0281-4)
- [Local Thyroid Hormone Action in Brain Development](https://doi.org/10.3390/ijms241512352)

## In vitro experiment with iPSC and neural models

* In vitro studies on the impact of high doses of T3 on the development of neural precursors cells to functional cortical neuronal circuits highlighted that this treatment leads to hyperactivation of neurons, as recorded by electrophysiology. Hyperactivation of neurons was related to an increase in action potential frequency in 6-week-old cortical neurons, but other timepoints (W3, W8) were less affected. The activity of T3-treated neurons were similar to W8 controls, suggesting an acceleration of the maturation caused by the treatment. Additionally, the same study identified an increase in basal Ca2+ levels in cortical neurons at 3 weeks but a decrease in 6-weeks old neurons. All timepoints showed increased growth of dendrites, soma cell bodies, abnormal morphology and incremented MAP2 stainings. Whole transcriptome analysis revealed an increased number of reads mapping to intergenic regions upon the addition of T3, pointing to a possible effect of the treatment on splicing and expression of non-coding RNAs. Transcriptionally, W6-treated cells showed more changes than treated cells at other timpoints when looking both at PCA and number of DEGs. Specifically, dysregulated functions included protein synthesis, RNA damage and repair, cell death and survival and RNA post-transcriptional modification while among the canonical pathways there were cellular growth, proliferation and development; cellular stress and intracellular and second messenger signaling, particularly the EIF2 signaling. Comparisons with external data highlighted similarities with studies on ADHD and ASD. :books: [The Potential Role of Thyroid Hormone Therapy in Neural Progenitor Cell Differentiation and Its Impact on Neurodevelopmental Disorders](https://doi.org/10.1007/s12035-023-03751-8).

* __Salas-Lucia et al.__ studied Allan-Herndon-Dudley syndrome through the use of brain organoids derived from hiPSCs from patients with MTC8 mutations. These organoids showed altered early neurodevelopment, with smaller neural rosettes and thinner cortical units; defective T3 transport in developing neural cells; decreased inducibility of TH-regulated genes and an overall downregulation of genes involved in cerebral cortex development. These phenotypes were rescued by the use of TH analogs: 3,5-diiodothyropropionic acid (DITPA) and 3,3′,5-triiodothyroacetic acid (TRIAC). Throught RNA-seq analysis, the authors identified some genes downregulated in MCT8-deficient organoids such as SATB2, KNCF1, CRYM, MMP2, CFB and HOXB3. GSEA revealed alterations in pathways rlevant for cerebral cortex development, including MAPK, cAMP, Wnt, mTOR, TNF, TGF-beta, Hippo and NFkB. :books: [Impaired T3 uptake and action in MCT8-deficient cerebral organoids underlie Allan-Herndon-Dudley syndrome](https://doi.org/10.1172/jci.insight.174645)

* Another study from __Gruffender et al.__ (https://doi.org/10.1038/s41598-024-59533-2) focuses on the spatiotemporal expression of the thyroid hormone transporter MCT8 (monocarboxylate transporter 8) and THRA (thyroid hormone receptor alpha) mRNA in human cerebral organoids (hCOs), differentiated with Lancaster protocol. hCOs were treated with T3 for 48 hours at different developmental stages, and the expression of known T3-responsive genes was analyzed via RT-qPCR to assess their responsiveness. Immunostaining showed MCT8 membrane expression in neuronal progenitor cell types, including early neuroepithelial cells, radial glia cells (RGCs), intermediate progenitors (IPCs), and outer RGCs (oRGCs). Enhanced MCT8 staining intensity was noted on the apical surface of the neuroepithelium lining the lumen. FISH detected THRA mRNA expression as early as in the neuroepithelium before the onset of neurogenesis. THRA mRNA expression remained low in the ventricular zone (VZ), increased in the subventricular zone (SVZ), and strong expression was observed in excitatory neurons of the cortical plate (CP). This suggests a marked increase in THRA mRNA expression along the neuronal differentiation trajectory. hCOs showed a robust up-regulation of known T3-responsive genes (KLF9, DBP, FLYWCH2, CADM2) following T3 treatment. The induction of these genes was lower in early-stage hCOs compared to later-stage hCOs, which had a more extensive neuronal population with higher THRA mRNA levels. :books: [Spatiotemporal expression of thyroid hormone transporter MCT8 and THRA mRNA in human cerebral organoids recapitulating first trimester cortex development](https://doi.org/10.1038/s41598-024-59533-2)

## Epidemiological data

* Many clinical and epidemiological studies highlights the importance of thyroid hormone signaling during fetal development.

* The human fetal thyroid gland only becomes fully functional in the second trimester of gestation, specifically between the 12th and 14th weeks. Before this, especially during the first half of pregnancy, the fetus is almost entirely dependent on maternal thyroid hormones, particularly thyroxine (T4). This maternal T4 is crucial for neuronal proliferation and migration as early as the 9th week of gestation. The maternal demand for TH increases significantly (by approximately 20% to 50%) during pregnancy to maintain a euthyroid state. This increase is influenced by human chorionic gonadotropin (hCG), which acts as a TSH receptor agonist, and rising levels of thyroxine-binding globulin (TBG). The requirement for levothyroxine (L-T4) can increase as early as the fifth week of gestation, with an average increase of 47% in the first half of pregnancy, plateauing by week 16. It is suggested that women with hypothyroidism should increase their levothyroxine dose by approximately 30% as soon as pregnancy is confirmed. 

* Severe developmental hypothyroidism is profoundly detrimental to neurodevelopment. Maternal hypothyroidism is associated with severe adverse pregnancy outcomes such as miscarriage, preterm delivery, and preeclampsia. It leads to significant fetal damage, including impaired nerve cell differentiation, inadequate CNS development, increased risk of perinatal defects, low birth weight, and impacts on motor and cognitive development. Even mild maternal thyroid insufficiency can be associated with indicators of intellectual disability in offspring. A meta-analysis (https://doi.org/10.1111/cen.13550) showed a significant association (SCH: odds ratio [OR] 2.14; hypothyroxinemia: OR 1.63). :books: 
However, the effect on Autism Spectrum Disorder is unclear, and no association has been found with attention deficit hyperactivity disorder (ADHD).

* SCH is also linked to obstetric complications such as placental abruption, increased risk of preterm birth, gestational hypertension, and severe preeclampsia. Preterm birth itself can lead to neuropsychological dysfunction, including reduced IQ and ADHD. Even mild thyroid dysfunction in preterm infants (consistently in the top decile of gestationally age-adjusted TSH or bottom decile of age-adjusted thyroxine) is associated with lower cognitive and motor scores at 2 years of age.

* The consequences of TH insufficiency depend on the precise developmental timing of the deficiency:
  + Early pregnancy: Can cause problems in visual attention, visual processing (acuity, strabismus), and gross motor skills.
  + Later in pregnancy: Can lead to an additional risk of poor visual and visuospatial skills, slower response times, and fine motor deficits.
  + After birth: Primarily affects language and memory skills.

__:books: References__

- [Developmental thyroid hormone disruption: Prevalence, environmental contaminants and neurodevelopmental consequences](https://doi.org/10.1016/j.neuro.2011.11.005)
- [Evaluation of maternal thyroid function during pregnancy: the importance of using gestational age-specific reference intervals](https://doi.org/10.1530/eje-07-0249)
- [Hypothyroidism and chronic autoimmune thyroiditis in the pregnant state: maternal aspects](https://doi.org/10.1016/j.beem.2004.03.006)
- [Is thyroid inadequacy during gestation a risk factor for adverse pregnancy and developmental outcomes?](https://doi.org/10.1089/thy.2005.15.60)
- [Maternal thyroid hormone insufficiency during pregnancy and risk of neurodevelopmental disorders in offspring: A systematic review and meta-analysis](https://doi.org/10.1111/cen.13550)
- [Clinical and subclinical maternal hypothyroidism and their effects on neurodevelopment, behavior and cognition](https://doi.org/10.20945/2359-3997000000201)
- [Thyroid function in preterm infants and neurodevelopment at 2 years](https://fn.bmj.com/content/105/5/504.abstract)
- [Timing and magnitude of increases in levothyroxine requirements during pregnancy in women with hypothyroidism](https://doi.org/10.1056/nejmoa040079)
- [Timing of thyroid hormone action in the developing brain: clinical observations and experimental findings](https://doi.org/10.1111/j.1365-2826.2004.01243.x)
- [Maternal Thyroid Function in Early Pregnancy and Child Neurodevelopmental Disorders: A Danish Nationwide Case-Cohort Study](https://doi.org/10.1089/thy.2017.0425)


## Evidence of EDC disruption of the pathway

* From Caporale et al, the thyroid pathway is the most impacted by the exposure to MIX N. Cell proliferation-related genes (CDC20B and histone-related genes) as well as NEUROG1, which was shown to act as a negative regulator of neocortical neurogenesis, were significantly up-regulated by MIX N and down-regulated by T3. A subset of genes related to neural differentiation (DCX, SYP, MAP2, and RBFOX3) were down-regulated, whereas cell proliferation genes (MKI67, CCNB1, CDC20, and HMGB2) were up-regulated by MIX N. Gene sets related to cell proliferation were among the top positively enriched sets in GSEA. When the cellular phenotypic analysis was applied to T3-treated samples, it displayed the opposite effect of MIX N, both for proliferation and maturation, whereas BPA partially recapitulated MIX N effects only for progenitor proliferation. In Xenopus MIX N decreases the expression of the TH dependent TFs klf9 at 1X concentration, whereas it decreased expression of the TH transporter oatp1c1 at 10X and increased the expression of the TH transporter mct8 at 1000X. In zebrafish exposure to the mix led to significantly decreased expression of the thyroid hormone receptors thra (P = 0.0006) and thrb (P = 0.0010) at 100X concentration. :books: [From cohorts to molecules: Adverse impacts of endocrine disrupting mixtures]( 10.1126/science.abe8244)

* Numerous environmental agents, including polychlorinated biphenyls (PCBs) and polybrominated diphenyl ethers (PBDEs), can disrupt TH signaling, contributing to negative effects on brain function by acting at the receptor level. :books: [Regulation of Brain Development by Thyroid Hormone and its Modulation by Environmental Chemicals](https://doi.org/10.1507/endocrj.KR-69); [Developmental thyroid hormone disruption: Prevalence, environmental contaminants and neurodevelopmental consequences](https://doi.org/10.1016/j.neuro.2011.11.005)

## Crosstalk with other pathways

* Retinoic acids share carrier proteins with THs and can increase MCT8 expression to increase their import. :books: [Retinoic Acid Induces Expression of the Thyroid Hormone Transporter, Monocarboxylate Transporter 8 (Mct8)](https://doi.org/10.1074/jbc.m110.123158)

* Nuclear liver X receptor beta interacts with TH signaling in regulating cortical layering by acting on the receptor ApoER2. :books: [Liver X receptor β and thyroid hormone receptor α in brain cortical layering](https://doi.org/10.1073/pnas.1006162107)
