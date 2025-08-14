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

## Compounds tested

- __Agonist__: GW3965 (0.2 $\mu M$)
- __Antagonist__: SR9238 (0.2 $\mu M$)

---------
  
## Results from reference compounds

### Summary scheme

[Hormonal Card - LiverX](Schemes/Cards_Hormones_LiverX.jpg)

### Result Description

* [Gene Signatures](1.GeneSignatures.md#Liver-X)
* [Differential Expression Analysis - Top Genes](2.DEATopGenes.md#Liver-X)
* [Differential Expression Analysis - Functional Enrichment](3.DEAFunctional.md#Liver-X)
* [WGCNA](4.WGCNA.md#Liver-X)

-----

## General pathway information

* Liver X Receptors (LXRs) function as heterodimers with Retinoid X Receptor (RXR) and bind to DNA at response elements known as DR4s (direct repeats of the half-site sequence 5′-G/AGGTCA-3′, separated by four nucleotides).
* These DR4 elements are also utilized by thyroid hormone receptors (TRs). Given that both TRs and LXRs form heterodimers with RXRs and bind to DR-4 elements with identical geometry and polarity, significant crosstalk between these two receptors has been reported, particularly concerning lipid metabolism-related genes.
* This highlights a strong relationship between thyroid hormone and LXR signaling.
* :books: [Crosstalk of thyroid hormone receptor and liver X receptor in lipid metabolism and beyond](https://doi.org/10.1507/endocrj.ej11-0114); [Liver X receptor β controls thyroid hormone feedback in the brain and regulates browning of subcutaneous white adipose tissue](https://doi.org/10.1073/pnas.1519358112)

------

## Physiological effect of the hormone

### Macrophages/Liver

* Increased levels of cholesterol in peripheral cells, such as macrophages, lead to LXR-dependent upregulation of genes encoding proteins like ABCA1 and ABCG1, which facilitate the transport of intracellular cholesterol out of cells to HDL particles. Concurrently, regulation of IDOL/*MYLIP* leads to downregulation of LDL receptors and decreased cholesterol uptake.

* In the liver, LXR-dependent regulation of genes encoding ABCG5 and ABCG8 promotes the excretion of hepatic cholesterol into the bile. Hepatic regulation of the gene encoding SREBP1c, as well as genes for other enzymes involved in fatty acid synthesis, increases triglyceride levels and promotes the secretion of VLDL.

### Liver X Role in Controlling Inflammation

* Liver X receptor activation promotes polyunsaturated fatty acid (PUFA) synthesis :books: [Liver X Receptor Activation Promotes Polyunsaturated Fatty Acid Synthesis in Macrophages: Relevance in the Context of Atherosclerosis](https://doi.org/10.1161/atvbaha.115.305539).
* PUFA deficiency has been identified as a risk factor for schizophrenia. Studies in mice have shown that gestational and early postnatal dietary deprivation of arachidonic acid (AA) and docosahexaenoic acid (DHA) elicited schizophrenia-like phenotypes in adult offspring :books: [Polyunsaturated fatty acid deficiency during neurodevelopment in mice models the prodromal state of schizophrenia through epigenetic changes in nuclear receptor genes](https://doi.org/10.1038/tp.2017.182).
* Long-chain PUFAs are precursors of various metabolites with diverse effects on inflammation and neuron outgrowth :books: [Perinatal Dietary Polyunsaturated Fatty Acids in Brain Development, Role in Neurodevelopmental Disorders](https://doi.org/10.3390/nu13041185).
* Human studies indicate a positive link between maternal intake of PUFAs during pregnancy and child neurodevelopment :books: [Polyunsaturated fatty acids and child neurodevelopment among a population exposed to DDT: a cohort study](https://doi.org/10.1186/s12940-019-0456-8).

### Liver X and Cholesterol Metabolism

* LXRs act as cholesterol sensors: when cellular oxysterols accumulate due to increasing cholesterol concentrations, LXRs induce the transcription of genes that protect cells from cholesterol overload. LXR activation regulates bile acid synthesis and metabolism/excretion, reverse cholesterol transport (RCT), cholesterol biosynthesis, and cholesterol absorption/excretion in the intestine.  In general, the LXR pathway is crucial for cholesterol efflux and metabolism; however, it is not the primary pathway leading to/regulating *de novo* cholesterol biosynthesis. The activity of fatty acid and cholesterol synthesis pathways is increased in LXR-deficient macrophages. :books: [Liver X receptor in cholesterol metabolism](https://doi.org/10.1677/JOE-09-0271).

* Regarding cholesterol biosynthesis, Wang et al. (2008) demonstrated that LXR $\alpha$ negatively regulates two genes, squalene synthase (FDFT1) and lanosterol 14 $\alpha$-demethylase (CYP51A1), which encode key enzymes in the cholesterol biosynthesis pathway. LXREs that confer LXR-mediated repression were identified in these genes. Based on these observations, it was proposed that LXR $\alpha$ plays an important role in the suppression of cholesterol biosynthesis. :books: [Regulation of Cholesterologenesis by the Oxysterol Receptor, LXRα](https://doi.org/10.1074/jbc.m804808200).

### Astrocytes/Neurons Cholesterol Shuttle

* Cholesterol is essential for proper neuronal function and integrity. It is synthesized by astrocytes and transported to neurons via vesicles, a process crucial for astrocyte-neuron crosstalk. When 24-hydroxycholesterol, a product of cholesterol catabolism, is released by neurons, it enters astrocytes to bind Liver X receptors.
* These receptors induce the expression of ApoE and ABCA1/G1, which in turn contribute to the formation of cholesterol-laden vesicles that are shuttled back to neurons. This process ensures the proper cholesterol supply to neurons. Without functional LXR pathways, cholesterol can accumulate in the brain, potentially contributing to neurological disorders such as Alzheimer's disease. :books: [Liver X receptors in the central nervous system: From lipid homeostasis to neuronal degeneration](https://doi.org/10.1073/pnas.172510899);  [Regulation of Cholesterol Homeostasis by the Liver X Receptors in the Central Nervous System](https://doi.org/10.1210/mend.16.6.0835)
* Overall, Liver X receptors induce the expression of genes important for cholesterol efflux from astrocytes to neurons (e.g., APOE, ABCA1, ABCG1).
* Comprehensive review on LXR effects on different cell types. :books: [Liver X Receptor Regulation of Glial Cell Functions in the CNS](https://doi.org/10.3390/biomedicines10092165000)

### Other Roles in Astrocytes

* The water channel aquaporin 4 (AQP4) is an LXR-regulated gene. AQP4 is expressed in astrocytic end feet and ependymal cells, where it regulates cerebrospinal fluid (CSF) homeostasis by facilitating fluid clearance from the parenchyma.
* In LXR $\alpha\beta^{-/-}$ mice, severe defects in CSF maintenance result in occlusion of the lateral ventricles and degeneration of the choroid plexus. :books: [Liver X Receptor Regulation of Glial Cell Functions in the CNS](https://doi.org/10.3390/biomedicines10092165)

### Liver X and Microglia

LXR agonists like GW3965 (the one we used in our study) inhibit the production of NO, IL-1 $\beta$, IL-6, and monocyte chemoattractant protein-1 (MCP-1) in microglia and astrocytes, thereby inhibiting stress-induced neuroinflammation that can potentially lead to cognitive disorders. This suggests a neuroprotective effect. :books: [Activation of liver X receptors prevents emotional and cognitive dysfunction by suppressing microglial M1-polarization and restoring synaptic plasticity in the hippocampus of mice](https://doi.org/10.1016/j.bbi.2021.02.026)

### Liver X and Oligodendrocytes

* LXRs differentially affect the mRNA levels of myelin genes in myelin-rich tissues such as the spinal cord, corpus callosum, optic nerve, and cerebellum. :books: [Liver X Receptor Regulation of Glial Cell Functions in the CNS](https://doi.org/10.3390/biomedicines10092165). 
* LXR ligands influence the mRNA levels of myelin-related genes, proteolipid protein (PLP), and myelin basic protein (MBP). :books: [Liver X Receptors differentially modulate central myelin gene mRNA levels in a region-, age- and isoform-specific manner](https://doi.org/10.1016/j.jsbmb.2016.02.032).
* Activation of LXRs also promotes oligodendrocyte maturation.

-----

## Impact on Neurodevelopment

Most studies on LXR's impact on neurodevelopment have been performed on mice.

1. Studies on LXR $\beta^{-/-}$ knockout mice have shown that the expression of Liver X Receptor $\beta$ is essential for the formation of superficial cortical layers and the migration of later-born neurons. :books: [Expression of liver X receptor β is essential for formation of superficial cortical layers and migration of later-born neurons](https://doi.org/10.1073/pnas.0806974105).

      * Key findings:
          * LXR $\beta$ contributes to cortical lamination, especially for upper layers.
          * The cortical plate was thinner in the LXR $\beta^{-/-}$ mice.
          * LXR $\beta$ is essential for the development and migration of upper layer neurons, while the formation of deep layers is largely independent of LXR $\beta$ signaling.
          * The processes of radial glial cells in LXR $\beta^{-/-}$ mice appeared truncated or less organized into radial formations, contributing to the abnormal radial migration of neurons in the upper cortical layers.
          * LXR $\beta$'s effect on radial glial cells is exerted through ApoE signaling. In ApoER2 mutant mice, early-generated layers form almost normally, but the formation of superficial, late-generated layers is severely altered, identically to LXR $\beta$ KO.
          * This phenotype (altered neuron migration) is rescued postnatally by thyroid hormone (TH). TH levels increase after mouse birth, inducing the expression of APOER2 (APOE and reelin receptor), compensating for the loss of LXR $\beta$. TH stimulates the retarded neurons to migrate to their correct position, rescuing the phenotype.

2. Liver X Receptor promotes neurogenesis in the hippocampus (dentate gyrus) and midbrain during neurodevelopment. Loss of LXR can induce autism-like behaviors in mice and a loss of progenitor cells. :books: [Liver X receptor β regulates the development of the dentate gyrus and autistic-like behavior in the mouse](https://doi.org/10.1073/pnas.1800184115); [Brain endogenous liver X receptor ligands selectively promote midbrain neurogenesis](https://doi.org/10.1038/nchembio.1156).

3. Liver X Receptors and oxysterols promote ventral midbrain neurogenesis *in vivo* and in human embryonic stem cells. :books: [Liver X Receptors and Oxysterols Promote Ventral Midbrain Neurogenesis In Vivo and in Human Embryonic Stem Cells](https://doi.org/10.1016/j.stem.2009.08.019).

4. After neurogenesis and migration are completed, most radial glial (RG) cells transform into astrocytes. Liver X receptor beta delays the transformation of radial glial cells into astrocytes during mouse cerebral cortical development. :books: [Liver X receptor β delays transformation of radial glial cells into astrocytes during mouse cerebral cortical development](https://doi.org/10.1016/j.neuint.2014.03.009).

5. LXRs also have effects on cerebellar development and function. :books: [Liver X receptor agonist treatment promotes the migration of granule neurons during cerebellar development](https://doi.org/10.1111/j.1471-4159.2010.07053.x)
   
-----

## *In vitro* experiment with iPSC and neural models

* Studies using human stem-cell models have uncovered a clear role for LXRs in neural differentiation. 

* In a human iPSC-derived cortical model of the 15q11.2 deletion (affecting the autism risk gene CYFIP1), researchers found that disrupted LXR signaling perturbed neurogenesis. Specifically, CYFIP1 loss caused dysregulated cholesterol metabolism and altered levels of 24S,25-epoxycholesterol (a potent endogenous LXR ligand), leading to premature cortical neuron differentiation. Application of 24S,25-epoxycholesterol to control neural progenitors mimicked this effect (more neurons, fewer progenitors), whereas deletion of LXR $\beta$ abolished it. Conversely, blocking LXR (by compound deletion of LXR $\beta$ in the CYFIP1-deficient cells) restored normal timing of neurogenesis. Thus, LXR $\beta$-mediated oxysterol signaling appears to fine-tune the switch from progenitor proliferation to neuronal differentiation in developing cortex. :books: [Impaired oxysterol-liver X receptor signaling underlies aberrant cortical neurogenesis in a stem cell model of neurodevelopmental disorder](https://doi.org/10.1016/j.celrep.2024.113946)

* LXR activation has been shown to drive dopaminergic neuron differentiation in human pluripotent cultures. LXRs and oxysterols were found to promote ventral midbrain neurogenesis *in vitro* in human embryonic stem cells. This highlights a direct role for LXRs in regulating the early stages of human neuronal development. :books: [Liver X Receptors and Oxysterols Promote Ventral Midbrain Neurogenesis In Vivo and in Human Embryonic Stem Cells](https://doi.org/10.1016/j.stem.2009.08.019)

* Other human iPSC studies similarly implicate LXR-driven cholesterol metabolism in neural development. For instance, human neural progenitors differentiated in the presence of retinoic acid (which activates LXR–ABCA1 cholesterol efflux pathways) showed boosted neuronal maturation. :books:  [Human iPSC-Based Models for the Development of Therapeutics Targeting Neurodegenerative Lysosomal Storage Diseases](https://doi.org/10.3389/fmolb.2020.00224).

* In neural stem cell models of Niemann-Pick disease, restoring LXR $\beta$ -regulated cholesterol homeostasis (by VPA or other agents) rescued defects in self-renewal and neurogenesis. :books: [Generation of patient specific human neural stem cells from Niemann-Pick disease type C patient-derived fibroblasts](https://doi.org/10.18632/oncotarget.19976).

* Recent comprehensive screen using fetal human neurospheres exposed to various hormone agonists found that LXR activation altered multiple neurodevelopmental processes. :books: [Nuclear hormone receptors control fundamental processes of human fetal neurodevelopment: Basis for endocrine disruption assessment](https://doi.org/10.1016/j.envint.2025.109400).
    
-----

## Epidemiological data

* Epidemiological studies in humans, while not directly measuring LXR activity, often reveal associations between factors that influence LXR pathways (e.g., lipid metabolism, exposure to certain chemicals) and neurodevelopmental outcomes or neurological diseases. Human studies increasingly link cholesterol/LXR pathway disturbances and endocrine exposures to neurodevelopmental outcomes.

* For Autism Spectrum Disorder, several clinical cohorts report altered lipid profiles. In one case-control study, children with ASD had significantly higher total cholesterol, LDL cholesterol and triglycerides than controls, while levels of 27-hydroxycholesterol (an LXR ligand) were markedly lower; these findings suggest dysregulated oxysterol–LXR signaling in ASD. :books: [Investigation of Liver X Receptor Gene Variants and Oxysterol Dysregulation in Autism Spectrum Disorder](https://doi.org/10.3390/children11050551). Intriguingly, other analyses find the opposite: a French-Canadian cohort showed hypocholesterolemia was about four times more prevalent in ASD patients than expected, and very low cholesterol (<25th percentile) strongly predicted ASD diagnosis. In that study, low cholesterol in ASD was associated with higher rates of intellectual disability and anxiety :books: [Implication of hypocholesterolemia in autism spectrum disorder and its associated comorbidities: A retrospective case–control study](https://doi.org/10.1002/aur.2183) Together, such data imply that deviations in cholesterol homeostasis, whether abnormally high or low, are linked to ASD risk, perhaps via LXR-mediated pathways. Genetic studies further hint at this link: polymorphisms in the LXR $\beta$ gene (NR1H2) have been assessed in ASD, and although SNP frequencies did not differ from controls, the ASD group showed altered LXR ligand (oxysterol) and lipid levels, consistent with a metabolic phenotype in ASD. :books: [Investigation of Liver X Receptor Gene Variants and Oxysterol Dysregulation in Autism Spectrum Disorder](https://doi.org/10.3390/children11050551)

* Lipid metabolism is also associated with ADHD. In a large Chinese pediatric sample, higher serum cholesterol and LDL were significantly associated with increased odds of ADHD, even after adjusting for age, sex, BMI and other factors. Hypertriglyceridemia and low HDL showed similar patterns. (These analyses, stratified by obesity status, remained significant in both obese and non-obese subgroups.) In other words, dyslipidemia, potentially reflecting perturbed LXR-regulated lipid handling, correlates with ADHD diagnoses in children. :books: [The Relationship Between Blood Lipid and Attention-Deficit/Hyperactivity Disorder (ADHD) in an Obese Population of Chinese Children: An Obesity-Stratified Cross-Sectional Study](https://doi.org/10.2147/IJGM.S333247)

* Regarding environmental disruptors, human epidemiology links prenatal and childhood chemical exposures to neurodevelopmental disorders, often via endocrine or metabolic disruption. Some flame retardants or phthalates are known or suspected to disrupt LXR activity. Several birth-cohort studies have found that prenatal phthalate exposure (e.g. Di-2-ethylhexyl phthalate (DEHP) and related metabolites) is associated with later ASD symptoms or traits. For example, higher maternal urinary DEHP during pregnancy predicted more ASD behaviors in the child at ages 2–4 years and it is known that DEHP can act via activating the LXR signaling pathway for what concerns lipid metabolism disorders . In adolescents, higher current levels of phthalate metabolites (and certain chlorinated phenols) were linked to a higher risk of ADHD-like behaviors; a two-fold increase in anti-androgenic phthalate metabolites raised the odds of ADHD-related behavior problems by around 34%. :books: [Prenatal environmental risk factors for autism spectrum disorder and their potential mechanisms](https://doi.org/10.1186/s12916-024-03617-3); [Di-2-ethylhexyl phthalate (DEHP) induced lipid metabolism disorder in liver via activating the LXR/SREBP-1c/PPARα/γ and NF-κB signaling pathway](https://doi.org/10.1016/j.fct.2022.113119); [Association of Exposure to Endocrine-Disrupting Chemicals During Adolescence With Attention-Deficit/Hyperactivity Disorder–Related Behaviors](https://doi.org/10.1001/jamanetworkopen.2020.15041)

* More broadly, systematic reviews report that prenatal exposure to various endocrine-disrupting chemicals (EDCs) – including phthalates, heavy metals, and persistent organic pollutants – is associated with deficits in cognitive, motor, and language development. This 2024 meta-analysis (26,000+ mother-child pairs) found significant negative impacts of prenatal EDC exposure on early neurobehavior: for instance, prenatal phthalates were linked to poorer motor development in infants, and other EDCs (PFAS, metals) were tied to language or cognitive delays. These effects are thought to involve hormonal and metabolic pathways (potentially including LXR/RXR disruption) during critical brain growth periods. :books: [Prenatal endocrine-disrupting chemicals exposure and impact on offspring neurodevelopment: A systematic review and meta-analysis](https://doi.org/10.1016/j.neuro.2024.07.006)

-----

## Evidence of EDC disruption of the pathway

Several studies have identified EDCs that disrupt the LXR pathway, with potential implications for metabolic and neurological health.

1. Bisphenol A (BPA):

      * In mice liver data the effect of BPA on the expression of *de novo* lipogenesis followed a nonmonotonic dose-response curve, with more important effects at lower doses than at the higher dose. In addition to lipogenic enzymes (*Acc*, *Fasn*, *Scd1*), the expression of transcription factors such as Liver X Receptor, the sterol regulatory element binding protein-1c, and the carbohydrate responsive element binding protein that govern the expression of lipogenic genes also followed a nonmonotonic dose-response curve in response to BPA. This suggests that even low doses of BPA can interfere with LXR-mediated metabolic regulation. :books: [Low Doses of Bisphenol A Induce Gene Expression Related to Lipid Synthesis and Trigger Triglyceride Accumulation in Adult Mouse Liver](https://doi.org/10.1002/hep.24685)
      * While direct human LXR binding data for BPA is still developing, the widespread human exposure to BPA and its established endocrine-disrupting properties raise concerns about its potential to perturb LXR-regulated metabolic and developmental processes in humans.

2. Flame Retardants:

      * The flame retardants triphenyl phosphate (TPHP) and 2-ethylhexyl diphenyl phosphate (EHDPP) were identified in house dust samples and demonstrated their ability to antagonize LXRs. The potency of TPHP was similar to that of the LXR-antagonist SR9238. TPHP could also inhibit cholesterol efflux and promote foam cell formation in RAW264.7 macrophages and mouse peritoneal macrophages and significantly promoted atherosclerotic lesion formation in the APOE $^{-/-}$ mouse model.  :books: [Screening of House Dust from Chinese Homes for Chemicals with Liver X Receptors Binding Activities and Characterization of Atherosclerotic Activity Using an in Vitro Macrophage Cell Line and ApoE−/− Mice](https://doi.org/10.1289/EHP5039).
      * Human exposure to these flame retardants is pervasive, raising concerns about their potential to disrupt LXR function and contribute to dyslipidemia and metabolic disorders in humans.

3. Environmental pollutants directly affect the liver X receptor alpha activity: Kinetic and thermodynamic characterization of binding described in :books: [Environmental pollutants directly affect the liver X receptor alpha activity: Kinetic and thermodynamic characterization of binding](https://doi.org/10.1016/j.jsbmb.2015.04.011).

4. Screening of chemicals with binding activities of liver X receptors from reclaimed waters described in :books: [Screening of chemicals with binding activities of liver X receptors from reclaimed waters](https://doi.org/10.1016/j.scitotenv.2020.136570).

-----

## Crosstalk with other pathways

* LXRs form obligate heterodimers with retinoid X receptors (RXRs) to regulate target genes. In differentiating stem cells, RXR partners orchestrate neuronal gene programs: for example, chromatin profiling shows LXR:RXR (and RAR:RXR) complexes drive distinct gene networks during ESC→neuron differentiation. Thus, LXR activation often overlaps with retinoic acid (RAR/RXR) signaling in neurogenesis. :books: [RXR heterodimers orchestrate transcriptional control of neurogenesis and cell fate specification](https://doi.org/10.1016/j.mce.2017.07.033)

* Beyond RXR, LXRs intersect with other nuclear receptors. In mouse cortex, LXR $\beta$ and thyroid hormone receptor $\alpha$ (TR $\alpha$) share DNA response elements and co-regulate genes critical for cortical layering. Loss of LXR $\beta$ delays neuronal migration, but rising thyroid hormone postnatally can compensate via TR $\alpha$ , suggesting dynamic LXR–TR crosstalk in cortical development. Several types of crosstalk between TRs and LXRs have been identified and crosstalk has also been observed in other physiological systems such as central nervous system rather than lipid metabolism. :books: [RXR heterodimers orchestrate transcriptional control of neurogenesis and cell fate specification](https://doi.org/10.1073/pnas.1006162107); [Crosstalk of thyroid hormone receptor and liver X receptor in lipid metabolism and beyond](https://doi.org/10.1507/endocrj.ej11-0114); [Liver X receptor β and thyroid hormone receptor α in brain cortical layering](https://doi.org/10.1073/pnas.1006162107)

* In metabolic tissues LXR activation has been shown to directly modulate glucocorticoid receptor (GR) signaling. LXRs bind to glucocorticoid response elements (GREs) and can inhibit GR-driven transcription :books: [Liver X Receptors Regulate the Transcriptional Activity of the Glucocorticoid Receptor: Implications for the Carbohydrate Metabolism](https://doi.org/10.1371/journal.pone.0026751). For example, GW3965 (our LXR agonist) attenuated dexamethasone-induced gene expression of gluconeogenic enzymes by recruiting LXR/RXR to GREs and displacing GR. Although this work was in liver cells, it shows that LXR activation can repress glucocorticoid pathways via competitive DNA binding, a mechanism that could extend to neural progenitors exposed to stress hormones.

* Estrogen signaling is also modulated by LXR. In rat sensory neurons, LXR agonism counteracted ER$\alpha$ -driven inflammatory gene expression :books: [Cholesterol-dependent LXR transcription factor activity represses pronociceptive effects of estrogen in sensory neurons and pain induced by myelin basic protein fragments](https://doi.org/10.1016/j.bbih.2024.100757). For instance, 24S,25-epoxycholesterol (an LXR ligand) suppressed ER$\alpha$ -mediated IL-6 and $Cav\alpha2\delta$ induction in DRG neurons, reducing estrogen’s pro-nociceptive effects. This implies that LXR activation can dampen estrogen receptor pathways in neurons, which may influence neuronal maturation or responses to estrogenic endocrine disruptors.

* Importantly, these hormone cross-regulations influence neural stem cell behavior. LXR agonists promote progenitor proliferation: in neural cultures, GW3965 and LXR623 significantly stimulate NPC expansion (via MEK/ERK signaling) in wild-type but not LXR-knockout cells :books: [LXR agonists promote the proliferation of neural progenitor cells through MEK-ERK pathway](https://doi.org/10.1016/j.bbrc.2016.12.163). 

* In human models, the recent comprehensive screen using fetal human neurospheres previously cited  show LXR (and RAR, RXR, GR) agonists changed NPC proliferation and differentiation, and RNA-seq revealed activation of conserved pathways (Notch, Wnt) and extensive receptor crosstalk. Together these data suggest that LXR signaling intermingles with other hormonal programs to shape gene expression, progenitor proliferation, and neuronal differentiation. This crosstalk also affects vulnerability to chemicals. This human neurosphere study highlights that perturbing LXR (or GR, RXR, etc.) at physiological concentrations disrupts key neurodevelopmental processes. In other words, endocrine disruptors that target one receptor may indirectly influence LXR-dependent pathways (and vice versa), altering cortical neurogenesis. :books: [Nuclear hormone receptors control fundamental processes of human fetal neurodevelopment: Basis for endocrine disruption assessment](https://doi.org/10.1016/j.envint.2025.109400)
