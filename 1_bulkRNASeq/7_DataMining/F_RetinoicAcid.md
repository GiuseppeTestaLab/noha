## Table of Content

* [Compounds tested](#Compounds-tested)
* [Results from reference compounds](#Results-from-reference-compounds)
* [General pathway information](#General-pathway-information)
* [Physiological effect of the hormone](#Physiological-effect-of-the-hormone)
* [Impact on neurodevelopment](#Impact-on-neurodevelopment)
* [In vitro experiment with iPSC and neural models](#in-vitro-experiment-with-iPSC-and-neural-models)
* [Epidemiological data](#Epidemiological-data)
* [Evidence of EDC disruption of the pathway](#Evidence-of-EDC-disruption-of-the-pathway)
* [Crosstalk with other pathways](#Crosstalk-with-other-pathways)  

------------------
## Compounds tested:

* __Agonist__: [ATRA](https://www.tocris.com/products/retinoic-acid_0695?gclid=Cj0KCQjwrfymBhCTARIsADXTabk7--arS5mFDKZf5EJp8XihM-CF8RASVlXG7yGC2zUdiZYomcyuz5gaAvdAEALw_wcB&gclsrc=aw.ds), selective pan-RAR agonist (RARα, β, or γ), does not bind RXR. 

* __Antagonist__: [AGN193109](https://www.tocris.com/products/agn-193109_5758), high affinity pan-retinoic acid receptor (RAR) antagonist. Exhibits no significant affinity for retinoic X receptors.

------------

## Results from reference compounds

### Summary scheme

[Hormonal Card - Retinoic Acid](Schemes/Cards_Hormones_Retinoids.jpg)

### Result Description
* [Gene Signatures](1.GeneSignatures.md#retinoic-acid)
* [Differential Expression Analysis - Top Genes](2.DEATopGenes.md#retinoic-acid)
* [Differential Expression Analysis - Functional Enrichment](3.DEAFunctional.md#retinoic-acid)
* [WGCNA](4.WGCNA.md#retinoic)


------------

## General pathway information

* Retinoic acid (RA) is derived from retinol (vitamin A) as a metabolic product. This small lipophilic molecule is produced in a two-step oxidative reaction involving RDH10, the main retinol dehydrogenase present in the embryo, and retinaldehyde dehydrogenases (RALDH1, 2, 3). RA exists in several isomeric forms including all-trans-RA, 9-cis-RA and 13-cis-RA; however, all-trans-RA (ATRA) is the primary ligand during development. :books: [Retinoic acid controls early neurogenesis in the developing mouse cerebral cortex](https://pubmed.ncbi.nlm.nih.gov/28790015/)

* RA acts as an agonist of two nuclear receptor families that bind DNA and directly regulate transcription. These families are (i) the RA receptors (RAR) and (ii) the retinoid X receptors, (RXR). The RARs are highly conserved in vertebrates and are primarily activated by all-trans-RA (ATRA). By contrast, the RXRs are activated only by 9-cis-RA, a stereoisomer of ATRA that is detected only when vitamin A is in excess. Of notice, the heterodimerization of RXR with RAR and other nuclear receptors such as LXR, FXR, or PPAR is mutually competitive, and atRA signaling not only triggers the transcription of its target genes, but also competitively suppresses the transcription of others.

* Upon ATRA binding to its receptor, the RAR will form a heterodimer with RXR, increasing   their joint affinity and selectivity for retinoic acid response elements (RAREs) that contain direct repeats (DR) of 5'-AGGTCA-3' separated by one to five base pairs (termed DR1-DR5). DRs (DR1-5) determine RA-activated RAR-RXR complex target gene expression. DR5 and DR2 activate transcription, whereas DR1-containing genes display transcriptional repression. :books: [LXR, a nuclear receptor that defines a distinct retinoid response pathway](https://genesdev.cshlp.org/content/9/9/1033); [Polarity-specific activities of retinoic acid receptors determined by a co-repressor](https://www.nature.com/articles/377451a0); [Transcriptional Factors Mediating Retinoic Acid Signals in the Control of Energy Metabolism](https://www.mdpi.com/1422-0067/16/6/14210).
 
* Moreover, RARE has been identified in the promoter regions of RARα and β, indicating that the expression of these RARs is also under the control of atRA. *This is consistent with our gene signatures*
 
* All-trans retinoic acid is degraded by CYP26 enzymes, which belong to cytochrome P450 hydroxylase family. A number of CYP26 family including CYP26A1, B1, C1, and D1 have been characterized and all of them possess the ability to degrade atRA into less bio-active retinoid. Since RALDH2 and CYP26A1 are both regulated by atRA itself, the metabolism of atRA therefore forms an auto-regulatory loop that regulates and balances atRA levels in embryos. Such regulation not only maintains the endogenous atRA level within a normal range, but also allows the organisms to respond to exogenous atRA fluctuation.

-------

## Physiological effect of the hormone

* The RA signaling pathway has been implicated in various developmental processes. During early embryonic development, retinoids act as an important morphogen across different species from invertebrate to metazoan including human.

* Indeed decades ago studies uncovered the roles of RA during embryogenesis by subjecting mammalian or avian embryos to vitamin A deficiency, revealing that retinol is essential for the development of many organs including the hindbrain, spinal cord, forelimb buds, skeleton, heart, eye, pancreas, lung and genitourinary tract. Interestingly, each of the RAR subtypes functions redundantly and most of the RXR subtypes are not critical for the embryonic development. :books: [Retinoic acid signaling pathways](https://journals.biologists.com/dev/article/146/13/dev167502/19797/Retinoic-acid-signaling-pathways)

* The RA metabolic enzymes show distinct differential expression pattern during early embryonic development, and interestingly their expression is regulated by the RA signaling.

* Among the most important and better-established effects of RA is anterior-posterior patterning of the non-axial mesoderm across the dorsal-ventral axis. The retinoic acid-metabolizing enzyme, CYP26A1, is essential for normal hindbrain patterning, vertebral identity, and development of posterior structures. Retinoic acid signaling modulation is ultimately able to profoundly alter patterning of the primary body axis in embryos of the frog Xenopus laevis. 

* In a study, researchers from the RIKEN Brain Science Institute in Japan report a new technique that allows them to visualize the distribution of retinoic acid in a live zebrafish embryo, in real-time. This technique enabled them to observe two concentration gradients going in opposing directions along the head-to-tail axis of the embryo, thus providing evidence that retinoic acid is a morphogen. :books: [Visualization of an endogenous retinoic acid gradient across embryonic development](https://www.nature.com/articles/nature12037)

### Interplay with other developmental pathways

* ATRA can regulate __HOX transcription factors__, and HOX transcription factors play a role in the specification of the embryonic anterior/posterior axis through Hox gene regulation, functioning as downstream effectors of retinoids. 
* ATRA antagonizes __BMP signaling__, which can also be inhibited by HOX transcription factors as part of a negative feedback loop.
* ATRA also downregulates Fgf8 expression during embryonic development and overall shows antagonizing effect with __FGF pathway__. :books: [Retinoic acid controls body axis extension by directly repressing Fgf8 transcription](https://pubmed.ncbi.nlm.nih.gov/25053430/) 
* RA signaling also plays a pivotal role in expressing __sonic hedgehog__ genes during early development, especially of endoderm.
* Also __TGFB__ Pathway is profoundly linked with retinoid signalling.

-----------

## Impact on neurodevelopment

* Retinoid signalling is intrinsically linked to neurogenesis in utero along with neuronal homeostasis in the mature brain. In the developing brain, Retinoic Acid is involved in neural tube patterning, neurogenesis, cell differentiation and synaptic function. :books:[Regulation of prefrontal patterning and connectivity by retinoic acid](https://www.nature.com/articles/s41586-021-03953-x)  

* Throughout neurogenesis, RA is readily available due to the presence of retinaldehyde dehydrogenase 2 (RALDH2), which synthesizes RA in the paraxial mesoderm of the developing embryo. RA then diffuses to the neural plate and spinal cord to promote the differentiation of neural progenitors. :books: [Opposing FGF and Retinoid Pathways Control Ventral Neural Pattern, Neuronal Differentiation, and Segmentation during Body Axis Extension](https://www.cell.com/neuron/fulltext/S0896-6273(03)00565-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0896627303005658%3Fshowall%3Dtrue). 

* The nested gene expression and combinatorial action of RDH10 and RALDH2 causes a posteriorward flow of retinal that generates an initial RA gradient in the trunk mesoderm with a peak at the level of the hindbrain-spinal cord boundary. Subsequent diffusion of RA creates two gradients across the hindbrain and spinal cord, which acquire their final shape through CYP26A1-mediated RA degradation at the anterior (ant.) and posterior (post.) ends of the neural plate. FB, forebrain; MB, midbrain; HB, hindbrain; SC, spinal cord. 🖼️ [Model for the establishment of retinoic acid morphogen gradients in the early embryo](https://www.researchgate.net/publication/23785181_Retinol_dehydrogenase_10_is_a_feedback_regulator_of_retinoic_acid_signalling_during_axis_formation_and_patterning_of_the_central_nervous_system/figures?lo=1)

* Previous evidence shows that RA signaling during brain development cooperating with FGF signaling to repress BMP and ZIC proteins.
In particular, FGF is a general repressor of differentiation, including ventral neural patterning, while RA attenuates Fgf in the neuroepithelium and paraxial mesoderm. RA is further required for neuronal differentiation and expression of key ventral neural patterning genes. FGF and RA pathways are mutually inhibitory and suggest that their opposing actions provide a global mechanism that controls differentiation during body axis extension and neurodevelopment. :books: [Retinoic acid signaling and neuronal differentiation](https://link.springer.com/article/10.1007/s00018-014-1815-9);  
[Opposing FGF and Retinoid Pathways Control Ventral Neural Pattern, Neuronal Differentiation, and Segmentation during Body Axis Extension](https://www.cell.com/neuron/fulltext/S0896-6273(03)00565-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0896627303005658%3Fshowall%3Dtrue)

#### Impact on cortex development

* Several __studies in mice__ highlight the importance of ATRA in regulating cortical neuron generation, primarily when produced from the meninges. RA also induces the differentiation of various types of neurons and glia, by activating the transcription of genes that encode various transcription factors, cell signaling molecules, structural proteins, enzymes and cell-surface receptors. :books: [Retinoic acid from the meninges regulates cortical neuron generation](https://pubmed.ncbi.nlm.nih.gov/19879845/)

* Moreover, __RA regulates the balance between direct and indirect neurogenesis__ in the developing mouse cerebral cortex. :books: [Retinoic acid controls early neurogenesis in the developing mouse cerebral cortex](https://pubmed.ncbi.nlm.nih.gov/28790015/).

* Particularly, __RA Regulates prefrontal cortex development, patterning and connectivity__, as ATRA froms an anterior-posterior gradient in the developing cortex. Since this brain region and its connections with the thalamus are crucial for cognitive flexibility and working memory, an alteration of RA signalling can be causative of NDDs such as autism spectrum disorder and schizophrenia. :books: [Regulation of prefrontal patterning and connectivity by retinoic acid](https://www.nature.com/articles/s41586-021-03953-x)

* In light of these evidences, alterations in RA signaling have been implicated in the __pathophysiology of ASD  and SCZ__. :books: [Impaired neurodevelopmental pathways in autism spectrum disorder: a review of signaling mechanisms and crosstalk](https://jneurodevdisorders.biomedcentral.com/articles/10.1186/s11689-019-9268-y) [Polygenic disruption of retinoid signalling in schizophrenia and a severe cognitive deficit subtype](https://www.nature.com/articles/s41380-018-0305-0). 

* RA pathway has also been implicated in __synaptic plasticity__, suggesting the impact on brain function and behavior persists in adults. :books: [Synaptic Signaling by All-Trans Retinoic Acid in Homeostatic Synaptic Plasticity](https://www.cell.com/neuron/fulltext/S0896-6273(08)00707-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0896627308007071%3Fshowall%3Dtrue)

* Overall the effects of ATRA are multiple in both embryonic and brain development, and peculiarly, in cortical development it __facilitates differentiation of neurons at the expense of proliferating progenitors__. This ability can be harnessed to __induce the differentiation of stem cells into neural cell types__. Indeed, due to its well-established role in promoting neuronal differentiation, RA is commonly added to ESC cultures to induce neural differentiation. :books: [Retinoic acid signaling and neuronal differentiation](https://link.springer.com/article/10.1007/s00018-014-1815-9); [Retinoic acid in the development, regeneration and maintenance of the nervous system](https://www.nature.com/articles/nrn2212)

* Interestingly, it is widely accepted that increased RA levels adversely affect the formation of the head and in particular, inhibit forebrain development by promoting hindbrain expansion resulting in a microcephalic phenotype. On the other hand, reduced RA signaling levels have been linked to an increasing number of developmental syndromes that exhibit microcephaly including FAS, DiGeorge, Smith-Magenis, Matthew-Wood, and Vitamin A Deficiency (VAD) Syndromes, and others. :books: [Reduced Retinoic Acid Signaling During Gastrulation Induces Developmental Microcephaly](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8967241/) 

#### Sex differences 

* Few data are available about the differences between RAR distribution in males and females, as well as about the possible changes in a vitamin A deficient state (VAD) response in both sexes.
* Different isoforms of both RXR and RAR are differentially expressed across male and female developing rat brains. Particularly, RARα and β isotypes were detected in practically all the male brain areas whereas immunostaining was weak or absent in the female brain except RARα. RXRγ was absent in the female brain, while it was observed in some regions in the male. This underlies that __the brain management of retinoic acid differs between males and females__, also leading to differences in their response to VAD diet in terms of receptor expression. :books: [Expression of retinoic acid receptors and retinoid X receptors in normal and vitamin A deficient adult rat brain](https://www.sciencedirect.com/science/article/abs/pii/S094096021200129X)

#### ATRA target genes

* The known __direct RAR target, HoxA1, is required for the differentiation of embryonic stem (ES) cells into neurons__. :books: [Hoxa1 is required for the retinoic acid-induced differentiation of embryonic stem cells into neurons](https://pubmed.ncbi.nlm.nih.gov/18512762/).
* Evidence from a variety of cell culture systems shows that RA directly and indirectly regulates the expression of many other genes such as Btg2, in addition to HoxA1, that facilitate cell cycle exit and differentiation. Importantly, Erf and Etv3 are upregulated by RAR agonists, and down-regulated by RAR antagonists, and knockdown of ERF or ETV3 results in paralysis, loss of primary neurons and increased proliferation of undifferentiated neural progenitors. Thus, these Ets-repressors are key effectors that inhibit neural progenitor identity and promote differentiation.   

#### Review on RA effect during neurogenesis
:books: [__Retinoic acid signaling and neuronal differentiation__](https://link.springer.com/article/10.1007/s00018-014-1815-9)
 
-----------

## In vitro experiment with iPSC and neural models

* Long-term all‑trans RA exposure in high-density cultures induces upregulation of SOX2, PAX6, NeuroD1, and results in neuronal differentiation of hESCs/hiPSCs; Wnt/GSK3β inhibition disrupts this effect. :books:  [Retinoids: Mechanisms of Action in Neuronal Cell Fate Acquisition](https://pubmed.ncbi.nlm.nih.gov/38137880/)

* Hua et al. (2021) used RA alongside Wnt and SHH activators to pattern hiPSC spheroids into cerebellar-layer-like structures. (I) generated Purkinje-like and granule-like cells with electrophysiological properties; (II) RA was critical for achieving correct layer identity. :books: [Cerebellar Differentiation from Human Stem Cells Through Retinoid, Wnt, and Sonic Hedgehog Pathways](https://pubmed.ncbi.nlm.nih.gov/32873223/) :

* In ocular organoid cultures, RA (500 nM–10 µM) modulated neuroretinal maturation, pigmentation, and corneal transparency from days 30–90. :books: [All-trans retinoic acid modulates pigmentation, neuroretinal maturation, and corneal transparency in human multiocular organoids](https://stemcellres.biomedcentral.com/articles/10.1186/s13287-022-03053-1)
  
* Human gastruloids exposed to RA (0–24 h) show segmentation and recapitulate neural tube development – dorsal-ventral axis, neural crest, and anterior-posterior HOX gene patterns :books: [Retinoic acid induces human gastruloids with posterior embryo-like structures](https://pubmed.ncbi.nlm.nih.gov/39164488/)

------------

## Epidemiological data 

#### :books: [Polygenic disruption of retinoid signalling in schizophrenia and a severe cognitive deficit subtype](https://www.nature.com/articles/s41380-018-0305-0)
Progress in identifying variants that may disrupt at-RA functionality in schizophrenia has recently been facilitated by collaborative genome-wide association studies (GWAS) assembled by the Psychiatric Genomics Consortium (PGC). In the largest PGC mega analysis of schizophrenia, over 100 loci were uncovered at rigorous genome-wide significance levels. Interestingly, five genes involved in retinoid biology are located in genome-wide significant loci. Common polygenic risk in retinoid genes implicated by GWAS is associated with schizophrenia, with the retinoic acid receptor beta gene RARB significantly associated with CD subtype. Increasing burden of rare RARB variation was also associated with decreased grey matter volume in a left posterior cerebellar region in the CD subgroup after voxel-wise correction, along with covariation amongst grey matter concentration (GMC) in several brain regions. Analysis of non-coding sequence in DR5-RARE, proximal to genes, revealed an excess of rare variation that could disrupt the expression of retinoid target transcripts. This, and other forms of retinoid dysregulation, was supported by differential expression of DR5-RARE proximal genes in two independent schizophrenia cohorts, which were enriched in functionally significant pathways. 

#### :books: [Epidemiological aspects of prenatal exposure to high doses of vitamin a in Spain](https://link.springer.com/article/10.1007/BF00145783)
Epidemiological study of prenatal exposure to high doses of vitamin A in Spain, using data from the Spanish hospital-based, case-control registry. Although it is difficult to reach conclusions with such a very low exposure level (1.3 per 1,000 livebirths), their results suggest that a teratogenic effect might exist for exposures to high doses of vitamin A. This is concordant with other sources but kind of counterintuitive, since RA is a pro-differentiative agents and RA-based treatment have been proposed in the past few decades to halt cancer cellular division, malignancy and progression. 

#### :books: [Retinoic acid-related orphan receptor alpha (RORA) variants are associated with autism spectrum disorder](https://pubmed.ncbi.nlm.nih.gov/28608249/)
The observed sex bias in ASD towards male has prompted investigators to propose sex-dependent mechanisms for ASD. Retinoic acid-related orphan receptor-alpha (RORA) is a new ASD candidate gene that has been shown to be differentially regulated by male and female hormones. Previous studies have shown deregulation of its expression in the prefrontal cortex and the cerebellum of ASD patients. In this study they looked for possible associations between two functional polymorphisms in the RORA gene (rs11639084 and rs4774388) and the risk of ASD in 518 Iranian ASD patients and 472 age, gender, and ethnic-matched healthy controls, and found that RORA participation is involved in the pathogenesis of ASD, further highlighting RA role in ASD development. 

#### :books: [Vitamin A deficiency increases the risk of gastrointestinal comorbidity and exacerbates core symptoms in children with autism spectrum disorder](https://pubmed.ncbi.nlm.nih.gov/32225174/)

---------------

## Evidence of EDC disruption of the pathway

* :books: [__Highlighting the gaps in hazard and risk assessment of unregulated Endocrine Active Substances in surface waters: retinoids as a European case study__](https://enveurope.springeropen.com/articles/10.1186/s12302-020-00428-0#:~:text=Besides%20retinoids%20themselves%2C%20other%20environmental,biomarkers%20and%20induce%20retinoid%2Dlike)

* :books: [__A screen for disruptors of the retinol (vitamin A) signaling pathway__](https://pubmed.ncbi.nlm.nih.gov/23696197/). Identified BPA as a disruptor of retinol signaling along with other chemicals. 

* :books: [__Bisphenol A modulates germ cell differentiation and retinoic acid signaling in mouse ES cells__](https://www.sciencedirect.com/science/article/abs/pii/S0890623812002274)

* :books: [__Disruption of Retinol (Vitamin A) Signaling by Phthalate Esters: SAR and Mechanism Studies__](https://pubmed.ncbi.nlm.nih.gov/27532513/) (among them there is DphP).

* :books: [__Nuclear receptors are the major targets of endocrine
disrupting chemicals__](https://hal.science/hal-03488727/document): Multiple studies have identified several chemicals among different classes of EDCs, chemicals that bind to and activate RARs in the micromolar range

* :books: [__Retinoid-xenobiotic interactions: the Ying and the Yang__](https://hbsn.amegroups.org/article/view/6877/html): a review on __retinoid signalling involved in xenobiotic metabolism__, finds BPA as a direct activator of RAR-RXR heterodimer, and highlights functional retinoic acid responsive elements (RAREs) have been identified in several genes encoding phase I and phase II xenobiotic-metabolizing enzymes and phase III transporters.

* :books: [__Activation of retinoic acid receptor-dependent transcription by organochlorine pesticides__](https://pubmed.ncbi.nlm.nih.gov/15589975/): organochlorine pesticides, including chlordane, dieldrin, aldrin, endrin, and endosulfan, can activate RARβ and RARγ, but not RXR isoforms, and are able to induce __RAR-mediated transcription of the gene encoding CYP26A1__.

----------

## Crosstalk with other pathways

LXRs form obligate heterodimers with retinoid X receptors (RXRs) to regulate target genes. In differentiating stem cells, RXR partners orchestrate neuronal gene programs: for example, chromatin profiling shows LXR:RXR (and RAR:RXR) complexes drive distinct gene networks during ESC→neuron differentiation. Thus LXR activation often overlaps with retinoic acid (RAR/RXR) signaling in neurogenesis.

:books: [RXR heterodimers orchestrate transcriptional control of neurogenesis and cell fate specification](https://doi.org/10.1016/j.mce.2017.07.033)
