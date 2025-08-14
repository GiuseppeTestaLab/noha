## Table of Content

* [Compounds tested](#Compounds-tested)
* [Results from reference compounds](Results-from-reference-compounds)
* [General pathway information](#General-pathway-information)
* [Physiological effect of the hormone](#Physiological-effects-of-estrogens)
* [Impact on neurodevelopment](#Impact-on-neurodevelopment)
* [In vitro experiment with iPSC and neural models](#In-vitro-experiment-with-iPSC-and-neural-models)
* [Epidemiology and clinical data](#Epidemiology-and-clinical-data)
* [Evidence of EDCs disruption of the pathway](#Evidence-of-EDCs-disruption-of-the-pathway)
* [Crosstalk with other pathways](#Crosstalk-with-other-pathways)

------------------
## Compounds tested: 

* Agonist: [E2 Estradiol](https://www.ncbi.nlm.nih.gov/books/NBK549797/)
* Antgonist: [ICI 182780 Fulvestrant](https://www.ncbi.nlm.nih.gov/books/NBK560854/) --> __complete agonist__: it acts by binding, blocking, and degradation of estrogen receptors, which in turn leads to the complete cessation of estrogen signaling through the estrogen receptors.

-----

## Results from reference compounds
 
### Summary scheme

[Hormonal Card - Estrogen](Schemes/Cards_Hormones_Estrogen.jpg)

### Result Description

* [Gene Signatures](1.GeneSignatures.md#Estrogen)
* [Differential Expression Analysis - Top Genes](2.DEATopGenes.md#Estrogen)
* [Differential Expression Analysis - Functional Enrichment](3.DEAFunctional.md#Estrogen)
* [WGCNA](4.WGCNA.md#Estrogen)

------------------

## General pathway information

#### :books: [Estrogen receptor signaling mechanisms](https://www.sciencedirect.com/science/article/abs/pii/S1876162319300112?via%3Dihub)
* __Estrogens: definition and hystory__
  * The term **“estrogens”** [derives from the Greek words oistros (frenzy, in heat) and gennan (to produce)] refers to a group of female hormones, including **estrone, estradiol, estriol, and estretrol**. However, the word estrogen is commonly used to refer to estradiol (or 17β-estradiol), due to its physiological relevance and predominance during reproductive years.
  * Chemically, estrogens belong to the family of organic compounds known as **steroids**. As such, their core structure is composed of *17 carboncarbon bonds* arranged as *four fused rings* (three cyclohexane rings and a cyclopentane ring). All four estrogens contain 18 carbons (C18H24O2) and are collectively known as *C18 steroids*. They consist of one benzene ring, a phenolic hydroxyl group, and a ketone group (estrone), or one (17β-estradiol), two (estriol), or three (estretrol) hydroxyl groups.
  * Estrogens are primarily synthesized in the *ovaries*, but also in the *adrenal glands* and *adipose tissue*. While females produce all estrogens throughout life, the hormones 16-hydroxyestradiol (estriol) and 15α-hydroxyestriol (estretrol) are predominantly found during pregnancy, and estrone is usually found at higher levels during menopause (Samavat & Kurzer, 2015).
    * Estradiol, the predominant circulating estrogen in humans, it is mainly secreted by the granulosa cells of the ovarian follicles, and the corpora lutea.
    * Estretrol is synthesized exclusively by the fetal liver and reaches maternal circulation through the placenta (Coelingh Bennink, Holinka, Visser, & Coelingh Bennink, 2008; Holinka, Diczfalusy, & Coelingh Bennink, 2008).
    * Estrone, which is produced by aromatization of androstenedione in extraglandular tissues, can be reversibly transformed to estradiol by the enzyme 17β-hydroxysteroid dehydrogenase in peripheral tissues (Bulun, Zeitoun, Sasano, & Simpson, 1999; RYAN, 1959).   
  *  All four estrogens are able to bind to both *nuclear and membrane estrogen receptors*, with different affinity and strength of the response (Watson, Jeng, & Kochukov, 2008).
  *  In **males**, estrogens exert pleiotropic effects by acting on several tissue and organs, including the male reproductive system. The action of estrogens is manifest from prenatal life during which the exposure to estrogen excess might influence the development of some structures of the male reproductive tract. The role of estrogen in regulating human male reproduction is less clear compared to mice. ([Estrogens, Male Reproduction and Beyond](https://www.ncbi.nlm.nih.gov/books/NBK278933/)) (Fig 1 of [Milestones in the advancement of research in the area of estrogens in men](https://www.ncbi.nlm.nih.gov/books/NBK278933/bin/estrogens-male-repro-Image001.jpg))

* __Estrogens biosynthesis__
  * The main substrate for steroid hormone biosynthesis is dietary cholesterol, specifically lowdensity lipoprotein (LDL)-cholesterol (Carr, MacDonald, & Simpson, 1982). Through a process called steroidogenesis, cholesterol is converted to the 21-carbon (pregnanes, progestogens), 19-carbon (androstanes), and 18-carbon (estranes) steroid hormones in gonads, adrenal cortex, and adipose tissue (Miller, 2017). The main site of estrogen synthesis is the ovaries, and specifically the granulosa cells  


#### :books: [Steroid Transport, Local Synthesis, and Signaling within the Brain: Roles in Neurogenesis, Neuroprotection, and Sexual Behaviors](https://www.frontiersin.org/articles/10.3389/fnins.2018.00084/full)

-----

## Physiological effects of estrogens

Estrogens are sex steroid hormones, and as such display a broad spectrum of physiological functions. These include regulation of the menstrual cycle and reproduction, bone density, brain function, cholesterol mobilization, development of breast tissue and sexual organs, and control of inflammation (Liang & Shang, 2013). While estrogens play diverse roles in normal male and female physiology, in certain physiological situations, they can play similar roles in both sexes (Rotstein). In **females**, estrogens are responsible for primary and secondary sexual characteristics. Estradiol promotes epithelial cell proliferation in the uterine endometrium and mammary glands starting in puberty (Gruber, Tschugguel, Schneeberger, & Huber, 2002; Koos, 2011; Simpson et al., 2005). During pregnancy, estrogens produced by the placenta help prepare the mammary gland for milk production (Voogt, 1978). On the other hand, lower levels of estrogens produced **in men** are essential for functions including sperm maturation, erectile function and maintenance of a healthy libido (Schulster, Bernie, & Ramasamy, 2016). All these roles are mediated by estrogen receptors.

#### The estrogen receptors
* **nuclear estrogen receptors** (intracellular and associated with HSP90 when inactive): ER-alpha, ER-beta
  * 6 different structural and functional domains: N-terminal (NTD, A/B domains, AF-1), DNA binding domain (DBD, C domain), the hinge (D domain), the C-terminal region containing the ligand binding domain (LBD, E/F domain, AF-2).
  * when estrogen binds a receptor, --> conformational changes and no more bound to HSP90 --> omo or heterodimerization --> DBD binds ERE (estrogen responsive elements, on DNA) --> expression of estrogen responsive genes
*  **membrane estrogen receptors** : GPER1, ERX, Gq-mER

#### Mechanisms of estrogen receptor signaling
* *via nuclear receptors*: as a steroid hormone, estrogen can enter the plasma membrane and interact with intracellular ERα and ERβ to exert direct effects by binding to DNA sequences.
*  *via membrane receptors*: estrogen can activate intracellular signaling cascades via interaction with the GPER1 and/or ERα and ERβ.
*  *estrogen-mediated signaling events* can be divided into **genomic** and **non-genomic**, due to differences in the cellular and molecular events leading to gene expression regulation in which estrogen-receptor complexes can either bind directly or indirectly to DNA.
*   <ins>NUCLEAR ESTROGEN RECEPTORS (DIRECT GENOMIC SIGNALING)</ins>:  Upon binding of estradiol to ERα or ERβ in the cytoplasm, a conformational change occurs, inducing receptor dimerization. This complex is then translocated to the nucleus, where it binds to the chromatin at ERE sequences, enhancer regions within or close to promoters, and/or 3’-untranslated regions of target genes. While EREs have been
identified in several gene promoters and regulatory regions, it has been reported that more than one third of human genes regulated by estrogen receptors do not contain ERE sequence elements (O’Lone, Frith, Karlsson, & Hansen, 2004)--> *see indirect genomic siganlling*. Over 70,000 EREs in the human and mouse genomes (Bourdeau et al., 2004). Interestingly, 17,000 of these EREs were located near mRNA transcriptional start sites, and only 660 were conserved sites. While these elements share a high degree of sequence similarity, it is important to recognize that the intrinsic sequence composition of the EREs can alter the affinity of the receptor to bind DNA.
* <ins>NUCLEAR ESTROGEN RECEPTORS (INDIRECT GENOMIC SIGNALING)</ins>: an estimated 35% of genes targeted by estrogen lack ERE-like sequences --> “indirect genomic signaling” or “transcriptional cross-talk”, and are based on activation of gene expression by estrogen receptors not binding DNA directly --> estrogen receptor complexes act through protein-protein interactions with other transcription factors and response elements. An important mediator of indirect genomic signaling is the stimulating protein-1 (**Sp-1**). *Examples of genes induced by estrogen via the Sp-1* mechanism are: low-density lipoprotein (LDL) receptor (C. Li, Briggs, Ahlborn, Kraemer, & Liu, 2001), progesterone receptor B (O’Lone et al., 2004), endothelial nitric oxide synthase (eNOS) (Chambliss & Shaul, 2002), GATA binding protein 1 (GATA1), signal transducer and activator of transcription 5 (STAT5) (Björnström & Sjöberg, 2005), and the retinoic acid receptor-1α genes (Sun, Porter, & Safe, 1998). The nuclear estrogen receptors also induce the expression of genes containing the activator protein-1 (**AP-1**) sites through protein-protein interactions (Gaub, Bellard, Scheuer, Chambon, & Sassone-Corsi, 1990). AP-1 is a transcription factor that *regulates key cellular processes such as cell differentiation, proliferation, and apoptosis* (ex. IGF1).
* <ins>MEMBRANE RECEPTOR (INDIRECT NON-GENOMIC SIGNALING)</ins>: GPER1: activation of signal-transduction mechanisms with the subsequent production of intracellular second messengers, cAMP regulation, and protein-kinase activation of signaling cascades that result in indirect changes in gene expression (Lösel & Wehling, 2003).  The protein-kinase cascades can be classified into *four major ones*:
    1. **phospholipase C (PLC)/protein kinase C (PKCs) pathway** (Marino, Pallottini, & Trentalance, 1998) 
    2. **Ras/Raf/MAPK cascade** (Dos Santos et al., 2002; Watters, Campbell, Cunningham, Krebs,& Dorsa, 1997)
    3. **phosphatidyl inositol 3 kinase (PI3K)/Akt kinase cascade** (Marino, Acconcia, & Trentalance, 2003)
    4. **cAMP/protein kinase A (PKA)** signaling pathway (Q. Gu & Moss, 1996; Picotto, Massheimer, & Boland, 1996).   
  *  Both ERα and ERβ are also targets for phosphorylation by protein kinases, including MAPKs, indicating that non-genomic actions of estrogens may also involve *self-regulation of receptor expression*.
  *  These mechanisms are cell-type specific and activated under certain physiological events, and by specific receptor variants

#### Genomic and non-genomic signaling crosstalk 
Two primary mechanisms of signaling "cross-talk" have been identified, both of which involve interactions between proteins from distinct estrogen signaling pathways. The first mechanism entails the binding of estrogen to nuclear estrogen receptors, leading to their dimerization and translocation into the nucleus. Once inside, these receptor complexes associate with transcription factors that have been phosphorylated through GPER1-dependent signaling pathways. The resulting complexes may then interact either with estrogen response elements (EREs) via the nuclear receptors or with specific DNA-binding sites recognized by transcription factors such as AP-1, STATs, ATF-2/c-Jun, Sp1, and NF-κB (Björnström & Sjöberg, 2005). In the second mechanism, GPER1 interacts with membrane-associated ERα and ERβ, initiating kinase cascades that phosphorylate a range of transcription factors—including AP-1, STATs, Elk-1, CREB, and NF-κB—as well as the estrogen receptors themselves. These phosphorylated molecules subsequently bind to DNA to modulate gene expression (Björnström & Sjöberg, 2005). Collectively, these integrated signaling pathways can potentiate transcriptional responses in a tissue-specific and functionally significant manner.
    
#### Estrogen receptor ligand independent signaling 
Mainly triggered by phosphorylation on specific residues (e.g., serine and tyrosine) in the receptors themselves, or their association with coregulators.

#### Estrogen receptor co-regulators and transcriptional control 
The cell also expresses a battery of coregulators that can either enhance or decrease transcriptional activity of steroid hormone receptors. Ex.   steroid receptor coactivator (SRC)/p160 group, the histone acetyltransferase cAMP responsive element binding protein (CREB)-binding protein (CBP)/p300, ATP-dependent chromatin remodeling complexes like SWI/SNF, E3 ubiquitin-protein ligases, and steroid RNA activator (SRA). *still topic of ongoing research*

#### :books: [The Impact of Estrogens and Their Receptors on Immunity and Inflammation during Infection](https://www.mdpi.com/2072-6694/14/4/909)

#### ERE
   *  :books: [Estrogen Receptor Interacting with DNA](https://www.ks.uiuc.edu/Research/pro_DNA/ster_horm_rec/dbd/) 
   *  :books: [Anatomy of the estrogen response element](https://www.cell.com/trends/endocrinology-metabolism/fulltext/S1043-2760%2804%2900009-8)


#### Endogenous and exogenous estrogen receptor ligands
   * There are five main classes of ER ligands: endoestrogens, phytoestrogens, xenoestrogens, selective estrogen receptor modulators (SERMs) and metalloestrogens.

📑 [Table: types of Estrogen Receptor Ligands](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6533072/table/T1/?report=objectonly)

---------

## Impact on neurodevelopment
* :books: [The Impact of Estrogen and Estrogen-Like Molecules in Neurogenesis and Neurodegeneration: Beneficial or Harmful?](https://www.frontiersin.org/articles/10.3389/fncel.2021.636176/full)
  * As estrogen receptors are ubiquitously distributed, they can modulate cell proliferation, differentiation, and survival in several tissues and organs, including the central nervous system (CNS). Estrogens can exert neuroprotective roles by acting as antioxidants, promoting DNA repair, inducing the expression of growth factors, and modulating cerebral blood flow. Additionally, estrogen-dependent signaling pathways are involved in regulating the balance between proliferation and differentiation of neural stem/progenitor cells (NSPCs), thus influencing neurogenic processes.
  * Estrogens are recognized as a significant *steroidal mitogen for epithelial cells*, commonly related to oncogenesis (Anderson et al., 1998). On the other hand, E2 is also involved in the *differentiation of diverse types of normal tissues*, (ex. lobular-alveolar cells and galactophorous ducts). Finally, the *pro-survival effects* of estrogens have been demonstrated in diverse cell types, including germ cells from human testes (Pentikäinen et al., 2000), breast cancer cells (Weldon et al., 2004; Yu et al., 2012), and neurons (Kishi et al., 2005)
  * **Estrogens and Their Role in Cell (Neuronal) Survival: Neuroprotection**: estradiol and metabolites are capable of promoting survival in neurons subjected to diverse stress conditions (Behl and Lezoualc’h, 1998; Robb and Stuart, 2010; Choi et al., 2020) by:
      *  **Anti-oxidant Activity (vs ROS)**:  Estrogens can produce reactive oxygen species by increasing mitochondrial activity and redox cycling of estrogen metabolites (Kumar et al., 2010). On the other hand, estrogens’ phenolic hydroxyl group can act as an antioxidant agent, being a protective factor against cardiovascular and neurodegenerative diseases (Kumar et al., 2010). estrogen *positively regulates* the expression of *proteins related to mitochondrial biogenesis*, such as the nuclear respiratory factor-1 (**NRF-1**) and the peroxisome proliferator-activated receptor-gamma coactivator 1 (**PGC-1**) (Kemper et al., 2013; Klinge, 2017).
      *  **DNA Repair**: estrogen supplementation positively regulates the transcription of **APE1 and NTH1** in the dorsal raphe of OVX old macaques (Bethea et al., 2016). Additionally, estrogen appears to enhance the transcription of **DNA repair enzymes** in the cerebral cortex after hypoxia, contributing to reducing oxidative stress (Dietrich et al., 2013)
      *  Estrogens can induce **BDNF** expression through direct binding to an estrogen-sensitive response element (ERE) in the BDNF gene (Scharfman and MacLusky, 2006). consequent upregulation of **NGF** and **NPY**.
      *  **Synaptic Plasticity**: estrogen receptors ERα and ERβ are located in the axons and synaptic terminals of pre-synaptic cells and also in the dendritic spines (structures crucial for synaptic transmission) of post-synaptic cells (McEwen and Milner, 2007), indicating a role in local regulation of synapses instead of exerting an effect at nuclear transcription (Zarate et al., 2017).
      *  **Cerebral Blood Flow**: By binding to ERα in endothelial cells, estrogens induce nitric oxide (**NO**) release, consequently producing vasodilation

  * **Estrogens and estrogen-like molecules in neurogenesis and gliogenesis**:
      * Experiments performed by Wang et al. (2008) using human-derived NSPCs revealed a significant increase in their proliferative activity after estrogen stimulation. In fact, 17β-estradiol (E2) induces time- and dose-dependent effects associated with an increase in DNA replication and up-regulation of cell cycle proteins. Remarkably, ERβ-, but not ERα-specific ligands, lead to increased phosphorylation of ERK 1/2, which initiates cell cycle entry (demonstrated by the up-regulation of mitotic markers such as CDK1/cdc2 and PCNA), DNA replication, and duplication of the centrosome -the major microtubule-organizing center in the cell- in dividing cells (Wang et al., 2008). Thus, in these cells, E2 effects appear to be predominantly mediated by ERβ. Similarly, primary cultures of embryonic rat-derived NSPCs respond to E,2 increasing their proliferation (Brannvall et al., 2002).
      *  A growing body of evidence suggests estrogens are relevant in oligodendrogenesis by promoting oligodendrocyte differentiation in a GPER-1-dependent manner. In this regard, it has been shown in a mouse model that the pharmacological treatment with the classical estrogen receptor antagonist Fulvestrant (ICI-182,780) modifies the proliferative activity of NSPCs without affecting the differentiation toward the oligodendroglial lineage (Okada et al., 2010).
  
  * **Estrogens and Neurodegenerative Diseases**: Experimental and clinical data suggest that estrogens have protective effects against neurodegenerative diseases (Kumar et al., 2010; Vail and Roepke, 2019).
      * **AD** post-menopausal women have a higher prevalence of this pathology than men
      *  **PD** Clinical data show that women with low estrogen exposure can develop PD earlier than women with high estrogen exposure (Ragonese et al., 2004), indicating a role for this sex hormone in PD’s pathogenesis.
      *  **Brain Ischemia and Stroke**  estrogens can act as anti-inflammatory agents in the blood vessel wall, protecting it from cytokines and free radicals. 

* :books: [Precise Therapeutic Targeting of Distinct NRXN1+/− Mutations](https://www.biorxiv.org/content/10.1101/2023.10.28.564543v1)
  * Treatment with β-estradiol increases **NRXN1** expression in **glutamatergic neurons**, while antisense oligonucleotides knockdown mutant isoform expression across both glutamatergic and GABAergic neurons. Direct or indirect manipulation of NRXN1 splicing isoforms provides a promising therapeutic strategy for treating humans with 2p16.3 deletions.

* :books: [Estradiol and the Development of the Cerebral Cortex: An Unexpected Role?](https://www.frontiersin.org/articles/10.3389/fnins.2018.00245/full)
   * The cerebral cortex undergoes rapid folding in an “inside-outside” manner during embryonic development, resulting in the establishment of six discrete cortical layers. This unique cytoarchitecture occurs via the coordinated processes of neurogenesis and cell migration. In addition, these processes are fine-tuned by a number of extracellular cues, which exert their effects by regulating intracellular signaling pathways. Interestingly, multiple brain regions have been shown to develop in a sexually dimorphic manner. In many cases, estrogens have been demonstrated to play an integral role in mediating these sexual dimorphisms in both males and females.  
      
* :books: [Steroid Transport, Local Synthesis, and Signaling within the Brain: Roles in Neurogenesis, Neuroprotection, and Sexual Behaviors](https://www.frontiersin.org/articles/10.3389/fnins.2018.00084/full)
   * 17β-E2 treatment of rat embryonic NSCs increases cell proliferation and neuronal differentiation (Brännvall et al., 2002), and promotes human embryonic NSC differentiation into dopaminergic neurons (Kishi et al., 2005).
   * The **brain is a cholesterol-rich organ**, accounting for about 25% of the total amount present in humans. In the CNS, *cholesterol* is mainly synthesized by astrocytes, oligodendrocytes, microglia and to a lesser extent by neurons, where it is essentially present in its unesterified form (Björkhem and Meaney, 2004; Dietschy and Turley, 2004; Zhang and Liu, 2015). Brain cholesterol is *involved in myelin sheath genesis, in synaptogenesis, and neurotransmission as well as in neurosteroidogenesis* (Mauch et al., 2001; Do Rego et al., 2009; Linetti et al., 2010; Liu et al., 2010; Zhang and Liu, 2015).

* :books: [Genes, Gender, Environment, and Novel Functions of Estrogen Receptor Beta in the Susceptibility to Neurodevelopmental Disorders](https://www.mdpi.com/2076-3425/7/3/24):
   * Estrogen receptor beta (ERβ) signaling during perinatal brain development and put it in the context of sex-specific differences in neurodevelopmental disorders.
   * ERβ directs DNA demethylation to specific sites, of which one such site may bear consequences for the susceptibility to the neurological reading disorder *dyslexia* (candidate genes such as CYP19A1 and DYX1C1 among others).
 

 * **Relevant works from [Tollkunh lab](https://www.tollkuhnlab.org/)**
   * :books: [Estrogen receptor alpha is required in GABAergic, but not glutamatergic, neurons to masculinize behavior](https://www.sciencedirect.com/science/article/abs/pii/S0018506X17300727?via%3Dihub)
   * :books: [Regulation of neural gene expression by estrogen receptor alpha](https://www.biorxiv.org/content/10.1101/2020.10.21.349290v1): We profiled gene expression and chromatin accessibility and show *ERα induces a neurodevelopmental gene program in adulthood*. We further demonstrate that ERα binds with Nuclear factor I X-type (**Nfix**) to regulate a *male-biased gene expression program* that initiates in early life. Our results reveal a neural strategy for ERα-mediated gene regulation and provide molecular targets that underlie estrogen’s effects on brain development, behavior, and disease.

* **Estrogen and sex differences**
  * :books: [Relating sex differences in cortical and hippocampal microstructure to sex hormones](https://www.biorxiv.org/content/10.1101/2023.11.01.565213v1)
    * Albeit correlative, this study underscores the importance of incorporating sex hormone variables into the investigation of brain structure and plasticity.
  * :books: [Genes, Gender, Environment, and Novel Functions of Estrogen Receptor Beta in the Susceptibility to Neurodevelopmental Disorders](https://www.mdpi.com/2076-3425/7/3/24):
    * Dysregulations in sex-hormone signaling, like those evoked by endocrine-disrupting chemicals, may affect this and other neurodevelopmental disorders in a sex-specific manner through ERβ.

* Other references
   * :books: [Estradiol and the Developing Brain](https://journals.physiology.org/doi/full/10.1152/physrev.00010.2007)
   * :books: [Estrogen promotes differentiation and survival of dopaminergic neurons derived from human neural stem cells](https://onlinelibrary.wiley.com/doi/10.1002/jnr.20362)

## In vitro experiment with iPSC and neural models
* :books: [Androgens increase excitatory neurogenic potential in human brain organoids](https://www.nature.com/articles/s41586-021-04330-4):
    * it is oestradiol, and not testosterone, which masculinizes the mouse brain.
    * In both male and female organoids, androgens (DHT or testosterone) led to an increase in intermediate progenitors, which are basally located cells that produce excitatory cortical neurons (Fig. 1b, c, Extended Data Fig. 2a). This was not seen after treatment with oestradiol.
    * Because inhibitory neurons are generated outside the cortex, in the ventral forebrain (Extended Data Fig. 9a), the authors ventralized organoids and treated them with sex steroid hormones; this led to no change or a slight decrease in ventral intermediate progenitors after treatment with DHT (Fig. 3a, b), whereas treatment with oestradiol had no effect.
    * Notably, they observed no expression of aromatase (CYP19A1), an enzyme that converts testosterone to oestrogen (Extended Data Fig. 9c), or oestrogen receptors (ESR1 and ESR2), supporting a nonresponsive state to oestrogen as they observed, and in contrast to the mouse. The differences between these observations in human tissue and previous studies in mice prompted them to explore the potential evolutionary divergence between these two species using organoids. The authors observed that in mouse organoids, it was oestradiol, and not DHT, that caused an increase in intermediate progenitors (Extended Data Fig. 9e–g), in contrast to human.
    * The authors encourage doing further studies of the role of oestrogen itself on other stages or processes during human brain development.

* :books: [Parallel in vivo analysis of large-effect autism genes implicates cortical neurogenesis and estrogen in risk and resilience](https://www.sciencedirect.com/science/article/pii/S0896627321000027)

   * **Estrogen confers resilience to multiple ASD genetic risks in human in vitro models:** The authors tested the interaction of estrogen with the convergent ASD phenotype in human in vitro model systems. To do this, they assayed whether treatment with 17-β-estradiol could rescue the relative increase in proliferative Ki67+ NPCs derived from human CRISPRi cell lines repressing three genes of apparently disparate cellular function (DYRK1A, NRXN1, or ADNP). Indeed, in all three ASD risk gene CRISPRi lines, treatment with **estrogen significantly moved the proportion of Ki67+ cells closer to that of nontargeting control lines** (Figure 8A). They validated this interaction in a 3D model of human brain development. They treated human cortical organoids with DMSO, 17-β-estradiol, harmine, or both and assayed proliferation by Ki67 staining. Consistent with Xenopus and human NPCs, **17-β-estradiol treatment decreased the proportion of proliferating cells**, DYRK1A inhibition by harmine increased it, and treatment together moved it closer to control levels (Figures 8B–8E). Thus, in both Xenopus and human model systems, **activating estrogen signaling reduces the number of proliferative progenitor cells and modulates the convergent phenotype resulting from ASD risk gene disruption**. They next tested the endogenous relevance of this finding by **inhibiting estrogen signaling with a stable CRISPRi** human dorsal forebrain NPC line repressing ERβ/ESR2 and assaying whether these cells showed NPC-relevant defects during neural differentiation. They observed a significant **increase in NPC markers such as EMX2, VIM, and HES1** after exposure to conditions promoting differentiation to excitatory neurons (Figure S5; Table S5), suggesting disruptions in neurogenesis. Overall, these experiments suggest that in human cells, **estrogen affects NPC biology and can suppress the convergent ASD-relevant phenotype**.


* :books: [Mitigating Effect of Estrogen in Alzheimer’s Disease-Mimicking Cerebral Organoid](https://www.frontiersin.org/articles/10.3389/fnins.2022.816174/full):
   * Cerebral organoids protocol adapted from Lancaster and Knoblich (2014)
   * Dose-Dependent Treatment With Estrogen Increased the Expression of Neuronal Markers in Induced Pluripotent Stem Cells-Derived Cerebral Organoids (co-stimulated with amyloid-beta (Aβ) to mimic AD)

-------------

## Epidemiology and clinical data 
#### :books: [Implications of Prenatal Steroid Perturbations for Neurodevelopment, Behavior, and Autism](https://academic.oup.com/edrv/article/35/6/961/2354721)
  
#### :books: [Estrogen Signaling as a Therapeutic Target in Neurodevelopmental Disorders](https://jpet.aspetjournals.org/content/360/1/48)
  * A number of studies from both clinical as well as preclinical research suggest a protective role of estrogen in neurodevelopmental disorders, including autism spectrum disorder (ASD) and schizophrenia. Alterations in the levels of estrogen receptors have been found in subjects with ASD or schizophrenia, and adjunctive estrogen therapy has been shown to be effective in enhancing the treatment of schizophrenia.
  *  The **neuroprotective effects** of estrogen occur through genomic and nongenomic signaling, antioxidant functions, and the maintenance of neuronal ATP through estrogen receptors (Brann et al., 2007). Estrogens can *upregulate antiapoptotic genes and downregulate proapoptotic genes* through diffusion of estradiol through the cell membrane and translocation and binding of estrogen receptors to genes in the nucleus, promoting or inhibiting transcription (Dubal et al., 1999; Choi et al., 2004; Scott et al., 2012). This regulation seems to be largely dependent on ERa and occurs on the scale of hours (Dubal et al., 2001; Scott et al., 2012).
  *  Estrogen is known to promote the synthesis of neurotrophins (ex. **BDNF**) and protects the brain against inflammation and stress.
  *  The fetal brain is exposed to maternal estradiol and that which is locally produced by the fetus (McCarthy, 2008). The developing brain shows high expression levels of estrogen receptors, which regulate gene expression and signal transduction (McCarthy, 2008). This expression is somewhat sex specific and varies over time as brain development progresses (Yokosuka et al., 1997; McCarthy, 2008).
  *  **ASD:** Many studies over the years have shown that increased testosterone exposure during pregnancy (Tordjman et al., 1997; Auyeung et al., 2009, 2012; Whitehouse et al., 2010; Palomba et al., 2012; Xu et al., 2015), decreased aromatase expression (Chakrabarti et al., 2009; Pfaff et al., 2011; Crider et al., 2014), and reduced estrogen or estrogen receptor expression (Chakrabarti et al., 2009; Crider et al., 2014) are significantly correlated with the development of ASD.
  *  **Schizofrenia:**  estrogen levels correlate with the symptoms of schizophrenia. Several clinical studies have reported a correlation between low plasma estrogen levels and an increase in the risk for schizophrenic symptoms in women (Mahé and Dumaine, 2001)
  *  **Bipolar disorder:** Studies have shown that reduced levels of estrogen, for example, postpartum, result in increased psychosis (Arnold, 2003).

#### :books: [Maternal hormonal milieu influence on fetal brain development](https://onlinelibrary.wiley.com/doi/10.1002/brb3.920)
   * An adverse maternal hormonal environment during pregnancy can be associated with abnormal brain growth. Subtle changes in fetal brain development have been observed even for maternal hormone levels within the currently accepted physiologic ranges.
   * Besides sexual behavior, both fetal androgens and estrogens have also been associated with differential growth in other sexually dimorphic brain areas, including prefrontal cortex, cerebellum, amygdala, and hippocampus (Neufang, Specht, & Hausmann, 2009; Peper et al., 2009), explaining sex differences found in cognitive (Berman et al., 1997) and affective skills later in life (Amin, Epperson, Constable, & Canli, 2006; Lombardo et al., 2012). Other factors contribute to sex differences in brain morphology, including enzymes involved in sex steroid biosynthesis and metabolism, such as aromatase (Biegon et al., 2010).
   * Specific structural effects of prenatal androgens, estrogens, and progesterone are neurite outgrowth and synaptogenesis, dendritic branching, and myelination (Garcia‐Segura & Melcangi, 2006; Haraguchi et al., 2012). The precise molecular mechanisms for these hormonal effects are still being elucidated but include activation of programmed cell death, differential expression of transcription factors and microRNAs, DNA methylation, and histone post‐translational modifications of genes, including steroid receptor genes (Gore, Martien, Gagnidze, & Pfaff, 2014).
   * Interactions between sex hormones and neurotransmitters, such as serotonin, dopamine, GABA, and glutamate, have also been described to be of relevance in animals and humans (Barth, Villringer, & Sacher, 2015).
   * In humans, there were observed sexually dimorphic effects of sex steroids on behavior, with increased cord sex hormones (androgens, estrogens, and progesterone) affecting mood and cognition (Jacklin, Wilcox, & Maccoby, 1988; Marcus, Maccoby, Jacklin, & Doering, 1985). Increased levels of sex hormones during fetal development might also predispose male children to autism spectrum disorders (Malkki, 2014).
   * Maternal functional polymorphisms in the sex steroid synthesis and metabolism pathways that are associated with higher estrogen levels were related to attention problems, hyperactivity, and poorer adaptive skills in male offspring (Miodovnik et al., 2012)

-----------
 
## EDC exposures and estrogens

* :books: [Implications of Prenatal Steroid Perturbations for Neurodevelopment, Behavior, and Autism](https://academic.oup.com/edrv/article/35/6/961/2354721) :
   * Hydroxylated Polychlorinated bisphenols (OH-PCBs) are major PCB metabolites and high-affinity inhibitors of human estrogen sulfotransferase (SULT1E1), which sulfonates estrogens and thus prevents them from binding to and activating their receptors [(Wang T. et al, J Biol Chem. 20219](https://pubmed.ncbi.nlm.nih.gov/33524392/)
   * Jacobson and Jacobson (128) studied a population of humans living near Lake Michigan whose mothers, although pregnant, consumed fish contaminated with relatively high levels of PCBs (Polychlorinated biphenyls). Children with the highest PCB concentrations had the lowest IQ scores and deficits in memory and attention.

* :books: [The Impact of Estrogen and Estrogen-Like Molecules in Neurogenesis and Neurodegeneration: Beneficial or Harmful?](https://www.frontiersin.org/articles/10.3389/fncel.2021.636176/full)
   * Bisphenol-A --> **Deregulation of Wnt pathway** --> Reduced proliferative activity and neuro/oligodendrogenesis in rat
   * **BPA** Bisphenol-A and the Proliferation/Differentiation Balance of NSPCs and neuronal maturation --> read there, but mainly in rat

 * :books: [Estrogen Signaling as a Therapeutic Target in Neurodevelopmental Disorders](https://jpet.aspetjournals.org/content/360/1/48) 
   * One study showed a significant association with increased BPA exposure and risk for the development of ADHD-hyperactivity symptoms at 4 years of age, although this effect disappeared by 7 years of age (Casas et al., 2015). Another study showed that gestational exposure to BPA resulted in an increased risk for the development of hyperactivity and aggressive behavior in children (Perera et al., 2012). The effect was stronger in girls than their male counterparts (Perera et al., 2012). Also, an aromatase gene (CYP19) variant is correlated with increased hyperactivity in boys (Miodovnik et al., 2012), suggesting that alterations in the estrogen-testosterone balance have significant behavioral outcomes.

* :books: [Exposure to phthalates impaired neurodevelopment through estrogenic effects and induced DNA damage in neurons](https://www.sciencedirect.com/science/article/pii/S0166445X19310100)
   * Induced pluripotent stem cells (iPSC)-derived human neurons exposed to phthalates triggered double-strand DNA breaks in vitro.
   * Exposure to phthalates (i.e. DBP, DINP, BBP) affects neurodevelopment in zebrafish embryos and induces neurotoxicity in human neurons partly through disrupting the expression of estrogen receptors (esr1, esr2a, esr2b).

 * :books: [Maternal hormonal milieu influence on fetal brain development](https://onlinelibrary.wiley.com/doi/10.1002/brb3.920)
   * Developmental exposure to estrogenic endocrine disruptors in rodents was associated with sexually dimorphic social and anxiety behaviors and learning difficulties in adolescence and adulthood (Carbone et al., 2013; Kundakovic et al., 2013)

 * :books: [Developmental estrogen exposure in mice disrupts uterine epithelial cell differentiation and causes adenocarcinoma via Wnt/β-catenin and PI3K/AKT signaling](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002334)
   *   Mouse model of endometrial adenocarcinoma that results from brief developmental exposure to an estrogenic chemical, **diethylstilbestrol (DES)**, to determine causative factors. We defined a mechanistic pathway from developmental exposure to an endocrine-disrupting chemical to the *development of adult-onset cancer*. These findings provide an explanation for how human cancers, which are often associated with abnormal activation of PI3K/AKT signaling, could result from exposure to environmental insults during development.
  
 * :books: [Cellular and molecular features of EDC exposure: consequences for the GnRH network](https://www.nature.com/articles/s41574-020-00436-3)
    * The epigenetic, molecular and cellular organization of the gonadotropin-releasing hormone network is most vulnerable to EDCs during early development.
    * Effects of EDCs are not limited to classic agonist or antagonist action on sex steroid receptors but also induce long-lasting gene expression and epigenetic changes in the developing brain.
 
 * :books: [Endocrine-Disrupting Chemicals and Hippocampal Development: The Role of Estrogen and Androgen Signaling](https://karger.com/nen/article/doi/10.1159/000531669/853476/Endocrine-Disrupting-Chemicals-and-Hippocampal)

----------

## Crosstalk with other pathways
* :books: [Stress, Sex, and Sugar: Glucocorticoids and Sex-Steroid Crosstalk in the Sex-Specific Misprogramming of Metabolism](https://pmc.ncbi.nlm.nih.gov/articles/PMC7382384/#s5): ER and GR have been shown to affect each other’s action both by altering receptor and ligand availability as well as by modulating each other’s genomic binding and transcriptional end points.

* :books: [Estradiol Rapidly Rescues Synaptic Transmission from Corticosterone-induced Suppression via Synaptic/Extranuclear Steroid Receptors in the Hippocampus](https://academic.oup.com/cercor/article/22/4/926/424336?login=true)

* :books: [Thyroid Hormone Causes Mitogen-Activated Protein Kinase-Dependent Phosphorylation of the Nuclear Estrogen Receptor](https://academic.oup.com/endo/article/145/7/3265/2500000)

* :books: [Thyroid hormone- and estrogen receptor interactions with natural ligands and endocrine disruptors in the cerebellum](https://www.sciencedirect.com/science/article/pii/S0091302217300602)

* :books: [AhR and ARNT modulate ER signaling](https://www.sciencedirect.com/science/article/pii/S0300483X09004697)

* :books: [Chemical activation of estrogen and aryl hydrocarbon receptor signaling pathways and their interaction in toxicology and metabolism](https://www.tandfonline.com/doi/10.1080/17425255.2019.1569627?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)

* :books: [Endocrine disrupting chemicals targeting estrogen receptor signaling: Identification and mechanisms of action](https://pmc.ncbi.nlm.nih.gov/articles/PMC3119362/)

* :books: [Limiting Effects of RIP140 in Estrogen Signaling: POTENTIAL MEDIATION OF ANTI-ESTROGENIC EFFECTS OF RETINOIC ACID](https://www.jbc.org/article/S0021-9258(19)30550-2/fulltext#FN1)

* :books: [An R package for generic modular response analysis and its application to estrogen and retinoic acid receptor crosstalk](https://www.nature.com/articles/s41598-021-86544-0#Sec9)

* :books: [LXRβ/estrogen receptor-α signaling in lipid rafts preserves endothelial integrity](https://pmc.ncbi.nlm.nih.gov/articles/PMC3726156/)
