---
title: Loci multiplicity, mutations and transcriptional profiling
created: 2023-03-31T07:09:14.787Z
modified: 2023-04-11T15:26:57.060Z
---

# Loci multiplicity, mutations and transcriptional profiling


Starting from a general assumption as direct extension of the central dogma of MB which if cells have multiple copies of genes than their expressions are proportional to their respective copy number, it follows that major cell transcriptional differences can be quantitatively explained by presence of cell populations with different karyotypes.

Hence, if a gene $g$ as an expression per unit of copy $\mu_g$ then expression of cell $c$  is: 
$$e_{c,g} = s_c \mu_g C_{c,g}$$ 

for all genes,where $C_{c,g}$ is the copy number of the cell in the chromosome macroregion containing the gene $g$ and $s_c$ is the cell sampling depth.

Clonealign, an R package, makes this assumption and, based on a GLMM, associate  transcript profiles of each cell to one of the defined copy number profiles. 

There are various reasons why this assumption does not hold true. For example, SNPs and indels in a genes may result in nonsense or missense mutations  causing exonic disfunction as loss or gain of stop codons.
Furthermore, mutations may occurs before or after the copy number alterations increasing the variability in gene expressions among samples and single cells derived from the same source.

I propose the following change to the relation between average gene expression e copy number:
$$e_{c,g} = s_c \mu_g \sum_jM^{\epsilon}_{c,g,j}.$$

--------------

# Clonealign model and notation

Gene expression transformation and inverse:

$\xi(x)=\log(1+e^x)$
$\xi^{-1}(y)=\log(e^y-1)$
$\xi_{\mathrm{soft}}^{-1}(y)=\log(1-e^{-|y|})+ \max(0,y)$

Karyotype selection function
$$\kappa_j(\Gamma_{i,j})=\begin{cases}
\argmax_i (\Gamma_{i,j}) & \text{if } \max_i(\Gamma_{i,j})>p_{min}\\
-1
\end{cases}$$


Inputs

$Y \in \mathbb{R}^{N \times G}$
$L \in \mathbb{R}^{G \times C}$

$K \geq 0$



$C_{\mathrm{allele}} \in \mathbb{R}^{V \times C}$
$V = \# \mathrm{variants} \in \mathrm{VCF}$

$M_\mathrm{Cov} \in \mathbb{R}^{V \times N}$
$M_\mathrm{Ref} \in \mathbb{R}^{V \times N}$
$M_\mathrm{Alt} = M_\mathrm{Cov}- M_\mathrm{Ref} \in \mathbb{R}^{V \times N}$

$$
[\mathbb{1}\otimes \log\rho_{\mathrm{var}}]_{c,v,n} =
1_c \times \log \begin{cases} 
\mathrm{B}_{\mathrm{Bin}}(M_\mathrm{Alt}, M_\mathrm{Cov},\alpha_{mid},\beta_{mid})_{v,n}
& [C_{\mathrm{allele}}^T]_{c,v}=2\\
\frac{1}{2} \mathrm{B}_{\mathrm{Bin}}(M_\mathrm{Alt}, M_\mathrm{Cov},\alpha_{low},\beta_{low})_{v,n}+
\frac{1}{2} \mathrm{B}_{\mathrm{Bin}}(M_\mathrm{Alt}, M_\mathrm{Cov},\alpha_{hi},\beta_{hi})_{v,n} & \text{otherwise}
\end{cases}
$$

$$
\mathcal{L}_{\mathrm{allele}}^T=\sum_v [\mathbb{1}\otimes \log\rho_{\mathrm{var}}]_{c,v,n}
$$

$\mathcal{L}_{\mathrm{allele}} \in \mathbb{R}^{N \times C}$



$X \in \mathbb{R}^{N \times P}$

$\mathrm{PCs}= \mathrm{PCA}_K(\xi(Y))+ \mathcal{N}(\mathbb{\mu}_0,\mathbb{\sigma}_{1/2}) \in \mathbb{R}^{N \times K}$

$\beta \in \mathbb{R}^{G \times P}$


$W \in \mathbb{R}^{G \times K}$


$\psi=\mathrm{PCs}$


$f=\psi W^T+ X\beta^t \in \mathbb{R}^{N \times G}$



$\chi \in \mathbb{R}^{K}$

$\chi_e=\exp(\chi) \in \mathbb{R}^{K}$

$\log_\alpha=\vec{\alpha}-\mathbb{1} \otimes log(\mathrm{tr}(e^{\vec{\alpha}})) \in \mathbb{R}^{C}$

-----------

$[\mu]_{g}= E_n(\frac{Y_{n,g}}{E_g(Y_{n,g})})\otimes \mathbb{1} \in \mathbb{R}^{G}$

$\sigma_i=\mathbb{1} \in \mathbb{R}^{G}$

$\mathcal{Q}^\mu=\mathcal{N}(\xi^{-1}(\mu),\sigma_i) \in \mathbb{R}^{G}$

-----------

$S = \text{\# samples}$

$\mu_s(S)=\mathbb{1} \otimes \mathcal{Q}^\mu \in \mathbb{R}^{S \times G}$


$\mu_{s,c,n,g}=\frac{\mu_s(S)L_{g,c} e^{f_{n,g}}} {\sum_g\mu_s(S)L_{g,c} e^{f_{n,g}}}$

-----------

$s_n=\sum_g Y_{n,g}$


${{y_\mathrm{pdf}}_{s,c,n}}^g \sim \mathrm{Multinom(s_n,\mu_{s,c,n,g})}$

$\sum_g{y_{pdf}}^g=1$

$\mathcal{L}_\mathrm{clone}(Y)= \log({y_{\mathrm{pdf}}}^g(Y)) \in ?$

-----------

### (i) $E_q[\log(\rho(y | z, \theta))]$

$\mathcal{L}={\mathcal{L}_\mathrm{clone}}_{s,c,n} + {\mathcal{L}_\mathrm{allele}^T}_{c,n}$

$\Gamma_{\mathrm{logit}} \in \mathrm{R}^{N \times C}$
$\Gamma_{n,c} = \frac{{\Gamma_{\mathrm{logit}}}_{n,c}}{\sum_n{\Gamma_{\mathrm{logit}}}_{n,c}}$


$\mathbb{E}_q[\rho(y | \theta)] = \sum_{n,c} \Gamma_{n,c}\ {\mathbb{E}_s(\mathcal{L})^T}_{n,c})$


$W_{\log\rho} = \sum \log (\mathcal{N}(\mu_0,\frac{1}{\sqrt{\chi_e}})(W)) \in \mathbb{R}^{C \times G \times K}$



$\chi_{\log \rho} = \sum \log(\mathrm{Gamma}(2,1)(\chi)) \in \mathbb{R}^{K}$

$\psi_{\log \rho} = \sum\log(\mathcal{N}(\mu_0,\sigma_1)(\psi)) \in \mathbb{R}^{N \times K}$





  # (ii) E_q[log p(theta)]
  E_log_p_p <- tf$reduce_sum(log_alpha * gamma) +
    tf$reduce_sum(tfd$Normal(loc = tf$zeros(shape(1), dtype = dtype), scale = tf$ones(shape(1), dtype = dtype))$log_prob(tf$log(mu_samples))) / tf$to_float(S) +
    tf$reduce_sum(tfd$Dirichlet(tf$constant(rep(1/C, C), dtype = dtype))$log_prob(tf$exp(log_alpha) + tf$constant(1e-3, dtype = dtype)))

  if(K > 0) {
    E_log_p_p <- E_log_p_p + W_log_prob + chi_log_prob + tf$reduce_sum(p_psi)
  }


  # (iii) E_q[log q]
  E_log_q <- tf$reduce_sum(tf$reduce_mean(qmu$log_prob(mu_samples), 0L)) +
    tf$reduce_sum(tf$where(gamma == 0, tf$zeros(shape = gamma$shape, dtype = dtype), gamma * tf$nn$log_softmax(gamma_logits)))


  elbo <- EE_p_y + E_log_p_p - E_log_q




--------------

--------------

--------------

--------------

--------------

%These aspects may cause heterogeneity in the cells expression compared to few karyotypes. 


# C_{o,c,i}

$C_{o,c,r}$
$C_{o,c,i}$
$M_{o,c,i,m}$
$\epsilon_{i,m}$
$\mu_{g}$
$e_{o,c,g}$
$s_{o,c}$


$$e_{o,c,g} = s_{o,c}\ \mu_{g(i)}\ \frac{\sum_mM_{o,c,i,m}\ \epsilon_{i,m}}{}$$



$\sum_mM_{o,c,i,m} \prod_i\epsilon_{i,m}$

For each cell $c$ there are different mutational profiles $J_c$. For each region $r$ with CN equal to $C_{c,r}$, the number of mutational profiles $J_{c,r}\leq J_c \leq C_{c,r}$.
This relation hold for all cells. Hence, $J_{r}\leq J \leq C_{r}$.

Let the mutational profile $j$ be defined as $\vec{m}_j$ where in each position there is the transition from the original base $A$ to the mutated base $a$, $A\rightarrow a$ (indels included).

A mutational profile projected in a region $r$ is give by:
$$\mathbb{P}_r(\vec{m}_j)=\vec{m}_{j,r}$$

Each mutational profile $\vec{m}_{j,r}$ has a multiplicity $M$ which depends on the region $r$, $M_{c,r,j}$, such that for all $r$ of any size 
$$\sum_{j=1}^{J_r} M_{c,r,j}=C_{c,r} \quad\quad \forall r.$$

The size of regions biologically relevant are those with constant CN or where genes are defined or the single alleles. 

Let a mutation in $i$ of a profile $\vec{m}_{j,r}$ have an effectiveness $\epsilon_{j,r,i} \in \{0,1\}$, then

$$\vec{\epsilon}_{j,r}=\epsilon(\mathbb{P}_r\vec{m}_{j})$$
and a total effectiveness on the whole region $r$
$$E_{j,r}=\mathrm{sp}(\mathrm{diag}(\vec{\epsilon}_{j,r})) = \prod_i\vec{\epsilon}_{j,r}.$$

Multiplicity is introduced not only to address the problem in terms of mutational and CN events arsing in clonal evolution, but also to quantitatively define the quantity of transcript in cells or organoids.  
Therefore, we define an effective multiplicity

$$M^\epsilon_{c,r,j}=\begin{cases}
M_{c,r,j} & \text{if } \prod_i\vec{\epsilon}_{j,r}[i]=1\\
0
\end{cases},$$
which can be rewritten as
$$M^\epsilon_{c,r,j}=M_{c,r,j} \prod_i\vec{\epsilon}_{j,r}[i]$$

Let $r$ correspond to a gene $g$ then:
$$e_{c,g} = s_c \mu_g \sum_jM^{\epsilon}_{c,g,j}.$$



As already said, with in a cell, there are $J$ haplotypes, and in a gene region within a cell there are $J_g$ haplotypes.

Some haplotypes and multiplicities are shared among cells such that if $n_{g,j}$ is the number of cell in the organoid $o$ with the same haplotype and multiplicity, then:
$$e_{o,g} = s_o \mu_g \sum_jn_{g,j}M^{\epsilon}_{g,j}.$$

Let us notice that we dropped the $c$ index in the multiplicity and there is no cell depth $s_c$ because, experimentally, cells are leased and cellular material gets mixed all together, and only after that,  cDNA is fragmented and undergoes PCR cycles. 

One should notice that if and only if the ISA holds then:
$$e_{o,g}=\sum_{c\in o}e_{c,g}.$$

Would be nice to write the rhs in terms of loci $i$ and set of mutations in a gene.
What happens when mutations are clonal? When mutations are not clonal?
Furthermore who does it change for organoids?





---------------

Write

$\vec{p}_{c,r,j}= \vec{p}_{c,r}(\vec{m}_j)$


$?\sum_{j| \vec{m}_{j}(i)=m} 1=M_{c,i,m}$
$\sum_m M_{c,i,m}=C_{c,r} \quad\quad \forall i\in r$


$\epsilon_{i,m}M_{c,i,m}=M^\epsilon_{c,i,m}$

$\sum_mM^\epsilon_{c,i,m}=M^\epsilon_{c,i}$

$$M^\epsilon_{c,g}=\begin{cases}
M^\epsilon_{c,i,m} & \text{if } \prod_i\epsilon_{i,m}=1\\
0
\end{cases}$$
