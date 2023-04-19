---
title: Loci multiplicity, mutations and transcriptional profiling
created: 2023-03-31T07:09:14.787Z
modified: 2023-04-19T14:36:56.880Z
---

# Loci multiplicity, mutations and transcriptional profiling

Starting from a general assumption as direct extension of the central dogma of MB which if cells have multiple copies of genes than their expressions are proportional to their respective copy number, it follows that major cell transcriptional differences can be quantitatively explained by presence of cell populations with different karyotypes.

Hence, if a gene $g$ has an expression per unit of copy $\mu_g$, then the expected normalized expression of cells with copy number $C_{c,g}$ in that chromosome region containing the gene $g$ is:
$$\mathbb{E}_n(\frac{e_{n,g,c}}{s_n}) =  \mu_g C_{c,g}$$

for all genes, where $e_{n,g,c}$ is the number of reads of cell $n$ for gene $g$ with copy number $c$ and $s_n$ is the cell's sampling depth.

Clonealign, an R package, makes this assumption and, based on a GLMM, associate  transcript profiles of each cell to one of the defined copy number profiles.

There are various reasons why this assumption does not hold true. For example, SNPs and indels in a genes may result in nonsense or missense mutations  causing exonic disfunction as loss or gain of stop codons.
Furthermore, mutations may occurs before or after the copy number alterations increasing the variability in gene expressions among samples and single cells derived from the same source.

I propose the following change to the relation between average gene expression e multiplicity:
$$\mathbb{E}_n(\frac{e_{n,g,c}}{s_n}) =  \mu_g \sum_jM^{\epsilon}_{c,g,j}.$$

--------------

# Clonealign model and notation

### Transformation of Expressions and Clone Assignment
Gene expression transformation and inverse:

$\xi(x)=\log(1+e^x)$
$\xi^{-1}(y)=\log(e^y-1)$
$\xi_{\mathrm{soft}}^{-1}(y)=\log(1-e^{-|y|})+ \max(0,y)$

where $x$ is the number of reads in a given region/bin and $y$ is the log transformed expression.


Karyotype is assigned based on the following function:
$$\kappa_j(\Gamma_{i,j})=\begin{cases}
\argmax_i (\Gamma_{i,j}) & \text{if } \max_i(\Gamma_{i,j})>p_{min}\\
-1 & \text{otherwise unassigned}
\end{cases}$$

where $i$ is the clone of cell $j$

-------------

### Inputs

$Y$ is the expression of $N$ cells over $G$ genes:
$Y \in \mathbb{R}^{N \times G}$

$L$ is the copy number matrix  over the gene regions $G$ for $C$ clones
$L \in \mathbb{R}^{G \times C}$


$K \geq 0$ random components


The CN of $C$ clones per $V$ variants
$C_{\mathrm{allele}} \in \mathbb{R}^{V \times C}$

From variant calling analysis follows the reference and alternative count matrices as well as the coverage
$M_\mathrm{Cov} \in \mathbb{R}^{V \times N}$
$M_\mathrm{Ref} \in \mathbb{R}^{V \times N}$
$M_\mathrm{Alt} = M_\mathrm{Cov}- M_\mathrm{Ref} \in \mathbb{R}^{V \times N}$

----------

### Allele likelihood

The log-likelihood associated to the variant calling data for each cell $n$ given the allele specific copy number corresponding  to variant $v$ of clone $c$ is :  

$$
[\mathbb{1} \otimes \log\rho_{\mathrm{var}}]_{c,v,n} =
1_c \otimes \log \begin{cases}
\mathcal{B}_{\mathrm{Bin}}(M_\mathrm{Alt}, M_\mathrm{Cov},\alpha_\mathrm{mid},\beta_\mathrm{mid})_{v,n}
& [C_{\mathrm{allele}}^T]_{c,v}=2\\
\frac{1}{2} \mathcal{B}_{\mathrm{Bin}}(M_\mathrm{Alt}, M_\mathrm{Cov},\alpha_\mathrm{low},\beta_\mathrm{low})_{v,n}+
\frac{1}{2} \mathcal{B}_{\mathrm{Bin}}(M_\mathrm{Alt}, M_\mathrm{Cov},\alpha_\mathrm{hi},\beta_\mathrm{hi})_{v,n} & \text{otherwise}
\end{cases}
$$

and considering any variant the log-likelihood is:

$$
\mathcal{L}_{\mathrm{allele}}^T(M_{\_}| \alpha,\beta, C_\mathrm{allele})=\sum_v [\mathbb{1}\otimes \log\rho_{\mathrm{var}}]_{c,v,n}
$$

$\mathcal{L}_{\mathrm{allele}} \in \mathbb{R}^{N \times C}$

--------------
### Linear Mixed Model

Given the model matrix $X \in \mathbb{R}^{N \times P}$ together with the covariate coefficients $\beta \in \mathbb{R}^{G \times P}$,
the projection along the first $K$ principal components of the $N$ cells
$\mathrm{PCs}= \mathrm{PCA}_K(\xi(Y))+ \mathcal{N}(\mathbb{\mu}_0,\mathbb{\sigma}_{1/2}) \in \mathbb{R}^{N \times K}$ used in deriving the random effect matrix $W \in \mathbb{R}^{G \times K}$ and the random corresponding components $\psi=\mathrm{PCs}$,
then the  fixed and random effects can be described as:
$f=\psi W^T+ X\beta^t \in \mathbb{R}^{N \times G}$.

----------

### Clonal Likelihood

Given a library average expression of gene $g$ [per unit of CN?]:

$[\mu]_{g}= E_n(\frac{Y_{n,g}}{E_g(Y_{n,g})})\otimes \mathbb{1} \in \mathbb{R}^{G}$,

and an initial variance [why is it 1 and why does it not depend on CN?]

$\sigma_i=\mathbb{1} \in \mathbb{R}^{G}$,

then the distribution of expression from each copy of gene $g$ is:
$\mathcal{Q}^\mu=\mathcal{N}(\xi^{-1}(\mu),\sigma_i) \in \mathbb{R}^{G}$

Note^1^: there is no reasons why $Y$ should be back transformed.
Note^2^: $\mu$, or rather $\xi^{-1}(\mu)$, and $\sigma_i$ are part of the hidden variable/code $z$ (when $z$ is fixed, then we are describing the likelihood of data $Y$).

Given a set of samples:
$S = \text{\# samples}$

each with iid expressions per gene  
$\mu_s(S)=\mathbb{1} \otimes \mathcal{Q}^\mu \in \mathbb{R}^{S \times G}$


then, the probabilities associated to the events of reads falling over a gene $g$ specific for a cell $n$ associated to a clone $c$ and having fixed/random effects $f$ and copy number derived from the value $L_{g,c}$ is given by:

$\mu_{s,c,n,g}=\frac{\mu_s(S)L_{g,c} e^{f_{n,g}}} {\sum_g\mu_s(S)L_{g,c} e^{f_{n,g}}}$.


Considering the total reads for a cell $n$
$s_n=\sum_g Y_{n,g}$,

the distribution of having an expression follows a multinomial distribution:
${{y_\mathrm{pdf}}_{s,c,n}}^g \sim \mathrm{Multinom(s_n,\mu_{s,c,n,g})}$

which satisfies the relation:
$\sum_g{y_{pdf}}^g=1$.

Eventually, the log-likelihood of having an expression $Y$ given the clones $L$ is defined as:

$\mathcal{L}_\mathrm{clone}(Y)= \log({y_{\mathrm{pdf}}}^g(Y)) \in ?$

-----------

### (i) Expected log-likelihood $\mathbb{E}_q[\log(\rho(y | z, \theta))]$

$\mathcal{L}={\mathcal{L}_\mathrm{clone}}_{s,c,n} + {\mathcal{L}_\mathrm{allele}^T}_{c,n}$

$\Gamma_{\mathrm{logit}} \in \mathrm{R}^{N \times C}$
$\Gamma_{n,c} = \frac{e^{{\Gamma_{\mathrm{logit}}}_{n,c}}}{\sum_c{e^{{\Gamma_{\mathrm{logit}}}_{n,c}}}}$


$\mathbb{E}_q[\log(\rho(y | \theta))] = \sum_{n,c} \Gamma_{n,c}\ {\mathbb{E}_s(\mathcal{L})^T}_{n,c}$







### BhO?

$\chi \in \mathbb{R}^{K}$

$\chi_e=\exp(\chi) \in \mathbb{R}^{K}$

$\log\text{-}\alpha=\vec{\alpha}-\mathbb{1} \otimes log(\mathrm{tr}(e^{\vec{\alpha}})) \in \mathbb{R}^{C}$

-----------

### (ii) $\mathbb{E}_q[\log(\rho(\theta))]$



$\log\text{-}\rho_W = \sum \log [\mathcal{N}(\mu_0,\frac{1}{\chi_e})](W) \in \mathbb{R}^{C \times G \times K}$



$\log\text{-}\rho_\chi = \sum \log[\mathcal{G}(2,1)](\chi) \in \mathbb{R}^{K}$

$\log\text{-}\rho_\psi = \sum\log[\mathcal{N}(\mu_0,\sigma_1)](\psi) \in \mathbb{R}^{N \times K}$

Given data association to any of the clones $C$ is a categorically distributed likelihood, then a Dirichlet prior implies the posterior of of clones given the data is a Dirichlet distribution.

$$
\mathbb{E}_q[\log(\rho(\theta)]=\\
\sum_{n,c} \log\text{-}\alpha_c \ \ \Gamma_{n,c}+\\
\sum_{s,g} \frac{1}{S} \log[\mathcal{N}(\mu_1, \sigma_1)](\log(\mu_s(S)))+\\
\sum_c \log\mathcal{D}(\bigotimes_c^C \frac{1}{C})(e^{\log\text{-}\alpha})+\\
10^{-3}
$$

$$
\mathbb{E}_q(\log[\rho(\theta)]=
\mathbb{E}_q[\log(\rho(\theta)]+
W_{\log\rho}+
\chi_{\log \rho}+
\psi_{\log \rho}
$$




### (iii) $\mathbb{E}_q[\log(q)]$
$$
\mathbb{E}_q[\log(q)] =\\
\sum \mathbb{E}_s\log[\mathcal{Q}](\mu_s(S))+\\
\sum \Gamma_\mathrm{logits}\log\text{-softmax}(\Gamma_\mathrm{logits})
$$

$$
q(z_n=c)=\frac{e^{{\Gamma_\mathrm{logits}}_{n,c}}}{\sum_c e^{\Gamma_\mathrm{logits}}}=\mathrm{softmax}({\Gamma_\mathrm{logits}}_{n})_c
$$
$$
q(\vec{z}=\vec{c})=\prod_n\frac{e^{{\Gamma_\mathrm{logits}}_{n,c}}}{\sum_c e^{{\Gamma_\mathrm{logits}}_{n,c}}}
$$

$$
\log[q(\vec{z} =\vec{c})] = \sum_n \{{\Gamma_\mathrm{logits}}_{n,c}-\log[{\sum_c e^{{\Gamma_\mathrm{logits}}_{n,c}}}]\}=\sum_n\log\text{-}\mathrm{softmax}({\Gamma_\mathrm{logits}}_{n})
$$



$$
\mathbb{E}_q[\log(q(\vec{z}))] = \sum_{n,c} \frac{e^{{\Gamma_\mathrm{logits}}_{n,c}}}{\sum_ce^{{\Gamma_\mathrm{logits}}_{n,c}}}\log\text{-softmax}({\Gamma_\mathrm{logits}}_{n})_c =\sum_{n,c} \Gamma_{n,c}\log\text{-softmax}({\Gamma_\mathrm{logits}}_{n})_c
$$

$$
q(\vec{\mu})=\mu_s(S)
$$


$$
\log[q(\vec{\mu})]=\log[\mu_s(S)]=\sum_g \log[\mathcal{N}(\xi^{-1}(\mu_g),\sigma_1)]
$$


$$
\mathbb{E}_q\log[q(\vec{\mu})]=\frac{1}{S}\sum_{g,s} \mathcal{N}(\xi^{-1}(\mu_g),\sigma_1) \log[\mathcal{N}(\xi^{-1}(\mu_g),\sigma_1)]
$$

elbo <- EE_p_y + E_log_p_p - E_log_q




--------------

# VAE



given the the data $x$ the variables $z$ and the generative parameters $\theta$:
$p(z;\theta)$ is the prior on the latent variable $z$
$p(x|z;\theta)$ is the likelihood iid sampled data $x$ given $z$
$p(z|x;\theta)$ is the posterior
$p(x)=\int p(x|z)p(z)dz$ is the marginalization

$q(z|x;\phi)$ is a recognition model or a probabilistic encoder where $z$ is the code and $\phi$ are the model variational parameters.
$p(x|z;\theta)$ is the probabilistic decoder which from a code $z$ reproduce the distribution of possible $x$.


$$
KL(q(z|x)||p(z,x)) = \int dz \ q(z|x)\ og \frac{q(z|x)}{p(z|x)}\\
=\int dz \ q(z|x)\ \log \frac{q(z|x) p(xx|)}{p(x,z)}\\
=\int dz\ q(z|x)\ [\log \frac{q(z|x)}{p(x,z)}+\log{p(x)}]\\
=\log{p(x)} + \mathbb{E}_{q(z|x)}[\log\frac{q(z|x)}{p(x,z)}]\\
=\log{p(x)}-\mathcal{L}
$$

$$
\mathcal{L} = \mathbb{E}_{q(z|x)}[\log\frac{p(x,z)}{q(z|x)}]=\\
\mathbb{E}_{q(z|x)}[\log\frac{p(z)}{q(z|x)}]+\mathbb{E}_{q(z|x)}[\log{p(x|z)}]=\\
{}\\
\mathbb{E}_{q(z|x)}[\log p(z)] - \mathbb{E}_{q(z|x)}[\log{q(z|x)}] + \mathbb{E}_{q(z|x)}[\log{p(x|z)}] = \\
{}\\
-\underbrace{KL(q(z|x)||p(z))}_\text{regularizer}+\underbrace{\mathbb{E}_{q(z|x)}[\log{p(x|z)}]}_\text{expected negative error}
$$


$$
\nabla_\phi E_{q_\phi(z)}(f(z)) =
\nabla_\phi \int dz q_\phi(z) f(z) =\\
\int dz \nabla_\phi q_\phi(z) f(z) =\\
\int dz \  f(z) \ q_\phi(z) \nabla_\phi \log q_\phi(z) =\\
\int dz \  f(z) \ q_\phi(z) \nabla_\phi q_\phi(z) / q_\phi(z)
$$

$$
z=g(\epsilon, x, \phi)\sim q(z|x;\phi)\\
z_{i,l}=g(\epsilon_l, x_i, \phi)\sim q(z|x_i;\phi)
$$

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


## Model

For each cell $c$ there are different mutational profiles $J_c$. For each region $r$ with CN equal to $C_{c,r}$, the number of mutational profiles $J_{c,r}\leq J_c \leq C_{c,r}$.
This relation hold for all cells. Hence, $J_{r}\leq J \leq C_{r}$.




Let the mutational profile $j$ be defined as $\vec{m}_j$ where in each position there is the transition from the original base $A$ to the mutated base $a$, $A\rightarrow a$ (indels included).

A mutational profile projected in a region $r$ is give by:
$$\mathbb{P}_r(\vec{m}_j)=\vec{m}_{j,r}$$

Each mutational profile $\vec{m}_{j,r}$ has a multiplicity $M$ which depends on the region $r$, $M_{c,r,j}$, such that for all $r$ of any size
$$\sum_{j=1}^{J_r} M_{c,r,j}=C_{c,r} \quad\quad \forall r.$$

The size of regions biologically relevant are those with constant CN or where genes are defined or the single alleles.

Let a mutation in $i$ of a profile $\vec{m}_{j,r}$ have an effectiveness $\epsilon_{j,r,i} \in \{0,1\}$, then

$$\vec{\epsilon}_{j,r}=\epsilon[\mathbb{P}_r(\vec{m}_{j})]$$
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
