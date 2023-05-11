---
title: polyadenylation
created: 2023-05-04T07:54:51.125Z
modified: 2023-05-10T09:50:44.299Z
---

# polyadenylation

The expression is proportional to the Multiplicity times the factor
$$f=\frac{c\mathrm{UTR}_{len}}{a\mathrm{UTR}_{len}+c\mathrm{UTR}_{len}}\geq0.$$
$f$ is zero mainly due to IS.




$$\mathbb{E}_n(\frac{e_{n,g,c}}{s_n}) =  \mu_g \sum_jM^{\epsilon}_{c,g,j}.$$


$u$ is the new polyadenylation position in respect to the gene intial position

If $u_{max}$ is intended as the length of the gene or the closer p(A)
$$\epsilon = \begin{cases}
u/u_{max} & u\leq u_{max}\\
0 & \text{degradation}
\end{cases}$$

If $u_{max}$ is the length of the gene till the furthest p(A)
$$\epsilon = \begin{cases}
u/u_{max} & u\leq u_{max}\\
u_{max} & \text{no degradation}\\
0 & \text{degradation}
\end{cases}$$


## parameters for the model

| k=genotype ID | g=gene | S=start     | P=PA pos    | A=APA pos   | stop-codon  | $$\epsilon_{3'}=1-\\\theta(A-P)$$ | $$\epsilon_{full}=\epsilon_{3'}\\\tfrac{A-S}{P-S}$$ |
| ------------- | ------ | ----------- | ----------- | ----------- | ----------- | --------------------------------- | --------------------------------------------------- |
| 1             | P53    | chr1:445566 | chr1:449999 | chr1:448999 | chr1:849999 | abcd                              | abcd                                                |
| 1             | AKT1   | chr1:845566 | chr1:849999 | chr1:848999 |             |                                   |                                                     |
| 2             | AKT1   | chr1:845566 | chr1:849999 | chr1:999999 |             |                                   |                                                     |



