# Functional $\alpha$-transformations for the statistical analysis of density functions

## Abstract
The nature of probability density functions can cause problems when the classical methods of Functional Data Analysis are adopted, since square-integrability may fail and the constraints may not be taken into consideration. 
In that sense, the compositional log-ratio approach and multivariate Aitchison geometry have been generalized to the functional case with the definition of the Bayes space, that embeds densities in a suitable geometry.
The goal of this work is to study the possibility of generalizing the multivariate $\alpha$-transformations to the functional case. 
The $\alpha$-transformations have been recently introduced as an alternative to the log-ratio approach in Compositional Data Analysis, and their main strengths are to include zeros in compositions and to control the degree of data modification through the parameter $\alpha$. 
$\alpha$-transformations allow the classical Functional Data Analysis methods to work on densities by transforming them into square-integrable functions and selecting the parameter in order to optimize a performance metric depending on the problem at hand. 
Moreover, the log-ratio transformation can be retrieved as a limit case of the $\alpha$-transformation, hence the parameter selection process includes it as well.
A real case study is presented to test the functioning of the procedure in a regression and time-series modelling scenario, together with the comparison with the log-ratio approach.

## Author
Simone Giacomello      simone.giacomello@mail.polimi.it

MSc Thesis. Mathematical Engineering - Statistical Learning @Politecnico di Milano. Academic year 2023/2024

Advisor: Prof.ssa Alessandra Menafoglio

## Code 
In this repository you can find the R scripts that produced the analysis. The complete set of data can be found at https://polimi365-my.sharepoint.com/:f:/g/personal/10620752_polimi_it/ElGqmXWeVb9Eml8nekSPr1oBLD-81l-YEt4uDgvcbqCWcQ?e=mEUXGi

## References
Scimone R, Menafoglio A, Sangalli LM, Secchi P. A look at the spatio-temporal mortality patterns in Italy during the COVID-19 pandemic through the lens of mortality densities. Spat Stat. 2022 Jun;49:100541. doi: 10.1016/j.spasta.2021.100541. Epub 2021 Sep 28. PMID: 34631399; PMCID: PMC8486968.

The original study code can be found at https://github.com/RiccardoScimone/Mortality-densities-italy-analysis
