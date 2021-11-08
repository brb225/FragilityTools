
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FragilityTools

The goal of FragilityTools is to provide a toolbox for researchers
interested in the fragility index. The toolbox covers methodologies
developed across three articles (Baer, Gaudino, Fremes, et al. 2021;
Baer, Gaudino, Charlson, et al. 2021; Baer, Fremes, et al. 2021).

## Installation

You can install the package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("brb225/FragilityTools")
```

## Examples

The first step of any analysis is to load the package with:

``` r
library(FragilityTools)
```

The fragility index was created by Walsh et al. (2014) to study the
robustness of statistical decisions. They defined the fragility index to
be the smallest count of patient outcome modifications needed to reverse
statistical significance. At first, the definition was only for
initially statistical significant decisions, but it was quickly extended
to initially insignificant decisions as well. Suppose that we observe
the data below from a two group clinical trial with a binary (or
dichotomous)
endpoint:

``` r
mat <- matrix(c(40, 100, 20, 110), nrow=2, dimnames=list(c('treatment', 'control'), c('event', 'nonevent')))
mat
#>           event nonevent
#> treatment    40       20
#> control     100      110
```

On this contingency table, Fisher’s exact test returns a p value equal
to 0.0124 for testing for a treatment effect. With the default cutoff
0.05, the treatment effect here is statistically significant. This p
value can be difficult for a general audience to interpret, and the
fragility index (being a patient count) is much more interpretable. We
can find the fragility index with:

``` r
out <- bin.fi(mat, test='fisher')
out$mat
#>           event nonevent
#> treatment    37       23
#> control     100      110
```

The function `bin.fi`, short for “binary fragility index”, calculates
the fragility index by finding an alternative contingency table which
reverses statistical significance and has the fewest possible rowwise
modifications. Above we see that the returned contingency table has 3
fewer events in the treatment group. Since that contingency table
reverses statistical significance (with p = 0.0585605) and no
contingency table with fewer modifications does, the fragility index
equals 3. Therefore there exists three patients for which their outcome
being different would have led to a different statistical conclusion.
Out of a trial with 270 patients, this may be a surprisingly small
number, suggesting that the trial’s statistical conclusion of a
treatment effect could be fragile.

We now review the core functions in the package by separately
considering the methodological development in each of the above
articles. We will only briefly review the functions: more exentensive
documentation and examples are available in the help files for each
function and the executable `.R` files.

### Sample size calculation

The fragility index based sample size calculations in Baer, Gaudino,
Fremes, et al. (2021) can be reproduced using the code in
`FragilityTools/exec/`. The files `/sscalc_FAMOUS_analysis.R` and
`/sscalc_FAME_analysis.R` contain specifications of the statistical
parameters needed to perform the sample size calculations for the FAMOUS
and FAME studies. Each rely on the script
`/sscalc_example_simulation_script.R` to perform the core calculations.

The function `general.fi.samplesize` performs the sample size
calculation for given parameters, while the function
`genreal.fi.sscurve` conveniently loops over `general.fi.samplesize` to
perform sample sizes at varying fragility index tolerances. The
functions `get.rejection.rates` and `get.traditional.ss.params` to
understand the statistical properties of the designed tests.

### The incidence and generalized fragility indices

The incidence fragility index examples in Baer, Gaudino, Charlson, et
al. (2021) can be reproduced using the code in `FragilityTools/exec/`.
The file `/incidencefi_examples.R` is a script which provides all the
output concerning incidence fragility indices that appears in the
article.

The function `bin.fi` can calculate the incidence fragility index given
a constrain on the per-patient sufficiently likely threshold through the
`q` argument. The incidence fragility indices for each possible
threshold `q` can be calculated using the function `bin.fi.incidence`.
The output can be conveniently visualized using `incidence.plot`.

The generalized fragility index examples in Baer, Gaudino, Charlson, et
al. (2021) can be reproduced using the code in `FragilityTools/exec/`.
The file `/generalizedfi_examples.R` is a script which provides all the
output concerning generalized fragility indices that appears in the
article and several other examples.

The function `greedy.fi` is the workhorse which is under the hood of
each of the following functions. It efficiently approximates fragility
indices with a greedy algorithm. The function `surv.fi` allows for
calculating fragility indices with time-to-event outcomes. The function
`ttest.fi` allows for calculating fragility indices corresponding to
one-sample t tests. The function `binmeta.fi` allows for calculating
fragility indices corresponding to meta analyses with 2 x 2 contingency
tables. The function `glm.fi` returns the coefficient table from
`stats::glm` augmented with fragility indices for each coefficient test.
The function `glm.gaussian.covariate.fi` allows for calculating
fragility indices when a gaussian distributed covariate is modified in a
`glm` such a logistic regression.

### The LTFU-aware fragility indices

The LTFU-aware fragility index examples in Baer, Fremes, et al. (2021)
can be reproduced using the code in `FragilityTools/exec/`. The file
`/ltfufi_examples.R` is a script which provides all the output that
appears in the article.

The function `ltfu.fi` calculates the LTFU-aware fragility indices and
outputs a convenient visualization.

## References

<div id="refs" class="references">

<div id="ref-baer2021ltfu">

Baer, Benjamin R., Stephen E. Fremes, Mario Gaudino, Mary Charlson, and
Martin T. Wells. 2021. “On Clinical Trial Fragility Due to Patients Lost
to Follow Up.” *BMC Medical Research Methodology, in Press*.

</div>

<div id="ref-baer2021likely">

Baer, Benjamin R., Mario Gaudino, Mary Charlson, Stephen E. Fremes, and
Martin T. Wells. 2021. “Fragility Indices for Only Sufficiently Likely
Modifications.” *Proceedings of the National Academy of Sciences, in
Press*.

</div>

<div id="ref-baer2021samplesize">

Baer, Benjamin R., Mario Gaudino, Stephen E. Fremes, Mary Charlson, and
Martin T. Wells. 2021. “The Fragility Index Can Be Used for Sample Size
Calculations in Clinical Trials.” *Journal of Clinical Epidemiology*
139: 199–209.

</div>

<div id="ref-walsh2014fragility">

Walsh, Srinathan, McAuley, Mrkobrada, Levine, Ribic, Molnar, et al.
2014. “The Statistical Significance of Randomized Controlled Trial
Results Is Frequently Fragile: A Case for a Fragility Index.” *Journal
of Clinical Epidemiology*.

</div>

</div>
