capwire
=======

[![Build Status](https://travis-ci.org/mwpennell/capwire.png?branch=master)](https://travis-ci.org/mwpennell/capwire)

### Estimate population census size from non-invasive genetic sampling

This R package implements the Capwire population size estimators developed by Craig Miller, Paul Joyce and Lisette Waits [(paper here)](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2005.02577.x/abstract). They were designed for non-invasive genetic samplings which is different than standard mark-recapture surveys in that the same individual can be detected multiple times (often many times) in the same sampling session. Details on the estimators can be found in the original paper.

This package is the collaborative work of [Matthew Pennell](http://mwpennell.github.io), Carisa Stansbury, [Lisette Waits](http://www.uidaho.edu/cnr/fishwild/lisettewaits), and [Craig Miller](http://www.uidaho.edu/sci/math/faculty/craigmiller). The package is described in [this paper](http://mwpennell.github.io/pdfs/pennell-mer-2012.pdf).

To install the package directly from github, the easiest way is to install using
[devtools](https://github.com/hadley/devtools). Install `devtools`, then type

```
library(devtools)
install_github("mwpennell/capwire")
```

We hope this package is useful. If you run into issues or have any suggestions, please feel free to contact us: mwpennell@gmal.com.

