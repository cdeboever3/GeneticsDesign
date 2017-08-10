# GeneticsDesign

This is a fork of the Bioconductor mirror. I've added support for unselected
controls in case/control studies that matches the output from
http://zzz.bwh.harvard.edu/gpc/. I've also added a small python wrapper `gdpy`
that relies on rpy2.

You can install this package into R using
```
library(devtools)
install_github("cdeboever3/GeneticsDesign")
```

If you want to use the python wrapper, you'll need to first install the R
package into whatever R your rpy2 uses, then run 
```
cd cdpy
python setup.py install
```
