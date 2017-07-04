import itertools
import os

import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects import r

# Source the GPC.R file.
_gpc = os.path.realpath(os.path.join(os.path.split(__file__)[0],
                                     '../../R/GPC.R'))
r('source("{}")'.format(_gpc))

def gpc_default(
    pA, 
    RRAa, 
    pD, 
    nCase, 
    ratio=1, 
    pB=None, 
    Dprime=1, 
    alpha=5e-8,
    unselected=False, 
    protective=False,
):
    """
    Wrapper around the GPC.default function that accounts for some of the bugs
    in the function.  ratio is nCase / nControl
    if pB is None, it will be set to equal pA
    I'm assuming additive traits.
    """
    gpc_default_fnc = robjects.r['GPC.default']

    if pB is None:
        pB = pA
    nControl = nCase / ratio

    if protective:
        # I use the case number here but invert the ratio.
        res = gpc_default_fnc(pA=pA, pD=pD, RRAa=RRAa, RRAA=RRAa * 2,
                              Dprime=Dprime, pB=pB, nCase=nCase, ratio=1 /
                              float(ratio), alpha=alpha, unselected=unselected,
                              quiet=True)
    else:
        # Due to how the control/case counts and ratio are implemented in
        # GPC.default, I have to provide nControl instead of nCase.
        res = gpc_default_fnc(pA=pA, pD=pD, RRAa=RRAa, RRAA=RRAa * 2,
                              Dprime=Dprime, pB=pB, nCase=nControl, ratio=ratio,
                              alpha=alpha, unselected=unselected, quiet=True)
    # Return the power estimate.
    return res[0][0]

def gpc_default_param_sweep(
    pA, 
    RRAa, 
    pD, 
    nCase, 
    ratio, 
    alpha=[5e-8],
    unselected=[False], 
    protective=[False],
    total_samples=None,
):
    """Arguments are similar to gpc_default except all arguments should be
    provided as lists. All parameter combinations will be calculated. If total
    samples is provided, ratio is ignored and will be calculated on the fly
    assuming nCase + nControl = total_samples."""
    itp = itertools.product(pA, RRAa, pD, nCase, ratio, alpha, unselected,
                            protective)
    results = []
    while True:
        try:
            [pA, RRAa, pD, nCase, ratio, alpha, unselected, protective] = \
                    list(itp.next())
            if total_samples is not None:
                ratio = float(nCase) / (total_samples - nCase)
            nControl = nCase / ratio
            power = gpc_default(pA, RRAa, pD, nCase, ratio, alpha=alpha, 
                                unselected=unselected, protective=protective)
            results.append([pA, RRAa, pD, nCase, nControl, ratio, alpha,
                            unselected, protective, power])
        except StopIteration:
            break
    columns=['allele_freq', 'odds_ratio', 'prevalence', 'num_cases',
             'num_controls', 'ratio', 'alpha', 'unselected', 'protective',
             'power']
    out = pd.DataFrame(results, columns=columns)
    return out
