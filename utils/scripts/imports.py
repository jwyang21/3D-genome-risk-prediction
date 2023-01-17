import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import os
import glob

from collections import Counter, defaultdict
from tqdm import tqdm
from scipy import stats

mpl.rcParams['figure.dpi'] = 150
plt.rc('font', family='FreeSans', size=7)

plt.rc('figure', figsize=(1.5, 1.5))

pd.set_option("display.max_columns", None)

def save_figures(f, exts=['png', 'pdf']):
    for ext in exts:
        plt.savefig(f + f'.{ext}', dpi=300, bbox_inches='tight', transparent=True)
