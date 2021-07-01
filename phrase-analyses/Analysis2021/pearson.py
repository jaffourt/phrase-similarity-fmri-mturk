from Data import Data
from Stats import Stats
import numpy as np


def pearson(a, b):
    """

    :param a: 1-N vector
    :param b: 1-N vector

    """
    n = len(a)
    cov = sum((a - np.mean(a)) * (b - np.mean(b))) / n  # n-1
    r = cov / (np.std(a) * np.std(b))
    return r


data = Data.as_numpy(expt='Expt', map_type='spmT', mask='language_ROIs.nii')

stats = Stats(input=data.dataframe, function=pearson, conditions=data.condition_names)
#
# stats.run(participant_average=True)
# stats.plot_corr_matrix(save=True, title='pearson')
#
# stats.run(canonical_average=True)
# stats.plot_corr_matrix(save=True, title='pearson_averaged')

stats.run(average_residuals=True)
stats.plot_bar_canonical(show=True, save=True, title='residuals_pearson_averaged_bar')