from Data import Data
from Stats import Stats
import numpy as np


def spearman(a, b):
    """

    :param a: 1-N vector
    :param b: 1-N vector

    """
    n = len(a)

    temp = a.argsort()
    rank_a = np.empty_like(temp)
    rank_a[temp] = np.arange(len(a))

    temp = b.argsort()
    rank_b = np.empty_like(temp)
    rank_b[temp] = np.arange(len(b))

    cov = sum((rank_a - np.mean(rank_a)) * (rank_b - np.mean(rank_b))) / n  # n-1
    r = cov / (np.std(rank_a) * np.std(rank_b))
    return r


data = Data.as_numpy(expt='Expt', map_type='spmT', mask='language_ROIs.nii')
stats = Stats(input=data.dataframe, function=spearman, conditions=data.condition_names)

# horizontal bar plot
# stats.run(method='canonical_average')
# stats.plot_bar_canonical(show=True, title='spearman_averaged_bar')

# correlation matrix
# stats.plot_corr_matrix(save=True, title='spearman_averaged')
# stats.save_csv(values='individual', title='spearman_averaged', identifiers=data.subjectIDs)

# within and between horizontal bar plot
# stats.plot_within_between_averages(show=True, save=True, title='within_v_between')

# stats.run(method='average_residuals')
# stats.plot_bar_canonical(show=True, save=True, title='residuals_spearman_averaged_bar')

stats.run(method='average_preproc')
stats.plot_bar_canonical(show=True, save=False, title='preprocessed_spearman_averaged_bar')