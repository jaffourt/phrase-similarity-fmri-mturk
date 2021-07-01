import sys

import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import numpy as np
from pathlib import Path
from tqdm import tqdm

np.seterr(divide='ignore')


class Stats:
    def __init__(self, **kwargs):
        self.plot_data = None
        self.condition_averages = None
        self.semantic_classes = None
        self.function = kwargs.get('function', None)
        self.conditions = kwargs.get('conditions', [])
        self.c_ind = self._cond_indices()
        self.input = kwargs.get('input', None)

        self.run() if kwargs.get('run') else None

    def run(self, **kwargs):
        def invalid_name():
            raise AttributeError("Invalid method name or None")

        run_object = {
            'canonical_average': self.average_participants_and_conditions,
            'participant_average': self.average_participants,
            'condition_average': self.average_conditions,
            'within_phrase': self.correlate_within_phrase,
            'between_phrase': self.correlate_between_phrase,
            'average_residuals': self.correlate_between_residuals,
            'average_preproc': self.correlate_between_preprocessed,
            'None': invalid_name
        }
        self.plot_data = run_object[kwargs.get('method', 'None')]()

    def average_participants_and_conditions(self):
        self.semantic_classes = self.c_ind.keys()
        average = []
        # participant average:
        for m in tqdm(self.input):
            # condition averages
            values = self._square((self._apply_method(m)))
            # 4,4 -> manually sort condition indices together
            # np.tanh(np.arctanh(x).mean())
            values = [x.mean() for x in self.split(values, 4, 4)]
            # 5,5 -> manually spread into 5 x 5 phrase type matrix
            out = list(np.array(values).reshape(5, 5))
            average.append(out)
        # np.array(average), np.tanh(np.arctanh(np.array(average)).mean(axis=0))
        self.condition_averages = np.array(average)
        return np.array(average).mean(axis=0)

    def average_participants(self):
        average = []
        # participant average:
        for m in tqdm(self.input):
            average.append(list(self._apply_method(m)))
        return np.tanh(np.arctanh(np.array(average)).mean(axis=0))

    def average_conditions(self):
        # todo - iter participants
        values = squareform(self._apply_method(self.input[0]))
        values = [np.tanh(np.arctanh(x.mean())) for x in self.split(values, 4, 4)]
        return np.array(values).reshape(5, 5)

    def correlate_within_phrase(self):
        d = {}
        for phrase in self.conditions:
            if phrase.split('_')[0] not in d:
                d[phrase.split('_')[0]] = [self.conditions.index(phrase)]
            else:
                d[phrase.split('_')[0]].append(self.conditions.index(phrase))
        average = []
        for m in tqdm(self.input):
            data = []
            for c in d:
                indices = d[c]
                data.extend(list(self._apply_method(m[indices])))
            average.append(np.array(data).mean())
        return average

    def correlate_between_phrase(self):
        average = []
        for m in tqdm(self.input):
            average.append(self._apply_method(m).mean())
        return average

    def correlate_between_preprocessed(self):
        from sklearn import preprocessing
        d = {}
        for phrase in self.conditions:
            if phrase.split('_')[0] not in d:
                d[phrase.split('_')[0]] = [self.conditions.index(phrase)]
            else:
                d[phrase.split('_')[0]].append(self.conditions.index(phrase))

        preproc = []
        for m in self.input:
            data = np.empty_like(m)
            for c in d:
                indices = d[c]
                phrase = np.array(m[indices])
                scaler = preprocessing.StandardScaler().fit(phrase)

                # self._debug_plot(phrase)
                data[indices] = list(scaler.transform(phrase))
            preproc.append(data)

        raw = self.input
        self.input = np.array(preproc)
        plot_data = self.average_participants_and_conditions()
        self.input = raw

        return plot_data

    def correlate_between_residuals(self):
        d = {}
        for phrase in self.conditions:
            if phrase.split('_')[0] not in d:
                d[phrase.split('_')[0]] = [self.conditions.index(phrase)]
            else:
                d[phrase.split('_')[0]].append(self.conditions.index(phrase))

        residuals = []
        for m in self.input:
            data = np.empty_like(m)
            for c in d:
                indices = d[c]
                phrase = np.array(m[indices])
                # self._debug_plot(phrase)
                data[indices] = phrase - phrase.mean(axis=0)
            residuals.append(data)

        raw = self.input
        self.input = np.array(residuals)
        plot_data = self.average_participants_and_conditions()
        self.input = raw

        return plot_data

    def _apply_method(self, data):
        return pdist(data, self.function)

    def _cond_indices(self):
        # condition average:
        d = {'base': [],
             'meaningpres': [],
             'prepchange': [],
             'ajchange': [],
             'nchange': []}
        for i in range(0, len(self.conditions)):
            if "BASE" in self.conditions[i]:
                d['base'].append(i)
            elif "MEANINGPRES" in self.conditions[i]:
                d['meaningpres'].append(i)
            elif "PREPCHANGE" in self.conditions[i]:
                d['prepchange'].append(i)
            elif "AJCHANGE" in self.conditions[i]:
                d['ajchange'].append(i)
            elif "NCHANGE" in self.conditions[i]:
                d['nchange'].append(i)
        return d

    """
        Here begin the methods for visualization. 
    """

    @staticmethod
    def _debug_plot(data):

        n = data.shape[0]

        avg = data.mean(axis=0)
        x_max = 300
        colors = plt.cm.bone(np.linspace(0, 0.75, n))  # stick to whitish-grey

        ax1 = plt.subplot(211)
        for i in range(n):
            plt.plot(data[i].T[0:x_max], color=colors[i], linewidth=0.5)

        plt.setp(ax1.get_xticklabels(), fontsize=6)
        # plt.plot(avg.T[0:x_max], 'r', linewidth=2.5)

        # # linear regression
        # x = np.array(list(range(data.shape[1]))).reshape(-1, 1)
        # lm = linear_model.LinearRegression()
        # avg = []
        # for i in range(n):
        #     lm.fit(x, data[i])
        #     avg.append(list(lm.predict(x)))

        average_fit = np.array(data).mean(axis=0)
        residuals = data - average_fit

        plt.plot(average_fit[0:x_max], 'r', linewidth=2.5)

        # share x only
        ax2 = plt.subplot(212, sharex=ax1)
        for i in range(n):
            plt.plot(residuals[i][0: x_max], linewidth=0.5)

        # make these tick labels invisible
        plt.setp(ax2.get_xticklabels(), visible=False)

        plt.show()

    def _plot_helper(self, arg):
        # reshape the matrix to be square
        if self.plot_data is None and (arg == 1 or arg == 2):
            raise ValueError("average matrix not available for plotting")
        elif arg == 1 or arg == 2:
            if not self.plot_data.ndim > 1:
                data = self._square(self.plot_data)
            elif self.plot_data.ndim > 1 and not self.plot_data.shape[0] == self.plot_data.shape[1]:
                raise ValueError("asymmetrical matrix has 2 or more dimensions")
            else:
                data = self.plot_data

        if self.condition_averages is not None and arg == 2:
            return self.condition_averages

        if arg == 3:
            self.run(method='within_phrase')
            data = {'within': self.plot_data}
            self.run(method='between_phrase')
            data['between'] = self.plot_data
            return data

        # Who is requesting this method?
        return data

    def plot_corr_matrix(self, **kwargs):
        title = kwargs.get('title', '')
        output = Path('figures/').joinpath(title)

        # reshape the matrix to be square
        data = self._plot_helper(1)

        fig, ax = plt.subplots()
        img = ax.imshow(data, cmap='Greys')

        if self.semantic_classes is None:
            self.semantic_classes = self.conditions

        ax.set_xticks(np.arange(len(self.semantic_classes)))
        ax.set_yticks(np.arange(len(self.semantic_classes)))

        ax.set_xticklabels(self.semantic_classes)
        ax.set_yticklabels(self.semantic_classes)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8)

        # Loop over data dimensions and create text annotations.
        if len(self.semantic_classes) > 5:
            fontsize = 3
        else:
            fontsize = 5
        for i in range(len(self.semantic_classes)):
            for j in range(len(self.semantic_classes)):
                ax.text(j, i, round(data[i, j], 3), ha="center", va="center", color="w", fontsize=fontsize)

        fig.colorbar(img)
        img.set_clim(vmin=0, vmax=1)
        ax.set_title(kwargs.get('title', ''))
        fig.tight_layout()
        if kwargs.get('show', False):
            plt.show()
        if kwargs.get('save', False):
            plt.savefig(output.parent / (output.name + '.png'), dpi=600)

    def plot_bar_canonical(self, **kwargs):
        if self.condition_averages is None:
            raise ValueError("Canonical bar plot requires average=True")

        title = kwargs.get('title', '')
        output = Path('figures/').joinpath(title)

        # we are only interested in the BASE comparison for this plot
        ind_data = self._plot_helper(2)
        base = []
        for ind in ind_data:
            base.append(list(ind[0][1:5]))

        base = np.array(base)
        mean = np.mean(base, axis=0)
        se = np.std(base, axis=0) / np.sqrt(base.shape[0])

        fig, ax = plt.subplots()

        # plot mean phrase correlation to BASE
        conditions = ['meaning preserve', 'prep change', 'adj change', 'noun change']
        y_pos = np.arange(len(conditions))
        ax.barh(y_pos, mean, xerr=se, capsize=2, align='center', color='dimgray', zorder=0)

        # plot each phrase correlation to BASE
        ax.scatter(base, np.array([y_pos] * len(base)), zorder=1, color='black', alpha=0.3, marker='.')

        # plot in the original G&T style
        ax.invert_yaxis()

        # format axis
        ax.set_yticks(y_pos)
        ax.set_yticklabels(conditions)
        ax.set_title(title)
        ax.set_xlabel('Correlation to BASE')
        ax.set_xlim(right=1)

        fig.tight_layout()
        if kwargs.get('show', False):
            plt.show()
        if kwargs.get('save', False):
            fig.savefig(output.parent / (output.name + '.png'), dpi=600)

    def plot_within_between_averages(self, **kwargs):
        title = kwargs.get('title', '')
        output = Path('figures/').joinpath(title)

        data = self._plot_helper(3)

        base = [data['within'], data['between']]
        base = np.array(base).T
        mean = np.mean(base, axis=0)
        se = np.std(base, axis=0) / np.sqrt(base.shape[0])

        fig, ax = plt.subplots()

        # plot each set
        conditions = ['Within-phrase', 'Between-phrase']
        y_pos = np.arange(len(conditions))
        ax.barh(y_pos, mean, xerr=se, capsize=2, align='center', color='dimgray', zorder=0)
        ax.scatter(base, np.array([y_pos] * len(base)), zorder=1, color='black', alpha=0.3, marker='.')

        # plot in the original G&T style
        ax.invert_yaxis()

        # format axis
        ax.set_yticks(y_pos)
        ax.set_yticklabels(conditions)
        ax.set_title(title)
        ax.set_xlabel('Average spearman correlation')
        ax.set_xlim(right=1)

        fig.tight_layout()
        if kwargs.get('show', False):
            plt.show()
        if kwargs.get('save', False):
            fig.savefig(output.parent / (output.name + '.png'), dpi=600)

    def save_csv(self, **kwargs):
        title = kwargs.get('title', '')
        output = Path('figures/').joinpath(title)
        values = kwargs.get('values', None)
        if values == 'individual':
            if self.condition_averages is None:
                raise ValueError("Individual values requires average=True")
            a = self.condition_averages
            data = np.reshape(a, (a.shape[0] * a.shape[1], a.shape[2]))
            # print(kwargs.get('identifiers'))
            # TODO: make this output file more accessible and map data to subjectIDs
            np.savetxt(output.parent / (output.name + '.csv'), data, delimiter=",")

    """
        Here begins the static methods.
    """

    @staticmethod
    def _square(data):
        """

        :param data: pdist output
        :return: a square matrix with 1s in the diagonal
        """
        data = squareform(data)
        ind = np.diag_indices(data.shape[0])
        data[ind] = 1
        return data

    @staticmethod
    def split(array, nrows, ncols):
        """Split a matrix into sub-matrices."""

        r, h = array.shape
        return (array.reshape(h // nrows, nrows, -1, ncols)
                .swapaxes(1, 2)
                .reshape(-1, nrows, ncols))
