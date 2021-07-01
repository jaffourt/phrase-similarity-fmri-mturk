from pathlib import Path
import nibabel as nib
import numpy as np
import pandas as pd


class Data:
    """

    Required inputs:
        - expt - the exact string used for the experiment directory
        - map_type - 'con' or 'spmT', i.e. the beta weight contrasts or t-values

    """

    def __init__(self, **kwargs):
        self.map_type = kwargs.get('map_type', None)
        # the name of the directory containing neuroimaging data + masks
        self.expt = kwargs.get('expt', None)
        if self.expt:
            self.path = Path(self.expt)
        self._load_mask(kwargs.get('mask', None))
        self.condition_names = None
        self.dataframe = []
        self.subjectIDs = []
        self.__participant_data = self.__load_nii__()

    def _load_mask(self, roi):
        if roi is not None:
            self.mask = nib.load(self.path.joinpath(roi)).get_data() > 0

    def __load_nii__(self):
        files = [file for file in self.__parse_dirs__()]
        if files:
            self.condition_names = list(files[0].maps.keys())
        return files

    def __parse_dirs__(self):
        if self.expt is None:
            raise AttributeError("No value provided for experiment name")
        if self.map_type is None:
            raise AttributeError("No value provided for map type")
        for f in self.path.glob('*/'):
            if f.is_dir():
                yield Participant(subjectID=str(f).split('/')[1], files=f.glob('**/' + self.map_type + '/*nii'))

    @classmethod
    def as_numpy(cls, **kwargs):
        data = cls(expt=kwargs.get('expt', None), map_type=kwargs.get('map_type', None), mask=kwargs.get('mask', None))
        rois = data.mask > 0 if data.mask is not None else None
        for p in data.__participant_data:
            data.subjectIDs.append(p.subjectID)
            if rois is not None:
                data.dataframe.append([p.maps[condition].get_data()[rois] for condition in p.maps])
            else:
                data.dataframe.append([p.maps[condition].get_data().reshape(-1) for condition in p.maps])
        data.dataframe = np.array(data.dataframe)
        return data

    @classmethod
    def as_pandas(cls, **kwargs):
        data = cls(expt=kwargs.get('expt', None), map_type=kwargs.get('map_type', None), mask=kwargs.get('mask', None))
        rois = data.mask > 0 if data.mask is not None else None
        data.dataframe = {}
        for p in data.__participant_data:
            if rois is not None:
                data.dataframe[p.subjectID] = pd.DataFrame(
                    np.transpose([p.maps[condition].get_data()[rois] for condition in p.maps]),
                    columns=data.condition_names)
            else:
                data.dataframe[p.subjectID] = pd.DataFrame(
                    np.transpose([p.maps[condition].get_data().reshape(-1) for condition in p.maps]),
                    columns=data.condition_names)
        return data


class Participant:
    def __init__(self, **kwargs):
        self.subjectID = kwargs.get('subjectID', None)
        self.maps = {
            'MAN_BASE': [],
            'GLASS_BASE': [],
            'PEN_BASE': [],
            'CAT_BASE': [],

            'MAN_MEANINGPRES': [],
            'GLASS_MEANINGPRES': [],
            'PEN_MEANINGPRES': [],
            'CAT_MEANINGPRES': [],

            'MAN_PREPCHANGE': [],
            'GLASS_PREPCHANGE': [],
            'PEN_PREPCHANGE': [],
            'CAT_PREPCHANGE': [],

            'MAN_AJCHANGE': [],
            'GLASS_AJCHANGE': [],
            'PEN_AJCHANGE': [],
            'CAT_AJCHANGE': [],

            'MAN_NCHANGE': [],
            'GLASS_NCHANGE': [],
            'PEN_NCHANGE': [],
            'CAT_NCHANGE': []
        }

        for f in kwargs.get('files', None):
            self.maps[f.stem] = nib.load(f)
