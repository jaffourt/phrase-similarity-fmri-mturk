#brute force is best force
import warnings
import os
import sys
import pandas as pd
import h5py
from shutil import copyfile

def copy_img(*args):
    out=os.path.join(os.getcwd(),args[0],args[1])
    for img in args[3:]:
        s_out=os.path.join(out,img.split('/')[-1].split('_')[0])
        if not os.path.exists(s_out):
            os.makedirs(s_out)
        fout=os.path.join(s_out,'%s.nii'%args[2])
        print(fout)
        sys.exit()
        copyfile(img,fout)

def load_spio(*args):
    import scipy.io as spio
    f=spio.loadmat(args[0],squeeze_me=True,struct_as_record=False)
    return f

assert len(sys.argv)>1, "python3 cp_activation_maps expt.csv"

#pseudonecessary hardcoded stuff - robustness is not the goal!
m='/mindhive/evlab/u/Shared/SUBJECTS'
df=pd.read_csv(sys.argv[1])
exp=sys.argv[1].split('.')[0]
df['Contrasts']=df.drop('Subject',axis=1).values.tolist()
dat=df[['Subject','Contrasts']]
#data now ready for iteration

#we have to look at the spm.mat file for each index.
for ind, row in dat.iterrows():
    d_dir=os.path.join(m,row['Subject'],'firstlevel_'+exp)
    spm=os.path.join(d_dir,'SPM.mat')
    if os.path.exists(spm):
        copy_d={}
        try:
            f=h5py.File(spm,'r')
            for ind,ref in enumerate(f['SPM']['xCon']['name'][()]):
                contrast=''.join(chr(i[0]) for i in f[ref[0]])
                if contrast in row['Contrasts']:
                    copy_d[contrast]=ind+1
        except:
            warnings.warn('SPM.mat not in HDF5 format. Trying scipy.io load.')
            try: 
                f=load_spio(spm)
                for ind, contrast in enumerate(f['SPM'].xCon):
                    if contrast.name in row['Contrasts']:
                        copy_d[contrast.name]=ind+1
            except:
                warnings.warn("Multiple loads attempted, SPM.mat file may be corrupted. Skipping %s." % row['Subject'])
                continue
        if not len(copy_d)==len(row['Contrasts']):
            missing=','.join(list(set(row['Contrasts']).difference(set(copy_d.keys()))))
            warnings.warn("Contrasts not found in SPM.mat: %s" % missing)
        #but continue to copy, either way
        for key in copy_d:
            #this script is built around PL-2017, so just check .nii
            spmT=os.path.join(d_dir,"spmT_%04d.nii"%copy_d[key])
            con=os.path.join(d_dir,"con_%04d.nii"%copy_d[key])
            if not os.path.exists(spmT) or not os.path.exists(con):
                warnings.warn('An image was not found in %s, skipping.' % d_dir)
            else:
                copy_img(exp,row['Subject'],key,spmT,con)
    else:
        warnings.warn('SPM.mat file not found in %s'%d_dir)
