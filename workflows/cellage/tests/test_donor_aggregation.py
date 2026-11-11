import numpy as np
import pandas as pd
from scripts.common import aggregate_donors, adjusted_fit, attach_donors


def test_cell_numbers_do_not_change_equal_donor_summary():
    c=pd.DataFrame({'lineage':['fibro']*102,'donor_key':['D1']*100+['D2']*2,
                    'tissue':['lung']*102,'celltype':['Type_A']*102,'age':[20]*100+[80]*2,
                    'gse':['G1']*102,'hUSI':[0.]*100+[1.]*2})
    d=aggregate_donors(c,['hUSI'],min_cells=2)
    assert len(d)==2 and d.score_mean.mean()==.5
    assert len(aggregate_donors(c,['hUSI'],min_cells=20))==1


def test_donors_require_evidence():
    c=pd.DataFrame({'sample':['S1'],'gse':['G1']})
    assert attach_donors(c,{},'fibro').donor_key.isna().all()
    assert attach_donors(c,{'donor_column':'sample','donor_identity_evidence':'user confirmation'},'fibro').donor_key.iloc[0]=='S1'


def test_fully_confounded_age_is_not_reported_as_adjusted_effect():
    c=pd.DataFrame({'donor_key':[f'D{i}' for i in range(24)],'age':np.repeat([10,30,50,70],6),
                    'gse':np.repeat(['G1','G2','G3','G4'],6),'score_mean':np.arange(24)/24})
    out=adjusted_fit(c)
    assert out['status'].startswith('not_identifiable') and np.isnan(out['p'])
