import pandas as pd
import pytest
from scripts.common import strict_map


def test_exact_join_ignores_row_order():
    a=pd.DataFrame({'sample':['D1','D2']},index=['C1','C2'])
    s=pd.DataFrame({'cellid':['C2','C1'],'sample':['D2','D1'],'score':[20,10]})
    joined,qc=strict_map(a,s)
    assert qc['passed'] and joined.set_index('cell_key').loc['C1','score']==10


def test_duplicates_fail_and_metadata_mismatch_blocks():
    a=pd.DataFrame({'sample':['D1','D2']},index=['C1','C2'])
    with pytest.raises(ValueError,match='not one-to-one'):
        strict_map(a,pd.DataFrame({'cellid':['C1','C1'],'sample':['D1','D1']}))
    joined,qc=strict_map(a,pd.DataFrame({'cellid':['C1'],'sample':['D9']}))
    assert joined is None and qc['metadata_conflicts']==1


def test_repeated_barcode_requires_real_composite_key():
    a=pd.DataFrame({'sample':['D1','D2']},index=['BC','BC'])
    s=pd.DataFrame({'cellid':['BC','BC'],'sample':['D2','D1'],'score':[2,1]})
    joined,qc=strict_map(a,s)
    assert qc['passed'] and qc['keys']=='cell_key|sample'


def test_suffixes_never_stripped_and_low_mapping_blocked():
    a=pd.DataFrame(index=['ABC-1'])
    s=pd.DataFrame({'cellid':['ABC.1']})
    joined,qc=strict_map(a,s)
    assert joined is None and qc['match_rate']==0
