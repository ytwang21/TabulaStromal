import pandas as pd
from scripts.common import score_columns, fdr
from scripts.step06_gene_pathway_analysis import ora


def test_recommendations_and_gmm_are_not_senescence_intensity():
    d=pd.DataFrame({'hUSI':[.2,.8],'SenCID_rec_SID1':[5.,-1.],'husi_GMM':[1,10]})
    assert score_columns(d)==['hUSI']


def test_fdr_accounts_for_all_tested_terms_including_zero_hits():
    terms={'hit':set('ABCDE'),'nohit':set('FGHIJ')}
    out=ora(set('ABC'),set('ABCDEFGHIJKLMNO'),terms)
    d=fdr(pd.DataFrame(out),[])
    assert len(d)==2 and d.loc[d.term.eq('nohit'),'p'].iloc[0]==1
    assert d.q.ge(d.p).all()
