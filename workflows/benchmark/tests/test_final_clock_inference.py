"""SYNTHETIC inference-interface checks; no simulated biological result files."""
from pathlib import Path
from types import SimpleNamespace
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from final_clock_common import *
from sklearn.linear_model import ElasticNet

class FinalInferenceTests(unittest.TestCase):
    def setUp(self):
        self.schema=dict(celltypes=['cell_A','cell_B'],genes=['GENE_WITH_UNDERSCORE','G2'],
            features=['cell_A_GENE_WITH_UNDERSCORE','cell_A_G2','cell_B_GENE_WITH_UNDERSCORE','cell_B_G2'])
        self.raw=np.array([[1,2,0,0],[3,4,5,6],[2,1,1,3],[4,5,3,1]],np.float32)
        scale=StandardScaler().fit(self.raw);model=ElasticNet(alpha=.1,l1_ratio=.5,max_iter=10000).fit(scale.transform(self.raw),[20,50,30,60])
        ctx=SimpleNamespace(lineage='SYNTHETIC_TEST_ONLY',K=4,source_signature='synthetic',
            cfg=dict(fields=dict(sample='sample',gse='GSE',celltype='Celltype',age='Age'),features=dict(min_cells_per_sample=5,chunk_cells=3)))
        self.clock=make_bundle(ctx,model,scale,np.arange(4),self.schema,'SYNTHETIC_TEST_ONLY')
        self.frame=pd.DataFrame(self.raw,columns=self.schema['features'],index=['GSE|s0','GSE|s1','GSE|s2','GSE|s3'])
    def test_named_columns_are_reordered(self):
        expected=predict_age(self.frame,self.clock).age_pred
        actual=predict_age(self.frame[self.frame.columns[::-1]],self.clock).age_pred
        np.testing.assert_array_equal(expected,actual)
    def test_memory_layout_does_not_change_prediction(self):
        expected=bundle_predict(self.clock,np.array(self.raw,order='C'))
        storage=np.empty((len(self.raw),self.raw.shape[1]*2),dtype=np.float32)
        storage[:,::2]=self.raw
        for values in [np.array(self.raw,order='F'),storage[:,::2]]:
            np.testing.assert_array_equal(values,self.raw)
            np.testing.assert_array_equal(bundle_predict(self.clock,values),expected)
    def test_missing_features_use_zero_before_frozen_scaler(self):
        name=self.frame.columns[1];missing=self.frame.drop(columns=name)
        explicit=self.frame.copy();explicit[name]=0
        np.testing.assert_array_equal(predict_age(missing,self.clock).age_pred,predict_age(explicit,self.clock).age_pred)
    def test_no_celltype_ratios_enter_prediction(self):
        changed=self.frame.assign(cell_A_ratio=999.)
        np.testing.assert_array_equal(predict_age(changed,self.clock).age_pred,predict_age(self.frame,self.clock).age_pred)
    def test_duplicate_and_nonfinite_inputs_rejected(self):
        duplicate=self.frame.copy();duplicate.columns=['same']*4
        with self.assertRaises(ValueError):predict_age(duplicate,self.clock)
        invalid=self.frame.copy();invalid.iloc[0,0]=np.nan
        with self.assertRaises(ValueError):predict_age(invalid,self.clock)
    def test_reload_bundle_and_unlabeled_h5ad(self):
        import anndata as ad
        with tempfile.TemporaryDirectory(prefix='final_clock_inference_',dir='/tmp') as tmp:
            tmp=Path(tmp);modelpath=tmp/'clock.joblib';dump(modelpath,self.clock)
            np.testing.assert_array_equal(predict_age(self.frame,load_final_clock(modelpath)).age_pred,predict_age(self.frame,self.clock).age_pred)
            # Six cells/sample; Celltype/Gene names contain underscores and
            # the second sample has two measured subtype blocks.
            X=np.vstack([np.tile([1,2],(6,1)),np.tile([3,4],(3,1)),np.tile([5,6],(3,1))]).astype(np.float32)
            obs=pd.DataFrame(dict(sample=['s0']*6+['s1']*6,GSE=['GSE']*12,Celltype=['cell_A']*9+['cell_B']*3),index=[f'cell{i}' for i in range(12)])
            a=ad.AnnData(sp.csr_matrix(X),obs=obs,var=pd.DataFrame(index=self.schema['genes']));p=tmp/'unlabeled.h5ad';a.write_h5ad(p)
            result=predict_age(p,self.clock)
            expected=predict_age(self.frame.iloc[:2],self.clock)
            self.assertEqual(result.sample_id.tolist(),expected.sample_id.tolist())
            np.testing.assert_allclose(result.age_pred,expected.age_pred,atol=1e-7)
            self.assertNotIn('age_true',result)

if __name__=='__main__':unittest.main()
