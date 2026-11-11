"""Thirteen fresh estimators, train-only preprocessing, balanced cell-bag IDEAL variants."""
from pathlib import Path
import copy, json, logging, warnings
import joblib
import numpy as np
import scipy.sparse as sp
from scipy.stats import pearsonr
from sklearn.base import BaseEstimator,RegressorMixin
from sklearn.preprocessing import StandardScaler,QuantileTransformer
from sklearn.feature_selection import mutual_info_regression
from sklearn.linear_model import ElasticNet,Ridge,Lasso
from sklearn.ensemble import ExtraTreesRegressor,RandomForestRegressor,GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.decomposition import PCA
from sklearn.svm import SVR
from sklearn.exceptions import ConvergenceWarning
from .storage import dump,write_json,atomic_path
LOG=logging.getLogger(__name__)
NAMES=['XGBoost','ExtraTrees','ElasticNet','MLP','Ridge','Lasso','RandomForest','GradientBoosting',
       'Supervised_PCA_Ridge','SVR','DeepSets','Pseudobulk','Ensemble']


class Preprocessor:
    """No target/validation information enters transform. fit_ids is an audit trail."""
    def __init__(self, indices, names, scaler, fit_ids, K):
        self.indices=np.asarray(indices,int);self.names=list(names);self.scaler=scaler
        self.fit_ids=list(fit_ids);self.requested_K=int(K)

    @classmethod
    def fit(cls,X,y,names,ids,K,cfg,rank=None):
        rank=rank_features(X,y,cfg) if rank is None else rank
        selected=np.sort(rank[:min(K,len(rank))])
        if not len(selected):raise ValueError('No features pass train-only variance filter')
        dense=X[:,selected].toarray().astype(np.float32)
        scaler=StandardScaler().fit(dense)
        return cls(selected,np.asarray(names)[selected],scaler,ids,K)

    def transform(self,X):return self.scaler.transform(X[:,self.indices].toarray().astype(np.float32))


def rank_features(X,y,cfg):
    mean=np.asarray(X.mean(0)).ravel();v=np.asarray(X.power(2).mean(0)).ravel()-mean**2
    eligible=np.flatnonzero(v>cfg['features']['variance_threshold'])
    if cfg['features']['selection']=='variance':scores=v[eligible]
    else:
        scores=np.empty(len(eligible));batch=256
        for start in range(0,len(eligible),batch):
            ix=eligible[start:start+batch];dense=X[:,ix].toarray()
            scores[start:start+len(ix)]=mutual_info_regression(dense,y,random_state=cfg['seed'],n_neighbors=min(3,len(y)-1))
            if start%(batch*100)==0:LOG.info('MI ranking: %d/%d eligible features',start,len(eligible))
    return eligible[np.lexsort((eligible,-scores))]


class PCARidge(BaseEstimator,RegressorMixin):
    """Preserves legacy named model: supervised topK, then PCA and Ridge."""
    def __init__(self,alpha=10.,seed=42):self.alpha=alpha;self.seed=seed
    def fit(self,X,y):
        self.pca=PCA(n_components=min(50,len(X)-1,X.shape[1]),random_state=self.seed)
        self.model=Ridge(alpha=self.alpha).fit(self.pca.fit_transform(X),y);return self
    def predict(self,X):return self.model.predict(self.pca.transform(X))


class BulkRegressor:
    """Composition-free IDEAL pseudobulk: quantile transform and four-regressor mean."""
    def __init__(self,cfg,seed):self.cfg=cfg;self.seed=seed
    def fit(self,X,y):
        self.qt=QuantileTransformer(n_quantiles=min(100,len(X)),output_distribution='normal',random_state=self.seed)
        z=self.qt.fit_transform(X)
        self.models=[ElasticNet(alpha=.1,l1_ratio=.5,max_iter=self.cfg['models']['max_iter']),Ridge(alpha=10),
                     Lasso(alpha=.1,max_iter=self.cfg['models']['max_iter']),
                     GradientBoostingRegressor(n_estimators=min(200,self.cfg['models']['boost_estimators']),max_depth=3,learning_rate=.05,subsample=.8,random_state=self.seed)]
        for m in self.models:fit_checked(m,z,y)
        return self
    def predict(self,X):return np.mean([m.predict(self.qt.transform(X)) for m in self.models],axis=0)


def make_estimator(name,cfg,seed,params=None):
    c=cfg['models'];p=params or {};nj=cfg['resources']['threads']
    if name=='ElasticNet':return ElasticNet(alpha=p.get('alpha',.1),l1_ratio=p.get('l1_ratio',.5),max_iter=c['max_iter'],random_state=seed)
    if name=='Ridge':return Ridge(alpha=p.get('alpha',10.))
    if name=='Lasso':return Lasso(alpha=p.get('alpha',.1),max_iter=c['max_iter'],random_state=seed)
    if name=='XGBoost':
        from xgboost import XGBRegressor
        return XGBRegressor(n_estimators=c['xgb_estimators'],max_depth=3,learning_rate=.02,subsample=.8,
                            colsample_bytree=.7,reg_lambda=10,reg_alpha=1,random_state=seed,n_jobs=nj,tree_method='hist')
    if name=='ExtraTrees':return ExtraTreesRegressor(n_estimators=c['tree_estimators'],max_features=.3,min_samples_leaf=3,random_state=seed,n_jobs=nj)
    if name=='RandomForest':return RandomForestRegressor(n_estimators=c['tree_estimators'],max_depth=8,max_features=.3,min_samples_leaf=5,random_state=seed,n_jobs=nj)
    if name=='GradientBoosting':return GradientBoostingRegressor(n_estimators=c['boost_estimators'],learning_rate=.02,max_depth=2,loss='huber',random_state=seed)
    if name=='MLP':return MLPRegressor(hidden_layer_sizes=(128,64),alpha=1,max_iter=c['mlp_max_iter'],early_stopping=False,random_state=seed)
    if name=='Supervised_PCA_Ridge':return PCARidge(seed=seed)
    if name=='SVR':return SVR(kernel='rbf',C=10,epsilon=.2)
    if name=='Pseudobulk':return BulkRegressor(cfg,seed)
    raise KeyError(name)


def fit_checked(model,X,y):
    """Record warnings; retry linear coordinate descent, fail if it stays unconverged."""
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always');model.fit(X,y)
    notes=[str(w.message) for w in ws]
    bad=any(issubclass(w.category,ConvergenceWarning) for w in ws)
    if bad and isinstance(model,(ElasticNet,Lasso)):
        model.set_params(warm_start=True)
        original_budget=model.max_iter
        for factor in [5,25]:
            model.set_params(max_iter=int(original_budget*factor))
            with warnings.catch_warnings(record=True) as retry:
                warnings.simplefilter('always');model.fit(X,y)
            notes += [str(w.message) for w in retry]
            if not any(issubclass(w.category,ConvergenceWarning) for w in retry):bad=False;break
        if bad:raise RuntimeError('Selected linear model remains unconverged after retries: '+str(notes[-1:]))
    model.fit_warnings_=notes
    return model


def bags_for(ds,pre,rows):
    """Sparse selected cell features, scaled only within active subtype blocks.

    Centering zero-expanded vectors would turn absent-subtype genes into nonzero values;
    sparse cell inputs therefore use the training sample-level scale without centering.
    The ordinary sample-level regressors use the complete centered transform.
    """
    all_cells=ds['cells'][:,pre.indices].multiply(1/np.maximum(pre.scaler.scale_,1e-12)).tocsr()
    ix=ds['bag_index'];out=[]
    for row in rows:
        cellrows=np.flatnonzero(ix['sample']==row)
        out.append((all_cells[cellrows],ix['celltype'][cellrows]))
    return out


def torch_sparse(x,device):
    import torch
    co=x.tocoo();indices=np.vstack([co.row,co.col])
    return torch.sparse_coo_tensor(torch.as_tensor(indices,dtype=torch.long,device=device),
              torch.as_tensor(co.data,dtype=torch.float32,device=device),co.shape).coalesce()


def build_network(K,hidden):
    import torch
    from torch import nn
    class BalancedAttention(nn.Module):
        def __init__(self):
            super().__init__();self.weight=nn.Parameter(torch.empty(K,hidden));nn.init.xavier_uniform_(self.weight)
            self.bias=nn.Parameter(torch.zeros(hidden));self.encoder=nn.Sequential(nn.ReLU(),nn.LayerNorm(hidden),nn.Linear(hidden,hidden),nn.ReLU())
            self.attention=nn.Sequential(nn.Linear(hidden,hidden),nn.Tanh(),nn.Linear(hidden,1));self.head=nn.Linear(hidden,1)
        def forward(self,bags):
            predictions=[]
            for x,ct in bags:
                if x.shape[0]:
                    e=self.encoder(torch.sparse.mm(x,self.weight)+self.bias)
                    # Equal-weight mean within subtype removes the observed subtype frequency.
                    subtype=torch.stack([e[ct==v].mean(0) for v in torch.unique(ct)])
                    w=torch.softmax(self.attention(subtype).ravel(),dim=0)
                    pooled=(w[:,None]*subtype).sum(0)
                else:pooled=self.bias*0
                predictions.append(self.head(pooled).reshape(()))
            return torch.stack(predictions)
    return BalancedAttention()


class DeepRegressor:
    def __init__(self,K,cfg,seed):
        self.K=K;self.cfg=cfg;self.seed=seed;self.hidden=cfg['models']['deep_hidden']
        self.device=cfg['resources']['device']
    def fit(self,bags,y,checkpoint=None,resume=True,fine_tune=False):
        import torch
        torch.manual_seed(self.seed);torch.set_num_threads(self.cfg['resources']['threads'])
        if self.device=='cuda' and not torch.cuda.is_available():raise RuntimeError('CUDA requested but unavailable')
        if not fine_tune:
            self.net=build_network(self.K,self.hidden).to(self.device);self.ymean=float(np.mean(y));self.yscale=max(float(np.std(y)),1.)
        optimizer=torch.optim.Adam(self.net.parameters(),lr=self.cfg['models']['deep_lr']*(.1 if fine_tune else 1),weight_decay=1e-4)
        start=0;epochs=self.cfg['models']['deep_epochs']
        if checkpoint and resume and Path(checkpoint).exists():
            state=torch.load(checkpoint,map_location=self.device,weights_only=False)
            self.net.load_state_dict(state['state_dict']);optimizer.load_state_dict(state['optimizer'])
            self.ymean=state['ymean'];self.yscale=state['yscale'];start=state['epoch']+1
        ready=[(torch_sparse(x,self.device),torch.as_tensor(ct,dtype=torch.long,device=self.device)) for x,ct in bags]
        target=torch.as_tensor((np.asarray(y)-self.ymean)/self.yscale,dtype=torch.float32,device=self.device)
        batch=self.cfg['models']['deep_batch_samples'];self.history=[]
        for epoch in range(start,epochs):
            order=np.random.default_rng(self.seed+epoch).permutation(len(y));self.net.train();losses=[]
            for off in range(0,len(y),batch):
                ix=order[off:off+batch];optimizer.zero_grad();pred=self.net([ready[i] for i in ix])
                loss=((pred-target[ix])**2).mean()
                if not torch.isfinite(loss):raise RuntimeError('Nonfinite DeepSets loss')
                loss.backward();torch.nn.utils.clip_grad_norm_(self.net.parameters(),5.);optimizer.step();losses.append(float(loss.detach()))
            self.history.append(dict(epoch=epoch,loss=float(np.mean(losses))))
            if checkpoint:
                with atomic_path(checkpoint) as tmp:torch.save(dict(state_dict=self.net.state_dict(),optimizer=optimizer.state_dict(),epoch=epoch,ymean=self.ymean,yscale=self.yscale),tmp)
            if epoch%10==0:LOG.info('DeepSets epoch %d/%d train loss %.4f',epoch+1,epochs,np.mean(losses))
        return self
    def predict(self,bags):
        import torch
        self.net.eval();out=[]
        with torch.no_grad():
            for x,ct in bags:
                pred=self.net([(torch_sparse(x,self.device),torch.as_tensor(ct,dtype=torch.long,device=self.device))])
                out.append(float(pred.cpu()[0])*self.yscale+self.ymean)
        return np.asarray(out)
    def save(self,path):
        import torch
        with atomic_path(path) as tmp:torch.save(dict(K=self.K,hidden=self.hidden,seed=self.seed,ymean=self.ymean,yscale=self.yscale,state_dict=self.net.state_dict()),tmp)
    @classmethod
    def load(cls,path,cfg):
        import torch
        d=torch.load(path,map_location=cfg['resources']['device'],weights_only=False)
        obj=cls(d['K'],cfg,d['seed']);obj.hidden=d['hidden'];obj.ymean=d['ymean'];obj.yscale=d['yscale']
        obj.net=build_network(obj.K,obj.hidden).to(obj.device);obj.net.load_state_dict(d['state_dict']);return obj


class BlendRegressor:
    def __init__(self,deep,bulk):self.deep=deep;self.bulk=bulk
    def predict(self,X,bags):return .5*self.deep.predict(bags)+.5*self.bulk.predict(X)


def fit_model(name,X,y,cfg,seed,params=None,bags=None,checkpoint=None,resume=True):
    if name=='DeepSets':return DeepRegressor(X.shape[1],cfg,seed).fit(bags,y,checkpoint,resume)
    if name=='Ensemble':
        deep=DeepRegressor(X.shape[1],cfg,seed).fit(bags,y,checkpoint,resume)
        bulk=fit_checked(BulkRegressor(cfg,seed),X,y);return BlendRegressor(deep,bulk)
    return fit_checked(make_estimator(name,cfg,seed,params),X,y)


def predict_model(model,name,X,bags=None):
    if name=='DeepSets':return model.predict(bags)
    if name=='Ensemble':return model.predict(X,bags)
    return model.predict(X)


def save_bundle(path,name,model,pre,info):
    path=Path(path);path.mkdir(parents=True,exist_ok=True)
    dump(path/'preprocessor.joblib',pre);write_json(path/'model.json',dict(model=name,**info,selected_features=pre.names,
                                                                         preprocessing_fit_ids=pre.fit_ids,requested_K=pre.requested_K,effective_K=len(pre.indices)))
    outputs=[path/'preprocessor.joblib',path/'model.json']
    if name in ['DeepSets','Ensemble']:
        deep=model if name=='DeepSets' else model.deep;deep.save(path/'state_dict.pt');outputs.append(path/'state_dict.pt')
        if name=='Ensemble':dump(path/'bulk.joblib',model.bulk);outputs.append(path/'bulk.joblib')
    else:dump(path/'estimator.joblib',model);outputs.append(path/'estimator.joblib')
    return outputs


def load_bundle(path,cfg):
    path=Path(path);info=json.loads((path/'model.json').read_text());name=info['model'];pre=joblib.load(path/'preprocessor.joblib')
    if name=='DeepSets':model=DeepRegressor.load(path/'state_dict.pt',cfg)
    elif name=='Ensemble':model=BlendRegressor(DeepRegressor.load(path/'state_dict.pt',cfg),joblib.load(path/'bulk.joblib'))
    else:model=joblib.load(path/'estimator.joblib')
    return name,model,pre,info
