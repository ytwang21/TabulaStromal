"""MAD intentionally preserves the audited old benchmark definition."""
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
METRICS = ['R2', 'MAE', 'MAD', 'Age_MAD', 'MedAE', 'RMSE', 'PCC']


def evaluate(y, pred):
    y, p = np.asarray(y, float), np.asarray(pred, float)
    if y.shape != p.shape: raise ValueError('Prediction/target shape mismatch')
    if not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError('Nonfinite ages/predictions; no posthoc prediction filtering is allowed')
    n = len(y); varying = n > 1 and np.ptp(y) > 0
    err = np.abs(y-p)
    age_mad = float(np.mean(np.abs(y-y.mean()))) if n else np.nan
    return dict(n=n, R2=float(r2_score(y,p)) if varying else np.nan,
                MAE=float(err.mean()) if n else np.nan, MAD=age_mad, Age_MAD=age_mad,
                MedAE=float(np.median(err)) if n else np.nan,
                RMSE=float(np.sqrt(np.mean(err**2))) if n else np.nan,
                PCC=float(pearsonr(y,p)[0]) if varying and np.ptp(p)>0 else np.nan,
                metric_note=('insufficient_n_or_constant_age' if not varying else
                             'constant_prediction' if np.ptp(p)==0 else ''))
