"""Detect mixed 70/100 scoring, pairing errors and misleading R2 comparisons."""
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from audit_external_metric_scope import independent_metrics, verified_external_metrics, METRICS
from tabulastromal.metrics import evaluate


class ExternalMetricScopeTests(unittest.TestCase):
    def fixture(self):
        y = np.array([20., 35., 45., 60., 80.])
        p = np.array([22., 34., 43., 62., 75.])
        split = pd.DataFrame(dict(sample_id=list('abcde'), Age=y, set=['adaptation']*2+['holdout']*3))
        predictions = split[['sample_id', 'Age']].assign(prediction=p, model='ElasticNet', K=8000, fold=1,
            source='external_all_post_adaptation')
        score = independent_metrics(y, p)
        metrics = pd.DataFrame([dict(model='ElasticNet', K=8000, fold=1, n_evaluation=5,
            source='external_all_post_adaptation', **{k:score[k] for k in ['n', *METRICS]})])
        return metrics, predictions, split

    def test_direct_formulas_match_existing_scorer(self):
        rng = np.random.default_rng(46)
        y = rng.normal(50, 18, 93)
        p = .8*y + rng.normal(10, 8, 93)
        actual = independent_metrics(y, p)
        reference = evaluate(y, p)
        for key in METRICS:
            self.assertAlmostEqual(actual[key], reference[key], places=12)
        self.assertAlmostEqual(actual['R2'], 1-actual['RMSE']**2/actual['Age_SD']**2)

    def test_all_metrics_have_identical_counts_and_ids(self):
        m, p, split = self.fixture()
        _, evidence = verified_external_metrics(m, p.sample(frac=1, random_state=6), split)
        self.assertTrue(evidence.same_pairs_for_all_metrics.all())
        self.assertEqual(evidence[['n_for_R2','n_for_MAE','n_for_MAD','n_for_PCC']].iloc[0].tolist(), [5]*4)

    def test_mixed_70_r2_and_100_mae_is_rejected(self):
        m, p, split = self.fixture()
        m['R2'] = independent_metrics(p.Age.iloc[2:], p.prediction.iloc[2:])['R2']
        with self.assertRaises(AssertionError):
            verified_external_metrics(m, p, split)

    def test_missing_adaptation_samples_is_rejected(self):
        m, p, split = self.fixture()
        with self.assertRaises(ValueError):
            verified_external_metrics(m, p.iloc[2:], split)

    def test_duplicate_or_wrong_age_pair_is_rejected(self):
        m, p, split = self.fixture()
        with self.assertRaises(ValueError):
            verified_external_metrics(m, pd.concat([p, p.iloc[:1]]), split)
        p.loc[0, 'Age'] += 1
        with self.assertRaises(AssertionError):
            verified_external_metrics(m, p, split)

    def test_nonfinite_predictions_are_not_silently_filtered(self):
        m, p, split = self.fixture()
        p.loc[0, 'prediction'] = np.nan
        with self.assertRaises(ValueError):
            verified_external_metrics(m, p, split)

    def test_lower_absolute_error_can_have_lower_r2(self):
        valid = independent_metrics([0., 25., 50., 75., 100.], [8., 33., 58., 83., 108.])
        external = independent_metrics([45., 48., 50., 52., 55.], [49., 52., 54., 56., 59.])
        self.assertLess(external['MAE'], valid['MAE'])
        self.assertLess(external['R2'], valid['R2'])
        self.assertLess(external['RMSE'], valid['RMSE'])


if __name__ == '__main__':
    unittest.main()
