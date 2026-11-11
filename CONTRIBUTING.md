# Contributing

Report a concrete problem with the command, configuration, package versions and a small synthetic reproducer. Do not attach patient metadata, credentials, private expression matrices, or model bundles containing sample identifiers.

For code changes, run `python scripts/validate_repository.py` and the affected module's tests in [docs/testing.md](docs/testing.md). Keep original inputs read-only and write outputs to a new local run directory. Do not relax checkpoint integrity or group-isolation checks to make a run pass.

Describe any change to a feature universe, split rule, model-selection criterion, expression layer, age unit, gene-set version, or statistical correction. New results need a source manifest, explicit evaluation scope, and a clear distinction between synthetic verification and study evidence.

Please discuss large refactors and resource redistribution with the maintainer first. A project-wide open-source license is pending; see NOTICE.
