# Publishing this source release

The local Git repository is prepared independently of GitHub authentication. Publication is complete only after the remote repository exists and its main-branch commit matches the reviewed local commit. `provenance/validation.json` records local checks; it is not evidence of remote publication.

Use the official GitHub CLI for interactive login. Never paste access tokens into issues, chat messages or repository files.

```bash
gh auth login --hostname github.com --git-protocol https --web --scopes repo,workflow
python scripts/publish_github.py --repo ytwang21/TabulaStromal --visibility public
```

The publisher checks the clean committed tree, verifies the authenticated owner, creates an empty repository if absent, uploads main, and verifies the remote commit. It stops if a nonempty remote exists, visibility differs, or origin points elsewhere. It never force-pushes. The `workflow` OAuth scope is needed to upload the included GitHub Actions workflow.

This source release does not automatically publish original study data, model weights, a manuscript, a DOI release or a paper. Those have separate resource and provenance requirements.
