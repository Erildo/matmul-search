# Run logs

Raw stdout of each run (JSON lines, one per instance or step), kept so every outcome can be traced, one folder per attempt.
The summary of each run is also in `results/runs.jsonl`, and the interpretation is in `attempts/<id>.md`.
Large artefacts (CNFs, DRAT proofs, scheme pools) stay in git-ignored `scratch/`. They can all be regenerated from the code.

Not kept: the raw logs of 003, 005 and 006. 003's are not in `scratch/`, and 005 and 006 ran in a worktree that has since been removed.
Their per-instance outcomes are in `results/runs.jsonl`.
