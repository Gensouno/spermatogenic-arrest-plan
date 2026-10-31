# Research and manuscript working agreement

## Start and finish every project task

Read `research/progress.md`, `research/paper-todo.md`, `paper/metadata.yaml`, and the relevant analysis code/results when those local research files exist. Complete the user's project task and appropriate verification first. Decide whether it changes methods, experimental design, observations, interpretations, limitations, or reproducibility. For a substantive contribution, read and apply [.agents/skills/paper-sync/SKILL.md](.agents/skills/paper-sync/SKILL.md); otherwise leave manuscript prose alone. Maintain necessary research state in either case. Run `bash envs/paper_python.sh scripts/check_paper.py` and a Git synchronization readiness check before reporting. An infrastructure-only check is not a manuscript check.

Maintain ONE English research manuscript, incrementally, in `paper/chapters/`. Markdown is the daily prose source; Word is an export. Use formal, clear academic English for the title, abstract, keywords, prose, captions and academic declarations. Chinese communication and management notes are fine. Do not turn development logs into manuscript paragraphs. Export Word only when requested or at an agreed milestone, apart from an ignored initialization smoke test.

## Scientific boundaries

At task completion, run `bash envs/python.sh scripts/check_all_snapshots.py` to check every available frozen branch. The legacy `check_project.py` default checks a historical subset; it is not a complete current snapshot inventory. This integrity check does not replace numerical reproduction or scientific review.

Preserve original data, earlier snapshots and their bound sources. Never silently rewrite frozen input hashes to accept changed files. Recheck or explicitly downgrade affected evidence when code/config/data change. Distinguish implemented, executed, numerically reproduced, file-integrity checked, and biologically validated work. `verified` in the registry requires an explicit verification scope; it never means a biological mechanism has been established.

Preserve zero formal discoveries, multiplicity families, missing denominators, non-independent ages, original controls, count/TPM distinctions, source-ID ambiguities and reference-transfer failures. A stable descriptive list is not FDR-controlled; fitted RNA weights are not validated cell fractions. Do not assert novelty or priority without an adequate literature comparison. Record real citations and their publication/verification status; unresolved references cannot support confirmed literature claims.

Important manuscript updates require a Codex semantic review against actual code and evidence, recorded under `research/reviews/` with reviewed-file hashes. Structural checks cannot certify scientific truth, human review, animal independence or submission readiness. Keep research questions, contributions, baselines, metrics and relevant sensitivities explicit; do not add irrelevant benchmark/ablation sections mechanically.

## Privacy, Git and recurring checks

The existing `Gensouno/spermatogenic-arrest-plan` remote was public when checked on 2026-10-26. Do not change visibility or create another remote without authorization. Earlier authorization to push code does not authorize public release of an unpublished manuscript. `paper/` and `research/` are ignored by default, along with exports, tools and private configs. Keep their release policy explicit; use an authorized private destination before enabling academic-material synchronization. Do not force-add these paths to a public commit. Local snapshots/checkpoints are not Git commits.

Prepare meaningful Conventional Commit suggestions after verified milestones, run appropriate code tests and paper checks, inspect staged/proposed files for private material and secrets, and record code/data/config/environment hashes. Execute commits or pushes only within the user's actual authorized operation and repository scope; a readiness report is not permission. No force push, submission, preprint upload, external manuscript delivery or author contact without explicit authorization. Do not discard existing uncommitted work.

Weekly review is a real read-only scheduled check when configured, not this document or an end-of-task status check. It may inspect local changes, run checks and prepare a summary; it must not commit, push, upload private material or change repository settings. Do not claim a schedule is active until the scheduler confirms it. CI in the public repository checks infrastructure without requiring or publishing the private manuscript.

## Scientific skill selection

The user prefers effective use of relevant modules from [K-Dense Scientific Agent Skills](https://github.com/K-Dense-AI/scientific-agent-skills). Check the available skill catalog for each scientific task and use a matching installed module when it contributes a concrete workflow, source or validation procedure. Load only its `SKILL.md` and the references needed for the current step; do not load or install the full collection by default. `database-lookup` is the established route for reproducible public-database retrieval. Check availability before using other modules, and inspect their instructions and dependencies before any task-justified selective installation. A skill's presence does not establish that its software dependencies, database accounts or external services are installed or connected. Keep the existing evidence, environment, proxy and private-release rules in force; skill guidance is not scientific validation.

## Environment and continuity

Use `bash envs/python.sh` / `bash envs/rscript.sh` for existing analyses and `bash envs/paper_python.sh` for paper tools. Preserve the shared `.venv` / `.Rlib` symlink targets. Network downloads must inherit the ignored local proxy configuration from the first attempt, including uppercase/lowercase HTTP(S) variables and loopback exclusions; report failures without a direct fallback. Never publish machine-specific configuration or the parent workspace note.

Use one manuscript writer at a time. No subagents unless the user or an applicable instruction explicitly requests delegation. At task end briefly report project outcome, English-paper changes or no contribution, evidence and remaining uncertainty, checks, Git preparation/status and next work. User instructions and established conversation constraints take precedence over generic templates.
