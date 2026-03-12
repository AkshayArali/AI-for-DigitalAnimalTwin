What did you accomplish towards your research over the past week?

Finalized the Research Plan document, including the hypothesis, annotated bibliography (5 sources), Confirmed on two datasets


What (if anything) is blocking your progress?

Dataset selection is still confirmed as of now Tox21 — 12,000+ compounds screened against 12 toxicity assays; has both in vitro and some in vivo linkable endpoints. Freely available from the EPA/NIEHS. Likely your best starting point.
ToxCast — EPA's broader screening program with 1,000+ assays across thousands of compounds. Pairs well with Tox21 and gives you rich in vitro readouts for the multi-modal feature set your plan emphasizes.
ChemIDplus / DSSTox — EPA's curated chemical database; useful for linking compounds to rodent LD50 values, which gives you a clean binary toxicity label. This was the main blocker before any coding can begin.


What is your plan for the next week?

Identify and download the primary dataset; log version and access date in Research Notebook.
Set up Git repo with folder structure (e.g., /data, /notebooks, /src).
Begin data cleaning and compound-level train/val/test splitting.
Start RDKit fingerprint pipeline for the chemical-only baseline.
