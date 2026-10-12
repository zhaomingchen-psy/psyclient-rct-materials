# PsyCLIENT classroom RCT: analysis code, prompts, and summary data

Companion materials for *Practicing Counseling Skills with AI-Simulated Clients: A Pilot Classroom RCT* (under review, International Journal of Human-Computer Studies).

The simulated-client framework (PsyCLIENT-CP profiles, conversational trajectory modeling) is released separately at https://github.com/qiuhuachuan/PsyCLIENT (Qiu et al., 2026, arXiv:2601.07312).

## Layout

- `prompts/` — the three prompt templates used by the feedback module (Appendix A of the paper).
- `analysis/R/` — R scripts for the Study 2 analyses: HLM trajectories, two-occasion regressions, style-confidence gap, copying analyses, lagged models.
- `analysis/python/` — script that builds the analysis tables from the raw exports (needs the restricted raw data).
- `figures/` — `make_figures.py` and `style.py` reproduce Figures 3-6; vector PDFs of the figures are included.
- `data/aggregate/` — group-level summaries behind Figures 4-5 and Tables 5-8.
- `data/study1/` — Study 1 summaries (Tables 2-3) and the values behind Figure 3.

See `data/README.md` for column definitions.

## Participant-level data

Weekly self-efficacy scores, behavioral coding of sessions, practice logs, and session transcripts were collected under a consent form that limits access to the research team and allows publication of aggregate statistics only, so they are not in this repository. De-identified participant-level data (without transcripts) can be shared with qualified researchers under a data use agreement.

## Contact

- Huachuan Qiu (first author): qiuhuachuan@westlake.edu.cn
- Zhaoming Chen (classroom study and analyses): zhaoming.chen@utah.edu
- Zhenzhong Lan (corresponding author): lanzhenzhong@westlake.edu.cn

Two participant identifiers that were hard-coded in the R scripts for record linkage have been replaced with `<STUDENT_ID_xxx>` placeholders. The scripts document the analyses but will not run without the raw data.

## Figures

```
pip install -r requirements.txt
python figures/make_figures.py
```

`make_figures.py` expects the raw `.sav` and CSV inputs described at the top of the script.

## License

Code: MIT (see `LICENSE`). Data and prompts: CC BY 4.0 (see `LICENSE-DATA.md`).

## Citation

Please cite the paper once published; until then cite the PsyCLIENT framework paper (Qiu et al., 2026) and this repository.
