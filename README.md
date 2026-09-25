# LINCS Gene Query Tool

An interactive **Streamlit** app for exploring the **LINCS L1000** reverse-search API: given a gene, it finds the chemical perturbagens (drugs/compounds) that most strongly **up- or down-regulate** it, then visualizes and exports the results.

<p align="left">
  <img src="https://img.shields.io/badge/Python-3776AB?style=flat-square&logo=python&logoColor=white" alt="Python" />
  <img src="https://img.shields.io/badge/Streamlit-FF4B4B?style=flat-square&logo=streamlit&logoColor=white" alt="Streamlit" />
  <img src="https://img.shields.io/badge/pandas-150458?style=flat-square&logo=pandas&logoColor=white" alt="pandas" />
  <img src="https://img.shields.io/badge/Plotly-3F4F75?style=flat-square&logo=plotly&logoColor=white" alt="Plotly" />
</p>

## What it does

- Queries the **LINCS reverse-search** API for chemical perturbagens by gene and direction (up/down).
- Ranks results by **CD Coefficient** (effect strength).
- Presents interactive **Plotly** charts and sortable tables.
- Exports query results to CSV for downstream analysis.

## Use case

Useful in **drug-discovery and bioinformatics** workflows for connectivity-map style analysis — finding compounds that reverse or reinforce a gene-expression signature.

## Tech stack

Python · Streamlit · pandas · NumPy · Plotly · Requests

## Getting started

```bash
git clone https://github.com/omaaarsh/LINCS-Gene-Query-Tool.git
cd LINCS-Gene-Query-Tool
pip install streamlit requests pandas numpy plotly
streamlit run app.py
```

Then open the local URL Streamlit prints, enter a gene symbol, and explore.

## Author

**Omar Sherif Elghamry** — [LinkedIn](https://www.linkedin.com/in/omar-elghamry-3a7256248/) · [GitHub](https://github.com/omaaarsh)
