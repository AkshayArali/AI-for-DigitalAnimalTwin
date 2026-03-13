### Building AI “Digital Animal Twins” as an Ethical Alternative to In Vivo Toxicity Testing

**Author:** Akshay Aralikatti  
**GitHub:** [@AkshayArali](https://github.com/AkshayArali)

---

### Setup (Python environment)

- **Python:** 3.10 or 3.11 recommended (tested with 3.10+).
- **Virtual environment (optional but recommended):**

  ```bash
  python3 -m venv venv
  source venv/bin/activate   # Windows: venv\Scripts\activate
  ```

- **Install dependencies:**

  ```bash
  pip install -r requirements.txt
  ```

The virtual environment is not included in the repo; `requirements.txt` is the source of truth for dependencies.

---

### Abstract

Traditional toxicity testing relies heavily on animal (e.g., rodent) experiments, which are slow, expensive, ethically contentious, and often of limited relevance to human safety. Regulators and the scientific community are seeking New Approach Methodologies (NAMs) that can reduce or replace animal use while still providing reliable safety information. This project investigates whether AI-based “digital animal twin” models—trained on historical *in vivo* and *in vitro* toxicology data—can reliably predict rodent toxicity outcomes and thereby replace a substantial fraction of new animal tests for small molecules. The primary research question is: *To what extent can such models predict rodent toxicity well enough to prioritize which compounds truly require *in vivo* follow-up, and how many hypothetical animal tests could be avoided while maintaining acceptable sensitivity for truly toxic compounds?*
