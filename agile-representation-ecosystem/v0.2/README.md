# Agile Representation Ecosystem v0.2

A **runnable, artifact-first scientific agent experiment**. No LLM credentials required.

## Run

```bash
pip install numpy
python ecosystem.py --mode need --out artifacts.json
python ecosystem.py --mode fixed --out fixed_artifacts.json
python -m unittest -v
```

## Mechanism

Three independent scientific specialists fit competing hypotheses: linear, polynomial and harmonic. A verifier scores predictions on a held-out dataset. If all fail a predefined threshold, a need-driven coordinator recruits a fourth *hybrid* representation combining their feature spaces. Proposals, verification reports and decisions become SHA-256 content-addressed artifacts with parent lineage.

## Measured on seed 9

- Fixed ontology (best polynomial): held-out MSE 0.288565
- Need-driven expanded ontology (hybrid): held-out MSE 0.005781
- 3 tests passed locally.

**Caveat:** This benchmark deliberately contains a signal requiring both polynomial and harmonic features, and the hybrid model is supplied in the feature library. The experiment demonstrates the mechanics of recruitment and verification, **not** autonomous scientific creativity or general superiority of multi-agent systems. The verifier has direct access to held-out labels; a production scientific system would require stronger independent replication and data leakage controls.

## Roadmap

1. Replace these feature-fitting specialists with independent LLM/tool agents.
2. Build a pressure-ranked need index and typed artifact schema.
3. Introduce genuinely new feature synthesis, independent experiments and replication.
4. Evaluate against strong single-agent and context-only baselines with matched compute.
5. Test sparse burst recruitment, representational agility, and novelty under distribution shift.

[Research spine in Notion](https://app.notion.com/p/3fc633681625813790bfc2208e5fea70)
