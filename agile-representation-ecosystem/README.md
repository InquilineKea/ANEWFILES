# Agile Representation Ecosystem: runnable toy prototype

A small, dependency-light simulator for a **need-driven, modular, burst-activated scientific agent ecosystem**. It is a proof of mechanics, not a validated scientific discovery system.

## Run

```bash
pip install numpy
python prototype.py
```

## Mechanisms

- Four specialists preserve separate representations (causal, spectral, experimental, structural).
- A shifting need vector selects the most useful perspective.
- Five schedulers: fixed, round-robin, random, deterministic contextual, and stochastic contextual with phase rotation, novelty bonus and saturation penalty.
- Each activation produces a provenance-bearing observation artifact, with a toy verification flag.
- `coverage` measures weighted fraction of currently needed representations that have ever been activated. It is **not** discovery quality.

## What is deliberately not implemented

No LLM calls, no genuine external experiments, no true representational invention, no structured pseudospectral computations, no multi-agent deliberation, and no persistent artifact registry. The toy verifier compares against synthetic ground truth, so its `verified` flag should never be presented as independent scientific validation.

## Next engineering milestones

1. Replace synthetic specialists with independent tool-backed agents and typed artifact schemas.
2. Implement a pressure-ranked need index with unresolved questions, evidence requirements, and marginal-value estimates.
3. Add artifact deduplication, content hashes, dependency lineage and independent replication.
4. Add a task where every initial ontology fails, then allow one agent to propose a new descriptor.
5. Benchmark against simpler contextual gating with matched tokens, time and verification budgets.

## Research spine

[Controlled Access to Possibility: Pseudospectra, Dimensionality, Aging & Intelligence](https://app.notion.com/p/3fc633681625813790bfc2208e5fea70)

The motivating hypothesis is that productive scientific agent ecosystems preserve independently useful representations, recruit them selectively, coordinate in bursts, and turn collaboration into verifiable reusable artifacts. Resonance, modularity, rotating dominance, and cognitive control are design inspirations, not demonstrated biological equivalences.
