# Fluxes, Repertoires, Disagreement

**A review and rebuild of the "Computational Matter Zoo" ChatGPT threads, with a hackathon-sized core.**

- Sources: [thread 1](https://chatgpt.com/share/6ad18d80-0df0-83ea-8682-7c274ef045a2) (18 turns) and [thread 2](https://chatgpt.com/share/6ad197e7-cde4-83e9-8828-7f49bad97836). Thread 2 is thread 1 plus two new turns: "why decentralized?" and "how do agents communicate?"
- Every numerical claim marked **[checked]** below is reproduced by `checks.py` in this folder. All 17 checks pass, and the script only needs numpy.
- Written 2026-10-16. Confidence tags: *fairly sure*, *think so*, *speculating*.

---

## 0. The short version

1. **About 60% sensible.** The last third (QuotientLab, plus thread 2's event-triggered communication) is good and close to buildable. The first two-thirds stack metaphors and almost never say which parts are metaphor and which are mechanism.
2. **The best idea in both threads** is checkable and small. Two lattices can have the same spectrum and route signals oppositely. Two matrices can differ entry by entry and behave identically. The quantity that tells these cases apart is a **loop flux, the gauge-invariant part of a 1-form on a graph**.
3. **That same object shows up in your collaborator's attention idea** once you set it up right. Softmax logits are a 1-form on the query–key bipartite graph, softmax's row normalization is the gauge, and the 4-cycle fluxes are exactly the log-odds ratios of the attention matrix [checked]. So the two of you are working on *the same mathematical object*. That's your pitch.
4. **Three bugs worth fixing before anyone sees this:**
   - The flagship "isospectral, opposite routing" example is a mirror pair, and the isospectrality comes from time-reversal by construction [checked].
   - Thread 2's Laplacian-entropy "agent diversity" metric is inverted: six identical agents score about 5 out of a maximum of 5 [checked]. Use the kernel (Vendi score) instead.
   - The proposed ablation has no control for *how much* the agents communicate, as opposed to *when*.
5. **A lot was reinvented without being named.** That includes Spivak's functorial data migration (Buehler co-wrote the 2011 olog paper with Spivak), resolvent analysis, Kane–Lubensky topological kagome, MacNet, query-by-committee, island-model GAs and Lazer–Friedman. Naming them makes the project look *more* serious, not less novel.
6. **Deadline:** thread 1 says the ScienceClaw early-application deadline is Oct 16, 2026 (today), with the event Oct 30–Nov 1 at the MIT Media Lab. I could not confirm this from search. Check the official page.

---

## 1. Map of the threads, turn by turn

| # | You asked | It answered | Verdict |
|---|---|---|---|
| 1 | How would Ghrist, Papillon, Peyré, Levin, topology, higher-order and curvature people classify test-time compute (TTC)? | An 8-tuple 𝔗 = (K_t, F, ρ, U, M, Π, B, τ); a lens table; "universality classes of TTC" | **Mixed.** The sheaf idea is good (but discourse sheaves already exist). Curvature, persistent homology and Hodge percentages on 5–20 agents are mostly noise. The topology-sweep experiment was already run (MacNet). |
| 2 | Condensed-matter diversity: tensegrity, tensors, Weyl, magnonics… | "Modal diversity", a carrier zoo, phase regimes, "Computational Matter Zoo" | **Mostly metaphor.** ker R vs ker Rᵀ is correct. The Weyl metaphor is backwards (§3.3). The 4-regime table (mode diversity × coupling sparsity) is the useful part. |
| 3 | Andrea Liu; impurities in steel | Jamming and soft spots, allostery by pruning ~1% of bonds; "spectral doping"; QD over mode phenotypes | **Solid core.** The steel correction is right. It missed that Liu's isostaticity leads straight to Kane–Lubensky kagome (§4). |
| 4 | How would swarms help? Would Buehler care? Clune QD | Swarm ⊂ QD ⊂ physics; QDHF/QDAIF; agents invent descriptors | **Good.** Uncited precedents: ELM, AURORA, OMNI. |
| 5 | Kagome, MXenes, the OEA gauge-field paper | Geometry vs material vs coupling-matrix layers; "physical adjacency ≠ computational adjacency" | **Good framing.** Paper details not verified. |
| 6 | Agent roles for kagome | 20-row role table; gauge-invariant loop descriptors; multiple archives | **Good on the invariant-descriptor point, inflated on roles.** Most "agents" should be functions. |
| 7–8 | Copresheaves (Wang & Buehler); how many DOF? | Typed ecology, schema revision, Kan transport, residual; then admits the demos are ~2–3 DOF | **Correct once pushed.** It's Spivak's Σ-migration (§4). |
| 9 | Pseudospectra? | No. The Kirchhoff matrix is normal, so pseudospectra add nothing | **Correct.** |
| 10 | Resolvent, moiré, carrier steering | Transfer function H(ω) = C(iω−A)⁻¹B, SVD channels, effective rank, Landauer transport, moiré valley/layer | **Correct physics.** The SVD-of-resolvent part is McKeon–Sharma resolvent analysis. "Diversity of the resolvent repertoire" is too vague (fix in §5.2). |
| 11 | Why don't CPUs fit? | They do: an extreme projection onto {0,1}; a continuum from CMOS to quasiparticle computing | **Good, but misses Keyes** (gain, fan-out, restoration), the real reason exotic substrates lose. |
| 12–13 | Moiré/spintronic routing | Wire address → internal-state address; static, control-dependent and state-dependent routing; T(E, V, s) | **Nice taxonomy.** Experimental claims not verified, and most of it is cryogenic. |
| 14 | How do we use swarms? | Executable sandbox, fidelity ladder, mechanism-based mutations, humans as taste | **Good.** |
| 15 | Notion page | Created a private draft | — |
| 16 | Where does arXiv 2610.18174 (QRR) fit? Collaborator's "loops orthogonal to linear attention" | QRR's transport is pure gauge (zero holonomy); three meanings of "loop"; Hodge is linear; triangle S-matrix; QuotientLab | **Best turn.** All numbers reproduce. Misses the mirror-pair issue and the attention/bipartite reformulation. |
| 17 | Why a decentralized collective? (thread 2) | Two coupled searches (materials and representations); disagreements become experiments, not votes | **Good draft.** Rewritten in §7. |
| 18 | How do they communicate, strategically? (thread 2) | Private representations; disagreement trigger Δ_rs; Laplacian diversity; VOI router; "scientific synapses"; 3 timescales; stratified memory; schema forks; 3 coupled archives; ablation A–D | **Strong ideas with one inverted metric, one missing control, and too many archives.** Rebuilt in §6. |

---

## 2. What I checked

- **Triangle S-matrix [checked].** With H(Φ) on three sites, unit hoppings and flux Φ, at E = 0 with contact broadening γ = 2: Φ = 0 gives 20/40/40; Φ = +π/2 gives 0/100/0; Φ = −π/2 gives 0/0/100. The spectra for ±π/2 are both {−√3, 0, √3}. A random local rephasing changes the matrix by ~3 (Frobenius norm) and changes |S|² by ~1e-16.
- **QRR transport is pure gauge [checked].** U_ij = γ_i γ_j\* gives holonomy 1 around every loop.
- **arXiv 2610.18174 exists.** "Looped Reasoning Finishes Earlier Than You Think," Sumin Park and T. Konstantin Rusch, submitted 13 Oct 2026. QRR is a gauge-equivariant recurrent model with exact per-token phase symmetry, tested on Sudoku, Maze and Mini-ARC.
- **The Wang & Buehler paper exists.** "Self-Revising Discovery Systems for Science: A Categorical Framework for Agentic Artificial Intelligence," roughly June 2026, listed as arXiv 2606.01444 (ID not confirmed on arxiv.org).
- **Not checked:** the OEA artificial-gauge-field splitter numbers, the 2025–26 moiré experiments (gate-defined moiré collimators, nonlinear valley Hall, ferroelectric nonlinear Hall, trilayer hBN tessellations), XPd₅ flat bands, the SparksMatter description, and the hackathon details. ChatGPT was searching live, so they're probably roughly right, but **open each paper before citing it in an application.**

---

## 3. Corrections, worst first

### 3.1 The flagship counterexample is a mirror pair

H(−π/2) is H(+π/2) with ports 2 and 3 swapped (up to gauge) [checked]. More fundamentally, **H(−Φ) = H(Φ)\*** for every Φ [checked], so the two are *always* isospectral because complex conjugation (time reversal) preserves eigenvalues.

Restate it like this: **any spectrum-only descriptor is blind to the sign of every loop flux, by construction.** That's a stronger and more honest claim than "here's one counterexample."

Also, the ±π/2 three-mode loop is the known coupled-mode **circulator condition** (Ranzani & Aumentado 2015, *think so*). Don't present it as a discovery. Present it as a **rediscovery benchmark**: can the swarm, starting from a spectrum-only archive, find the circulator and invent "signed loop flux" as the missing descriptor without being told? That makes a good unit test for the whole system.

### 3.2 Thread 2's agent-diversity metric is upside down

The thread proposes building an agent-similarity graph W, taking the Laplacian L = D − W, and using exp(entropy of normalized Laplacian eigenvalues) as "agent spectral diversity."

For six **identical** agents (complete graph), L has eigenvalues {0, 6, 6, 6, 6, 6}, so the score is 5.0, essentially the maximum [checked]. It's worse than that: with uniform similarity s, the Laplacian eigenvalues are {0, 6s, …, 6s}, so the score is 5.0 for *every* s > 0 and undefined at s = 0. **The metric rewards homogeneity at best and ignores similarity at worst.**

Fix: take the entropy of the **similarity kernel** K, not the Laplacian. That is the **Vendi score** (Friedman & Dieng 2023). Identical agents score 1, and N orthogonal agents score N [checked]. If you want subspace-level comparisons between agents' representations, use CKA (Kornblith et al. 2019) or principal angles. Thread 2 mentions these in passing; make them the main tool.

### 3.3 The Weyl metaphor runs the wrong way

Thread 1 says rare hypotheses live "where two conceptual manifolds touch," like a Weyl crossing. But a Weyl point is where bands touch **without hybridizing**, protected by topology. "Two modes meet and make something new" describes an **avoided crossing**, where the modes do mix.

The correct and more useful version: topology lets two channels share a point in parameter space without contaminating each other. That's an argument for *protected specialists*, which is the point the thread wanted to make about decoupling anyway.

### 3.4 MAP-Elites does not scale to a "50–500+ dimensional ontology"

Grid archives break down past ~5–6 descriptor dimensions, and CVT-MAP-Elites (Vassiliades et al. 2018) only stretches that so far. Above that you need learned low-dimensional descriptors, such as AURORA (Grillotti & Cully), which already does the "descriptor inventor" job without an LLM, or QDHF. Alternatively use kNN novelty in an embedding.

The thread's own better instinct, "count independent mechanisms, not DOF," is right. Keep the archive at 2–4 descriptors and let schema revision swap them, not pile them up.

### 3.5 Graph diagnostics at agent scale are mostly noise

Ollivier–Ricci curvature, persistent homology and "72% gradient / 18% curl" on 5–20 agents won't be stable. The Hodge split also needs edge-valued data, such as pairwise preferences or agreement scores (HodgeRank), not raw messages. Keep the *idea* (decompose agreement into a global ranking plus cycles) and apply it to pairwise judgments, where HodgeRank already works.

### 3.6 Thread 2's ablation is missing the decisive control

Comparing (A) one agent, (B) everyone shares everything, (C) isolated agents and (D) event-triggered sharing cannot tell whether D wins because of the *trigger* or simply because there is *less communication*. Add:

- **(E) random sparse communication at the same rate as D.**

If D ≈ E, the trigger is decoration. That's the single most important experiment in the plan.

### 3.7 Budget size can flip the expected result

Lazer & Friedman (2007) found that densely connected networks win **early**, and sparser networks win **late**, because they keep more diversity alive. A hackathon budget is the *early* regime. So event-triggered sparsity might lose at small budgets and win at large ones. **Report performance-vs-budget curves, not one number.** A crossover would be a more interesting result than a clean win.

### 3.8 CPUs: the argument the thread skipped

Digital won not only because of noise margins but because of Keyes's (1985) criteria: **gain, fan-out, signal restoration, and isolation between input and output.** Most magnonic, skyrmion, moiré and photonic "computing" proposals fail on cascadability. Any substrate claim in the pitch should say which Keyes criterion it gives up and why that's acceptable for the target task.

### 3.9 Smaller items

- "Linear attention" is not "sparse attention." Thread 1 already caught this, so keep it.
- Excess heat phonons in a CPU are high-entropy waste, not an untapped channel.
- "Levin-style ratio of goal-space sizes" is not measurable as stated. Drop it or operationalize it.
- Sycophancy pattern: every reply opens "Yes." and escalates. Pushback only appeared when you forced it.

---

## 4. Reinvented prior work, and why naming it helps

| What the thread called it | What it already is | Why it matters for the pitch |
|---|---|---|
| Left Kan extension carrying old data into a richer schema | **Spivak's functorial data migration** (Σ_F), 2012; implemented in CQL (fast left Kan extension via the chase) and AlgebraicJulia/Catlab **ACSets** | You can *build* the copresheaf layer in a day with Catlab ACSets. Buehler co-wrote the **2011 olog paper with Spivak** on hierarchical protein materials, so this is his own lineage. |
| SVD of the resolvent to find the most amplified forcing/response | **Resolvent analysis** (McKeon & Sharma 2010); Hankel singular values and Gramians in control theory | Mature numerics, and it makes "transfer operator" language credible to engineers. |
| Kagome + zero modes + topology + Liu, as four topics | **Kane & Lubensky (2014)** topological polarization of isostatic (Maxwell) lattices, built on *deformed kagome* | Liu's jamming is isostatic (z = 2d). Isostatic lattices are Maxwell lattices. Kane–Lubensky kagome puts floppy modes on chosen edges. One literature, and the cleanest physical anchor you have. |
| Substrate that computes | Liu group: **tuning by pruning** (Goodrich, Liu & Nagel 2015), allostery (Rocks et al. 2017), **coupled learning** (Stern et al. 2021), physical learning machines (Dillavou et al. 2022) | Real precedents for "the material learns locally," not metaphors. |
| Sheaf-based agent consistency | **Discourse sheaves** (Hansen & Ghrist 2021) | Gives you private opinion spaces plus translation maps: the precise version of "agents speak different languages." |
| Hold the budget fixed and vary the agent communication topology | **MacNet** (Qian et al. 2024); **GPTSwarm** (Zhuge et al. 2024) | Already run, with a "collaborative scaling law." Cite it and differentiate: you vary *when* agents talk, not *the graph*. |
| LLM proposes mutations inside QD | **ELM** (Lehman et al. 2022) | |
| Learned "interestingness" | **OMNI** (Zhang, Lehman, Stanley & Clune 2024), QDHF, QDAIF | |
| Search over agent architectures | **ADAS** (Hu, Lu & Clune 2024) | Connects back to the test-time-compute question in thread 1. |
| Communicate where representations disagree | **Query-by-committee** (Seung, Opper & Sompolinsky 1992) | The core trigger has 30 years of theory behind it, which helps. |
| Isolated agents with sparse, periodic integration | **Island-model GAs** (migration interval and rate); **Lazer & Friedman 2007** | Known sweet spot, and the known budget-dependence in §3.7. |
| Effective rank of channel singular values | Roy & Vetterli effective rank; **Clements mesh** (2016) as a universal-unitary baseline | Gives "repertoire" a ceiling (§5.2). |

---

## 5. The rebuilt core: three precise objects

### 5.1 Fluxes mod gauge, one object for both of you

**Lattice side.** Hopping phases φ_ij form a 1-form on the lattice graph. Local rephasing (φ_ij → φ_ij + θ_i − θ_j) is gauge and changes nothing physical [checked]. Loop sums Φ_γ = Σ_γ φ_ij are gauge-invariant and change routing [checked]. Spectra cannot see the *sign* of Φ [checked].

**Attention side** (*my reconstruction of the collaborator's idea; maybe 35–45% that it's what they meant, but the math is solid*). Logits S_ij = q_i·k_j live on the edges of the complete bipartite query–key graph K_{m,n}. Then:

- Softmax over j is invariant to **row shifts** S_ij → S_ij + c_i. That is a gauge [checked].
- Key biases S_ij → S_ij + b_j are **not** gauge. They are a potential: a global "salience" per key [checked].
- What remains after removing additive row and column effects is the **double-centered matrix R**. Its space has dimension (m−1)(n−1), exactly the cycle space of K_{m,n} [checked].
- R's plaquette fluxes are S_ij − S_il − S_kj + S_kl = **log odds ratio of the attention matrix**, log(A_ij A_kl / A_il A_kj), which is independent of the softmax normalizers [checked].

So "nonlinear relations as loops orthogonal to the linear part" becomes: **the non-additive (interaction) part of attention equals the cycle-space component, the fluxes are log-odds ratios, and the SVD of R is the spectral representation.** Linear attention also has well-defined odds ratios, so the real question is which flux patterns each attention family can *reach*, not "loops vs linearity."

QRR fits as the degenerate case: its node-phase transport has zero flux everywhere [checked]. It removes gauge redundancy but carries no relational content.

**One-line pitch:** in both halves the data is a 1-form on a graph modulo gauge, and the content is the fluxes. One library handles both lattice hopping phases and attention logits.

*Anti-mode:* someone could fairly say the attention half is a relabeling of two-way ANOVA. True. The test of whether it's more than a relabeling is whether flux features (top singular modes of R) predict anything a key-salience baseline misses on a task where relational structure matters. Run that before claiming it.

*Janusian aside (speculating):* for a row-stochastic attention matrix read as a Markov chain, Schnakenberg's cycle affinities are another gauge-invariant loop quantity. They measure irreversibility rather than interaction. Interaction (odds ratio) and irreversibility (cycle affinity) are different loop quantities on different graphs. It may be worth checking which one the collaborator means by "indirect relation."

### 5.2 Repertoire dimension instead of "resolvent diversity"

Thread 1 proposes maximizing "diversity of {C R(E, λ) B}" without defining it. Define it:

> **Repertoire dimension** r(x, λ) = rank of the Jacobian ∂|S|²/∂λ at operating point λ, after quotienting out gauge directions.

- It is the local dimension of the set of routing behaviors your controls can reach.
- It has **hard ceilings**: (N−1)² for routing probabilities (doubly stochastic matrices, by Birkhoff) and N² for full unitaries.
- It has a free **baseline**: a Clements MZI mesh reaches full U(N) with N(N−1)/2 knobs.
- Examples [checked]: flux-only triangle, rank 1. Triangle with three hoppings plus a flux, rank 4, which is the Birkhoff ceiling for N = 3.

Use it as one QD descriptor ("how much can this substrate do") alongside one *type* descriptor ("what kind of thing does it do", e.g. circulating, splitting, localizing).

### 5.3 Representational disagreement as the only communication trigger

This is thread 2's best idea, tightened. Each agent r keeps its own representation and distance d_r. For candidate pairs (i, j):

  Δ_rs(i, j) = | rank_r(i, j) − rank_s(i, j) |

Use rank-normalized distances so scales don't matter, computed on a sampled set of pairs rather than all n² of them.

A large Δ means two scientific vocabularies disagree about whether two designs differ. Only the top few disagreements per round open a **temporary 2–3 agent coalition**. The coalition must return exactly one of:

- a new descriptor (as code),
- a discriminating simulation, or
- a proof of equivalence (a symmetry or gauge map).

This is query-by-committee applied to representations. In the rediscovery benchmark it fires exactly where it should: spectral distance ≈ 0 while routing distance is large, for the ±Φ pair.

---

## 6. The swarm, rebuilt

### 6.1 Agents vs code

Most of the 20 roles in thread 1 are **functions**: spectral ecologist, localization agent and transport agent all compute descriptors deterministically. Make those code. Use LLM agents only where judgment or invention is needed:

| Agent | Has its own | Emits |
|---|---|---|
| **Geometry/phase explorer** | design-space prior, mutation history | mutation operators (as code) |
| **Spectral reader** | spectrum-based representation | anomalies, disagreement reports |
| **Transport reader** | S-matrix and resolvent representation | anomalies, disagreement reports |
| **Descriptor inventor** | candidate descriptor library | new descriptor functions |
| **Quotient verifier / skeptic** | symmetry and gauge test suite | equivalence proofs, artifact kills |

Five agents. Thread 2 suggested six; merge "QD archivist" into deterministic code.

### 6.2 Communication protocol

- **Private by default.** Each agent has its own context, representation and candidate ranking. The spectral reader never sees the transport reader's interpretation, only published artifacts.
- **Published evidence only**, as typed records: `{claim, candidate ids, simulation hash, uncertainty, counterexample}`. Store them as ACSets/copresheaf instances so provenance comes free.
- **Trigger = top-k Δ_rs per round** (§5.3). Nothing else opens a channel, except one random cross-pair probe per round for exploration.
- **"Scientific synapses" are a bandit.** Thread 2's coupling weights J_ab updated by yield, plus an anti-echo penalty, are a reinvented multi-armed bandit over agent pairs. Use Thompson sampling: exploration comes built in and the ad-hoc concentration penalty goes away.
- **Drop VOI as stated.** Expected information gain isn't computable here. Use predicted reduction in Δ for that pair as the proxy, or let the bandit learn it.
- **Three timescales** (local ≪ coalition ≪ global) are an island-model migration schedule. Keep them, and sweep the migration interval as an experimental variable.
- **Schema forks** ("representational speciation") are worth keeping as a concept, but for the hackathon allow at most two coexisting descriptor sets.
- **Drop the third archive.** Thread 2's "three coupled archives" (matter, representation, swarm) is too much. Keep matter + representation and leave the swarm archive as future work.

---

## 7. Experiment design (hackathon-sized)

**Sandbox:** 3- to 12-site tight-binding graphs (triangle → small kagome patches) with complex hoppings, 2–4 ports, and exact S(E) from the formula in `checks.py`. Thousands of evaluations per second on a laptop. Use Kwant later for larger systems.

**Rediscovery benchmark (main result):** start the archive with spectrum-only descriptors. Success means the swarm (a) finds circulating designs, (b) flags the ±Φ disagreement, (c) proposes signed loop flux or an equivalent descriptor, and (d) the quotient verifier rejects gauge copies. Record the number of solver calls each step takes.

**Conditions** (matched LLM tokens and matched solver calls):

| | Condition | What it isolates |
|---|---|---|
| A | single agent with all tools | does division of labor help at all |
| B | 5 agents, share everything continuously | cost of consensus |
| C | 5 agents, isolated, merge only at end | cost of isolation |
| D | 5 agents, disagreement-triggered coalitions | the proposal |
| **E** | 5 agents, random sparse coalitions at D's rate | **trigger vs. mere sparsity** |

**Metrics**, on a *fixed* held-out evaluation set the swarm cannot redefine:

- validated distinct behaviors
- false novelty (gauge or mirror copies counted as new)
- false merging (distinct behaviors collapsed)
- solver calls per validated discovery
- Vendi score of agents' representations over time

Plot everything **against budget** (§3.7).

**Stretch:** the attention module. Extract logits from a small model on a relational task, compute R and its fluxes, and test whether flux features beat a key-salience baseline. Owned by your collaborator.

**48-hour split:**

| Window | Work |
|---|---|
| Day 0 (now) | sandbox plus `checks.py` → unit tests; descriptor library v0 |
| Hours 0–12 | agents A and D running end to end |
| Hours 12–30 | conditions B, C, E; budget sweeps |
| Hours 30–42 | rediscovery benchmark; plots |
| Hours 42–48 | write-up; attention stretch only if ahead |

---

## 8. Hypotheses about how it turns out

| | Hypothesis | Confidence | What would move it |
|---|---|---|---|
| H1 | D beats B, C and E on validated discoveries per solver call | plausible (~0.4) | D ≫ E at matched rate supports H1. D ≈ E sends mass to H2. |
| H2 | Sparsity helps but the trigger doesn't (D ≈ E > B) | plausible (~0.3) | Same test. Also check whether Δ concentrates on a few pairs regardless of content. |
| H3 | At hackathon budgets, B (dense) wins and D only catches up late | possible (~0.2) | Budget curves cross. That's still a publishable "Lazer–Friedman for LLM swarms" result. |
| H4 | The single agent A wins outright because the sandbox is too easy to need division of labor | can't rule out (~0.1) | A solves the rediscovery benchmark in few calls. Then make the sandbox harder (bigger kagome patches, noise, more ports). |
| — | something I haven't considered | ~0.1 | |

**The mode, questioned:** my default would be "D wins." That's partly because it's the story everyone in both threads wants. The evidence for it is the analogy to modal decoupling plus query-by-committee theory. The analogy is not evidence, and QBC's theory is about label efficiency, not discovery. Treat H2 and H3 as live.

---

## 9. Application text, rewritten

**Project (≈120 words).**
*QuotientLab* is a decentralized agent collective that learns which differences matter. In programmable wave and lattice systems, many distinct-looking designs are the same device in different coordinates (gauge copies), and some designs that look identical to standard descriptors behave oppositely (time-reversed partners share a spectrum but route signals in opposite directions). Agents with deliberately different representations (spectral, transport, geometric) work privately and communicate only when their representations disagree about whether two designs differ. Each disagreement must end in a new descriptor, a discriminating simulation, or a proof of equivalence. We benchmark against single-agent, always-communicating, isolated, and randomly-sparse baselines at matched compute. Our first test is whether the collective rediscovers the coupled-mode circulator and invents "signed loop flux" without being told.

**Why decentralized (≈90 words).**
The hard part is that the vocabulary is wrong at the start, not that there are too many candidates. A single optimizer only explores the descriptor space it is given. Agents holding different representations act as a committee: where they disagree is exactly where the current vocabulary is wrong. Keeping them private preserves that disagreement; triggering short coalitions on it turns it into experiments rather than votes. We test this claim directly, including against random sparse communication, so it can fail.

**The bridge (one line).** We study the same object on two substrates: a 1-form on a graph modulo gauge, whose loop fluxes carry the content, in lattice hopping phases and in attention logits.

---

## 10. References

*From memory unless marked verified. Check before citing.*

- Park & Rusch, *Looped Reasoning Finishes Earlier Than You Think*, arXiv:2610.18174 (2026). **verified**
- Wang & Buehler, *Self-Revising Discovery Systems for Science: A Categorical Framework for Agentic AI*, ~June 2026, arXiv 2606.01444 (ID unconfirmed). **verified exists**
- Spivak, *Functorial data migration*, Information and Computation (2012).
- Spivak, Giesa, Wood & Buehler, *Category theoretic analysis of hierarchical protein materials and social networks*, PLoS ONE (2011).
- Patterson, Lynch & Fairbanks, *Categorical data structures for technical computing*, Compositionality (2022). (ACSets, Catlab.jl)
- *Fast Left Kan Extensions Using the Chase*, arXiv:2205.02425.
- Hansen & Ghrist, *Opinion dynamics on discourse sheaves*, SIAM J. Appl. Math. (2021).
- Jiang, Lim, Yao & Ye, *Statistical ranking and combinatorial Hodge theory*, Math. Programming (2011).
- Kane & Lubensky, *Topological boundary modes in isostatic lattices*, Nature Physics (2014).
- Goodrich, Liu & Nagel, *The principle of independent bond-level response: tuning by pruning*, PRL (2015).
- Rocks et al., *Designing allostery-inspired response in mechanical networks*, PNAS (2017).
- Stern, Hexner, Rocks & Liu, *Supervised learning in physical networks*, PRX (2021).
- Dillavou et al., *Demonstration of decentralized physics-driven learning*, Phys. Rev. Applied (2022).
- Ranzani & Aumentado, *Graph-based analysis of nonreciprocity in coupled-mode systems*, New J. Phys. (2015).
- McKeon & Sharma, *A critical-layer framework for turbulent pipe flow*, JFM (2010).
- Clements et al., *Optimal design for universal multiport interferometers*, Optica (2016).
- Keyes, *What makes a good computer device?*, Science (1985).
- Katharopoulos et al., *Transformers are RNNs*, ICML (2020).
- Qian et al., *Scaling Large-Language-Model-based Multi-Agent Collaboration* (MacNet), 2024.
- Zhuge et al., *GPTSwarm: Language Agents as Optimizable Graphs*, ICML (2024).
- Lehman et al., *Evolution through Large Models* (2022).
- Zhang, Lehman, Stanley & Clune, *OMNI*, ICLR (2024).
- Ding et al., *Quality Diversity through Human Feedback* (2024); Bradley et al., *QDAIF* (2024).
- Hu, Lu & Clune, *Automated Design of Agentic Systems* (2024).
- Grillotti & Cully, *Unsupervised behavior discovery with QD optimization* (AURORA), IEEE TEVC (2022).
- Vassiliades, Chatzilygeroudis & Mouret, *CVT-MAP-Elites*, IEEE TEVC (2018).
- Friedman & Dieng, *The Vendi Score*, TMLR (2023).
- Kornblith et al., *Similarity of Neural Network Representations Revisited*, ICML (2019).
- Seung, Opper & Sompolinsky, *Query by committee*, COLT (1992).
- Lazer & Friedman, *The network structure of exploration and exploitation*, Admin. Sci. Quarterly (2007).
