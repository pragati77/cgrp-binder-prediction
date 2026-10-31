# CGRP Binder Prediction

**Using chemistry transformers to rank molecules that might block the CGRP (migraine) receptor.**

> **Honest headline:** this project built and carefully tested a working pipeline, but it did **not** discover a new migraine drug. Screening ~500 approved drugs found no credible candidate. The sections below explain what was done, what worked, what failed, and why, in plain language, from both the biology and the technical side.

---

## Table of contents

0. [TL;DR](#0-tldr)
1. [Biology from scratch](#1-biology-from-scratch)
2. [The data](#2-the-data)
3. [The modeling, step by step](#3-the-modeling-step-by-step)
4. [Results](#4-results)
5. [What went wrong first (and what we learned)](#5-what-went-wrong-first-and-what-we-learned)
6. [Limits and honest caveats](#6-limits-and-honest-caveats)
7. [How MT-DTI fits in](#7-how-mt-dti-fits-in)
8. [How to run it](#8-how-to-run-it)
9. [Repository layout](#9-repository-layout)
10. [Roadmap](#10-roadmap)
11. [FAQ: every question asked while building this](#11-faq-every-question-asked-while-building-this)
12. [Glossary](#12-glossary)
13. [References](#13-references)

---

## 0. TL;DR

| Question | Short answer |
|---|---|
| What is the goal? | Teach a computer to guess how strongly a chemical blocks the CGRP receptor, then use it to look for new candidates. |
| What is the target? | The **CGRP receptor** (a "lock" made of two proteins, CALCRL + RAMP1). Gepant migraine drugs block it. |
| Where did the training data come from? | ChEMBL, a public database of published lab experiments: ~900 chemicals with measured scores. |
| What does the model read? | Only the chemical's structure written as text (a SMILES string). |
| Did it work? | **Partly.** With the 6 real gepants hidden, the decoy-trained model still ranked 5 of the 6 in the top 10% of 487 chemicals. |
| Did it find a new drug? | **No.** None of 495 approved drugs scored as a credible blocker. |
| Is this proof of anything? | No. A prediction is a lead for lab testing, never proof. |

```
   THE BIG PICTURE

   CGRP (the KEY) ──fits into──▶ CGRP receptor (the LOCK) ──▶ pain signalling ON
                                          ▲
        gepant drugs (BLOCKERS) ──────────┘  jam the lock so CGRP can't get in

   Our model: chemical structure ──▶ [ model ] ──▶ predicted blocking strength (pActivity)
```

---

## 1. Biology from scratch

### 1.1 Headaches and migraine

Most headaches are either **primary** (the headache is the condition: tension-type, cluster, migraine) or **secondary** (caused by something else: dehydration, poor sleep, sinus infection, caffeine withdrawal, medication overuse, and, rarely, more serious problems).

**Migraine** is not just a bad headache. It is a neurological condition, typically with:

- moderate to severe throbbing pain, often on one side
- nausea, and sensitivity to light, sound, or smells
- attacks that last roughly 4 to 72 hours
- sometimes an **aura** (visual zigzags, tingling) beforehand

*(When to get urgent care: a sudden "worst headache of my life," or a headache with fever and stiff neck, confusion, weakness, numbness, vision loss, or trouble speaking. This README is a data-science project, not medical advice.)*

### 1.2 Genes, proteins, amino acids

- A **gene** is a stretch of DNA containing the instructions for building one protein.
- A **protein** is a chain built from **amino acids**, like beads on a string. There are 20 types, each written as a one-letter code (A, C, D, E, ...).
- A protein's **sequence** is the order of those letters. Its **length** is how many there are ("aa" = amino acids).
- Genes and their proteins often share a name, and scientists are loose about it. Italic (*CALCA*) usually means the gene; plain capitals (CALCA) often mean the protein.

Almost every protein sequence starts with **M** (methionine), the "start building here" signal.

### 1.3 The CGRP story: a key and a lock

During a migraine, nerve cells release a small protein called **CGRP** (calcitonin gene-related peptide). It travels a tiny distance and binds to a **receptor** on nearby cells, such as cells in the walls of blood vessels around the brain. Binding triggers signals that widen blood vessels and amplify pain.

Think lock and key:

| Role | What it is | Made from gene | Protein used here (UniProt ID) | Length |
|---|---|---|---|---|
| **Key** | α-CGRP (nerve version) | *CALCA* | P06881 | 37 aa active peptide (128 aa precursor) |
| **Key** | β-CGRP (close relative) | *CALCB* | P10092 | 37 aa active peptide (127 aa precursor) |
| **Lock, half 1** | CLR (calcitonin receptor-like receptor) | *CALCRL* | Q16602 | 461 aa |
| **Lock, half 2** | RAMP1 (receptor activity-modifying protein 1) | *RAMP1* | O60894 | 148 aa |

- The key is **released** by nerve cells and floats between cells.
- The lock is **anchored in the membrane** (outer skin) of the receiving cell and never moves.
- The lock needs **both** halves. CALCRL is a long protein that threads through the membrane several times. RAMP1 is a small helper that sits beside it and reshapes it so it recognizes CGRP. (CALCRL paired with other RAMP helpers responds to different signals.)

### 1.4 One gene, two products (splicing)

The *CALCA* gene produces different results in different tissues, through **alternative splicing** (the cell edits the gene's RNA message differently):

- in the **thyroid** → **calcitonin** (a hormone for calcium balance). UniProt entry **P01258** (141 aa) holds this version: calcitonin (32 aa) and katacalcin (21 aa), **no CGRP**.
- in **nerve cells** → **α-CGRP**. UniProt entry **P06881** (128 aa) holds this version.

> A lesson from this project: our first UniProt search asked for one result per gene and returned the thyroid entry for *CALCA*. Always check that a database search returned the entry you meant.

After a precursor is made, the cell cuts out the active 37-amino-acid CGRP peptide (a signal piece at the front and extra pieces around it are removed).

### 1.5 What the sequences actually show

- **α-CGRP and β-CGRP differ at only 3 of 37 positions** (position 3: D vs N; 22: V vs M; 25: N vs S). One receptor responds to both, so a drug that blocks the receptor handles both.
- The first seven letters of CGRP include two **cysteines** (C, at positions 2 and 7) that link to each other and form a small ring that helps CGRP fit its receptor. A model that only reads letters cannot see this ring.
- **CALCRL is rich in L, I, V** (oily amino acids). Runs of them (for example `LFYLTIIGHGLSIASLLISLGIFFYF`) mark stretches that cross the fatty cell membrane.
- **RAMP1 has just one such stretch, near its end**, so it crosses the membrane once, and most of the protein sits outside the cell. That fits its role as a small helper.
- Both receptor proteins begin with an oily stretch that acts as an "address label" and is cut off later.

### 1.6 Medicines that act on this system

| Class | Examples | How it works |
|---|---|---|
| **Gepants** (small-molecule pills/sprays) | rimegepant, ubrogepant, atogepant, zavegepant (marketed); olcegepant, telcagepant (earlier compounds, not marketed) | Sit in the receptor and **block** it. This is what our project models. |
| **Antibody drugs** (injections/infusions) | erenumab (blocks the receptor); fremanezumab, galcanezumab, eptinezumab (grab CGRP itself) | Used mainly to prevent attacks. |
| **Other migraine drugs** | triptans, ergots | Act mainly on **serotonin** receptors, not the CGRP receptor, as far as we know. |

Small molecules block the **receptor** (a lock with a pocket) much more easily than they can grab a tiny floppy peptide floating around. That is why the receptor, not CGRP, is our target.

### 1.7 Is migraine "in the DNA"?

Partly. Genes set how likely you are to get certain headaches; they don't directly cause the pain.

- Migraine runs in families; twin studies suggest genetics explain roughly **40 to 60%** of risk.
- Most migraine is **polygenic**: well over 100 gene variants each nudge risk up slightly.
- A few rare forms come from single genes. Familial hemiplegic migraine is linked to mutations in *CACNA1A*, *ATP1A2*, and *SCN1A*, which encode ion channels (gates that control nerve-cell excitability).
- Triggers still matter (sleep, meals, stress, hormones).
- There is no routine genetic test for common migraine.

CGRP is the best-known **protein** driver: its levels rise in migraine attacks, and blocking it prevents or stops attacks in many people. It is one important link in a longer chain, not the whole story.

---

## 2. The data

### 2.1 Protein sequences (UniProt)

[UniProt](https://www.uniprot.org) is a public protein database. Every protein has a unique ID (like a library card number), so everyone agrees which protein they mean. We used it for the sequences in Section 1. **Notebook 01** does this.

### 2.2 Chemicals and measurements (ChEMBL)

[ChEMBL](https://www.ebi.ac.uk/chembl/) is a public database run by the European Bioinformatics Institute. Scientists publish papers testing chemicals against proteins; ChEMBL staff read the papers and enter the results.

- In ChEMBL, a **target** is *the protein that was tested*. The entry "Calcitonin gene-related peptide type 1 receptor" means "the protein CALCRL." **It is the lock, not the medicine and not the key.**
- Each row is basically: *"chemical X was tested against this protein and gave this binding strength."*
- The chemicals are **candidate blockers**. A few are approved drugs (the gepants). Many are experimental compounds from research programs that never became medicines.

We searched by UniProt ID and found these "folders" of experiments:

| ChEMBL ID | What it is | Raw measurements | Used? |
|---|---|---|---|
| `CHEMBL3798` | CALCRL alone ("CGRP type 1 receptor") | 1,626 | Yes: the biggest pool |
| `CHEMBL2107838` | CALCRL + RAMP1 complex (the true CGRP receptor) | 242 | Yes: most accurate label |
| `CHEMBL2109232` | CALCRL + RAMP2 (adrenomedullin AM1 receptor) | 1 | No: different receptor, too little data |
| `CHEMBL2111191` | CALCRL + RAMP3 (adrenomedullin AM2 receptor) | 9 | No |
| `CHEMBL2111189` | Amylin receptor AMY1 (CALCR + RAMP1) | 1 | No |

We are **not discovering anything new in ChEMBL.** It is a record of past experiments. We use it as **training data**, like a deck of flashcards: front = a chemical, back = how strongly it blocked the receptor.

### 2.3 Where a "measurement" comes from

1. Scientists grow cells that carry the CGRP receptor in a lab dish.
2. They add the natural key, CGRP, and the cells respond.
3. They add a test chemical in different amounts and see whether it stops the response.
4. A strong blocker works even in a tiny amount; a weak one needs a lot.
5. They record the amount needed (for IC50: the amount that cuts the response in half).

Nobody ran these experiments for this project. We reuse the results of lab work done by many groups over about two decades. Because several labs may test the same chemical, one chemical can have several measurements. We take the **median** (middle value) so a single odd result doesn't dominate.

### 2.4 The score: from IC50 to pActivity

| Term | Full name | Meaning |
|---|---|---|
| **IC50** | Half-maximal inhibitory concentration | Amount of chemical that cuts the response in half. Lower = stronger. Depends on how the experiment was set up. |
| **Ki** | Inhibition constant | How tightly a blocker sticks, a property of the chemical and receptor. Lower = tighter. |
| **Kd** | Dissociation constant | Concentration at which half the receptors have a chemical stuck on. Lower = tighter. |

These are measured in nanomolar (nM) and run from 0.01 to 10,000+, with **lower = better**, which is awkward for a model. ChEMBL converts them to **pChEMBL** (our "pActivity"):

> **pActivity = 9 − log10(strength in nM)**

Every +1 means **10× stronger**, and higher is now better:

| Strength | pActivity | In words |
|---|---|---|
| 10,000 nM | 5 | weak |
| 100 nM | 7 | okay |
| 1 nM | 9 | strong |
| 0.1 nM | 10 | very strong |

So the score is *derived from* IC50/Ki/Kd, not equal to IC50. We pool all three measurement types, which is a common shortcut (they aren't perfectly comparable; a stricter project could separate them).

### 2.5 SMILES: a chemical written as text

A **SMILES** string writes a molecule's structure in plain text, e.g. `CCO` is ethanol. It lets a model "read" chemistry the way a language model reads words, and it's how the chemical goes into every model here.

### 2.6 Cleaning, in numbers

| Step | Count |
|---|---|
| Usable measurements (clean IC50/Ki/Kd values) | 1,564 + 229 = **1,793** |
| Unique chemicals (after taking the median per chemical) | **989** |
| Large molecules set aside (molecular weight > 1000, probably peptides) | 73 |
| Small drug-like molecules kept | **916** |

`n_measurements` in the chemical table means *how many times that one chemical was tested*. A chemical with three agreeing measurements is more trustworthy than one with a single test.

### 2.7 What the scores look like

| Group | Chemicals | Share |
|---|---|---|
| weak (below 6) | 90 | 10% |
| medium (6 to 8) | 334 | 36% |
| strong (8 to 10) | 329 | 36% |
| very strong (10+) | 163 | 18% |

Scores range from about 4 to 11. The data leans somewhat strong, because scientists publish compounds that worked; flops are often never reported.

### 2.8 The known-drug check

All six gepants appear in the data with very strong scores (atogepant 10.8, zavegepant 10.6, rimegepant 10.6, olcegepant 10.5, ubrogepant 10.2, telcagepant 8.7). This does two jobs:

1. **Sanity check:** the famous drugs are there with sensible numbers, so we downloaded the right thing.
2. **A built-in exam:** we know their true scores. We hide them from the model while it learns, then see whether it ranks them near the top.

---

## 3. The modeling, step by step

### 3.1 What the model predicts

**Input:** a chemical's structure (SMILES). **Output:** a guess at its pActivity, i.e. how strongly it would block the CGRP receptor in a lab test.

It does **not** compare chemicals to "proven medicines" or measure similarity to anything. Analogy: a house-price predictor learns from past sales and then guesses the price of a new house; it isn't asked whether the house resembles a famous house. The gepants are used only as a *test*.

All our data is about **one receptor**, so the protein sequence is identical in every row. The protein half of an MT-DTI-style model would carry no information, so our models use the **chemical only**.

### 3.2 Splitting the data fairly

If a model studies and is tested on the same chemicals it just memorizes the answers (like a student who sees the exam questions in advance). So we split:

| Group | Job | Analogy |
|---|---|---|
| **Train** (~80%) | the model learns from these | homework |
| **Validation** (~10%) | check progress | practice quiz |
| **Test** (~10%) | final score, used once | final exam |
| **6 gepants** | held out completely | bonus question |

**Scaffold split.** Chemists often make many variations of one core ring skeleton (a *scaffold*). If some variations are in training and others in testing, the exam is too easy. So whole scaffold families stay on one side. The test chemicals then have skeletons the model has never seen. *(We verified: zero scaffolds shared between train and test.)*

**Group sizes after the final split:** train 640, validation 91, test 91, gepants 6, removed 88. Average scores are similar across train/val/test (8.12, 7.97, 8.17), so no group is unfairly easy or hard.

**Leakage and the look-alike problem.** Drug companies make many near-copies of a promising compound, so each gepant had close relatives in the training set. Similarity is measured with the **Tanimoto score** (0 = nothing alike, 1 = identical) on molecular fingerprints:

| Drug | Similarity to nearest training chemical (before removal) |
|---|---|
| atogepant | **0.99** (almost certainly the same molecule recorded twice) |
| rimegepant | 0.85 |
| ubrogepant | 0.84 |
| telcagepant | 0.84 |
| olcegepant | 0.76 |
| zavegepant | 0.68 |

Leaving those in would let the model "rediscover" the drugs by copying a twin. So we moved **88 training chemicals with similarity ≥ 0.6** to a "removed" group. Afterward, every drug's nearest remaining training chemical was below 0.6. Cost: ~12% of training data. Benefit: an honest test.

### 3.3 The baseline: fingerprints + random forest

Before any fancy model, we built a plain one to set the bar:

- **Fingerprint:** each chemical becomes a checklist of 2,048 yes/no boxes ("contains this small fragment?").
- **Random forest:** hundreds of simple decision trees each guess a score, and the answers are averaged.

If a transformer can't beat this, it isn't earning its keep.

### 3.4 Pretrained chemistry transformers

Instead of training a transformer from scratch on only ~640 chemicals, we borrowed models already trained on enormous molecule collections:

| Model | Pretrained on | Notes |
|---|---|---|
| **ChemBERTa-77M-MLM** (`DeepChem/ChemBERTa-77M-MLM`) | PubChem SMILES (ChemBERTa-2 work used up to 77M compounds) | BERT-style; 384-number summary per molecule |
| **MoLFormer-XL** (`ibm/MoLFormer-XL-both-10pct`) | 10% of ZINC + PubChem (the full MoLFormer-XL paper used 1.1B molecules) | 768-number summary per molecule |

How we use them:

1. The transformer reads a SMILES and outputs a list of numbers (an **embedding**: a summary of what it "understands" about the molecule). The transformer stays frozen.
2. A **ridge regression** (a very simple model: a weighted sum of those numbers, with built-in protection against overfitting) predicts the score from the embedding.

Harmless warnings you may see when loading: `UNEXPECTED lm_head` (the unused fill-in-the-blank head) and `MISSING pooler` (a part we don't use).

### 3.5 How we score the models

| Measure | Plain meaning | Good |
|---|---|---|
| **Typical miss** (RMSE) | usual size of the error in score units; 1.0 ≈ 10× off in strength | lower |
| **Ranking skill** (Spearman) | does sorting by the model's guesses match the true order? 0 = none, 1 = perfect | higher |
| **Rank** | position on a leaderboard sorted from strongest (1) to weakest | the model's rank should be close to the true rank |
| **AUC** | how well it separates real blockers from decoys; 0.5 = coin flip, 1.0 = perfect | higher |

**Rank matters most for our goal.** We'll never trust the exact numbers; we'll sort a big list and test the top few in a lab. A model whose guesses are a bit off but whose *order* is right is still useful.

Reading a line such as `ATOGEPANT predicted 9.8 (actual 10.8) | rank 14 / 2`: the model guessed 9.8, the lab value is 10.8, and the model placed it 14th while it is truly 2nd.

### 3.6 Repeating the exam

One 91-chemical test set can flatter or punish a model by luck, and the gap between two models might be noise. So we repeat the exam **five times** with different scaffold splits, and report an average and a spread (`±`). If a gap is smaller than the spread, we call it a tie. We compare models **within the same round** (paired), since some rounds are easier for everyone.

### 3.7 Decoys: teaching the model what a non-blocker looks like

Every chemical in our data was made or tested because someone hoped it would block this receptor. The model therefore never saw a clear "nothing like a blocker" example, so it could only say "stronger" or "weaker," never "no." Our first screen showed the consequence (Section 5).

The fix: add **decoys**, ~4,700 random small molecules from ChEMBL that were never tested on this receptor, labeled with the floor score of **4.0**. Random chemicals very rarely block a specific receptor, though the assumption is not proven for each one.

Because real blockers are heavier than random chemicals (median molecular weight 546 vs 409), molecular weight alone separates them with AUC 0.90. So we also trained on **size-matched decoys** (molecular weight 440 to 665) and tested everything on a size-matched test.

### 3.8 The screen of approved drugs

We scored ~2,400 approved small-molecule drugs (ChEMBL `max_phase = 4`), excluding anything already tested on the CGRP receptor, with safeguards:

- **Size range:** only drugs between 440 and 665 in molecular weight (495 drugs), since that is the range the size-matched model has seen.
- **Consensus:** a drug must score well under **both** the all-decoys model and the size-matched model (we take the lower of the two).
- **Closeness:** report each drug's similarity to the nearest known blocker, a measure of how much the prediction is extrapolation.
- **Clipping:** scores are limited to the 4 to 11 range seen in the data.

---

## 4. Results

### 4.1 Models compared (5-round average ± spread)

| Model | Typical miss | Ranking skill | The 6 medicines' average rank (lower = better; ~19 perfect, ~90 random) |
|---|---|---|---|
| fingerprint + random forest | 0.94 ± 0.25 | 0.74 ± 0.19 | 76.7 |
| ChemBERTa + ridge | 0.93 ± 0.07 | 0.79 ± 0.06 | **56.9** |
| MoLFormer + ridge | **0.88 ± 0.06** | **0.82 ± 0.04** | 60.4 |

- On ordinary accuracy the three are close, but the fingerprint model is much less consistent (it dropped to ranking skill 0.42 in one round).
- **Both transformers beat the fingerprint model at ranking the hidden medicines in all 5 rounds.**
- ChemBERTa vs MoLFormer is essentially a tie.
- Absolute performance is modest: a hidden medicine lands around position 57 to 60 in a list of ~170, not near the top.

### 4.2 Hidden-drug test, single split (ChemBERTa vs baseline; rank out of 97, 1 = strongest)

| Drug | True rank | Baseline | ChemBERTa |
|---|---|---|---|
| atogepant | 2 | 50 | **14** |
| zavegepant | 4 | **9** | 12 |
| rimegepant | 6 | 78 | 90 |
| olcegepant | 7 | **48** | 61 |
| ubrogepant | 19 | 54 | **29** |
| telcagepant | 45 | 75 | **43** |

Rimegepant is the stubborn miss for every model (predicted ~6 vs actual 10.6). The reason is unknown.

### 4.3 Effect of decoys (ranks of the 6 medicines among 487 chemicals: 91 test blockers + 390 size-matched decoys + 6 drugs)

| Training | Blocker-vs-decoy AUC | Ranking among real blockers | Decoys predicted ≥ 8 | The 6 medicines' ranks |
|---|---|---|---|---|
| no decoys | 0.83 | 0.79 | 46 of 390 | 18, 20, 44, 78, 135, 275 |
| all decoys | 1.00 | 0.62 | 0 of 390 | 12, 22, 29, 52, 54, 71 |
| size-matched decoys | 1.00 | 0.67 | 0 of 390 | 7, 11, 30, 45, 47, 91 |

- Decoys eliminate absurdly high scores and move the medicines up (5 of 6 in the top 10% with size-matched decoys).
- **Trade-off:** ranking *among real blockers* drops from 0.79 to 0.62–0.67. The model got better at "is it a blocker?" and somewhat worse at "how strong?"
- AUC 1.00 should be read cautiously: blockers and random ChEMBL chemicals are very different kinds of chemistry, so this is an easy separation. It shows the model learned what known blocker families look like, not that it can recognize a brand-new kind.

### 4.4 Screen of 495 approved drugs

| Check | Result |
|---|---|
| Predicted ≥ 8 by the all-decoys model | 0 of 495 |
| Predicted ≥ 8 by the size-matched model | 1 of 495 |
| Predicted ≥ 8 by **both** | **0** |
| Best consensus score | ~7.2 ("medium") |
| Closeness of top drugs to known blockers | 0.11 to 0.26 (low: mostly guesswork) |

The top of the list is a mix of antivirals, cancer drugs, antibiotics, and blood-pressure drugs (e.g. danicopan, boceprevir, fostemsavir, nirmatrelvir). They have in common that they are large, many-ringed, amide-rich molecules, like the gepants. That is a resemblance in shape and size, **not** evidence of binding.

**Conclusion:** no credible repurposing candidate. The other ~1,900 approved drugs fall outside the size range the models can handle, so nothing can be said about them.

---

## 5. What went wrong first (and what we learned)

**The first screen was nonsense.** Scoring ~2,330 approved drugs with the original model gave **338 predicted ≥ 8**, and the top of the list included digoxin, cetyl alcohol (a waxy fatty alcohol), docosanol, and several antibiotics, with scores above 11, higher than any chemical in the whole dataset. Cetyl alcohol was predicted to out-block the real gepants.

**Why:**

1. The model had only seen CGRP-receptor blockers, so it had no concept of "not a blocker."
2. Our earlier tests couldn't catch this, because every test chemical was also blocker-like. They measured "which blocker is stronger?" not "is this a blocker at all?"
3. Ridge regression extrapolates in straight lines, so unfamiliar molecules can produce wild numbers.
4. All those drugs had a closeness of 0.09 to 0.26, nowhere near anything the model studied.

**Lessons:**

- A test must match the way the model will be used.
- Always report how far a prediction is from the training data.
- A hidden-drug test is only meaningful once look-alikes are removed (atogepant had a 0.99 twin in training).
- A single small test set can mislead, so repeat it.
- Check that a database search returned the entry you meant (the *CALCA* thyroid-vs-nerve mix-up).
- Catching a failure early with a cheap test is better than announcing a discovery.

---

## 6. Limits and honest caveats

- **Chemical structure only.** The models don't see the 3D shape of the receptor or the binding pocket, so they cannot say *where* or *how* something binds. The CGRP ring structure and the cysteine bond are invisible to them.
- **Blocking vs activating.** The data are inhibition measurements; the model cannot tell a blocker from an activator.
- **Binding is only one step.** A real drug must also be safe, absorbed, reach the right place, stay in the body the right time, and avoid side effects. We predict none of that.
- **Small data.** ~900 chemicals is small for deep learning, and the data are dominated by a few gepant-like chemical series.
- **Publication bias.** Failed compounds are often unpublished, so weak binders are under-represented.
- **Mixed measurement types.** IC50, Ki, and Kd are pooled, which isn't perfectly apples-to-apples.
- **Decoys are assumed inactive**, not proven.
- **Tiny drug test.** Six drugs is a very small sample; one lucky or unlucky drug swings the picture.
- **"Safe enough to reach patients"** applies to a drug's existing use, not automatically to a new use, dose, or long-term migraine treatment.
- **Not medical advice.** Nothing here should guide treatment decisions.
- Results can shift slightly if ChEMBL adds data or library versions change.

---

## 7. How MT-DTI fits in

[**MT-DTI**](https://github.com/deargen/mt-dti) (Molecule Transformer Drug-Target Interaction; Shin et al., 2019) predicts how strongly a drug binds a protein using (1) a molecule transformer pretrained on a huge set of PubChem SMILES, and (2) a protein encoder, combined into a binding prediction.

**Why we didn't run the original repo:**

- It requires **TensorFlow 1.12** (2018), which doesn't run on modern Colab without heavy repair.
- The pretrained and fine-tuned weights are in a `data.tar.gz` on Google Drive, and the pretraining checkpoint is on a Google Cloud Storage path.
- It was fine-tuned on **KIBA**, a kinase-inhibitor dataset, a different protein family from ours.

**What we kept of its idea:** a pretrained molecule transformer plus a simple prediction layer. ChemBERTa/MoLFormer are modern equivalents of the pretrained transformer, and ridge regression is the prediction layer.

**What we left out:** the protein branch. With one receptor, the protein input never varies. To use MT-DTI as designed, we'd need binding data for many related proteins (see the roadmap).

---

## 8. How to run it

Everything is built for **Google Colab** (free tier is enough; a GPU speeds up the transformer steps).

1. **Notebook 01** (`01_cgrp_fetch_proteins.ipynb`): downloads the protein sequences from UniProt and saves `targets_clean.csv`. *Not needed for the later notebooks, which use chemicals only.* Needs internet access to `rest.uniprot.org`; if you get `403 Forbidden`, the environment or server is refusing the request (see FAQ).
2. **Notebook 02** (`02_cgrp_chembl_and_models.ipynb`): **Runtime → Run all.** It will:
   - download chemicals and scores from ChEMBL
   - build the fair scaffold split and remove gepant look-alikes
   - run the baseline, ChemBERTa, the 5-round test, and MoLFormer
   - save progress to **Google Drive** (`MyDrive/cgrp_project/`) so a Colab restart doesn't wipe your work; saved files are reloaded automatically
3. **Notebook 03** *(in progress)*: decoys (Section 3.7) and the approved-drug screen (Section 3.8).

Files saved to Drive: `cgrp_compounds.csv`, `cgrp_split_final.csv`, `decoys.csv`, `approved_drugs.csv`, and results tables.

**Requirements:** Python 3, RDKit, pandas, numpy, scikit-learn, scipy, requests, PyTorch, Hugging Face `transformers` (all present or one-line installs in Colab).

---

## 9. Repository layout

```
cgrp-binder-prediction/
├── README.md
├── 01_cgrp_fetch_proteins.ipynb        # UniProt: key and lock sequences
├── 02_cgrp_chembl_and_models.ipynb     # ChEMBL data, fair split, baseline, ChemBERTa, MoLFormer
├── 03_decoys_and_screening.ipynb       # (in progress) decoys + approved-drug screen
└── results/                            # saved tables (optional)
```

---

## 10. Roadmap

1. **Package notebook 03** (decoys + screening) so the repo reproduces every result in this README.
2. **Better screening libraries.** Approved drugs were a long shot. Chemical libraries of drug-like, purchasable compounds designed for similar targets are a fairer hunt.
3. **Fine-tune the transformer** end to end on our data (closest to MT-DTI). With only ~900 chemicals this risks overfitting, so it needs the repeated-split test and decoys.
4. **Multi-receptor, MT-DTI-style model.** Add binding data for related receptors (for example other class B GPCRs that share parts with CALCRL) so the protein sequence actually varies and the protein branch matters. This may help it generalize beyond gepant-like chemistry.
5. **Add 3D structure.** Dock the model's top picks into a published CGRP-receptor structure to bring in the information that sequence-only methods lack.
6. **Separate "is it a blocker?" from "how strong?"** with two models, since squeezing both into one score cost ranking accuracy (Section 4.3).
7. **A blind prediction test on commercially available drugs** (details below).

### Planned blind test: predicting on drugs "without knowing the answer"

The aim is a genuinely prospective test: make predictions first, check afterward. To keep it honest:

- **Freeze everything first:** model, features, decoys, thresholds. No tuning after looking at results.
- **Write the ranked predictions to a dated file** before checking any literature.
- **Compare afterward** against known pharmacology and published CGRP-receptor testing.
- **Be careful about what "answer" means.** Many drugs that treat migraine (triptans, ergots, some preventives) work through *other* targets, so "treats migraine" is not the same as "blocks the CGRP receptor." Our model predicts only the latter, so judging it against a list of migraine drugs would mix up two different questions.
- Hide known CGRP-receptor drugs (and their close relatives) during training, as already done for the gepants.

---

## 11. FAQ: every question asked while building this

### Biology

**What causes headaches, and how is migraine different?**
See Section 1.1. Tension headaches are dull pressure from muscle tension, stress, posture, or eye strain; migraine is a neurological condition with throbbing pain, nausea, light/sound sensitivity, and sometimes aura.

**Is some code in the DNA that causes this kind of pain?**
Partly. Migraine has a strong genetic component (roughly 40 to 60% of risk) but is mostly polygenic, with a few rare single-gene forms. See Section 1.7. Genes raise or lower likelihood; they don't directly cause pain, and triggers still matter.

**Is there a specific protein that causes this?**
For migraine, **CGRP** is the best-established one: its levels rise during attacks, and drugs that block it (or its receptor) prevent or stop attacks. It's a key link, not the whole story. Other molecules, such as PACAP and serotonin, are also involved.

**Are CALCA and CALCB genes?**
Yes. They are names of two stretches of DNA. Each gene has a protein it makes. Italic usually means the gene, plain capitals often mean the protein (Section 1.2).

**Is CGRP the protein "released" by CALCA and CALCB too?**
Yes. *CALCA* makes α-CGRP and *CALCB* makes β-CGRP (they differ at only 3 of 37 amino acids). Strictly, genes don't release anything: the cell builds the protein from the gene, and nerve cells release it.

**How are CALCRL and RAMP1 different from CALCA and CALCB?**
CALCA/CALCB make the **key** (the CGRP signal, small, released, floats around). CALCRL/RAMP1 make the **lock** (the receptor, anchored in the cell membrane, receives the signal). CALCRL is the main part; RAMP1 is a helper that reshapes it so it recognizes CGRP. Both are needed (Section 1.3).

**Does CALCA make CGRP or calcitonin?**
Both, depending on tissue: calcitonin in the thyroid, CGRP in nerve cells, through alternative splicing (Section 1.4). In UniProt these are two entries (P01258 and P06881).

**What are the "key" and "lock" sequences?**
The key is a 37-letter peptide (α- or β-CGRP). The lock is two longer sequences (CALCRL, 461 aa; RAMP1, 148 aa). See Section 1.5 for what the letters reveal.

**Why target the receptor and not CGRP?**
Small drug-like molecules block a receptor's binding pocket much more easily than they can grab a tiny floppy peptide. Gepants work this way.

### The tables from UniProt (Notebook 01)

**What do `gene`, `role`, `uniprot_id`, and `length_aa` mean?**
`gene`: the gene's name (a code). `role`: my own label for the protein's job in the story. `uniprot_id`: the protein's unique database ID. `length_aa`: number of amino acids in the chain.

**Why is CALCA 141 amino acids if CGRP is only 37?**
The cell first makes a long precursor and then cuts out the active piece. In addition, the UniProt entry first returned for *CALCA* was the thyroid version, which contains calcitonin and no CGRP (Section 1.4).

### Data and ChEMBL

**What are we doing in ChEMBL?**
Collecting **training data**: past lab results of chemicals tested against the CGRP receptor, like a deck of flashcards. We aren't discovering anything there (Section 2.2).

**What is meant by "CGRP receptor" in ChEMBL? Is it the medicine?**
No. It is the **lock** (the receptor protein). The medicines are chemicals tested *against* it. The data are candidate blockers with measured strengths, a few approved drugs among many experimental compounds.

**Are we finding chemicals similar to the lock?**
No. The chemicals aren't similar to the lock; they're objects that **fit into** it. They also usually don't resemble the natural key (CGRP is a peptide; gepants are small drug-like molecules that fit the same pocket).

**Where do the measurements come from?**
Published lab experiments, entered into ChEMBL by curators (Section 2.3). No one ran new experiments for this project.

**What is `n_measurements`?**
In the target table: how many test results exist for that protein. In the chemical table: how many times *that one chemical* was tested (the median of those values is its score).

**Why did the 1,626 measurements become 989 chemicals?**
Many chemicals were tested more than once, and we kept only clean standardized values; repeats are combined by taking the median.

**Is the score IC50?**
Not exactly. It's pActivity (pChEMBL) = 9 − log10(nM), derived from IC50, Ki, or Kd. Higher is stronger, and each +1 is 10× stronger (Section 2.4).

**What does Kd stand for?**
**Dissociation constant**: the concentration at which half the receptors are occupied. Low Kd = grips tightly. Related: Ki (inhibition constant) and IC50 (half-maximal inhibitory concentration).

**What is RDKit?**
A free Python toolkit for chemistry. Plain Python doesn't know what a molecule is; RDKit reads SMILES, computes properties like molecular weight, finds scaffolds, draws structures, and makes fingerprints.

**What is the "do we have the known migraine drugs?" code for?**
A sanity check (are the famous drugs in our data?) and the setup for our hidden-drug exam (Section 2.8).

**What do the numbers 8.17, 7.97, 8.17 mean?**
The average pActivity of the train, validation, and test groups. They're nearly equal, which means the split is balanced: no group is unfairly easy or hard.

### The models

**What is the model actually predicting? Is it similarity to the original medicines?**
It predicts the **pActivity** (blocking strength) of a chemical from its structure. It doesn't measure similarity to medicines. The medicines are used only as a hidden test (Section 3.1).

**What does "predicted" mean in the output?**
The model's guess at the score; "actual" is the real lab value.

**What is "rank"?**
Position on a leaderboard sorted strongest-first. The "model's rank / true rank" pair shows how close the model's ordering is to reality (Section 3.5).

**Why repeat the testing?**
A single small test set can mislead: a 0.05 gap between two models could be luck. Five different splits give an average and a spread, so we can tell a real difference from noise (Section 3.6).

**Why not use the original MT-DTI code, or "already existing transformers"?**
We did use existing pretrained transformers (ChemBERTa, MoLFormer). The original MT-DTI repo needs TensorFlow 1.12 and Google-Drive-hosted weights, and was tuned on kinase data (Section 7).

**What does "doesn't beat the baseline" mean for a fancy model?**
If a quick, simple method does as well, the fancy model isn't adding value. Our transformers did beat the baseline at ranking the hidden medicines in all five rounds.

### Results

**What was the conclusion?**
A working, carefully validated pipeline, but no new drug. The best model recognizes gepant-like molecules well; screening ~500 approved drugs produced no credible candidate (Section 4 and 5).

**Does danicopan help with migraines?**
No, as far as we know. It's a medicine for a rare blood disorder (paroxysmal nocturnal hemoglobinuria) that blocks a protein in the immune system's complement pathway. It ranked first in our screen only because it's a large, ring-rich molecule like the gepants, and even then its score was just "medium" with a closeness of 0.21 to anything the model had seen. A high rank is a reason to look closer, not a sign the drug works.

**Ergotamine and dihydroergotamine ranked near the top, and they're migraine drugs. Is that a validation?**
Probably not. As far as we know they act mainly on serotonin receptors, and the model never saw any serotonin data. They're large many-ringed molecules, which score well here regardless.

**Why did the first screen give such absurd answers?**
See Section 5: the model had never seen a non-blocker, our tests couldn't reveal that, and ridge regression extrapolates wildly.

### Technical hiccups

**`ModuleNotFoundError: No module named 'rdkit'`**
RDKit isn't installed in the session. Run `%pip install -q rdkit` first and wait for it to finish (restart the session if needed).

**`NameError: name 'dev' is not defined`**
A variable from an earlier cell is missing, usually because Colab restarted and cleared everything. That is why Notebook 02 saves progress to Google Drive.

**"Whole progress is gone"**
Colab runtimes reset after being idle. Notebook 02 rebuilds everything and caches files on Drive so a restart takes seconds to recover from.

**`403 Forbidden` from `rest.uniprot.org`**
The server or a network layer refused the request, either because the environment restricts which sites it can reach or because UniProt temporarily rejected it. Only Notebook 01 needs UniProt, and the models don't need the protein sequences, so you can skip it.

**Warnings like `HF_TOKEN`, `UNEXPECTED lm_head`, `MISSING pooler`**
Harmless: an unauthenticated-download notice and unused parts of the pretrained model files.

---

## 12. Glossary

| Term | Plain meaning |
|---|---|
| **Amino acid** | One of 20 building blocks of proteins, each written as a letter |
| **AUC** | Score for separating two groups (0.5 coin flip, 1.0 perfect) |
| **Baseline** | A simple model used as the bar a fancier one must beat |
| **CGRP** | Calcitonin gene-related peptide; the "key" that rises during migraine |
| **ChEMBL** | Public database of published chemical-vs-protein experiments |
| **Consensus** | Requiring two models to both score a chemical highly |
| **Decoy** | A chemical assumed not to block the receptor, used to teach "no" |
| **Embedding** | A list of numbers summarizing a molecule, produced by a transformer |
| **Fingerprint** | A yes/no checklist of small chemical fragments in a molecule |
| **Gene** | A stretch of DNA with instructions for one protein |
| **Gepant** | A class of small-molecule drugs that block the CGRP receptor |
| **IC50 / Ki / Kd** | Measures of how much chemical is needed to block or bind (lower = stronger) |
| **Leakage** | Test answers sneaking into training, making a test too easy |
| **Peptide** | A very short protein |
| **pActivity / pChEMBL** | 9 − log10(strength in nM); higher = stronger, +1 = 10× |
| **Precursor** | A long protein made first and cut down to the active piece |
| **Random forest** | Many simple decision trees voting |
| **Receptor** | A protein on a cell's surface that receives signals (the lock) |
| **Ridge regression** | A simple weighted-sum model that resists overfitting |
| **SMILES** | A chemical structure written as text |
| **Scaffold** | The core ring skeleton of a molecule |
| **Splicing** | Editing a gene's RNA message differently in different tissues |
| **Tanimoto similarity** | Fingerprint-based closeness of two molecules (0 to 1) |
| **Transformer** | A neural network that reads sequences (text, SMILES, proteins) |
| **UniProt** | Public protein database with unique IDs |

---

## 13. References

- Shin, B., Park, S., Kang, K., Ho, J.C. (2019). *Self-Attention Based Molecule Representation for Predicting Drug-Target Interaction.* PMLR 106:230–248. Code: <https://github.com/deargen/mt-dti>
- Ahmad, W., Simon, E., Chithrananda, S., Grand, G., Ramsundar, B. (2022). *ChemBERTa-2: Towards Chemical Foundation Models.* <https://arxiv.org/abs/2209.01712>
- Chithrananda, S., Grand, G., Ramsundar, B. (2020). *ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction.* <https://arxiv.org/abs/2010.09885>
- Ross, J. et al. *Large-Scale Chemical Language Representations Capture Molecular Structure and Properties* (MoLFormer). <https://arxiv.org/abs/2106.09553>
- Models on Hugging Face: `DeepChem/ChemBERTa-77M-MLM`, `ibm/MoLFormer-XL-both-10pct`
- UniProt: <https://www.uniprot.org> · ChEMBL: <https://www.ebi.ac.uk/chembl/>
- RDKit: <https://www.rdkit.org>

---

*This project is for learning and research. Model predictions are leads for laboratory testing, not evidence of drug activity, and nothing here is medical advice.*
