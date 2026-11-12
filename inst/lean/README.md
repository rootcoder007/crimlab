# Proofs

The Lean 4 theorems `csls::proof()` names. Each file imports Mathlib and
nothing else.

To check them (the first run downloads a prebuilt Mathlib, about 5 GB):

```sh
lake exe cache get
lake build
```

A clean build with no `sorry` is the check; CI runs it on every change here.
