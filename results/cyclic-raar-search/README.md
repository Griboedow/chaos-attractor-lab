# Cyclic three-circle RAAR screen

See the [report](../../docs/cyclic-raar-search.md). Four completed directories retain all screen outcomes and all five follow-ups. `inventory.json` lists 97 files; `omissions.json` records six optional arrays totaling 23,820,473 bytes. Large arrays remain outside Git, with exact input/metric subsets in separately named compact archives. No complete trajectory is replaced by a subset under its original name.

The first screen predates the addition of explicit beta-range arguments. Its executed source is preserved. The second uses the extended source with the same formula and derivative implementation. The protocol files bind each run to its actual sources.

```text
python research/search_cyclic_raar.py --out D:/Codex/new-cyclic-scout01 --power 9 --beta-min 0.6 --beta-max 0.999
python research/follow_cyclic_raar.py --source D:/Codex/new-cyclic-scout01/screen.npz --out D:/Codex/new-cyclic-follow01
python research/search_cyclic_raar.py --out D:/Codex/new-cyclic-scout02 --power 12 --beta-min 0.93 --beta-max 0.9999
python research/follow_cyclic_raar.py --source D:/Codex/new-cyclic-scout02/screen.npz --out D:/Codex/new-cyclic-follow02
```

All outputs must be new directories. The five selected clouds retain positive finite-time growth but were deprioritized after inspecting their thin native-coordinate geometry. The algorithm is known; an exhaustive equivalence search for these configurations has not been completed. This archive records a search outcome, not a new-attractor claim.
