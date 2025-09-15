# EdgesCGALWeightedDelaunayND

Calcul des **arêtes du 1‑squelette** de la **Delaunay pondérée (regular triangulation) en dimension d quelconque** à partir de deux `.npy` : `points (N,d)` et `weights (N,)`. Sortie : `.npy` `(M,2)` en `uint64` avec des **paires triées (i<j)**.

## Dépendances (Ubuntu ≥ 22.04)
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libcgal-dev libeigen3-dev
```

## Build
```bash
cmake -S . -B build
cmake --build build -j
```

## Usage
```bash
./build/EdgesCGALWeightedDelaunayND points.npy weights.npy out_edges.npy
```

- `points.npy` : tableau `float32` ou `float64` de forme `(N,d)`, C‑contigu. Endianness supportées : `<` (little), `=` (native), `>` (big). En big‑endian, les données sont **byteswap** automatiquement.
- `weights.npy` : tableau `float32`/`float64` de forme `(N,)` ou `(N,1)`.
- `out_edges.npy` : tableau `uint64` `(M,2)` trié.

## Détails d’implémentation
- **CGAL ND** : `Epick_d<Dynamic_dimension_tag>` comme traits; TDS avec `Triangulation_vertex<Adapter, uint64_t>` pour embarquer l’index d’origine et `Triangulation_full_cell<Adapter>`.
- **Insertion** : `K::Weighted_point_d` construits depuis `K::Point_d` et le poids.
- **Extraction des arêtes** : pas d’itérateur d’arêtes en dD, on parcourt les `finite_full_cells` et on prend toutes les paires de sommets de chaque cellule, avec déduplication.

## Performance
- Compilation en `-O3 -DNDEBUG`.
- Pour de très grands N, l’extraction via full cells peut être coûteuse. C’est la voie standard en dD avec la TDS CGAL.

## Tests rapides
Des exemples sont fournis dans `data/`. Par exemple :
```bash
./build/EdgesCGALWeightedDelaunayND data/example_points.npy data/example_weights.npy out_edges.npy
python3 - <<'PY'
import numpy as np; E=np.load('out_edges.npy'); print(E.shape, E[:10])
PY
```

## Notes
- Fichiers `.npy` Fortran‑order non supportés (comme le reste des projets précédents).
- Si vous manipulez des `points` 2D `(N,2)`, le programme fonctionne, la dimension est lue dans la shape et traitée comme d=2.
- Le nom du projet **n’a pas changé**.