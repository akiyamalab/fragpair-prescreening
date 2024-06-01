# FragPair-PreScreening
---

有望なフラグメント空間配置対に基づく化合物プレスクリーニング

## Installation

Dockerfile is at `.devcontainer/Dockerfile`. It explains the requirement and what packages are needed.

### Requirement
- Boost (latest)
- Open Babel 2.4.1
  - **OpenBabel 3.x (latest) is not supported.**

### Build

```sh
make all
```

## Usage

### Execution Sample
```sh
./atomgrid-gen dude/testdata/testconf.in   # make atom grid
./fragment-query dude/testdata/testconf.in # create fragpair query
```

### Output
フラグメントベースの化合物立体配座検索システム（非公開）に渡すクエリ
```
FRAGMENT_1 FRAGMENT_2 Distance_g Distance_r
```
- `FRAGMENT_1`,  `FRAGMENT_2`: フラグメント対のSMILES
- `Distance`: フラグメント対間の相対距離
  - 相対距離範囲は $[$ `Distance_g`-`Distance_r`, `Distance_g`+`Distance_r` $)$
