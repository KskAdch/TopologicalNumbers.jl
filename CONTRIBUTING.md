# TopologicalNumbers.jl への貢献

不具合報告、文書改善、数値アルゴリズムや模型の追加を歓迎します。大きな変更に着手する前に Issue で目的と設計を共有すると、重複作業やAPIの不整合を避けやすくなります。

## 開発環境

Julia のパッケージ環境を有効にし、依存関係を準備します。

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

描画テストと Pfaffian の参照計算では、Python 環境を初回実行時に準備するため時間がかかることがあります。

## 実装方針

- 変更は一つの目的に絞り、既存APIとの互換性を保ってください。
- Julia や固体物性を学び始めた利用者にも意図が伝わるよう、公開APIには日本語の docstring を付けてください。
- コメントは数式上の前提や実装上の注意を補足するために使い、自明な処理の説明は避けてください。
- ループ内の一時配列、不要なコピー、型不安定なコンテナを避けてください。
- 計算結果や入力契約を変更した場合は、docstring、`docs/`、README のうち影響する文書も更新してください。

## 数値変更の検証

数値アルゴリズムや Hamiltonian を変更する場合は、変更箇所に応じて次をテストしてください。

- 既知の模型で期待するトポロジカル数が得られること
- メッシュ数を変えたときに期待値へ収束すること
- `rounds=false` の未丸め値と返値型が妥当であること
- 縮退、ギャップ閉鎖、周期境界、最小メッシュなどの境界条件
- MPI 対応経路を変更した場合は、逐次計算と複数rank計算が一致すること

乱数を使うテストは再現可能な seed を固定し、浮動小数点比較には物理的に妥当な許容誤差を明示してください。

## 品質確認

まず変更箇所に近いテストを実行し、提出前に全テストを実行します。

```powershell
julia --project=. -e 'using Pkg; Pkg.test()'
```

性能試験も実行し、実行時間と割り当て量に意図しない退行がないことを確認します。

```powershell
julia --project=. test/performance/runtests.jl
```

ソースとテストを整形します。

```powershell
julia --project=test -e 'using Pkg; Pkg.instantiate(); using JuliaFormatter; format(".")'
```

公開APIや数式を変更した場合は、ローカルパッケージを使って文書をビルドします。

```powershell
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

生成された `Manifest.toml` や描画ファイルは、変更目的に必要でない限りコミットしないでください。

## Pull Request

Pull Request には次を記載してください。

- 解決する問題と変更方針
- API、計算結果、性能への影響
- 追加または更新したテスト
- 実行した品質確認とその結果
- 関連する Issue

レビューしやすい規模を保ち、無関係な整形やリファクタリングを同じ Pull Request に混ぜないでください。すべてのCIが成功し、指摘へ対応した後にマージ対象となります。
