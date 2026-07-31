# テスト依存関係

Julia パッケージは `test/Project.toml`、Python パッケージは
`test/conda/CondaPkg.toml` で宣言します。テスト実行中に依存関係を追加してはいけません。
`runtests.jl` は Python 依存宣言専用の環境を `LOAD_PATH` に追加し、
`Pkg.test()` の一時環境でもこの宣言が参照されるようにします。

複数の Julia バージョンで互換性を検証できるよう、`test/Manifest.toml` は管理しません。
依存関係を更新するときは互換範囲を見直し、対応する全テストを実行してください。
