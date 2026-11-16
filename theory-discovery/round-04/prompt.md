# TDF v3.0 — Round 4: Failure-Driven Theory Discovery

Reconstructed executable summary of the prior Round 4 prompt, not a verbatim copy. Round 4 has not been executed.

あなたは理論産業組織論、ミクロ経済理論、ゲーム理論、情報経済学、契約理論、メカニズムデザイン、経済地理学、金融経済学、経済学説史の研究チームである。目的は既存研究を救済することではなく、数学的成果・反例・負の結果から経済現象の説明構造を変える原理を発見することである。

正規の research-paper-workflow (https://github.com/ryotamatsuki/research-paper-workflow) の THEORY_PAPER_RESEARCH_PIPELINE.md、GOVERNANCE.md、theory-discovery 内の TDF v3.0 マスタープロンプト、Framework、Pattern Catalog、Database、User Guide を参照する。取得できない資料・版・読解深度を明示。データベースは着想資料であり新規性の網羅的審査資料ではない。このラボの research-log、stage-reports、computations、literature、archive を確認する。

## 出発点
Round 2 Geography of Comparison、Round 3 Geography of Informational Rents の成果と Stage 7.5 NO-GO を尊重する。個別の情報劣化と非効率均衡の断絶、分母を除去しても主要結果が残る ablation、協調の失敗への還元、診断比較モデルの契約基礎欠如、数学的完成度と経済的貢献の乖離を失敗資源として扱う。q、x、u の分解、均衡効率性、弱非循環性、純粋立地 SPNE の限定モデル、厚生比較、情報レント誘因逆転、gamma の比較モデル間の符号逆転を再確認する。未証明の FIP、新規性、制度的実装を確立済みと扱わない。

## 探索の方向
A: Economics of Information Rights（価値と収益帰属）。B: Endogenous Decision Architecture（比較単位・配分制度の内生化）。C: Limits of Market Self-Correction（外部性と持続する均衡歪みを区別）。D: Beyond Agglomeration Economics（地理・集積・資本市場を取り除いた抽象構造）。E: 自由探索。現行モデル・対象分野・数式を放棄してもよい。少なくとも一部は企業集積・資本配分を前提としない。

## 実行手順
1. Theme Map; 2. Baseline Theory Audit; 3. Pattern Selection; 4. Candidate Generation; 5. Prior-art Falsification; 6. Value Screen; 7. Dialogue Loop。
12のパターン P01 Endogenization、P02 Unit Redefinition、P03 Concept Decomposition、P04 Hidden Function、P05 Selection Feedback、P06 Dynamic State、P07 Rights Reinterpretation、P08 Strategic Interaction、P09 Robust Bounds/Invariance、P10 Institutional Boundary、P11 Equivalence Breaking、P12 Policy-Principle Domain Test を検討（正式名称と版は正規資料を優先）。特に P02/03/07/09/10/11。各失敗→原因となる仮定→説明できる親理論→説明できない現象→反転した問い→パターン→最短反証、の Failure-to-Discovery Matrix を作成。

構造的に異なる8〜15候補（数より差異を重視、類似案は統合）。候補ごとに ID、英和タイトル、現象、標準説明、説明限界、新原理、変更する所与、因果連鎖、主体・行動・制約・利得・タイミング・均衡、条件付き命題、反実仮想、経済学的意義、成立・破綻条件、パターン、危険な親理論、最短反証を明記。Mechanism Fingerprint: players/actions/state/information/timing/payoff channel/equilibrium/proposition/policy ranking。候補と最接近先行研究について同一の指紋を作り、分野名を除去し同型・特殊例・既知定理の系を攻撃する。

反証担当は生成担当と分ける（不可能なら限界明記）。一次文献、本文・証明、前後引用、後続研究、最新WP、著者の限界を確認。Competitive Information Design, Bayesian Persuasion, Correlation Games, Strategic Information Acquisition, Information Rights, Incomplete Contracts, Relative Performance, Organization, Endogenous Information/Choice Architecture, Coalition/Hedonic/Potential/Weakly Acyclic Games, Market Design, Financial Intermediation, Economic Geography を探索。取得不能資料は記録しその候補の新規性は UNRESOLVED。先行研究のタイトルだけで差異を判断しない。

上位候補のみ最小予備分析（全数値計算 Python、必要なら SymPy）。反例・符号・比較静学・最小均衡仮説を検討。Ablation Test を必須とする。固有メカニズムを除いて主要結果が不変なら降格・棄却。説明対象・厚生定義・人口・制約・権利・情報・外部選択肢を揃える。驚きを利得に直接仮定しない。新名称、変数追加、既存定理の系、既知の外部性や協調の失敗の再記述、ミクロ基礎の欠如、経済的学習の乏しさは棄却要因。全件棄却可。

正規 TDF の評価（説明構造の変更、独立メカニズム、非自明反実仮想、最小モデル化、経済的重要性、検証可能性、各0〜4点）と Round 4 固有評価（不可欠性、還元耐性、ミクロ基礎、個別誘因・均衡接続、分野横断、編集者価値、追加学習）を分け、点数を新規性証明としない。Stage 7.5 先取り審査として、数式なしに言える原理、既存理論が説明できない理由、除去時の結果、接続、掲載価値、最強のリジェクト理由と反証結果を明記。旧テーマの双方を比較対象とし、見込みの改善がなければ保存のみ。

## 出力
I Research Retrospective / II Theory Discovery / III Prior-Art Falsification / IV Mathematical Screening / V Research Investment Decision / VI Next Research Project Proposal。最終候補2〜3件または全件棄却。最有力案について仮題、問い、原理、先行との差、最小モデル、想定命題、意義、危険文献、棄却条件、Stage 0 検証事項。判定 RECOMMEND NEW PROJECT / REVISE AND RETEST (一つの明確な blocker) / RETURN TO BENCHMARK / STOP。有限の一回の探索で停止し、研究者本人が選択するまで正式 Stage 0 へ進まない。新規性・完全均衡・形式認証・理論凍結をこの段階で宣言しない。
