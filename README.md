# 3dith
## Scores
- score1: count total number entries whose value > cohort-dependent-threshold
- score2: inverse exponential of binned difference matrix -> compute PC1 -> distance between PC1 of each chromosome of each sample from normal PC1 (e.g., FIRE)
- score3: PC1 그래프에서 절댓값 y 적분 -> quantify how much the PC1 graph is perturbated 
  - PC1 그래프가 얼마나 fluctuation이 심한지를 측정. 
- score4: froms score2, replace reference PC1 from normal PC1 to pcbc PC1
- score4_v2: pcbc ref를 어떤 버전 쓸지 parameterize (all / SC / nonSC)
- score5: harmonic mean of score2, score3 and score4
- score6: score2 * score3 / score4 
  - intuition:
      - score2: sample과 normal과의 거리
      - score4: sample과 pcbc (pure stem cell)과의 거리
      - normal과 거리가 멀 수록, PC1 graph의 perturbation이 심할 수록, stem cell 에 가까울 수록 cancer cell-like 할 라고 생각
- score7:
  - cosine(theta)
    - theta: 각 샘플을 (x, y) = (score2, score4) 로 cartesian 좌표평면에 나타냈을 때, 원점으로부터 각 sample의 line과 x축간의 각도
    - theta 값 구하기 전에 score2와 score4는 각 cohort별로 minmax scaling해서 [0,1]로 맞춰줌.
- stem_closeness:
  - score7 계산 과정에서, 아래를 변경
    - reference 값 만들 때, all available sample들 중 절반을 randomly pick해서 사용. (info leak 방지)
    - sklearn의 minmax_scaler 쓰지 말고, 각 cohort의 각 score 별 min & max 값 저장한 table 만든 후, 이 table에 있는 값을 사용
      - 새로운 test sample이 input되었다고 할때, 그 sample 자체의 값은 minmax_scaling 과정에서 제외되어야 하므로.
## Benchmark
- https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0741-y
