# 3dith

- score1: count total number entries whose value > cohort-dependent-threshold
- score2: inverse exponential of binned difference matrix -> compute PC1 -> distance between PC1 of each chromosome of each sample from normal PC1 (e.g., FIRE)
- score3: PC1 그래프에서 절댓값 y 적분 -> quantify how much the PC1 graph is perturbated 
  - PC1 그래프가 얼마나 요동치는지를 측정. 
- score4: froms score2, replace reference PC1 from normal PC1 to pcbc PC1
- score4_v2: pcbc ref를 어떤 버전 쓸지 parameterize (all / SC / nonSC)
- score5: harmonic mean of score2, score3 and score4
- score6: score2 * score3 / score4 
  - intuition:
      - score2: sample과 normal과의 거리
      - score4: sample과 pcbc (pure stem cell)과의 거리
      - normal과 거리가 멀 수록, PC1 graph의 perturbation이 심할 수록, stem cell 에 가까울 수록 cancer cell-like 할 라고 생각
