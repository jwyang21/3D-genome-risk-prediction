cpg_type=opensea
score_f='/data/project/3dith/data/cohort-1-best-score-km.csv'

echo "python3 4_pcc-sc-si.py --cpg_type $cpg_type --si_fname /data/project/3dith/data/stemness-index/1-s2.0-S0092867418303581-mmc1.xlsx --score_fname $score_f"
python3 4_pcc-sc-si.py --cpg_type $cpg_type --si_fname /data/project/3dith/data/stemness-index/1-s2.0-S0092867418303581-mmc1.xlsx --score_fname $score_f