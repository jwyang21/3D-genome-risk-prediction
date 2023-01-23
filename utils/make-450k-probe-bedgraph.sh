for cpg_type in opensea island shelf_shore
do
        echo "python3 make-450k-probe-bedgraph.py --cpg_type $cpg_type"
        python3 make-450k-probe-bedgraph.py --cpg_type $cpg_type
done
