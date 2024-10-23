rsat matrix-clustering \
-v 2 \
-matrix_file_table  /path/to/best_auc.tab \
-hclust_method average \
-calc sum \
-title "Selected motifs" \
-label_motif name -metric_build_tree Ncor \
-lth w 5 -lth cor 0.6 -lth Ncor 0.4 \
-label_in_tree name -quick -return json \
--cores 20 \
-o /path/to/motifs_only
