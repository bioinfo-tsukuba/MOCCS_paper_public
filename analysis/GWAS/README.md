# このディレクトリは
GWAS-SNPに対応するΔMOCCS2scoreを計算する時に用いたスクリプトを置く場所

- `apply_get_dMOCCS2score_GWAS_peaks_rand_phenotype_XX.sh`は、Rのスクリプトをsingularityで実行する時のシェルスクリプト

  遺伝研にて`qsub -cwd -o out_o/ -e out_e/ -l s_vmem=24G -l mem_req=24G apply_get_dMOCCS2score_GWAS_peaks_rand_phenotype_XX.sh`で実行した

- `job_XX.R`は、上記のシェルスクリプトが指定するRファイル

- `get_dMOCCS2score_GWAS_peaks_rand_phenotype_XX.R`は、上記のRファイルで読み込むΔMOCCS2score計算のための関数
