# Feature_selection_candidate
Model Candidate (GSFSJNE, BPSO &amp; IBPSO, Co-ABC, mRMR, LARS, CCM, Fast-OSFS, DART, Chi-Squared Feature selection, Lasso)

Model Candidates Description
모델 후보군에 코드와 설명

GSFSJNE
논문 요약 : https://productive-screen-f5c.notion.site/2021-07-28-Joint-neighborhood-entropy-based-gene-selection-method-with-fisher-score-for-tumor-class-f2cc088138cf46deb97a99f42ae3f863
•	GSFSJNE.R : main code
•	fisher_score.R : function code
•	GSFSJNE_fn.R : function code

BPSO & iBPSO
논문 요약 : https://productive-screen-f5c.notion.site/2021-07-01-Correlation-feature-selection-based-improved-Binary-Particle-Swarm-Optimization-for-gene-b439586c9d894b9ab2c405f9333f7d5e
•	BPSO.R : main code
•	iBPSO.R : function code
•	CFS_0702.R : function code
Co-ABC
논문 요약 : https://productive-screen-f5c.notion.site/2021-04-05-Co-ABC-Correlation-artificial-bee-colony-algorithm-for-biomarker-gene-discovery-using-g-56c3520a19a04c789521ec3043318f4b
•	Co_ABC.R : main code
•	CFS : function code (CFS_0702.R 구버전)

mRMR
논문 요약 : https://productive-screen-f5c.notion.site/2020-11-20-A-Machine-Learning-Approach-for-Identifying-Gene-Biomarkers-Guiding-the-Treatment-of-Bre-a03f8bebc8bb41d19b5818b255f4639b

LARS
모델 요약 : https://productive-screen-f5c.notion.site/2021-03-03-High-dimensional-feature-selection-for-genomic-datasets-878eec52b7ae4db6a3ded4779196613f
다른 모델로 1차 feature selection 진행 후 사용 권장

CCM
모델 요약 : https://productive-screen-f5c.notion.site/2021-03-03-High-dimensional-feature-selection-for-genomic-datasets-878eec52b7ae4db6a3ded4779196613f
https://www.notion.so/2021-03-03-High-dimensional-feature-selection-for-genomic-datasets-878eec52b7ae4db6a3ded4779196613f#f1addb2164de4a58866184935d019375
 !git clone https://github.com/Jianbo-Lab/CCM
 !cd CCM
 !python /content/CCM/core/run_synthetic.py

Fast OSFS
모델 요약 : https://productive-screen-f5c.notion.site/2021-03-03-High-dimensional-feature-selection-for-genomic-datasets-878eec52b7ae4db6a3ded4779196613f
•	data_analysis_0304.R : main code
•	partial_corr_coef.R : function
•	my_cond_indep_fisher_z.R : function
•	computer_dep_2.R : function
•	optimal_compter_dep_2.R : function
•	fast_osfs_z.R : function

DART
논문 요약: https://productive-screen-f5c.notion.site/2021-03-03-High-dimensional-feature-selection-for-genomic-datasets-878eec52b7ae4db6a3ded4779196613f
Information Gain
모델 요약: https://productive-screen-f5c.notion.site/2021-01-29-Validation-of-miRNAs-as-Breast-Cancer-Biomarkers-with-a-Machine-Learning-Approach-6c921313215842b8ad5ba309afc44c01
Chi-Squared Feature Selection
모델 요약: https://productive-screen-f5c.notion.site/2021-01-29-Validation-of-miRNAs-as-Breast-Cancer-Biomarkers-with-a-Machine-Learning-Approach-6c921313215842b8ad5ba309afc44c01

Lasso 회귀분석
