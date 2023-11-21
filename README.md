# Subcortical_Heritability
Code of article "Heritability of functional gradients in the human subcortex"

Author list: Xinyu Wu, Yu Zhang, Mufan Xue, Jinlong Li, Xuesong Li, Zaixu Cui, Jia-Hong Gao, Guoyuan Yang

Contact: bit_arims_xinyuwu@bit.edu.cn/x.w.wdveloce@gmail.com (Xinyu Wu), yanggy@bit.edu.cn (Guoyuan Yang)

## Dependencies

AFNI,<br/>
BrainSpace: https://github.com/MICA-MNI/brainspace,<br/>
APACE: https://github.com/NISOx-BDI/APACE,<br/>
SUIT: https://github.com/jdiedrichsen/suit,<br/>
LIBSVM: https://github.com/cjlin1/libsvm,<br/>
plot_fig_subcortex: https://github.com/wd-veloce/plot_fig_subcortex,<br/>

## Code descriptions

```
└─ src
    │   main_grad.m: Code involved computing functional gradients.
    │   main_svm.m: Code involved computing SVM classifier.
    │   main_t1wt2w.m: Code involved computing T1w/T2w ratio.
    │
    ├─ functions
    │  ├─ compute_similarity
    │  │       compute_similarity_cifti.m: Constructing individual-level FC matrix.
    │  │       preprocess_data_cifti.m: Preprocessing individual fMRI data.
    │  │
    │  ├─ heritability
    │  │       APACE.m: Computing heritability.
    │  │       APACE_bootstrap.m: Permutation of computing heritability, using Bootstrap strategy.
    │  │       APACE_jackknife.m: Permutation of computing heritability, using Jack-knife strategy.
    │  │
    │  ├─ roc
    │  │       AUC_TP_FP_Calculate.m: Computing TP, FP and AUC in SVM classification.
    │  │
    │  └─ tools
    │          fdr.m: FDR correction of the heritability.
    │          load_file_spmd.m: Loading files if using spmd to parallel.
    │          procrustes_alignment.m: Derived from BrainSpace Toolbox, to align individual-level functional gradients.
    │          readlines.m: Loading subject lists.
    │          save_file_spmd.m: Saving files if using spmd to parallel.
    │          show_progress.m: Showing progress in for-loop.
    │
 ```
