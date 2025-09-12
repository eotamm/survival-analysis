# Use the Python env where TabPFN is installed
reticulate::use_condaenv("r-tabpfn", required = TRUE)

# Verify TabPFN import
if (!reticulate::py_module_available("tabpfn")) {
  stop("Python module 'tabpfn' not found in conda env 'r-tabpfn'. Install there:\n",
       "  conda activate r-tabpfn\n",
       "  pip install \"git+https://github.com/PriorLabs/TabPFN.git\"")
}
reticulate::py_run_string("
import tabpfn, importlib.metadata as md
print('tabpfn file   :', getattr(tabpfn, '__file__', '<unknown>'))
print('tabpfn module :', getattr(tabpfn, '__version__', '<no __version__ attr>'))
try:
    print('pkg version   :', md.version('tabpfn'))
except Exception as e:
    print('pkg version   : <unknown>', e)
")

# Python helper to build classifier across versions and ignore the 500-feature limit
reticulate::py_run_string("
import importlib
def _make_tabpfn_clf(device='cpu', ensemble=16, seed=1):
    tabpfn = importlib.import_module('tabpfn')
    clf = tabpfn.TabPFNClassifier(device=device, random_state=int(seed))
    if hasattr(clf, 'set_ensemble_size'):
        clf.set_ensemble_size(int(ensemble))
    elif hasattr(clf, 'N_ensemble_configurations'):
        clf.N_ensemble_configurations = int(ensemble)
    elif hasattr(clf, 'n_ensemble_configurations'):
        clf.n_ensemble_configurations = int(ensemble)
    if hasattr(clf, 'ignore_pretraining_limits'):
        clf.ignore_pretraining_limits = True
    return clf
")


# TabPFN wrappers
one_hot_mm <- function(df) {
  X <- stats::model.matrix(~ . - 1, data = df)
  storage.mode(X) <- "double"
  X
}

tabpfn_classify <- function(X_train_df, y_train, X_test_df,
                            device = "cpu", ensemble = 32, seed = 1) {
  Xtr <- one_hot_mm(X_train_df)
  Xte <- one_hot_mm(X_test_df)
  miss <- setdiff(colnames(Xtr), colnames(Xte))
  if (length(miss))
    Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
  Xte <- Xte[, colnames(Xtr), drop = FALSE]
  clf <- reticulate::py$`_make_tabpfn_clf`(
    device = device, ensemble = as.integer(ensemble), seed = as.integer(seed)
  )
  clf$fit(Xtr, as.integer(y_train))
  p <- clf$predict_proba(Xte)
  reticulate::py_to_r(p)
}

tabpfn_bin_cindex_cv5 <- function(
    df, method_name,
    seed = 1, K = 5, fold_id = NULL,
    device = "cpu", ensemble = 32
) {
  message(sprintf("[TabPFN Bin→C-index] %s: %d-fold CV", method_name, K))
  risk_fun <- function(train, test) {
    xtr <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    ytr <- as.integer(train$Event)
    pmat <- tabpfn_classify(xtr, ytr, xte, device = device, ensemble = ensemble, seed = seed)
    if (is.matrix(pmat) && ncol(pmat) >= 2) as.numeric(pmat[, 2]) else as.numeric(pmat)
  }
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  summarize_cv(method_name, sprintf("TabPFN_Binomial_%dCV_toC", K), res)
}
