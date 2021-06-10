MakeForest <- function(hypers, opts) {
  mf <- Module(module = "mod_forest", PACKAGE = "MultiskewBART")
  return(new(mf$Forest, hypers, opts))
}
