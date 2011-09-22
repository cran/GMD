.First.lib <- function(libname, pkgname) {
  packageStartupMessage( sprintf("libname: %s",libname) )
  packageStartupMessage( sprintf("pkgname: %s",pkgname) )
  library.dynam(pkgname, pkgname, lib.loc=libname)
}
