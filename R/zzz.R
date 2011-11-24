.First.lib <- function(libname, pkgname) {
  packageStartupMessage( sprintf("libname: %s",libname) )
  packageStartupMessage( sprintf("pkgname: %s",pkgname) )
  packageStartupMessage( sprintf("Version: 0.3") )
  packageStartupMessage( sprintf("Packaged date: 2011-11-24") )
  library.dynam(pkgname, pkgname, lib.loc=libname)
}
