.First.lib <-
function(lib, pkg) {
# call Fortran source code
    library.dynam("nbpMatching", pkg, lib)
}
