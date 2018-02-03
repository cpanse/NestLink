#R

.onAttach <- function(lib, pkg){
  if(interactive()){
    version <- packageVersion('NestLink')
    packageStartupMessage("Package 'NestLink' version ", version)
    
    invisible()
  }
}
