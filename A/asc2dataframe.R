
## we can then convert the matrices into ASCII and then into dataframe. The dataframe can then be named accordingly
## asc2dataframe
l_file<-list.files(path = "", pattern = "\\.asc$")
tdata=asc2dataframe(l_file)