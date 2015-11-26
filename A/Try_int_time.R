##
A<-data.frame(start = factor(c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 336)), end = factor(c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)), blah = 1:12, blah2 = 13:24, blahx = 24:35, blah3 = rep(x = 1, times = 12), blah4 = rep(x = 2, times = 12), blah5 = rep(x = 3, times = 12), k = runif(n = 12, min = 1, max = 3))

int <- function(v){
  k <- v[9]
  for (j in k)
  {
    B <- v[6:8]
    t_b <- max(B[B<j])
    t_a <- min(B[B>=j])
    Ta_b <- which(B == t_b, arr.ind = T)[2]  #get column number within B
    Ta_a <- which(B == t_a, arr.ind = T)[2]
    C <- v[3:5]
    Tair_b <- C[,Ta_b]  #Air temp value before j  #something wrong here
    Tair_a <- C[,Ta_a]  #Air temp value after j
    Tair_int <- Tair_b + (Tair_a - Tair_b)*((j - t_b)/(t_a - t_b))
    return(data.frame(Tair1 = Tair_int, j_col = j))   
  }
}
x <- adply(.data = A, .margins = 1, .fun = int)    #adply faster
x1 <- ddply(.data = A, .variables = .(start, end), .fun = int)
system.time(adply(.data = A, .margins = 1, .fun = int))
system.time(ddply(.data = A, .variables = .(start, end), .fun = int))
