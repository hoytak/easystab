test_ext<-function(){
  a<-c(1,2,3,4)
  a
}

test_internal<-function(size=5){
  if(!size) size=4
  a<-1:size
  a
}

test_call<- function(){
  val<-.Call("pretty_matrix", matrix(sample(1:99, 99), nrow=9), PACKAGE="clustering")
  val
}

test_more<- function(s=TRUE){
  val<-.Call("test_function_wrapper", 1:5, 2:6, as.integer(5), s, PACKAGE="clustering")
  val
}

print_matrix<-function(){
  s <- matrix(1:6, nrow = 2, byrow = TRUE)
  str(s)
  .Call("print_matrix", s, PACKAGE="clustering")
}
