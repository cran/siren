congru = function(A,B){

  B=transpose(B)

  C = solve(sqrt((transpose(A)%*%A))) * (transpose(A) %*% B) * solve(sqrt(transpose(B)%*%B))

}
