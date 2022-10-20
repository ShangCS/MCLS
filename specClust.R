
specClust <- function (data, centers = NULL, nn = 7, method = "symmetric", 
                       gmax = NULL, max.iter = 10000, ...) 
{
  call = match.call()
  if (is.data.frame(data)) 
    data = as.matrix(data)
  da = apply(data, 1, paste, collapse = "#")
  indUnique = which(!duplicated(da))
  indAll = match(da, da[indUnique])
  data2 = data
  data = data[indUnique, ]
  n <- nrow(data)
  data = scale(data, FALSE, TRUE)
  if (is.null(gmax)) {
    if (!is.null(centers)) 
      gmax = centers - 1L
    else gmax = 1L
  }
  test = TRUE
  while (test) {
    DC = mydist(data, nn)
    sif <- rbind(1:n, as.vector(DC[[2]]))
    g <- graph(sif, directed = FALSE)
    g <- decompose(g, min.vertices = 4)
    if (length(g) > 1) {
      if (length(g) >= gmax) 
        nn = nn + 2
      else test = FALSE
    }
    else test = FALSE
  }
  W <- DC[[1]]
  n <- nrow(data)
  wi <- W[, nn]
  SC <- matrix(1, nrow(W), nn)
  SC[] <- wi[DC[[2]]] * wi
  W = W^2/SC
  alpha = 1/(2 * (nn + 1))
  qua = abs(qnorm(alpha))
  W = W * qua
  W = dnorm(W, sd = 1)
  DC[[1]] = W
  L = Laplacian(DC, nn, method)
  f <- function(x, extra) as.vector(extra %*% x)
  if (is.null(centers)) 
    kmax = 25
  else kmax = max(centers)
  ###
  #add the maxiter parameter to the arpack call, below
  ###    
  U <- arpack(f, extra = L, options = list(n = n, which = "SM", 
                                           nev = kmax, ncv = 2 * kmax, mode = 1, maxiter=max.iter), sym = TRUE)
  ind <- order(U[[1]])
  
  U[[2]] = U[[2]][indAll, ind]
  U[[1]] = U[[1]][ind]
  if (is.null(centers)) {
    tmp = which.max(diff(U[[1]])) + 1
    centers = which.min(AUC(U[[1]][1:tmp]))
  }
  if (method == "symmetric") {
    rs = sqrt(rowSums(U[[2]]^2))
    U[[2]] = U[[2]]/rs
  }
  result = kmeans(U[[2]], centers = centers, nstart = 20, ...)
  archeType = getClosest(U[[2]][indAll, ], result$centers)
  result$eigenvalue = U[[1]]
  result$eigenvector = U[[2]]
  result$data = data2
  result$indAll = indAll
  result$indUnique = indUnique
  result$L = L
  result$archetype = archeType
  result$call = call
  class(result) = c("specClust", "kmeans")
  result
}