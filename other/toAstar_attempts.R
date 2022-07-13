set.seed(1234)
Qin <- orthogmatwith111vec()

Astar <- Qin %*% diag(c(1, 2, 0)) %*% t(Qin)
u <- c(0.1, 0.2, 0.7)

eigenspace <- eigen(Astar, symmetric = TRUE)
Q <- eigenspace$vectors
Q %*% diag(eigenspace$values) %*% t(Q)

# now try to decompose Astar
v1n <- Q[1:2, 1] - Q[3, 1]
v2n <- Q[1:2, 2] - Q[3, 2]
AL <- cbind(v1n, v2n)
t(u[1:2]) %*% AL %*% diag(c(1, 2)) %*% t(AL) %*% u[1:2] +
  1 * Q[3, 1]^2 + 2 * Q[3, 2]^2 +
  (2 * 1 * Q[3,1] * v1n + 2 * 2 *Q[3,2] * v2n) %*% u[1:2]

u %*% Astar %*% u
