test_that("Rotation matrix generation on high dimensions", {
  u <- c(1, 2, 3, 4, 5, 6)
  Rmat <- vec2northpole(u)
  newu <- Rmat %*% u
  expect_equal(newu[1], sqrt(sum(u^2)))
  expect_equal(newu[-1], rep(0, length(u) - 1))

  u2 <- c(10, 100, 3, -10, 7, 0)
  expect_equal(sqrt(sum((Rmat %*% u2)^2)), sqrt(sum(u2^2)))

  u3 <- c(2, 0, 0, 0, 0, 0)
  expect_equal(sqrt(sum((Rmat %*% u3)^2)), 2)

  u4 <- c(6, 5, 4, 3, 2, 1)
  expect_equal(sqrt(sum((Rmat %*% u4)^2)), sqrt(sum(u4^2)))
})

test_that("Rotation matrix with inverse gets to any direction", {
  u <- c(1, 2, 3, 4, 5, 6)
  uend <- c(3, 4, 1, 2, 5, 6)
  Rmat1 <- vec2northpole(u)
  Rmat2 <- vec2northpole(uend)
  Rmat <- solve(Rmat2) %*% Rmat1

  expect_equal(uend, Rmat %*% u, ignore_attr = TRUE) #because size of u == size of uend
})

test_that("Handedness of coord system is preserved", {
  vector.cross <- function(a, b) { #code from https://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
    if(length(a)!=3 || length(b)!=3){
      stop("Cross product is only defined for 3D vectors.");
    }
    i1 <- c(2,3,1)
    i2 <- c(3,1,2)
    return (a[i1]*b[i2] - a[i2]*b[i1])
  }

  e1 <- c(1, 0, 0)
  e2 <- c(0, 1, 0)
  e3 <- c(0, 0, 1)

  stopifnot(isTRUE(all.equal(vector.cross(e1, e2), e3)))

  #the vector cross product gives direction in the right-hand rule
  #so a transformation R of e1, e2, e3 that preserves sizes, angles and handedness should be s.t.
  # Re1 x Re2 = Re3
  u <- c(1, 2, 3)
  Rmat <- vec2northpole(u)

  expect_equal(vector.cross(Rmat %*% e1, Rmat %*% e2), Rmat %*% e3,
               ignore_attr = TRUE)

  RmatB <- Directional::rotation(u/sqrt(sum(u^2)), c(1, 0, 0))

  RmatB %*% e1
  RmatB %*% e2
  RmatB %*% e3
  expect_equal(vector.cross(RmatB %*% e1, RmatB %*% e2), RmatB %*% e3,
               ignore_attr = TRUE)
})

test_that("rotation is the simplest for e1, e2, e3", {
  # rotating v1 to v2 can be acheived in two ways.
  # Suppose the angle between v1 and v2 is theta.
  # A rotation of theta in the axis perpendicular to v1 and v2
  # OR
  # rotation of pi around v2 axis, then rotation of 2*pi - theta in the axis perpendicular to v1 and v2
  # the former seems SIMPLER and better
  # for z to x that means y remains unchanged, x goes to -z
  e1 <- c(1, 0, 0)
  e2 <- c(0, 1, 0)
  e3 <- c(0, 0, 1)
  Rmat <- vec2northpole(e3)
  # Rmat <- Directional::rotation(e3, e1)
  expect_equal(Rmat %*% e3, e1, ignore_attr = TRUE)
  expect_equal(Rmat %*% e2, e2, ignore_attr = TRUE)
  expect_equal(Rmat %*% e1, -e3, ignore_attr = TRUE)
})

test_that("rotationmatrix(a,b) = Directional::rotation(b,a)", {
  skip_on_cran() #so Directional doesn't need to be loaded
  u <- c(1, 2, 3, 4, 5, 6)
  u <- u/sqrt(u%*%u)[[1]]
  b <- c(6,5,1,2,3,4)
  b <- b/sqrt(b%*%b)[[1]]


  Q1 <- rotationmatrix(u, b)
  Q2 <- Directional::rotation(b, u)
  expect_equal(Q1, Q2)
})

test_that("Rotation matrix rotates opposite direction to one", {
  u <- c(0, 2, 0, 0, 0, 0)
  uend <- c(0, -1, 0, 0, 0, 0)
  Rmat <- rotationmatrix(uend, u)

  expect_equal(uend * 2, Rmat %*% u, ignore_attr = TRUE) 
})
