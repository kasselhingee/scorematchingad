tapeJacobian <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- pTapeJacobian(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("d", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

tapeHessian <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- pTapeHessian(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("d^2", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

tapeGradOffset <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- pTapeGradOffset(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("doffset:", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

tapeLogJacDet <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- ptapelogdetJ(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("logJdet:", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

tapeSwap <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- swapDynamic(tape$ptr, tape$dyntape, tape$xtape)
  ADFun$new(outptr, 
            name = paste0("d", tape$name), 
            xtape = tape$dyntape, 
            dyntape = tape$xtape, 
            usertheta = rep(NA_real_, length(tape$xtape)))
}

