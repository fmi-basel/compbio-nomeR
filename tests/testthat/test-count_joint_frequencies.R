test_that("count_joint_frequencies works", {
  
  ## check matrix input
  ## create dummy random data
  set.seed(3346)
  nc <- 50
  nr <- 50
  rmatr <- matrix(data = as.numeric(rnorm(nc * nr) >= 0.5),
                  ncol = nc,nrow=nr)
  rctbl <- count_joint_frequencies(data=rmatr,max_spacing = nc)
  expect_s3_class(rctbl,"data.frame")
  
  ## check list input
  rlist <- sapply(1:nr,function(rw){
    as.numeric(rnorm(nc) >= 0.5)
  },simplify = F,USE.NAMES = T)
  rlistctbl <- count_joint_frequencies(data=rlist,max_spacing = nc)
  expect_s3_class(rlistctbl,"data.frame")
  
  
})
