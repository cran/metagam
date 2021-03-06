library(mgcv)

set.seed(123)
# Generate some fits
ndat <- 3
n <- 100
fits <- lapply(1:ndat, function(x){
  dat <- gamSim(n = n, verbose = FALSE)
  b <- gam(y ~ s(x0, bs = "cr") + s(x1, bs = "cr"), data = dat)
  strip_rawdata(b)
})

# Meta-analyze
test_that("metagam works", {
  metafits <- metagam(fits, grid_size = 10)
  expect_s3_class(metafits, "metagam")
  expect_equal(metafits$terms, "s(x0)")
  # # Check that predictions are correct
  # # Vector to compare with has been generate by dput in a manual first run
  # expect_equal(round(metafits$meta_estimates$estimate, 8),
  #              c(-0.22194647, -0.13913644, -0.04638216, 0.04226372, 0.10112179,
  #                0.11856284, 0.09830259, 0.04701996, -0.03726664, -0.14627162))
  #
  # # Check that standard errors are correct
  # expect_equal(round(metafits$meta_estimates$se, 8),
  #              c(0.42231386, 0.32479295, 0.26294788, 0.23210114, 0.22074271,
  #                0.22037242, 0.23108217, 0.26216805, 0.32538673, 0.42494948))

  # Now with second term
  metafits <- metagam(fits, grid_size = 10, terms = "s(x1)")
  expect_s3_class(metafits, "metagam")
  expect_equal(metafits$terms, "s(x1)")

  # # Check that predictions are correct
  # # Vector to compare with has been generate by dput in a manual first run
  # expect_equal(round(metafits$meta_estimates$estimate, 8),
  #              c(-1.90254055, -1.54695893, -1.18290861, -0.74264595, -0.28591313,
  #                0.40713456, 0.76885555, 1.37551416, 2.08311411, 2.8655438))
  #
  # # Check that standard errors are correct
  # expect_equal(round(metafits$meta_estimates$se, 8),
  #              c(0.45711194, 0.34464661, 0.28070113, 0.25022086, 0.24067681,
  #                0.24936844, 0.27250813, 0.3185325, 0.3812979, 0.50396543))

  # Now with both terms
  grid <- data.frame(x0 = seq(0, 1, by = .1), x1 = seq(0, 1, by = .1))
  metafits <- metagam(fits, grid = grid, terms = c("s(x0)", "s(x1)"))
  expect_s3_class(metafits, "metagam")
  expect_equal(metafits$terms, c("s(x0)", "s(x1)"))

  # Check that predictions are correct
  # Vector to compare with has been generate by dput in a manual first run
  # expect_equal(round(metafits$meta_estimates$estimate, 8),
  #              c(-0.22238192, -1.90845939, -0.14829188, -1.54892801, -0.06563728,
  #                -1.31104563, 0.01759225, -0.82227955, 0.08236379, -0.537348,
  #                0.11504366, 0.07442868, 0.1144476, 0.59027953, 0.08619289, 0.93966274,
  #                0.03228429, 1.51758374, -0.04820967, 2.21008333, -0.14742953,
  #                2.88682888))
  #
  # # Check that standard errors are correct
  # expect_equal(round(metafits$meta_estimates$se, 8),
  #              c(0.42295604, 0.45839506, 0.33330303, 0.35355031, 0.27261384,
  #                0.29159583, 0.23866624, 0.25726103, 0.22347574, 0.24035071, 0.21934645,
  #                0.24135044, 0.2229785, 0.25485699, 0.23775134, 0.28412049, 0.27228278,
  #                0.32960658, 0.3344465, 0.38952799, 0.42611871, 0.51139217))

  # Use link function
  metafits <- metagam(fits, grid_size = 10, type = "link")
  expect_s3_class(metafits, "metagam")


  # Check that predictions are correct
  # Vector to compare with has been generate by dput in a manual first run
  # expect_equal(round(metafits$meta_estimates$estimate, 8),
  #              c(5.4670439, 5.46037915, 5.49255883, 5.53565579, 5.56065563,
  #                5.55408808, 5.51427511, 5.45682166, 5.4022208, 5.37047193, 5.77234916,
  #                5.83204971, 5.90167633, 5.96530988, 6.00213034, 6.00226181, 5.96604297,
  #                5.9079231, 5.84119815, 5.77202393, 6.10294759, 6.17831516, 6.25659141,
  #                6.322799, 6.35921531, 6.35965851, 6.32644533, 6.27139173, 6.20166903,
  #                6.1198474, 6.47675955, 6.55684504, 6.63719029, 6.70500065, 6.7440403,
  #                6.74761609, 6.7159831, 6.65709282, 6.57599408, 6.48020266, 6.81718625,
  #                6.92221445, 7.01954777, 7.09604684, 7.13919389, 7.14320217, 7.10459342,
  #                7.02466736, 6.90695279, 6.77011341, 7.75841574, 7.79355039, 7.83062972,
  #                7.86849486, 7.89089952, 7.88699169, 7.8562124, 7.81060755, 7.76039564,
  #                7.70799628, 7.99200066, 8.08583347, 8.17387255, 8.24139827, 8.27613891,
  #                8.27366052, 8.23374719, 8.16380394, 8.07067864, 7.96398062, 8.6347526,
  #                8.72592665, 8.80718192, 8.86597393, 8.89224555, 8.88268038, 8.83721235,
  #                8.76390866, 8.67098748, 8.56761759, 9.40907753, 9.48478576, 9.55602755,
  #                9.6132786, 9.64492716, 9.64456517, 9.61156831, 9.55484575, 9.48169577,
  #                9.39741113, 10.17843352, 10.21479917, 10.26564162, 10.31793464,
  #                10.3525207, 10.35735471, 10.32883243, 10.27644202, 10.21197469,
  #                10.1488514))

  # Check that standard errors are correct
  # expect_equal(round(metafits$meta_estimates$se, 8),
  #              c(0.62125251, 0.54734162, 0.50675527, 0.48851356, 0.48130605,
  #                0.47980643, 0.48460574, 0.5031818, 0.54790212, 0.62801304, 0.52279316,
  #                0.44423438, 0.39858858, 0.3766332, 0.36771751, 0.36638181, 0.37307217,
  #                0.39492924, 0.44316837, 0.5251978, 0.47845495, 0.39378464, 0.34320052,
  #                0.31845653, 0.30843931, 0.30700639, 0.3147795, 0.33973125, 0.3930118,
  #                0.4806462, 0.46770855, 0.37812971, 0.32271933, 0.29482284, 0.28315425,
  #                0.28105084, 0.28955724, 0.31763749, 0.37598225, 0.46848224, 0.47121043,
  #                0.37811471, 0.31893477, 0.28863955, 0.27610191, 0.27422384, 0.28415058,
  #                0.31507888, 0.37685648, 0.47184022, 0.46800644, 0.37795971, 0.32114599,
  #                0.2922408, 0.28074594, 0.27986466, 0.29070015, 0.3219654, 0.38333894,
  #                0.47790479, 0.47665762, 0.38992139, 0.33623719, 0.30915261, 0.29850066,
  #                0.29783363, 0.30787888, 0.33628698, 0.39257267, 0.48167837, 0.50397737,
  #                0.42347752, 0.37444496, 0.35000036, 0.34069043, 0.34056594, 0.35029184,
  #                0.37699855, 0.43017712, 0.51553861, 0.53747731, 0.46553916, 0.4243672,
  #                0.40528701, 0.3987934, 0.39934706, 0.40736265, 0.42890294, 0.47319647,
  #                0.5481428, 0.64774137, 0.58012839, 0.54134183, 0.52325335, 0.51702804,
  #                0.51756521, 0.5251356, 0.54501404, 0.58560428, 0.65469007))

})


