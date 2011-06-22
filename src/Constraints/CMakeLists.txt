add_library(constraints 
  BoundConstraint.C CompoundConstraint.C Constraint.C LinearConstraint.C
  LinearEquation.C LinearInequality.C NonLinearConstraint.C NonLinearEquation.C
  NonLinearInequality.C
)
target_link_libraries(constraints utils)
