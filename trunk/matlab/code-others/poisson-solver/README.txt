We provide a simple Geometric Multigrid Solver (MATLAB 7 m-files) for solving the Poisson Equation
Uxx+Uyy = -1 on a silhouette with the silhouette contours providing zero Dirichlet boundary conditions. 


To run the software use:
U = GMGmain(silhouetteMask);

To smooth the solution near the boundaries (for example before taking derivatives) use:
Urelaxed = relaxAfterSolvingU(U,4);

For more details see:

@ARTICLE{Gorelick&etal2006,
  AUTHOR =       "L. Gorelick and M. Galun and E. Sharon and A. Brandt and R. Basri",
  TITLE =        "Shape representation and classification using the Poisson Equation",
  JOURNAL =      "IEEE Transactions on Pattern Analysis and Machine Intelligence",
  YEAR =         "2006",
  volume =       "28",
  number =       "12",
  month =        "December" }.

