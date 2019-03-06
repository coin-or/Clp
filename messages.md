Messages
========

Some of the more common messages and codes passed by CLP are listed in
the tables below. This is list is not meant to exhaustive. The notation
is as for printf from \"C\":

-   `%s` is a string

-   `%d` is an integer

-   `%g` or `%f` is a floating point value

  Code   Area                                                                                                                                                                                                        Text and notes
  ------ -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- -- -----------------------------------------------------------------------
  1      MPSREAD                                                                                                                                                                                                     `At line %d %s`
         This just prints out NAME line, ROW line, etc                                                                                                                                                               
  2      MPSREAD                                                                                                                                                                                                     `Problem %s has %d rows, %d columns and %d elements
                                                                                                                                                                                                                           `
         This gives statistics after reading an MPS file                                                                                                                                                             
  8      MPSREAD                                                                                                                                                                                                     `%s read with %d errors
                                                                                                                                                                                                                           `
         This gives error statistics for file                                                                                                                                                                        
  505    PRESOLVE                                                                                                                                                                                                    `
                                                                                                                                                                                                                           Presolved poblem not optimal, resolve after postsolve
                                                                                                                                                                                                                           `
         This could be because it was not feasible or because of maximum iterations. If this message occurs then consider using primal clean up                                                                      
  506    PRESOLVE                                                                                                                                                                                                    `
                                                                                                                                                                                                                           Presolve %d (%d) rows, %d (%d) columns and %d (%d) elements
                                                                                                                                                                                                                           `
         The first number is the number after presolve and the number in parentheses is amount of reduction                                                                                                          
  510    PRESOLVE                                                                                                                                                                                                    `
                                                                                                                                                                                                                           Presolve is modifying %d integer bounds and re-presolving
                                                                                                                                                                                                                           `
         If presolve determines at the end that an integer variable have its bounds changed then it will repeat the entrire presolve                                                                                 
  511    PRESOLVE                                                                                                                                                                                                    `
                                                                                                                                                                                                                           After Postsolve, objective %g, infeasibilities - dual %g (%d), 
                                                                                                                                                                                                                           primal %g (%d)
                                                                                                                                                                                                                           `
         This gives the state after postsolve - this gives the objective value and the sum of dual and primal infeasibilities with the number of infeasibilities in parentheses. Hopefully these should be zero      
  512    PRESOLVE                                                                                                                                                                                                    `
                                                                                                                                                                                                                           Presolved model was optimal, full model needs cleaning up
                                                                                                                                                                                                                           `
         If the numbers in previous message (511) were large then maybe we need to know, if small then that\'s life                                                                                                  

  : COIN Messages passed at or above logging level 1

  Code   Area                                                                                                                                                                                                                                                             Text and notes
  ------ ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- -- ----------------------------------------------------------------
  1      SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Primal infeasible - objective value %g
                                                                                                                                                                                                                                                                                `
         You may need to look at previous messages or use methods. Such as sumPrimalInfeasibilities() to find cause                                                                                                                                                       
  2      SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Dual infeasible - objective value %g
                                                                                                                                                                                                                                                                                `
         You may need to look at previous messages or use methods. Such as sumDualInfeasibilities() to find cause                                                                                                                                                         
  3      SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Stopped - objective value %g
                                                                                                                                                                                                                                                                                `
         The algorithm stopped as requested by the user.                                                                                                                                                                                                                  
  4      SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Stopped due to errors - objective value %g
                                                                                                                                                                                                                                                                                `
         Switch on log level 2 to see information on size of elements etc. If they look reasonable then maybe we need to know.                                                                                                                                            
  5      SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                %d Obj %g Primal inf %g (%d) Dual inf %g (%d)
                                                                                                                                                                                                                                                                                `
         At each re-factorization this gives the number of iterations and the value of the objective function. If there are primal infeasibilities then the sum and number are given and similarly for dual infeasibilities. (This is a simplified form of message.)      
  14     SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Perturbing problem by %g % of %g
                                                                                                                                                                                                                                                                                `
         There is more to this message but if the user sees this then s/he has chosen to perturb the problem or the algorithm has decided to do so. If the numbers look too large the user may wish to think again.                                                       
  19     SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                %d variables/rows fixed as scaled bounds too close
                                                                                                                                                                                                                                                                                `
         If this occurs look carefully at your input data                                                                                                                                                                                                                 
  24     SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Matrix will be packed to eliminate small elements
                                                                                                                                                                                                                                                                                `
         If this occurs the user should look carefully at data.                                                                                                                                                                                                           
  26     SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Matrix will be packed to eliminate %d duplicate elements
                                                                                                                                                                                                                                                                                `
         If this occurs the user should look carefully at data.                                                                                                                                                                                                           
  28     SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                Crash put %d variables in basis, %d dual infeasibilities
                                                                                                                                                                                                                                                                                `
                                                                                                                                                                                                                                                                          
  29     SIMPLEX                                                                                                                                                                                                                                                          `
                                                                                                                                                                                                                                                                                End of values pass after %d iterations
                                                                                                                                                                                                                                                                                `
         ??? If primal(1) or dual(1) the a sweep through model is made and this signals end of pass.                                                                                                                                                                      

  : CLP Messages passed at or above logging level 1

  Code   Area                                                                                                                                Text and notes
  ------ -------------------------------------------------------------------------------------------------------------------------------- -- -----------------------------------------------------------------------------------
  3001   MPSREAD                                                                                                                             `
                                                                                                                                                   Illegal value for %s of %g
                                                                                                                                                   `
         String will be \"infinity\" if setInfinity passed bad value, or \"default integer bound\" if setDefaultBound passed bad value.      
  3002   MPSREAD                                                                                                                             `
                                                                                                                                                   Bad image at line %d < %s >
                                                                                                                                                   `
         This gives line number and the offending line                                                                                       
  3003   MPSREAD                                                                                                                             `
                                                                                                                                                   Duplicate objective at line %d < %s >
                                                                                                                                                   `
         An objective row appears twice in one column                                                                                        
  3004   MPSREAD                                                                                                                             `
                                                                                                                                                   Duplicate row %s at line %d %s
                                                                                                                                                   `
         The named row appears twice in one column.                                                                                          
  3005   MPSREAD                                                                                                                             `
                                                                                                                                                   No match for row %s at line %d < %s >
                                                                                                                                                   `
         The named row did not appear in ROWS section.                                                                                       
  3006   MPSREAD                                                                                                                             `
                                                                                                                                                   No match for column at line %d < %s >
                                                                                                                                                   `
         The named column (in BOUNDS section) did not appear in COLUMNS section.                                                             
  6001   MPSREAD                                                                                                                             `
                                                                                                                                                   Unable to open mps input file %s
                                                                                                                                                   `
                                                                                                                                             
  6002   MPSREAD                                                                                                                             `
                                                                                                                                                   Unknown image %s at line %d of file %s
                                                                                                                                                   `
         The Mps reader could not make sense of the image file specified.                                                                    
  6003   MPSREAD                                                                                                                             `
                                                                                                                                                   Consider the possibility of a compressed file which zlib is unable to read.
                                                                                                                                                   `
         Some .gz files can not be read by zlib. Using gunzip and then gzip normally cures problem.                                          
  6004   MPSREAD                                                                                                                             `
                                                                                                                                                   EOF on file %s
                                                                                                                                                   `
         The Mps reader did not find expected section marker.                                                                                
  6005   MPSREAD                                                                                                                             `
                                                                                                                                                   Returning as too many errors
                                                                                                                                                   `
         The reader has put out 100 messages and is giving up.                                                                               
  507    PRESOLVE                                                                                                                            `
                                                                                                                                                   Presolve determined that the problem is infeasible with tolerance of %g
                                                                                                                                                   `
         If you want you can try with a larger tolerance                                                                                     
  508    PRESOLVE                                                                                                                            `
                                                                                                                                                   Presolve thinks problem is unbounded
                                                                                                                                                   `
         Perhaps the user should maximize if initially minimizing or vice versa.                                                             
  509    PRESOLVE                                                                                                                            `
                                                                                                                                                   Presolve thinks problem is infeasible AND unbounded???
                                                                                                                                                   `
         If you get this message we want to know                                                                                             

  : COIN Messages passed at or above logging level 0

  Code   Area                                                                                                   Text and notes
  ------ --------------------------------------------------------------------------------------------------- -- -----------------------------------------------------------------------
  3002   SIMPLEX                                                                                                `
                                                                                                                      Not solving empty problem - %d rows, %d columns and %d elements
                                                                                                                      `
         Test problem size before solving.                                                                      
  6002   SIMPLEX                                                                                                `
                                                                                                                      %d bad bound pairs or bad objectives were found
                                                                                                                      `
         Either the value in the objective was too large or a lower bound was greater than an upper bound.      
  6003   SIMPLEX                                                                                                `
                                                                                                                      Matrix has %d large values, first at column %d, row %d is %g
                                                                                                                      `
         Some of the values in matrix are ridiculous.                                                           
  6004   SIMPLEX                                                                                                `
                                                                                                                      Can't get out of loop ...
                                                                                                                      `
                                                                                                                

  : CLP Messages passed at or above logging level 0

There are also messages available at log level 2 (the most likely useful
relate to scaling), and will be addressed in a future version of this
User Guide.
