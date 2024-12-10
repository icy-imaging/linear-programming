This library is designed to solve linear program for the Simplex algorithm.
By default the primal problem is solved using Gauss-Jordan elimination and the lexicographic order pivoting rule.

By definition, the optimization problem is of the form:
min_x c'*x or max_x c'*x
A*x <= b
with x >=0
where x, b, c are column vectors and c' is the transpose of c. * stands for the matrix multiplication.
Here each constraint A[i]*x can be either "equal" or "lower than or equal" constraints
based on the boolean value in the equalityConstraints vector:
A[i]*x == b if equalityConstraints[i]
A[i]*x <= b if equalityConstraints[i]

The library is part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org but can also be used as a stand-alone library.
The library can run user-defined problems defined in text files. Run java -jar linearProgrammingICY.jar -help
for the details about use and input options.
Most useful function have public access for a call within a Java program.

Author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
Version 1.0
Date 2014-04-21
License gpl v3.0

NB: the tooltip image is from wikimedia commons: http://commons.wikimedia.org

Content of help for stand-alone call:
Empty arguments to use the default example.
Enter an integer from 0 to 3 for different example scenarios.
Enter -1 for a user-defined scenario through text files.

For custom scenarios, enter filenames (text files) for different parameters of the problem preceded by the appropriate prefix:
-c for the objective function file.
-A for constraint matrix file.
-b for the constraint value file.
-e for the equality constraint file.
-max or -min to indicate a maximization or minimization problem, respectively. Default in minimization.
Example arguments: java -jar linearProgrammingICY.jar -c c.txt -A A.txt -b b.txt -e eq.txt -max -o solution.txt
Each text file must contain a series of double values separated by ',' in a single line, except for the constraint file which contains one line per constraint.
For the equality file '0' stands for 'false' and '1' for true.
