#INPUT: Q1, Q2, Q3, appropriate quadratic forms in x,y,z with rational coefficients, and a list of homogeneous degree 1 forms in x,y,z with rational coefficients
#OUTPUT: The 5 x 4 matrices formed by points of Deltatilde over the intersection
#of Delta := (Q1*Q3 - Q2^2=0) and each line in the list given (the code will ignore cases where the intersection contains real points).
#Checks whether the support of the points on Deltatilde spans a 3-plane in P^4.

#We note that in this (and all accompanying Sage code), we use the variables x,y,z (whereas in the paper we use the variables u,v,w).
import numpy as np

F=QQ

# Bound on Similar Roots
bound = 0.001;

T.<x,y,z>=PolynomialRing(F)

#-----------------------------------------
# Helper Functions:
#-----------------------------------------
# Given two expressions (x == x_1, y == x_2), (x == y_1, y == y_2), computes the
# distance between the points they represent
def difference(r, s):
    x1 = r[0].rhs();
    x2 = r[1].rhs();
    
    y1 = s[0].rhs();
    y2 = s[1].rhs();
    
    diff = (x1 - y1).abs()^2 + (x2 - y2).abs()^2
    
    return sqrt(diff);

# Given an expression and a list of expressions, checks if the expressions
# is close to any item of the list
def is_similar_in_list(elt, lst):
    for item in lst:
        if difference(elt, item) < bound:
            return True
    return False

# Due to Sage's implementations, using solve may produce multiple solutions that are close to
# to each other but are actually the same solution, this function is here to prune out these duplicate roots.
def prune_similar(lst):
    unique = [];
    
    for elt in lst:
        if not is_similar_in_list(elt, unique):
            unique.append(elt);
    return unique;

#check whether the discriminant curve Delta is smooth
def is_delta_smooth(f):
    Grad=ideal(f,diff(f,x),diff(f,y),diff(f,z))
    if Grad.dimension()==0:
        print("Delta is smooth.")
        return True
    else:
        print("Delta is not smooth.")
        return False

#check whether Deltatilde is non-singular
def is_double_cover_smooth(q1, q2, q3):
    R.<x,y,z,r,s> = ProjectiveSpace(F, 4);

    q1_t = q1 - r^2;
    q2_t = q2 - r*s;
    q3_t = q3 - s^2;
    
    delta_tilde = R.curve([q1_t, q2_t, q3_t]);

    if delta_tilde.is_singular():
        print("Deltatilde is not smooth.");
        return False;
    else:
        print("Deltatilde is also smooth.");
        return True;

# Given a list of expressions of the form (x == r, y == s), finds the real expressions in them
# We consider a value to be real if its imaginary component is less than a certain threshold
def find_real_roots(lst):
    real_roots = []
    for root in lst:
        c1 = (root[0] + I - I).rhs();
        c2 = (root[1] + I - I).rhs();
        
        if abs(c1.imag()) < 0.0001 and abs(c2.imag()) < 0.0001:
            real_roots.append(root)
    return real_roots
    
# Takes in the standard 5x4 Gal(C/R)-matrix (Double Array) and outputs the determinant of all 5 of its
# 4 x 4 submatrices
def process_matrix(M):
    gal_matrix = np.matrix(M, dtype=complex);
    gal_rk = np.linalg.matrix_rank(gal_matrix);
    # print("Rank of Matrix: " + str(gal_rk));
    
    # Note that M is full-rank if and only if at least one of its 4 x 4 submatrix has
    # a non-zero determinant
    for i in range(5):
        new_mat = np.delete(gal_matrix, i, 1);
        new_det = np.linalg.det(new_mat);
        print("Determinant with Column " + str(i) + " Removed: " + str(new_det));
        # print("Close to Zero: " + str(abs(new_det) < bound));
        # print("--------");

# Given Q1, Q2, Q3, and a line, and a z-value
# Produces the 5x4 matrix with the 4 points (typically) on Deltatilde given
# by the 4 intersections (typically) of Delta and the line.
# The matrix is checked to be full-rank or not
# Returns True if processed successfully, False otherwise
def galois_matrix(q1, q2, q3, line, z_val):
    x, y = var('x, y')
    Q1 = q1 + x - x;
    Q2 = q2 + x - x;
    Q3 = q3 + x - x;
    T = line + x - x;
    
    Delta = Q1*Q3 - Q2^2;

    print("Intersection of Delta and Line on z = 1");

    # Remove duplicate roots resulted from imprecision of Sage's solve function
    roots = solve([Delta == 0, T == 0], x, y)
    roots = prune_similar(roots);
    print(roots);
    print("-------------------------");
    
    if len(roots) < 4:
        print("The line does not have at least 4 intersections with Delta!");
        return False;
    
    real_roots = find_real_roots(roots);
    if len(real_roots) != 0:
        print("This line intersects with Delta on real points:");
        print(real_roots);
        return False;

    gal_matrix = [];
    
    # For each root in the list, construct the appropriate row in the matrix and add it in
    for r in roots:
        c1 = (r[0] + I - I).rhs();
        c2 = (r[1] + I - I).rhs();
        
        if abs(c1.imag()) < 0.0001:
            c1 = c1.real()
        
        if abs(c2.imag()) < 0.0001:
            c2 = c2.real()
        
        # precision
        k = 50;
        c1 = c1.n(k) # This is x
        c2 = c2.n(k) # This is y
        c3 = sqrt(Q1.subs(x == c1, y == c2).n(k)) # This is r
        c4 = sqrt(Q3.subs(x == c1, y == c2).n(k)) # This is s
        
        q2_val = Q2.subs(x == c1, y == c2);
        
        print("Q2 : " + str(q2_val));
        
        # Adjusting r and s based on Q2(x, y, z)
        if abs(q2_val - c3*c4) < abs(q2_val - -1*c3*c4):
            row = [c1, c2, z_val, c3, c4]
        else:
            row = [c1, c2, z_val, -c3, c4]

        gal_matrix.append(row);

    gal_mat = Matrix(gal_matrix);
    print("Matrix: \n" + str(gal_mat));
    process_matrix(gal_matrix);
    
    return True;

#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    # inputs:
    Q1 = -2*x^2 - 2*x*y + 4*x*z - 2*y^2 + 6*y*z - 5*z^2;
    Q2 = 10*x*y - 20*x*z + 5*y^2 - 20*y*z + 20*z^2;
    Q3 = -48*x^2 - 48*x*y + 96*x*z - 20*y^2 + 88*y*z - 92*z^2;
    lines = [y];
    
    Delta = Q1*Q3 - Q2^2;
    z_val = 1;
    
    print("Q1 := " + str(Q1) + ";")
    print("Q2 := " + str(Q2) + ";")
    print("Q3 := " + str(Q3) + ";")
    print("Delta: " + str(Delta))
    print("------------------------------------------")
    
    if is_delta_smooth(Delta) and is_double_cover_smooth(Q1, Q2, Q3):
        i = 1;
        for line in lines:
            print("Iteration " + str(i));
            print("Line: " + str(line));
            
            gq1 = Q1(x, y, z_val)
            gq2 = Q2(x, y, z_val)
            gq3 = Q3(x, y, z_val)
            gline = line(x, y, z_val)
            
            galois_matrix(gq1, gq2, gq3, gline, z_val);
            i = i + 1;
            print("------------------------------------------")
    

        
if __name__ == '__main__':
    main()
