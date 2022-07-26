#-----------------------------------------
# The following code is adapted from the Magma code accompanying ''Curve classes on conic bundle threefolds and applications to rationality'' by Frei, Ji, Sankar, Viray, and Vogt.
#-----------------------------------------
#INPUT: Q1, Q2, Q3, appropriate plane conics with Rational Coefficients
#OUTPUT: The real Weierstrass points of the associated hyperelliptic curve Gamma and the signatures in the interval
#We note that in this (and all accompanying Sage code), we use the variables x,y,z (whereas in the paper we use the variables u,v,w).

import warnings
warnings.filterwarnings("ignore")

F=QQ

T.<x,y,z>=PolynomialRing(F);
S.<u>=PolynomialRing(F);
G.<t>=PolynomialRing(RR);
#-----------------------------------------
# Helper Functions:
#-----------------------------------------
# Computes the number of positive and negative eigenvalues of a given matrix
def signature(M):
    eigen = M.eigenvalues();
    
    pos_count = 0;
    neg_count = 0;
    
    for e in eigen:
        if e > 0:
            pos_count += 1;
        elif e < 0:
            neg_count += 1;
    
    return (pos_count, neg_count);

# Returns the Matrix of the fiber at point p in P^1(R)
def make_matrix(M1, M2, M3, p):
    M = M1*p^2 + 2*M2*p + M3;
    A = matrix(QQ, 1, [-1])
    return block_diagonal_matrix(M, A);

# Given M1, M2, M3 representing the symmetric matrix of Q1, Q2, Q3
# Prints out the real roots of -det(t^2 M_1 + 2 t M_2 + M_3) and the signature in each interval
def print_signature(M1, M2, M3):
    signature_list = [];
    
    char_poly = -det(M1*t^2 + 2*M2*t + M3);
    print("Gamma: " + str(char_poly));
    root_list = list(map(lambda x: x[0], char_poly.roots()));
    length = len(root_list);
    print("Images of real Weierstrass Points (over the chart [t:1] of P^1): " + str(root_list));   
    
    if length == 0:
        print("No real Weierstrass points");
        num_sig = signature(make_matrix(M1, M2, M3, 1));
        signature_list.append(num_sig);
        print(num_sig);
    else:
        for i in range(0, length):
            pt = root_list[i] + 1 if i == length - 1 else (root_list[i] + root_list[i+1])/2
            print("Evaluated at " + str(pt));
            num_sig = signature(make_matrix(M1, M2, M3, pt));
            signature_list.append(num_sig);
            print(num_sig);
            
    return signature_list;

#check whether the discriminant curve Delta is smooth
def is_delta_smooth(f):
    Grad=ideal(f,diff(f,x),diff(f,y),diff(f,z))
    if Grad.dimension()==0:
        print("Delta is smooth.")
        return True
    else:
        print("Delta is not smooth.")
        return False

#check whether the curve Deltatilde is smooth
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

# Converts a given quadratic form in 3 variables to its associated symmetric matrix
def conic_to_matrix(qt):
    a = Rational(qt.coefficient(x^2))
    b = Rational(qt.coefficient(y^2))
    c = Rational(qt.coefficient(z^2))
    d = Rational(qt.coefficient(x*y))
    e = Rational(qt.coefficient(y*z))
    f = Rational(qt.coefficient(x*z))
    
    M = Matrix([
    [a, d/2, f/2],
    [d/2, b, e/2],
    [f/2, e/2, c]
    ]);
    
    return M;
    
#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    # specified inputs

    Q1 = -5*x^2 + 10*x*y + 5*y^2 + 10*x*z + 8*y*z;
    Q2 = x^2 - 10*x*y - 2*y^2 + 10*y*z + 4*z^2;
    Q3 = 2*x^2 - 4*x*y - 4*y^2 - 8*x*z + 8*y*z + 2*z^2;
    
    Delta = Q1*Q3 - Q2^2;
    
    print("Delta: " + str(Delta)); 
    
    M1 = conic_to_matrix(Q1);
    M2 = conic_to_matrix(Q2);
    M3 = conic_to_matrix(Q3);
    
    if is_delta_smooth(Delta) and is_double_cover_smooth(Q1, Q2, Q3):
        try:
            signature_list = print_signature(M1, M2, M3);
        except Exception as e:
            print(e)
    print("------------------------------------------")

if __name__ == '__main__':
    main()
