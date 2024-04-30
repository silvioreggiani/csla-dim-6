algebras_nums_dict = load('algebras_nums_dict')

algebras_nums_dict_bis = { 'h9_bis' : [ '0', '0' , '0', '0', '21', '15+23'],
                           'h9_tris' :  [ '0' , '0' , '0' , '0' , '12', '51+23' ],
                           'h18_bis' : [ '0', '0' , '0', '12', '13', '15+24'],
                           'h19+_bis' : [ '0', '0' , '0', '32', '12', '14+35'],
                           'h19+_tris' : [ '0', '0' , '0' , '23' , '21', '14+35'],
                           'h26-_bis' : [ '0', '0' , '12', '31', '32', '15+24']}
algebras_nums_dict.update(algebras_nums_dict_bis)

names_lower = load('Names_Lower_Algs_of_der_list.sobj')
names_bis = ['h9_bis','h18_bis','h19+_bis','h26-_bis']
names_tris = ['h9_bis','h9_tris', 'h18_bis','h19+_bis', 'h19+_tris', 'h26-_bis']
names = names_lower + names_tris
algebras_nums_l = {key: algebras_nums_dict[key] for key in names}

def coef_str(nums):
    
    coef_str_dict = dict()
    for k in range(6):
        st = nums[k].replace(" ", "")
        sg = -1 
        if len(st) == 2:
            coef_str_dict[("e"+st[0], "e"+st[1])] = {"e"+str(k+1) : sg }
        elif len(st) == 5:
            coef_str_dict[("e"+st[0], "e"+st[1])] = {"e"+str(k+1) : sg }
            if st[2] == "-":
                sg = 1
                coef_str_dict[("e"+st[3], "e"+st[4])] = {"e"+str(k+1) : sg }
            else:
                coef_str_dict[("e"+st[3], "e"+st[4])] = {"e"+str(k+1) : sg }
        if len(st) == 8:
            coef_str_dict[("e"+st[0], "e"+st[1])] = {"e"+str(k+1) : sg }
            coef_str_dict[("e"+st[3], "e"+st[4])] = {"e"+str(k+1) : sg }
            coef_str_dict[("e"+st[6], "e"+st[7])] = {"e"+str(k+1) : sg }
    return coef_str_dict

def is_square(M):
    return M.nrows() == M.ncols()
def is_lower_triang(M):
    a = 0
    if is_square(M) == True:
        n = M.ncols()
        for i in range(n):
            for j in range(i+1,n):
                if M[i,j] != 0:
                    a = a+1
        return a == 0
    else:
        print('The matrix is not square.')
def is_s_lower_triang(M):
    a = 0
    if is_square(M) == True:
        n = M.ncols()
        for i in range(n):
            for j in range(i,n):
                if M[i,j] != 0:
                    a = a+1
        return a == 0
    else:
        print('The matrix is not square.')
def is_upper_triang(M):
    a = 0
    if is_square(M) == True:
        n = M.ncols()
        for i in range(n):
            for j in range(i):
                if M[i,j] != 0:
                    a = a+1
        return a == 0
    else:
        print('The matrix is not square.')
def is_diagonal(M):
    a = False
    if is_lower_triang(M) == True:
        if is_upper_triang(M) == True:
            a = True
    return a

def construct_the_algebra(algebra_name):    
    global alg_name
    alg_name = algebra_name
    global label
    label = algebras_nums_dict[alg_name]
    global ce
    ce = coef_str(label)
    global ce_dict
    ce_dict = { alg_name : ce}
    global alg
    alg.<e1,e2,e3,e4,e5,e6> = LieAlgebra(SR, ce, nilpotent=True )
    global basis
    basis = alg.basis().values()
    global der_basis
    der_basis= alg.derivations_basis()
    global der_dim
    der_dim= len(der_basis)
    global der_vars
    der_vars = [var("d"+str(i)) for i in range(der_dim)]
    global der_gen
    der_gen = sum(der_vars[i] * der_basis[i] for i in range(der_dim))
    global Center
    Center = alg.center()
    global Center_dim
    Center_dim = Center.dimension()
    global Derivated
    Derivated = alg.derived_subalgebra()
    global Derivated_dim
    Derivated_dim = Derivated.dimension()

def info_alg_and_der(algebra):
    print("The Lie Algebra is: " + alg_name)
    print(alg)
#    print("Its basis is:")
#    print(basis)
    print("About the center of the Lie algebra Sagemath tells us: ")
    print(Center)
    print("About the commutator of the Lie algebra Sagemath tells us: ")
    print(Derivated)
    print("The algebra of derivations is of dimension: " + str(der_dim))
#    print (" and a basis given by: ")
#    display(der_basis)
    print("A generic derivation is of the form:")
    display(der_gen)
    if is_lower_triang(der_gen) == True:
        print("which is lower triangular.")
    else:
        print("which is not in lower triangular form.")

def generic_derivations(algebra_name):
    construct_the_algebra(algebra_name)
    global der_basis_diag
    der_basis_diag = []
    global index_der_basis_diag
    index_der_basis_diag = []
    global der_basis_slt
    der_basis_slt = []
    global index_der_basis_slt
    index_der_basis_slt = []
    for i in range(der_dim):
            deri = der_basis[i]
            if is_diagonal(deri) == True:
                der_basis_diag.append(deri)
                index_der_basis_diag.append(i)
            if is_s_lower_triang(deri) == True:
                der_basis_slt.append(deri)
                index_der_basis_slt.append(i)
    global der_gen_diag
    der_gen_diag = diagonal_matrix(der_gen.diagonal()) #generic diagonal derivation
    global der_gen_slt
    der_gen_slt = der_gen - der_gen_diag #generic strictly triangular derivation
    global der_vars_diag
    der_vars_diag = [der_vars[i] for i in index_der_basis_diag]
    global der_vars_slt
    der_vars_slt = [der_vars[i] for i in index_der_basis_slt]

def first_appearence(matriz, variables):
    coefficient_dict = {}
    for d in variables:
        coefficient_dict[d] = sorted(matriz.find(lambda coef: coef == d,indices=True))[0]
    return coefficient_dict

def generic_automorphisms(algebra_name):
#    construct_the_algebra(algebra_name)
    generic_derivations(algebra_name)
    global auto_gen
    auto_gen = exp(der_gen)
    global aut_gen_slt
    aut_gen_slt = exp(der_gen_slt)
    global aut_gen_diag
    aut_gen_diag = exp(der_gen_diag)
    global aut_vars
    aut_vars =  [var("a"+str(i)) for i in range(der_dim)]
    global aut_vars_diag
    aut_vars_diag = [aut_vars[i] for i in index_der_basis_diag]
    global aut_vars_slt
    aut_vars_slt = [aut_vars[i] for i in index_der_basis_slt]
    global der_vars_coeffs
    der_vars_coeffs = first_appearence(der_gen,der_vars)
    A = aut_gen_slt
    for k in range(len(der_vars_slt)):
        d = der_vars_slt[k]
        a = aut_vars_slt[k]
        i,j = der_vars_coeffs[d]
        A = A.subs(d==a-A[i,j]+d)
    global aut_gen_diag_list 
    aut_gen_diag_list = [exp(der_vars_diag[k] * der_basis_diag[k]).subs(der_vars_diag[k]==log(aut_vars_diag[k])) for k in range(len(der_vars_diag))]
    for A_diag in aut_gen_diag_list:
        A = A * A_diag
        for k in range(len(aut_vars_slt)):
            d = der_vars_slt[k]
            a = aut_vars_slt[k]
            i,j = der_vars_coeffs[d]
            A = A.subs(a==a/(A[i,j]/a))
    A.simplify_full()
    global aut_gen
    aut_gen = A.simplify_full()

def generic_metric(algebra_name):
#    construct_the_algebra(algebra_name)
#    generic_derivations(algebra_name)
    generic_automorphisms(algebra_name)
    global aut_vars_coeffs
    aut_vars_coeffs = first_appearence(aut_gen,aut_vars)
    global Sigma
    Sigma = identity_matrix(SR, 6)
    aut_vars_coeffs_values = aut_vars_coeffs.values()
    count = 0
    for i in range(6):
        for j in range(i+1):
            if (i,j) not in aut_vars_coeffs_values:
                Sigma[i,j]=var("s"+str(count))
                count += 1
    global g
    g = transpose(Sigma)*Sigma
