# group theory home work
# a toolbox to find class of each element, subgroups and inverent subgroups
# this script is developed for trivial calculations in homeworks, and is not tested or optimized
# for production propers please refer to packages such as pygroup
# contact kefeiwu@outlook.com for questions, bugs, suggestions and improvements

import numpy as np

element_table = ['E','A','B','C','D','F','I','J','K','L','M','N'] # sepsifing a element table, order matters
(E,A,B,C,D,F,I,J,K,L,M,N) = (0,1,2,3,4,5,6,7,8,9,10,11) # converting table into number, Identity must be 0
element_lens = len(element_table)


table = ([[E,A,B,C,D,F,I,J,K,L,M,N],
          [A,E,F,I,J,B,C,D,M,N,K,L],
          [B,F,A,K,L,E,M,N,I,J,C,D],
          [C,I,L,A,K,N,E,M,J,F,D,B],
          [D,J,K,L,A,M,N,E,F,I,B,C],
          [F,B,E,M,N,A,K,L,C,D,I,J],
          [I,C,N,E,M,L,A,K,D,B,J,F],
          [J,D,M,N,E,K,L,A,B,C,F,I],
          [K,M,J,F,I,D,B,C,N,E,L,A],
          [L,N,I,J,F,C,D,B,E,M,A,K],
          [M,K,D,B,C,J,F,I,L,A,N,E],
          [N,L,C,D,B,I,J,F,A,K,E,M]]) # input your multiplication table
table = np.array(table)

# global variables are element_table, table (as multiplication table) and element_lens



def group_multiplication(mtable, a, b):# find the multiplication of a*b (index only) using multiplication table
    return (mtable[a][b])

def group_multi_seq(mtable, slist):# find the result of a sequence of multiplication (mtable, [g0,g1, ..., gi]) = g0*g1*...*gi...
    e0 = 0
    for i in (slist):
        e0 = group_multiplication(mtable, e0, i)
    return e0

def element2index(etable,x):
    for i in range(element_lens):
        if etable[i]==x:
            
            return i

def find_inverse(mtable, x):#given an element (number or symble) and a multiplication table, find its inverse
    a0 = 0 # index of our element
    if type(x)==int:
        a0 = x
    else:
        a0 = element2index(element_table,x)

    tempt = mtable[a0,:]
    for i in range(element_lens):
        if tempt[i]==0:
            i_index = i

    if type(x)==int:
        return i_index
    else:
        return element_table[i_index]
    
def findclass(mtable,x): # find the class of element x, using multiplication table
    if type(x) == int:
        a0=x
    else:
        a0=element2index(element_table,x)
    
    b0 = a0
    classtable = [a0] 
    for i in range(element_lens):
        # try to calculate i^-1*a0*i
        p0 = i
        p1 = find_inverse(mtable,i)
        b0 = group_multi_seq(mtable, [p1,a0,p0])
        
        if b0 not in classtable:
            classtable.append(b0)
    if type(x)==int:
        return classtable
    else:
        index_table = []
        for i in classtable:
            index_table.append(element_table[i])
        return index_table

def tell_subgroup(mtable, elist): #is a collect (elist) a group, elist supports index only
    outsider_count = 0
    no_inverse = 0
    group_indicator=0
    # closure
    for i in elist:
        for j in elist:
            a0 = group_multiplication(mtable, i,j)
            if a0 not in elist:
                outsider_count += 1
    # inverse
    for i in elist:
        b0 = find_inverse(mtable,i)
        if b0 not in elist:
            no_inverse += 1
    
    judgement_crit = no_inverse+outsider_count
    if judgement_crit==0:
        return 1 # inverse of all elements is inside and closure rule is satisfied
    else:
        return 0


#print(tell_subgroup(table,[0,10,11]))



def find_all_subgroup(mtable, output_format=1): # find all subgroups of a given multi_table, output_format can be either index (output_format=1) or symble(output_format=0)
    #generate all possible collections of group
    # since E must be in one subgroup, the remining elements (n-1) have 2^(n-1) posible ways
    c_num = int(2**(element_lens-1))
    subgroup_list = []
    for i in range(c_num):
        bin_i = bin(i) #use the binary number of i to indicate a possible way of collection
        tempi = str(bin_i)
        tempi = tempi[2:len(tempi)]
        
        while len(tempi)<(element_lens-1):
            tempi='0'+tempi
        tempi = '1'+tempi #Every (sub)group must include Identity



        temp_collection = []

        for i in range(element_lens):
            if tempi[i]=='1':
                temp_collection.append((i))

        a = tell_subgroup(mtable, temp_collection)
        if a==1:
            subgroup_list.append(temp_collection)


    if output_format==1:
        return subgroup_list
    else:
        subgroup_list_symble = []
        for i in subgroup_list:
            b0 = []
            for j in i:
                b0.append(element_table[j])
            subgroup_list_symble.append(b0)
        
        return subgroup_list_symble

def tell_subgroup_invariant(mtable, elist): # tell if a subgroup is invariant, 
    isgroup = tell_subgroup(mtable,elist)
    if isgroup==0:
        print('WARNING! NOT A GROUP')

    outsider=0
    for i in range(element_lens):
        for j in elist:
            #calculate i^-1*a*i
            a0 = i
            a1 = find_inverse(mtable,i)
            p0 = group_multi_seq(mtable,[a1,j,a0])
            if p0 not in elist:
                #print(a0,a1,j,p0)
                outsider += 1
    
    if outsider==0:
        return 1
    else:
        return 0


def find_all_inverent_subgroup(mtable, output_format=1): # find all inverent subgroup of a table, output in the format of index(1) or symble(0)
    invariant_list = []
    all_subgroup = find_all_subgroup(mtable, output_format=1)
    for i in all_subgroup:
        if tell_subgroup_invariant(mtable,i)==1:
            invariant_list.append(i)
    
    if output_format==1:
        return invariant_list
    else:
        invariant_subgroup_list_symble = []
        for i in invariant_list:
            b0 = []
            for j in i:
                b0.append(element_table[j])
            invariant_subgroup_list_symble.append(b0)
        return invariant_subgroup_list_symble



print('All subgroup(s) \n', find_all_subgroup(table, output_format=0))
print('All invariant subgroup(s) \n', find_all_inverent_subgroup(table, output_format=0))

        
'''
for i in range(element_lens):
    e_symble = element_table[i]
    class_e = findclass(table,e_symble)
    print('the class of element '+e_symble+' is' + str(class_e))     

'''



#print(table[1], element_lens)
