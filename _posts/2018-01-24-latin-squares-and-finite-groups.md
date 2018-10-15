---
layout: post
title:  "Latin squares and finite groups"
date:   2018-01-24 23:45:00 -0500
categories: [math, python]
math: true
---

Last semester, for an algebra homework, I was trying to prove that there exist only 2 groups of order 6 (namely \\(\mathbb{Z}_6\\) and \\(S_3\\)). The usual argument uses the classification of groups with order \\(pq\\) (with \\(p\\) and \\(q\\) prime), which itself uses Sylow theorems, but I wondered if I could prove it computationally. Here's my attempt.

**Definition:** a **magma** is a pair \\( (G,\ast) \\) where \\(G\\) is a non-empty set and \\(\ast\\) is a binary operation on \\(G\\). A **quasigroup** is a magma in which the equations \\(gx=h\\) and \\(yg=h\\), with \\(g,h\in G\\), always have a unique solution for \\(x\\) and \\(y\\) in \\(G\\). A **group** is a magma in which \\(\ast\\) is associative, modulative and invertive.

Another way of defining a quasigroup is as a magma in which the left and right cancellation laws hold. Note that we're not asking \\(\ast\\) to be associative in a magma (let alone in a quasigroup).

All there is to know about a finite magma \\((G,\ast)\\) is encoded in its **Cayley table**: a matrix that shows how the elements in \\(M\\) operate among themselves. In general, if \\(G = \\{a_1, \dots, a_n\\}\\), the Cayley table looks like this:

$$\begin{array}{c\|ccc}
*&a_1&\dots &a_n\\
\hline
a_1&a_1a_1&\dots &a_1a_n\\
\vdots&\vdots&\dots&\vdots\\
a_n&a_na_1&\dots &a_na_n\\
\end{array}$$

That is, the element in the \\((i,j)\\)-th position is the result of multiplying \\(a_i\\) with \\(a_j\\) in that order.

You can tell a lot about a magma and its operation just by looking at its Cayley table. For example, consider **Klein's 4 group** \\(V\\):

$$\begin{array}{c\|cccc}
*&e&a&b&c\\
\hline
e&e&a&b&c\\
a&a&e&c&b\\
b&b&c&e&a\\
c&c&b&a&e\\
\end{array}$$

You can immediately tell its commutative (because the matrix is symmetric), you can tell there is an identity (namely \\(e\\)) and that each element is its own inverse. You could also tell if the binary operation is associative from its Cayley table using [Light's test](https://en.wikipedia.org/wiki/Light%27s_associativity_test), although it isn't any better, computationally speaking, that just verifying every case by hand.

We will use Cayley tables as the bridge between algebra and combinatorics.

## Latin squares and quasigroups

**Definition:** a \\(n\times n\\) matrix of \\(n\\) different entries is called a **latin square** if no element appears more than once in any row or column. This property is called the **latin square property**.

We will deal with latin squares of size \\(n\\) whose entries are the integers from \\(0\\) to \\(n-1\\). For example

$$\begin{bmatrix}
0&1&2\\
1&2&0\\
2&0&1\\
\end{bmatrix}$$

is a latin square of size \\(3\\). We will also start indexing by 0.

**Theorem 1:** if \\((G, \ast)\\) is a quasigroup, then its Cayley table is a latin square.

**Proof:** Suppose \\(G = \\{a_1, \dots, a_n\\}\\) and that an element \\(b\in G\\) appears twice in row \\(l\\) (say, in columns \\(j\\) and \\(k\\)), by the definition of the Cayley table, this means that

$$ a_la_j = b = a_la_k $$

and because \\(G\\) is a quasigroup, the left cancellation law implies that \\(a_j = a_k\\), which is absurd because we assumed that \\(j\\) and \\(k\\) were different. Analogously, the right cancellation law implies that no element appears twice in any column. **Q.E.D.** 

This theorem has a reciprocal of some sort:

**Theorem 2:** Given a latin square \\(L = (l_{ij})\\), one can construct a quasigroup whose Cayley table is \\(L\\).

**Proof:** Let \\(G = \\{l_{11}, \dots, l_{1n}\\}\\) and denote \\(g_i := l_{1i}\\). Define \\(\ast\\) by

$$g_i * g_j = g_{l_{ij}}$$

by definition, \\(\ast\\) is well defined (that is, it is closed in the set). We need to check that the equations \\(gx =h\\) and \\(yg = h\\) have unique solutions. Consider \\(g_lx = g_k\\), because \\(L\\) is a latin square, \\(g_k = l_{1k}\\) appears somewhere in row \\(l\\), call the column it appears in \\(m\\), then \\(x = g_m\\) is a solution to \\(g_lx = g_k\\). It is unique, because if there existed \\(g_{\widetilde{m}}\\) such that \\(g_lg_{\widetilde{m}} = g_k = g_lg_m\\), then \\(g_k\\) would appear twice in row \\(l\\), which contradicts the fact that \\(L\\) is a latin square. Analogously, now arguing with columns, \\(yg = h\\) has a unique solution in \\(G\\). **Q.E.D**.

So now we're set!, we only need to find all latin squares of size \\(n\\) and to verify if they represent a valid binary operation for a group. Moreover, we could force the existence of an identity by focusing on finding **normalized** (or **reduced**) latin squares (that is, latin squares where the first row and column are \\(0, 1, \dots, n-1\\)).

### The algorithm for finding normalized latin squares of size n.

I use a [depth-first-search](https://en.wikipedia.org/wiki/Depth-first_search) style algorithm, starting with a normalized \\(n\times n\\) matrix 

$$A = \begin{bmatrix}
0&1&\cdots&n-1\\
1&-1&\cdots&-1\\
\vdots&\vdots&\ddots&\vdots\\
n-1&-1&\cdots&-1\\
\end{bmatrix}$$

where the unvisited locations are labeled with a \\(-1\\). We also start with an empty [stack](https://en.wikipedia.org/wiki/Stack_(abstract_data_type)) \\(S\\). The algorithm goes as follows:

1. Put matrix \\(A\\) in the stack \\(S\\).
2. If the stack \\(S\\) is empty, stop; if it isn't, pop a matrix \\(B\\) from it.
3. Find the first unvisited position \\((i,j)\\) in \\(B\\) (i.e. the first \\(-1\\)), if there isn't any, it is finished, put it in a special list of finished latin squares and go to step 2.
4. Push into the stack the result of replacing this \\(-1\\) with every number from \\(0\\) to \\(n-1\\) that isn't already on its row or column.
5. Go to step 2.

Here's the algorithm implemented in python:

```python
def dfs_in_matrix(A):
    # First we create an empty stack and we put the initial matrix
    # in it.
    list_of_solutions = []
    stack = LifoQueue()
    stack.put(A)
    
    while not stack.empty():
        # We pop a matrix from the stack
        B = stack.get()
        
        # We check if it's finished.
        if is_finished(B):
            list_of_solutions.append(B)
            continue
        
        # We find an unvisited position
        position = find_first_unvisited_position(B)
        if position == None:
            continue
        
        span = span_of_position(position, B)
        for k in range(len(B)):
            if k not in span:
                C = B.copy()
                C[position] = k
                stack.put(C)

    return list_of_solutions

def find_normalized_latin_squares(number):
    A = np.zeros((number, number))
    for k in range(number):
        A[0, k] = k
        A[k, 0] = k
    for i in range(1, number):
        for j in range(1, number):
            A[i, j] = -1
    list_of_solutions = dfs_in_matrix(A)
    return list_of_solutions
```
(the functions `is_finished`, `find_first_unvisited_position` and `span_of_position` are auxiliary, check [this jupyter notebook](https://gist.github.com/miguelgondu/404619477e50ccec62db8c53b4901091) for all the code discussed in this post). It checks out with the literature on the topic[^1], saying that there are 9408 normalized latin squares of size 6.

### The `Magma` class

Once we have all the normalized latin squares, we can build up a `Magma` class in python and we can write a verification function to find which of these correspond to associative operations (and thus to groups).

```python
class Magma:
    def __init__(self, _matrix):
        self.cayley_table = _matrix
        self.set = set(range(0, len(_matrix[0,:])))
    
    def mult(self, a, b):
        return int(self.cayley_table[a, b])

def is_magma_associative(mag):
    '''
    This function verifies if magma `mag` is associative by brute force.
    '''
    n = len(mag.cayley_table)
    _flag = True
    for a in range(n):
        for b in range(n):
            for c in range(n):
                _flag = _flag and (mag.mult(a, mag.mult(b,c)) == mag.mult(mag.mult(a,b),c))
    return _flag

def find_groups(number):
    latin_squares = find_normalized_latin_squares(number)
    associative_magmas = [sol for sol in latin_squares if is_magma_associative(Magma(sol))]
    return associative_magmas
```

After running this `is_magma_associative` in all 9408 reduced latin squares of order 6 we're left with 80 reduced latin squares such that, when interpreted as quasigroups, are associative. That is, only 80 of the original 9408 reduced latin squares of size 6 can be interpreted as Cayley tables for groups.

## The main result

We're trying to prove the following theorem:

**Theorem 3:** There are only 2 groups of order 6, namely \\(S_3\\) and \\(\mathbb{Z}_6\\).

It is useful, then, to cleary state what we interpret as \\(S_3\\) and \\(\mathbb{Z}_6\\). \\(\mathbb{Z}_6\\) are the usual integers modulo 6 with sum modulo 6, but note that \\(\mathbb{Z}_6\\) can also be interpreted in the following way: its a group of six elements \\(\\{a_1, a_2, a_3, a_4, a_5, 0\\}\\) such that

- \\(a_1\\) and \\(a_5\\) have order 6.
- \\(a_2\\) and \\(a_4\\) have order 3.
- \\(a_3\\) has order 2.
- \\(a_2^2 = a_4\\)
- \\(a_1a_2 = a_3\\)

(note that we just changed \\(i\\) for \\(a_i\\)). We can use this information to find an isomorphism between a latin-square-generated group and \\(\mathbb{Z}_6\\).

\\(S_3 = \\{\sigma_1, \sigma_2, \sigma_3, \rho_1, \rho_2, \mbox{id}\\}\\) is usually interpreted as the group of symmetries of a triangle (where \\(\sigma_i\\) is the reflection that fixes vertex \\(i\\) and \\(\rho_j\\) is a rotation of \\(120*j\\) degrees, but we prefer the following presentation:

$$S_3 = \langle \sigma, \rho\,\vert\,\sigma^2 = \rho^3 = \mbox{id},\, \sigma\rho = \rho^2\sigma \rangle$$

In this presentation, the 6 different elements are \\(\mbox{id}, \sigma, \rho, \rho\sigma, \rho^2\sigma\\) and \\(\rho^2\\).

So, to prove theorem 3, we will follow this strategy: we will give an isomorphism from either \\(S_3\\) or \\(\mathbb{Z}_6\\) to each of the 80 groups found using latin squares:

**Proof (of Theorem 3):** Theorem 1 and 2 show that all possible groups of a certain order are restricted by the amount of normalized latin squares of said order. After filtering the normalized latin squares of size 6 by verifying which represent an associative binary operation, we're left with 80 Cayley tables for groups. In [this jupyter notebook](https://gist.github.com/miguelgondu/404619477e50ccec62db8c53b4901091) we show an explicit isomorphism between each of these 80 latin square generated groups and either \\(\mathbb{Z}_6\\) or \\(S_3\\), but for the sake of completeness we show how these isomorphisms were constructed with explicit examples for \\(\mathbb{Z}_6\\) and \\(S_3\\). Consider the following normalized latin square:

$$\begin{bmatrix}
0 & 1 & 2 & 3 & 4 & 5\\
1 & 5 & 4 & 2 & 3 & 0\\
2 & 3 & 0 & 1 & 5 & 4\\
3 & 4 & 5 & 0 & 1 & 2\\
4 & 2 & 1 & 5 & 0 & 3\\
5 & 0 & 3 & 4 & 2 & 1\end{bmatrix}$$

and call the group it induces \\(G\\). After inspecting it, we can tell that the orders of their elements are either \\(2\\) or \\(3\\), so it is a candidate for being isomorphic to \\(S_3\\). Choose \\(4\mapsto \sigma\\) and \\(5\mapsto \rho\\), and note that
$$(5*5)*4 = 1*4 = 3 = 4*5$$
that is, this group obeys the presentation given for \\(S_3\\). The isomorphism would then be given by

$$\begin{matrix}
S_3 & & G\\
\hline
\mbox{id}&\mapsto&0\\
\sigma_1&\mapsto&4\\
\sigma_2&\mapsto&3\\
\sigma_3&\mapsto&2\\
\rho&\mapsto&5\\
\rho^2&\mapsto&1
\end{matrix}$$

Now consider the group \\(H\\) given by the following reduced latin square:

$$\begin{bmatrix}
0 & 1 & 2 & 3 & 4 & 5\\
1 & 5 & 4 & 2 & 3 & 0\\
2 & 4 & 5 & 1 & 0 & 3\\
3 & 2 & 1 & 0 & 5 & 4\\
4 & 3 & 0 & 5 & 1 & 2\\
5 & 0 & 3 & 4 & 2 & 1\end{bmatrix}$$

because the elements of \\(H\\) have order either 2, 3 or 6, we will construct an isomorphism between \\(H\\) and \\(\mathbb{Z}_6\\) using the identification \\(\mathbb{Z}_6 = \\{a_1, a_2, a_3, a_4, a_5, 0\\} \\) stated before. First note that \\(3\in H\\) is an element of order 2, \\(2, 4\in H\\) have order 6 and \\(1, 5\in H\\) have order 3. Because \\(1\ast1 = 5\\) and \\(4\ast1 = 3\\), we construct the following isomorphism

$$\begin{matrix}
\mathbb{Z}_6 & & H\\
\hline
0&\mapsto&0\\
1&\mapsto&4\\
2&\mapsto&1\\
3&\mapsto&3\\
4&\mapsto&5\\
5&\mapsto&2\\
\end{matrix}$$

 **Q.E.D**

---

[^1]: The results the algorithm gave were in par with what's said in *[Small Latin Squares, Quasigroups and Loops](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.71.1011&rep=rep1&type=pdf)*, an article by Brendan D. Mackay, Alison Meynert and Wendy Myrvold. Check *Table 1* of their article for more details. 