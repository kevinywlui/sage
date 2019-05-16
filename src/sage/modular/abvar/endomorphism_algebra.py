r"""
Space of algebras obtained by tensoring an endomorphism ring of a modular
abelian variety by QQ.
"""
from sage.modules.free_module_element import vector
from sage.rings.ring import Ring
from sage.structure.element import RingElement
from sage.rings.all import QQ, ZZ
from sage.rings.number_field.number_field import NumberField
from sage.matrix.constructor import matrix
from sage.structure.unique_representation import UniqueRepresentation


class Morphism(RingElement):
    r"""
    Morphisms in the isogeny category of modular abelian varieties.
    """

    def __init__(self, parent, matrix):
        self._matrix = matrix
        RingElement.__init__(self, parent)
        self._parent = parent

    def parent(self):
        return self._parent

    def matrix(self):
        return self._matrix

    def _im_gens_(self, codomain, im_gens):
        r"""
        Return the image of ``self`` in codomain under the map that sends
        the images of the generators of the parent of ``self`` to the
        tuple of elements of im_gens.

        We assume that im_gens forms a rational basis for the endomorphism
        algebra.
        """
        Bs = self.parent().gens()
        coeffs = self._linear_combination_coefficients(Bs)
        return codomain(sum(x * y for x, y in zip(coeffs, im_gens)))

    def __repr__(self):
        return "Element of {} given by \n{}".format(self.parent(),
                                                    self.matrix())

    def _linear_combination_coefficients(self, Bs):
        r"""
        Return the coefficients needed to write self as a linear combination of
        elements of Bs.

        INPUT:

        - ``Bs`` -- a list of morphisms containing self in its span.

        OUTPUT:

        A list of numbers of length equal to ``Bs``.
        """
        Bmatrix = matrix(QQ, [b.list() for b in Bs])
        return Bmatrix.solve_left(vector(self.list()))

    def matrix_space(self):
        return self.parent().matrix_space()

    def list(self):
        r"""
        Return a list of elements in the matrix of self.
        """
        return self.matrix().list()

    def __mul__(self, other):
        M = self.matrix_space()
        try:
            return self.parent()(self.matrix() * M(other))
        except:
            pass
        return self.parent()(self.matrix() * other.matrix())

    def __add__(self, other):
        M = self.matrix_space()
        try:
            return self.parent()(self.matrix() + M(other))
        except:
            pass
        return self.parent()(self.matrix() + other.matrix())

    def __eq__(self, other):
        M = self.matrix_space()
        try:
            return self.matrix() == M(other)
        except:
            pass
        return self.matrix() == other.matrix()


class EndomorphismAlgebra(Ring, UniqueRepresentation):
    r"""
    Class of a algebras obtained by tensoring a endomorphism ring by QQ.
    """
    Element = Morphism

    def __init__(self, endomorphism_ring, category=None):
        self._E = endomorphism_ring
        self._E_free = self._E.free_module()
        self._A = self._E.abelian_variety()
        self._ambient_space = self._E_free.ambient_vector_space()

        A = self._A
        if category is None:
            homset_cat = A.category()
        else:
            homset_cat = category
        # Remark: Ring.__init__ will automatically form the join
        # of the category of rings and of homset_cat
        Ring.__init__(self, A.base_ring(), category=homset_cat.Endsets())

    def _repr_(self):
        return "Endomorphism ALGEBRA of {}".format(self.abelian_variety())

    def _element_constructor_(self, x):
        r"""
        This constructors element of this endomorphism algebra. We can
        construct from matrices or rationals.
        """
        M = self.matrix_space()
        if x in M:
            return self.element_class(self, M(x))
        elif x in self.endomorphism_ring():
            return self.element_class(self, M(x.matrix()))
        elif isinstance(x, Morphism):
            return self.element_class(self, M(x.matrix()))
        try:
            return self.element_class(self, M(x))
        except TypeError:
            pass
        raise ValueError('x must be either a Morphism'
                         ' or a matrix of the suitable size')

    def an_element(self):
        r"""
        Return an element of self by returning the last generator.

        OUTPUT: an element
        """
        return self.gens()[-1]

    def some_elements(self):
        r"""
        Return some elements of self by returning the generators.

        OUTPUT: a tuple of elements
        """
        return self.gens()

    def matrix_space(self):
        r"""
        Return the underlying matrix space of this endomorphism algebra.

        OUTPUT: a matrix space.

        EXAMPLES::

            sage: J = J0(23); J.dimension()
            2
            sage: J.endomorphism_algebra().matrix_space()
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field
        """
        return self.endomorphism_ring().matrix_space().change_ring(QQ)

    def identity(self):
        r"""
        Return the identity element of this endomorphism algebra.

        OUTPUT: the identity morphism of this endomorphism algebra.

        EXAMPLES::

            sage: J = J0(33)[2];
            sage: A = J.endomorphism_algebra()
            sage: one = A.identity()
            sage: all([one*x == x for x in A.gens()])
            True
        """
        M_one = self.matrix_space().one()
        return self.element_class(self, M_one)

    def _coerce_map_from_(self, S):
        r"""
        Can coerce from spaces that can be coerce into rational square matrices
        of dimension 2*d by 2*d, where d is the dimension of the abelian
        variety.
        """
        M = self.matrix_space()
        if M.coerce_map_from(S):
            return True
        elif S == self:
            return True

    def characteristic(self):
        r"""
        Return the characteristic of this endomorphism algebra which is zero.

        OUTPUT: the integer 0.
        """
        return ZZ(0)

    def is_field(self, proof=True):
        r"""
        Return whether this endomorphism algebra is a field.

        OUTPUT: a boolean
        """
        return self.abelian_variety().is_simple()

    def basis(self):
        r"""
        Return a $\QQ$-basis of this endomorphism algebra.

        OUTPUT: a tuple of elements.
        """
        F = self.free_module()
        return tuple(self(x) for x in F.basis())

    def gens(self):
        r"""
        Return a set of $\QQ$-module generators.

        OUTPUT: a tuple consisting of elements.
        """
        return self.basis()

    def gen(self, i):
        return self.gens()[i]

    def ngens(self):
        return len(self.gens())

    def random_element(self):
        r"""
        Return a random element of self by returning a random element in the
        endomorphism ring.
        """
        E = self.endomorphism_ring()
        return self(E.random_element())

    def endomorphism_ring(self):
        r"""
        Return the endomorphism ring this algebra was constructed from.
        """
        return self._E

    def abelian_variety(self):
        r"""
        Return the abelian variety attached to self.
        """
        return self.endomorphism_ring().abelian_variety()

    def power_basis(self):
        r"""
        Return a power basis for self when self is a field.

        We use a random algorithm. So we cached both for performance and
        consistency reasons.

        OUTPUT: a list of morphisms in this endomorphism algebra.
        """
        try:
            return self._a_power_basis
        except AttributeError:
            pass

        d = self.abelian_variety().dimension()

        # this will very likely work
        for i in range(100):
            phi = self.random_element()
            f = phi.matrix().minpoly()
            if f.degree() == d:
                break

        # self._a_power_basis = tuple(self(phi**i) for i in range(d))
        self._a_power_basis = tuple([self.one()] +
                                    [self(phi**i) for i in range(1, d)])
        return self._a_power_basis

    def isomorphic_field(self, both_maps=True, check=True):
        r"""
        Return maps to and from the endomorphism algebra as a number field when
        self is a field.

        We use a random algorithm. So we cached both for performance and
        consistency reasons.

        INPUT:

        - ``both_maps`` (default: True) -- boolean determining whether to
          return isomorphisms to and from the field.
        - ``check`` (default: True) -- boolean determining whether to check
          that the composition yields the identity.

        OUTPUT:

        - a tuple consisting of ``(K, K_to_EA, EA_to_K)`` or just ``K``, where
            - ``K`` is a number field isomorphic to self.
            - ``K_to_EA`` is an isomorphism from ``K`` to self, or ``None``.
            - ``EA_to_K`` is an isomorphism from self to ``K``, or ``None``.

        EXAMPLES::

            sage: J = J0(29)
            sage: A = J.endomorphism_algebra()
            sage: K, K_to_EA, EA_to_K = A.isomorphic_field()
            sage: K.disc()
            8
            sage: alpha = K.gens()[0]
            sage: EA_to_K(K_to_EA(alpha)) == alpha
            True
            sage: EA_to_K(K_to_EA(0)) == 0
            True
            sage: K_to_EA(EA_to_K(0)) == 0
            True

            sage: J = J0(11); J.dimension()
            1
            sage: J.endomorphism_algebra().isomorphic_field()
            (Number Field in alpha with defining polynomial x - 1, Ring morphism:
               From: Number Field in alpha with defining polynomial x - 1
               To:   Endomorphism ALGEBRA of Abelian variety J0(11) of dimension 1
               Defn: 1 |--> Element of Endomorphism ALGEBRA of Abelian variety J0(11) of dimension 1 given by
                     [1 0]
                     [0 1], Ring morphism:
               From: Endomorphism ALGEBRA of Abelian variety J0(11) of dimension 1
               To:   Number Field in alpha with defining polynomial x - 1
               Defn: Element of Endomorphism ALGEBRA of Abelian variety J0(11) of dimension 1 given by
                     [1 0]
                     [0 1] |--> 1)
        """
        if not self.is_field():
            raise ValueError("self must be a field")

        power_basis = self.power_basis()

        if self.abelian_variety().dimension() == 1:
            phi = power_basis[0]
        else:
            phi = power_basis[1]
        f = phi.matrix().minpoly()

        K = NumberField(f, names='alpha')
        alpha = K.gens()[0]

        # here EA represents the endomorphism algebra
        K_to_EA = K.hom([phi])

        # we know phi**i -> alpha**i so write elements of self.gens() in terms
        # of phi.
        P = matrix([p.list() for p in power_basis])
        G = matrix([g.list() for g in self.gens()])

        C = P.solve_left(G)

        d = len(power_basis)
        EA_to_K = self.hom(
            [sum(C[i][j] * alpha**j for j in range(d)) for i in range(d)],
            check=False)

        if check:
            assert (K_to_EA(EA_to_K(x)) == x for x in self.gens())
            assert (EA_to_K(K_to_EA(alpha)) == alpha)

        if both_maps:
            return K, K_to_EA, EA_to_K
        else:
            return K

    def free_module(self):
        r"""
        Return this endomorphism ring as a submodule of the free module
        `\QQ^{4d^2}`, where `d` is the dimension of the associated abelian
        variety.
        """
        return self.endomorphism_ring().free_module().change_ring(QQ)
