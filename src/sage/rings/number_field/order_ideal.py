"""
Ideals of Number Field Orders
"""

from __future__ import absolute_import
from sage.misc.cachefunc import cached_method
from sage.rings.ideal import Ideal_generic
from sage.rings.all import ZZ


class OrderIdeal(Ideal_generic):
    """
    An ideal of a number field order.
    """

    def __init__(self, order, gens, coerce=True):
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        self._order = order
        field = order.number_field()
        Ideal_generic.__init__(self, field, gens, coerce)

    # This should be changed
    def _repr_short(self):
        return '(%s)' % (', '.join(map(str, self.gens())))

    def __repr__(self):
        return "Fractional ideal %s of non-maximal order" % self._repr_short(
        )

    def is_integral(self):
        """
        Return True if this ideal is integral.

        EXAMPLES::

           sage: R.<x> = PolynomialRing(QQ)
           sage: K.<a> = NumberField(x^5-x+1)
           sage: K.ideal(a).is_integral()
           True
           sage: (K.ideal(1) / (3*a+1)).is_integral()
           False
        """
        try:
            return self.__is_integral
        except AttributeError:
            self.__is_integral = all(a in self._order for a in self.gens())
            return self.__is_integral

    def order(self):
        """
        Return the number field order containing this ideal.

        EXAMPLES::

           sage: K.<a> = NumberField(x^2-5)
           sage: O = K.order(a)
           sage: I = O.ideal(5)
           sage: I.order() == O
           True
        """

        return self._order

    def number_field(self):
        """
        Return the number field that this is a fractional ideal in.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 2); K
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(3).number_field()
            Number Field in a with defining polynomial x^2 + 2
            sage: K.ideal(0).number_field() # not tested (not implemented)
            Number Field in a with defining polynomial x^2 + 2
            sage: O = K.order(2*a); O.is_maximal()
            False
            sage: O.ideal(3).number_field()
            Number Field in a with defining polynomial x^2 + 2
        """
        return self.ring()


    # @cached_method
    # def basis(self):
    #     """
    #     Return a basis for self as a ZZ-module
    #     """
    #     K = self.number_field()
    #     order_basis = self.ring().basis()
    #     from itertools import product
    #     ZZ_gens = [x[0] * x[1] for x in product(order_basis, self.gens())]

    #     from sage.matrix.constructor import Matrix
    #     M = Matrix([x.vector() for x in ZZ_gens])
    #     d = M.denominator()
    #     dM = (d * M).change_ring(ZZ)
    #     H = dM.hermite_form(include_zero_rows=False)
    #     Hd = H / d
    #     return [K(x) for x in Hd.rows()]

    # @cached_method
    # def free_module(self):
    #     return basis_to_module(self.basis(), self.number_field())

    # def prime_factors(self):
    #     O = self.order()
    #     OK = O.integral_closure()

    #     aOK = OK.ideal(self.gens())
    #     above = aOK.prime_factors()

    #     return [O.intersection(p) for p in above]

def basis_to_module(B, K):
    r"""
    Given a basis `B` of elements for a `\ZZ`-submodule of a number
    field `K`, return the corresponding `\ZZ`-submodule.

    EXAMPLES::

        sage: K.<w> = NumberField(x^4 + 1)
        sage: from sage.rings.number_field.order_ideal import basis_to_module
        sage: basis_to_module([K.0, K.0^2 + 3], K)
        Free module of degree 4 and rank 2 over Integer Ring
        User basis matrix:
        [0 1 0 0]
        [3 0 1 0]
    """
    V, from_V, to_V = K.absolute_vector_space()
    M = ZZ**(V.dimension())
    C = [to_V(K(b)) for b in B]
    return M.span_of_basis(C)
