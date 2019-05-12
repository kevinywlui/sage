r"""
Space of algebras obtained by tensoring an endomorphism ring of a modular
abelian variety by QQ.
"""
from sage.categories.homset import HomsetWithBase
from sage.modules.matrix_morphism import MatrixMorphism
from sage.modules.free_module import FreeModule_submodule_field
from sage.rings.ring import Ring
from sage.structure.element import RingElement
from sage.rings.all import QQ, ZZ
from sage.rings.number_field.number_field import NumberField
from sage.structure.unique_representation import UniqueRepresentation


class Morphism(MatrixMorphism, RingElement):
    r"""
    Morphisms in the isogeny category of modular abelian varieties.
    """

    def __init__(self, parent, A):
        MatrixMorphism.__init__(self, parent, A)
        RingElement.__init__(self, parent)
        self._parent = parent

    def parent(self):
        return self._parent

class EndomorphismAlgebra(UniqueRepresentation, HomsetWithBase, Ring):
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
        HomsetWithBase.__init__(self, A, A, category=homset_cat)

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
        else:
            raise ValueError('x must be either a Morphism'
                             'or a matrix of the suitable size')

    def matrix_space(self):
        r"""
        Return the underlying matrix space of this endomorphism algebra.
        """
        return self.endomorphism_ring().matrix_space().change_ring(QQ)

    def identity(self):
        r"""
        Return the identity element of this endomorphism algebra.
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
        if M._coerce_map_from_(S):
            return True
        elif S == self:
            return True

    def characteristic(self):
        return ZZ(0)

    def gens(self):
        return tuple(self(x) for x in self.endomorphism_ring().gens())

    def gen(self, i):
        return self.gens()[i]

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

    def _maps_to_field(self):
        r"""
        Return the maps to and from the hecke_eigenvalue_field when simple and
        new.
        """
        A = self.abelian_variety()
        if not A.is_simple():
            return ValueError("self must be simple")

        d = A.dimension()

        # this will very likely work
        for i in range(100):
            phi = self.random_element()
            f = phi.matrix().minpoly()
            if f.degree() == d:
                K = NumberField(f, names='a')
                break
        return f, phi
        K_to_A = K.hom([phi])
        return K, K_to_A

    def free_module(self):
        r"""
        Return this endomorphism ring as a submodule of the free module
        `\QQ^{4d^2}`, where `d` is the dimension of the associated abelian
        variety.
        """
        dim = self.abelian_variety().dimension()
        V = QQ**(4 * dim * dim)
        return V.submodule([V(m.matrix().list()) for m in self.gens()])
