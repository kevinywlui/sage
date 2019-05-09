r"""
Space of algebras obtained by tensoring an endomorphism ring of a modular
abelian variety by QQ.
"""
from sage.categories.homset import HomsetWithBase
from sage.modules.matrix_morphism import MatrixMorphism
from sage.modules.free_module import FreeModule_submodule_field
from sage.rings.ring import Ring
from sage.rings.all import QQ

from . import morphism


class Morphism(MatrixMorphism):
    r"""
    Morphisms in the isogeny category of modular abelian varieties.
    """
    pass


class EndomorphismAlgebra(HomsetWithBase):
    r"""
    Class of a algebras obtained by tensoring a endomorphism ring by QQ.
    """
    Element = Morphism

    def __init__(self, endomorphism_ring, category=None):
        self._E = endomorphism_ring
        self._E_free = self._E.free_module()
        self._A = self._E.abelian_variety()
        self._ambient_space = self._E_free.ambient_vector_space()

        domain = self._A
        codomain = self._A
        HomsetWithBase.__init__(self, domain, codomain, category=category)

    def __call__(self, M):
        r"""
        Coerce something into this space assuming M is a matrix.
        """


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

        E = self.endomorphism_ring()
        d = A.dimension()

        # this will very likely work
        for i in range(100):
            phi = E.random_element()
            f = phi.matrix().minpoly()
            if f.degree() == d:
                K = NumberField(f, names='a')
                break
        K_to_A = K.hom([phi])
        return K, K_to_A
