# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import pytest

from oiseau.dg.nodal import (
    Element,
    RefElement,
    RefElementType,
    RefHexahedron,
    RefLine,
    RefQuadrilateral,
    RefTetrahedron,
    RefTriangle,
    get_ref_element,
)


@pytest.mark.parametrize(
    "element_class, kind, node_shape, face_nodes",
    [
        (RefLine, RefElementType.Line, (3,), 1),
        (RefTriangle, RefElementType.Triangle, (6, 2), 3),
        (RefQuadrilateral, RefElementType.Quadrilateral, (9, 2), 3),
        (RefTetrahedron, RefElementType.Tetrahedron, (10, 3), 6),
        (RefHexahedron, RefElementType.Hexahedron, (27, 3), 9),
    ],
)
def test_reference_element_operators(element_class, kind, node_shape, face_nodes):
    ref = element_class(order=2)
    assert isinstance(ref, RefElement)
    assert ref.order == 2
    assert ref.number_of_nodes == node_shape[0]
    assert ref.number_of_face_nodes == face_nodes
    assert ref.r.shape == node_shape

    # A noncontiguous slice exercises evaluation at a subset of points.
    np.testing.assert_allclose(ref.vandermonde(ref.r[::2]), ref.v[::2], atol=1e-12)
    assert isinstance(get_ref_element(kind, order=2), element_class)
    if ref.r.ndim == 1:
        np.testing.assert_allclose(ref.d @ ref.r, 1.0, atol=1e-12)
        np.testing.assert_allclose(ref.d @ ref.v, ref.gv, atol=1e-12)
    else:
        np.testing.assert_allclose(ref.d[:, :, 0] @ ref.r[:, 0], 1.0, atol=1e-12)
        np.testing.assert_allclose(ref.d[:, :, 0] @ ref.v, ref.gv[:, :, 0], atol=1e-12)


def test_nodal_element_owns_its_reference_and_coordinates():
    ref = RefTriangle(order=1)
    nodes = ref.r + 2.0
    element = Element(ref, nodes)
    expected = nodes.copy()

    nodes.fill(0.0)
    del ref

    assert element.order == 1
    assert isinstance(element.reference, RefTriangle)
    result = element.nodes
    np.testing.assert_array_equal(result, expected)
    result.fill(1.0)
    np.testing.assert_array_equal(element.nodes, expected)


def test_nodal_elements_reject_invalid_inputs():
    with pytest.raises(ValueError):
        RefTriangle(order=0)

    ref = RefTriangle(order=1)
    with pytest.raises(ValueError, match="Points must have shape"):
        ref.vandermonde(np.zeros((3, 1)))
    with pytest.raises(ValueError, match="per reference node"):
        Element(ref, np.zeros((2, 2)))
